"""Joint GP plus direct measured FMCW relative range-rate factors.

Keeps every historical observation exactly once at the current linearization.
Old target clones are marginalized by covariance subblocks, not weak resets.
Linear Gaussian specialization is augmented Kalman smoothing; no contrary claim.
"""
import time
from dataclasses import asdict
import numpy as np
from scipy.linalg import block_diag,cho_factor,cho_solve
from physics import rotate,exp,skew,spherical
from target import cartesian,spherical_jac,wrap
from solver import Parameters,gp_matrices,right_jacobian


def extend(mean,P,k,nav,cfg,par,joint):
    """Append GP state, propagate nav with raw IMU, marginalize oldest clone."""
    N=cfg.nodes;V=16*N if joint else 0;oldU=len(mean)-V;F,Q=gp_matrices(cfg.dt,par.jerk_psd)
    if oldU//9>=par.window:
        mean=mean[9:];P=P[9:,9:];oldU-=9
    newU=oldU+9;D=newU+V
    out=np.empty(D);out[:oldU]=mean[:oldU];out[oldU:newU]=F@mean[oldU-9:oldU]
    C=np.zeros((D,D));C[:oldU,:oldU]=P[:oldU,:oldU]
    C[:oldU,oldU:newU]=P[:oldU,oldU-9:oldU]@F.T;C[oldU:newU,:oldU]=C[:oldU,oldU:newU].T
    C[oldU:newU,oldU:newU]=F@P[oldU-9:oldU,oldU-9:oldU]@F.T+Q
    if joint:
        Fn=nav['model_F'][k];Qn=nav['model_Q'][k];oldn=mean[oldU:].reshape(N,16)
        out[newU:]=(np.einsum('nij,nj->ni',Fn,oldn)+nav['model_offset'][k]).ravel()
        Pnn=P[oldU:,oldU:].reshape(N,16,N,16)
        propagated=np.einsum('nai,nimj,mbj->namb',Fn,Pnn,Fn,optimize=True).reshape(V,V)
        for i in range(N):propagated[16*i:16*i+16,16*i:16*i+16]+=Qn[i]
        C[newU:,newU:]=propagated
        Pnu=P[oldU:,:oldU].reshape(N,16,oldU)
        transformed=np.einsum('nij,nju->niu',Fn,Pnu).reshape(V,oldU)
        C[newU:,:oldU]=transformed;C[:oldU,newU:]=transformed.T
        C[newU:,oldU:newU]=transformed[:,-9:]@F.T;C[oldU:newU,newU:]=C[newU:,oldU:newU].T
    return out,.5*(C+C.T),newU


def base_observation(x,U,k,radar,nav,packets,cfg,joint):
    N=cfg.nodes;D=len(x);u=x[U-9:U]
    stdrad=np.array([cfg.range_std,*([np.deg2rad(cfg.angle_std_deg)]*2)])
    if joint:
        e=x[U:].reshape(N,16);p=nav['p'][k]+e[:,:3];R=nav['R'][k]@exp(e[:,6:9]);Jr=right_jacobian(e[:,6:9])
        rows=14 if k else 3
        Hn=np.zeros((N,rows,16));Hu=np.zeros((N,rows,9));residual=np.empty((N,rows));std=np.empty(rows)
        offset=11 if k else 0
        if k:
            q=p[:,None]-packets['satellites'];distance=np.linalg.norm(q,axis=-1)
            residual[:,:8]=packets['pseudoranges'][k]-distance-nav['clock'][k,:,None]-e[:,None,15]
            Hn[:,:8,:3]=q/distance[...,None];Hn[:,:8,15]=1.;std[:8]=cfg.gnss_range_std
            mag=rotate(R.swapaxes(-1,-2),np.asarray(cfg.magnetic_world))
            residual[:,8:11]=packets['mag'][k]-mag;Hn[:,8:11,6:9]=skew(mag)@Jr;std[8:11]=cfg.mag_std
        body=rotate(R.swapaxes(-1,-2),u[None,:3]-p);B=spherical_jac(body);Hw=B@R.swapaxes(-1,-2)
        residual[:,offset:]=radar[k]-spherical(body);residual[:,offset+1:]=wrap(residual[:,offset+1:])
        Hu[:,offset:,:3]=Hw;Hn[:,offset:,:3]=-Hw;Hn[:,offset:,6:9]=B@skew(body)@Jr;std[offset:]=stdrad
        Hn/=std[None,:,None];Hu/=std[None,:,None];residual/=std
        H=np.zeros((N,rows,D));H[:,:,U-9:U]=Hu
        for i in range(N):H[i,:,U+16*i:U+16*i+16]=Hn[i]
        return residual.ravel(),H.reshape(N*rows,D)
    body=rotate(nav['R'][k].swapaxes(-1,-2),u[None,:3]-nav['p'][k]);B=spherical_jac(body)
    Hw=B@nav['R'][k].swapaxes(-1,-2)
    Jn=np.zeros((N,3,16));Jn[:,:,:3]=-Hw;Jn[:,:,6:9]=B@skew(body)
    Rmeas=np.diag(stdrad**2)+Jn@nav['P'][k]@Jn.swapaxes(-1,-2)
    chol=np.linalg.cholesky(Rmeas)
    residual=radar[k]-spherical(body);residual[:,1:]=wrap(residual[:,1:])
    res=np.linalg.solve(chol,residual[...,None])[...,0]
    H=np.zeros((N,3,U));H[:,:,U-9:U-6]=np.linalg.solve(chol,Hw)
    return res.ravel(),H.reshape(N*3,U)


def observation(x,U,k,radar,nav,packets,cfg,joint,doppler=None,ignore_buoy=False):
    residual,H=base_observation(x,U,k,radar,nav,packets,cfg,joint)
    # Exact QR elimination of GNSS residual directions independent of the
    # state at THIS linearization. Eight pseudoranges have only four local
    # columns (position3 + clock). All useful information is preserved.
    N=cfg.nodes;m=14 if k else 3
    residual=residual.reshape(N,m);H=H.reshape(N,m,len(x))
    if k:
        local=np.empty((N,8,4));columns=np.array([0,1,2,15])
        for i in range(N):local[i]=H[i,:8,U+16*i:U+16*i+16][:,columns]
        Q,_=np.linalg.qr(local,mode='reduced')
        rg=np.einsum('nij,ni->nj',Q,residual[:,:8])
        Hg=np.einsum('nij,nik->njk',Q,H[:,:8])
        residual=np.concatenate([rg,residual[:,8:]],axis=1)
        H=np.concatenate([Hg,H[:,8:]],axis=1)
    if doppler is None:return residual.ravel(),H.reshape(-1,len(x))
    if not joint:raise ValueError('Doppler study requires joint navigation states')
    N=cfg.nodes;e=x[U:].reshape(N,16)
    target=x[U-9:U];bp=nav['p'][k]+e[:,:3]
    bv=np.zeros_like(bp) if ignore_buoy else nav['v'][k]+e[:,3:6]
    q=target[:3]-bp;distance=np.linalg.norm(q,axis=-1);los=q/distance[:,None]
    relative_velocity=target[3:6]-bv
    prediction=np.sum(los*relative_velocity,axis=-1)
    g=(relative_velocity-los*prediction[:,None])/distance[:,None]
    Hd=np.zeros((N,len(x)));Hd[:,U-9:U-6]=g;Hd[:,U-6:U-3]=los
    for i in range(N):
        Hd[i,U+16*i:U+16*i+3]=-g[i]
        if not ignore_buoy:Hd[i,U+16*i+3:U+16*i+6]=-los[i]
    std=doppler['std'][k]
    residual=np.concatenate([residual,((doppler['velocity'][k]-prediction)/std)[:,None]],axis=1)
    H=np.concatenate([H,(Hd/std[:,None])[:,None,:]],axis=1)
    return residual.ravel(),H.reshape(-1,len(x))


def gain(H,P,U,N):
    """Exact sparse products: each row touches target and its own buoy only.

    Covariance P remains FULL; no cross-node or cross-time terms dropped.
    """
    m=len(H)//N;groups=H.reshape(N,m,len(P))
    HP=H[:,U-9:U]@P[U-9:U]
    for i in range(N):
        rows=slice(i*m,(i+1)*m);block=slice(U+16*i,U+16*i+16)
        HP[rows]+=groups[i,:,block]@P[block]
    S=HP[:,U-9:U]@H[:,U-9:U].T
    for i in range(N):
        rows=slice(i*m,(i+1)*m);block=slice(U+16*i,U+16*i+16)
        S[:,rows]+=HP[:,block]@groups[i,:,block].T
    S+=np.eye(len(H));factor=cho_factor(.5*(S+S.T),lower=True,check_finite=False)
    return cho_solve(factor,HP,check_finite=False).T,HP


def fit(cfg,t,radar,nav,par=Parameters(),joint=True,packets=None,doppler=None,ignore_buoy=False):
    if packets is None:raise ValueError('Raw navigation observations required')
    T,N=radar.shape[:2];V=16*N if joint else 0
    mean=np.zeros(9+V);P=np.zeros((9+V,9+V))
    P[:9,:9]=np.diag(np.repeat(np.square([1000.,par.initial_velocity_std,par.initial_acceleration_std]),3))
    if joint:P[9:,9:]=block_diag(*nav['P'][0])
    states=[];covariances=[];nav_states=[];relative_cov=[];seconds=[];cost=[]
    for k in range(T):
        tic=time.perf_counter()
        if k:prior,Pprior,U=extend(mean,P,k,nav,cfg,par,joint)
        else:prior,Pprior,U=mean.copy(),P.copy(),9
        x=prior.copy()
        if not k:x[:3]=(nav['p'][0]+rotate(nav['R'][0],cartesian(radar[0]))).mean(0)
        first=None
        for _ in range(par.iterations):
            residual,H=observation(x,U,k,radar,nav,packets,cfg,joint,doppler,ignore_buoy)
            K,HP=gain(H,Pprior,U,N)
            adjusted=residual+H@(x-prior)
            if first is None:first=float(residual@residual)
            x=prior+K@adjusted
        residual,H=observation(x,U,k,radar,nav,packets,cfg,joint,doppler,ignore_buoy)
        K,HP=gain(H,Pprior,U,N)
        P=Pprior-K@HP;P=.5*(P+P.T);mean=x
        state=x[U-9:U].copy();Pu=P[U-9:U,U-9:U].copy();states.append(state);covariances.append(Pu)
        if joint:
            e=x[U:].reshape(N,16);p=nav['p'][k]+e[:,:3];v=nav['v'][k]+e[:,3:6];R=nav['R'][k]@exp(e[:,6:9]);ba=nav['ba'][k]+e[:,9:12]
            J=np.zeros((N,9,16));J[:,:3,:3]=np.eye(3);J[:,3:6,3:6]=np.eye(3)
            J[:,6:9,6:9]=-R@skew(nav['accel'][k]-ba)@right_jacobian(e[:,6:9]);J[:,6:9,9:12]=-R
            Pnn=P[U:,U:].reshape(N,16,N,16);Pnu=P[U:,U-9:U].reshape(N,16,9)
            B=np.einsum('nai,nimj,mbj->namb',J,Pnn,J,optimize=True)
            C=np.einsum('nai,nib->nab',J,Pnu)
            relative=B+Pu[None,:,None,:]-C[:,:,None,:]-C.transpose(2,0,1)[None]
            # Above transpose encodes Cov(U,B_j) at (axis_i,node_j,axis_j).
            for i in range(N):relative[i,6:9,i,6:9]+=cfg.accelerometer_density**2/cfg.imu_dt*np.eye(3)
            relative_cov.append(relative)
            theta=e[:,6:9]
        else:p=nav['p'][k];v=nav['v'][k];R=nav['R'][k];ba=nav['ba'][k];theta=np.zeros_like(p)
        nav_states.append(dict(p=p.copy(),v=v.copy(),R=R.copy(),ba=ba.copy(),accel=nav['accel'][k],theta=theta.copy()))
        cost.append([first,float(residual@residual)]);seconds.append(time.perf_counter()-tic)
    return dict(state=np.asarray(states),covariance=np.asarray(covariances),nav=nav_states,relative_covariance=np.asarray(relative_cov) if joint else None,seconds=np.asarray(seconds),cost=np.asarray(cost),parameters=asdict(par))


def predict(history,h,cfg,acceleration=True,covariance=True):
    x=history['state'];hs=np.broadcast_to(np.asarray(h),(len(x),));p=x[:,:3]+hs[:,None]*x[:,3:6]+.5*hs[:,None]**2*x[:,6:9]
    buoy=[];distcov=[]
    for k,n in enumerate(history['nav']):
        hh=float(hs[k]);a=rotate(n['R'],n['accel']-n['ba'])+np.array([0.,0.,-9.81]) if acceleration else np.zeros_like(n['p'])
        b=n['p']+hh*n['v']+.5*hh*hh*a;buoy.append(b)
        if covariance:
            u=(p[k]-b);u/=np.linalg.norm(u,axis=-1,keepdims=True)
            G=np.concatenate([u,hh*u,.5*hh*hh*u],axis=-1)
            C=np.einsum('ia,iajb,jb->ij',G,history['relative_covariance'][k],G,optimize=True)
            _,Q=gp_matrices(hh,history['parameters']['jerk_psd']);C+=u@Q[:3,:3]@u.T
            C+=cfg.accelerometer_density**2*hh**3/3*np.eye(cfg.nodes)
            distcov.append(.5*(C+C.T))
    buoy=np.asarray(buoy);cmd=np.linalg.norm(p[:,None]-buoy,axis=-1)
    return p,buoy,cmd,np.asarray(distcov)
