"""T-Trend-inspired Student-t acceleration-jump factors + persistent joint GP.

Aravkin/Burke/Pillonetto, SIAM J Control Optim 52(5), 2014,
10.1137/130918861. This is a causal jump+Gaussian-GP specialization, not a
claim to reproduce their full batch algorithm or convergence guarantees.
"""
import time
from dataclasses import dataclass,asdict,replace
import numpy as np
from scipy.linalg import block_diag,cho_factor,cho_solve
from physics import rotate,exp,skew
from target import cartesian
from solver import Parameters,gp_matrices,right_jacobian
from sequential import extend,observation


@dataclass(frozen=True)
class AdaptiveParameters:
    background_jerk_psd: float = .001
    jump_scale_std: float = .02
    degrees_of_freedom: float = 3.
    iterations: int = 3
    starts: tuple = (1.,.001)


def sparse_product(H,P,U,N):
    """Exploit raw observation sparsity without changing joint covariance."""
    m=len(H)//N;groups=H.reshape(N,m,len(P));HP=H[:,U-9:U]@P[U-9:U]
    for i in range(N):
        a=slice(i*m,(i+1)*m);b=slice(U+16*i,U+16*i+16)
        HP[a]+=groups[i,:,b]@P[b]
    S=HP[:,U-9:U]@H[:,U-9:U].T
    for i in range(N):
        a=slice(i*m,(i+1)*m);b=slice(U+16*i,U+16*i+16)
        S[:,a]+=HP[:,b]@groups[i,:,b].T
    S+=np.eye(len(H))
    factor=cho_factor(.5*(S+S.T),lower=True,check_finite=False)
    K=cho_solve(factor,HP,check_finite=False).T
    return K,HP


def sensor_rows(x,U,k,radar,nav,packets,cfg,which):
    r,H=observation(x,U,k,radar,nav,packets,cfg,True)
    m=14 if k else 3;r=r.reshape(cfg.nodes,m);H=H.reshape(cfg.nodes,m,len(x))
    chosen=slice(0,11) if which=='navigation' else slice(m-3,m)
    if which=='navigation':
        # Eight Gaussian pseudoranges carry four local linear directions
        # (position3+clock). QR drops only an x-independent residual component.
        # This is lossless for the once-linearized navigation likelihood.
        columns=np.array([0,1,2,15]);local=np.empty((cfg.nodes,8,4))
        for i in range(cfg.nodes):local[i]=H[i,:8,U+16*i:U+16*i+16][:,columns]
        Q,_=np.linalg.qr(local,mode='reduced')
        rg=np.einsum('nij,ni->nj',Q,r[:,:8]);Hg=np.einsum('nij,nik->njk',Q,H[:,:8])
        return np.concatenate([rg,r[:,8:11]],axis=1).ravel(),np.concatenate([Hg,H[:,8:11]],axis=1).reshape(-1,len(x))
    return r[:,chosen].ravel(),H[:,chosen].reshape(-1,len(x))


def radar_map(prior,Pprior,initial,U,k,radar,nav,packets,cfg,iterations,augmented=False):
    x=initial.copy();first=None
    for _ in range(iterations):
        base=x[:-3] if augmented else x
        r,H=sensor_rows(base,U,k,radar,nav,packets,cfg,'radar')
        if augmented:H=np.c_[H,np.zeros((len(H),3))]
        if first is None:first=float(r@r)
        K,HP=sparse_product(H,Pprior,U,cfg.nodes)
        x=prior+K@(r+H@(x-prior))
    base=x[:-3] if augmented else x
    r,H=sensor_rows(base,U,k,radar,nav,packets,cfg,'radar')
    if augmented:H=np.c_[H,np.zeros((len(H),3))]
    K,HP=sparse_product(H,Pprior,U,cfg.nodes)
    P=Pprior-K@HP
    return x,.5*(P+P.T),first,float(r@r)


def fit(cfg,t,radar,nav,packets,student=True,apar=AdaptiveParameters()):
    # Only previous/current target clones are needed by future likelihoods.
    # Marginalizing older Gaussian clones changes no current posterior; all
    # old information and full current target/nav cross blocks remain.
    T,N=radar.shape[:2];par=replace(Parameters(),window=2,jerk_psd=apar.background_jerk_psd)
    mean=np.zeros(9+N*16);P=np.zeros((len(mean),len(mean)))
    P[:9,:9]=np.diag(np.repeat(np.square([1000.,par.initial_velocity_std,par.initial_acceleration_std]),3));P[9:,9:]=block_diag(*nav['P'][0])
    F,_=gp_matrices(cfg.dt,par.jerk_psd)
    G=np.vstack([.5*cfg.dt**2*np.eye(3),cfg.dt*np.eye(3),np.eye(3)])
    states=[];covs=[];nav_states=[];relcovs=[];seconds=[];cost=[];scales=[]
    nu=apar.degrees_of_freedom;s2=apar.jump_scale_std**2
    for k in range(T):
        tic=time.perf_counter()
        if k:prior,Pbase,U=extend(mean,P,k,nav,cfg,par,True)
        else:prior,Pbase,U=mean.copy(),P.copy(),9
        # Navigation likelihood is linearized once; not suppressed by robust
        # radar weights. Current radar is assimilated only in radar_map.
        if k:
            r,H=sensor_rows(prior,U,k,radar,nav,packets,cfg,'navigation')
            K,HP=sparse_product(H,Pbase,U,N);prior=prior+K@r
            Pbase=Pbase-K@HP;Pbase=.5*(Pbase+Pbase.T)
        if not k:
            initial=prior.copy();initial[:3]=(nav['p'][0]+rotate(nav['R'][0],cartesian(radar[0]))).mean(0)
            mean,P,c0,c1=radar_map(prior,Pbase,initial,U,k,radar,nav,packets,cfg,apar.iterations)
            chosen_scale=1.
        else:
            base_factor=cho_factor(Pbase,lower=True,check_finite=False);candidates=[]
            starts=apar.starts if student else (1.,)
            for start in starts:
                lam=float(start);x=np.r_[prior,np.zeros(3)]
                for _ in range(apar.iterations if student else 1):
                    used=lam
                    variance=s2/lam if student else s2*nu/(nu-2)
                    Pg=np.zeros((len(prior)+3,len(prior)+3));Pg[:-3,:-3]=Pbase
                    Pg[U-9:U,U-9:U]+=variance*(G@G.T)
                    Pg[U-9:U,-3:]=variance*G;Pg[-3:,U-9:U]=variance*G.T;Pg[-3:,-3:]=variance*np.eye(3)
                    aug_prior=np.r_[prior,np.zeros(3)]
                    x,Pcandidate,c0,c1=radar_map(aug_prior,Pg,x,U,k,radar,nav,packets,cfg,1 if student else apar.iterations,True)
                    if student:lam=(nu+3)/(nu+np.sum(x[-3:]**2)/s2)
                shifted=x[:-3].copy();shifted[U-9:U]-=G@x[-3:]
                d=shifted-prior;quadratic=float(d@cho_solve(base_factor,d,check_finite=False))
                process=(nu+3)*np.log1p(np.sum(x[-3:]**2)/(nu*s2)) if student else np.sum(x[-3:]**2)/(s2*nu/(nu-2))
                objective=quadratic+float(process)+c1
                candidates.append((objective,x,Pcandidate,c0,c1,used))
            chosen=min(candidates,key=lambda z:z[0]);_,xx,PP,c0,c1,chosen_scale=chosen
            mean=xx[:-3];P=PP[:-3,:-3]
        state=mean[U-9:U];Pu=P[U-9:U,U-9:U];e=mean[U:].reshape(N,16)
        p=nav['p'][k]+e[:,:3];v=nav['v'][k]+e[:,3:6];R=nav['R'][k]@exp(e[:,6:9]);ba=nav['ba'][k]+e[:,9:12]
        J=np.zeros((N,9,16));J[:,:3,:3]=np.eye(3);J[:,3:6,3:6]=np.eye(3)
        J[:,6:9,6:9]=-R@skew(nav['accel'][k]-ba)@right_jacobian(e[:,6:9]);J[:,6:9,9:12]=-R
        B=np.einsum('nai,nimj,mbj->namb',J,P[U:,U:].reshape(N,16,N,16),J,optimize=True)
        C=np.einsum('nai,nib->nab',J,P[U:,U-9:U].reshape(N,16,9))
        relative=B+Pu[None,:,None,:]-C[:,:,None,:]-C.transpose(2,0,1)[None]
        for i in range(N):relative[i,6:9,i,6:9]+=cfg.accelerometer_density**2/cfg.imu_dt*np.eye(3)
        states.append(state.copy());covs.append(Pu.copy());relcovs.append(relative)
        nav_states.append(dict(p=p.copy(),v=v.copy(),R=R.copy(),ba=ba.copy(),accel=nav['accel'][k],theta=e[:,6:9].copy()))
        scales.append(chosen_scale);cost.append([c0,c1]);seconds.append(time.perf_counter()-tic)
    result=dict(state=np.asarray(states),covariance=np.asarray(covs),nav=nav_states,relative_covariance=np.asarray(relcovs),seconds=np.asarray(seconds),cost=np.asarray(cost),parameters=asdict(par),precision_scale=np.asarray(scales),adaptive_parameters=asdict(apar))
    return result
