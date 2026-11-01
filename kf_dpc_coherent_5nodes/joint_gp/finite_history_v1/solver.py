"""Causal finite-history joint MAP with WNOJ GP prior and raw radar factors.

Theory: Anderson et al., Autonomous Robots 2015; Tang et al., RA-L 2019;
joint sensor/target MAP and Schur elimination, Taylor et al., IPSN 2006.
No target posterior is reinserted as a prior over reused window measurements.
"""
from dataclasses import dataclass,asdict
import time
import numpy as np
from scipy.linalg import cho_factor,cho_solve,block_diag
from physics import exp,skew,rotate,spherical
from target import spherical_jac,cartesian,wrap
from navigation import window_prior


@dataclass(frozen=True)
class Parameters:
    window: int = 11
    iterations: int = 3
    jerk_psd: float = 1.
    initial_position_std: float = 10000.
    initial_velocity_std: float = 100.
    initial_acceleration_std: float = 10.


def gp_matrices(h,q):
    F=np.kron(np.array([[1.,h,.5*h*h],[0.,1.,h],[0.,0.,1.]]),np.eye(3))
    Q=q*np.kron(np.array([[h**5/20,h**4/8,h**3/6],[h**4/8,h**3/3,h*h/2],[h**3/6,h*h/2,h]]),np.eye(3))
    return F,Q


def motion_precision(W,dt,par):
    H=np.zeros((9*W,9*W));F,Q=gp_matrices(dt,par.jerk_psd);Y=np.linalg.inv(Q)
    H[:9,:9]=np.diag(np.repeat(1/np.square([par.initial_position_std,par.initial_velocity_std,par.initial_acceleration_std]),3))
    for k in range(W-1):
        a=slice(9*k,9*k+9);b=slice(9*k+9,9*k+18)
        H[a,a]+=F.T@Y@F;H[b,b]+=Y;H[a,b]-=F.T@Y;H[b,a]-=Y@F
    return H


def right_jacobian(theta):
    a=np.linalg.norm(theta,axis=-1);K=skew(theta)
    a2=a*a;den=np.where(a>1e-6,a2,1.)
    c1=np.where(a>1e-6,(1-np.cos(a))/den,.5-a2/24)
    c2=np.where(a>1e-6,(a-np.sin(a))/np.where(a>1e-6,a*a2,1.),1/6-a2/120)
    return np.eye(3)-c1[...,None,None]*K+c2[...,None,None]*(K@K)


def build_system(x,e,radar,nav_p,nav_R,nav_P,prior,Hmotion,Rdiag,joint):
    """Block arrow normal matrix, no truth or future observations."""
    W,N=radar.shape[:2];U=9*W;V=6*W
    if joint:
        pose=e.reshape(N,W,6).transpose(1,0,2)
        p=nav_p+pose[...,:3];R=nav_R@exp(pose[...,3:])
    else:p=nav_p;R=nav_R
    Rbw=R.swapaxes(-1,-2)
    local=rotate(Rbw,x[:,None,:3]-p);pred=spherical(local)
    B=spherical_jac(local);Hw=B@Rbw
    residual=radar-pred;residual[...,1:]=wrap(residual[...,1:])
    Hpose=np.concatenate([-Hw,B@skew(local)],axis=-1)
    if joint:Hpose[...,3:]=Hpose[...,3:]@right_jacobian(pose[...,3:])
    if joint:Wmeas=np.broadcast_to(np.diag(1/Rdiag),(W,N,3,3))
    else:
        posecov=nav_P[:,:,[0,1,2,6,7,8]][:,:,:,[0,1,2,6,7,8]]
        effective=np.diag(Rdiag)+Hpose@posecov@Hpose.swapaxes(-1,-2)
        Wmeas=np.linalg.inv(effective)
    Huu=Hmotion.copy();gu=-Hmotion@x.ravel()
    Hu=Hw.swapaxes(-1,-2)@Wmeas
    Hpp=np.sum(Hu@Hw,axis=1);gp=np.sum(np.einsum('tnij,tnj->tni',Hu,residual),axis=1)
    for k in range(W):
        s=slice(9*k,9*k+3);Huu[s,s]+=Hpp[k];gu[s]+=gp[k]
    if not joint:return Huu,gu,None,None,None,float(x.ravel()@Hmotion@x.ravel()+np.einsum('tni,tnij,tnj->',residual,Wmeas,residual))
    Hnn=prior['precision'].copy();gn=np.einsum('nij,nj->ni',Hnn,prior['mean']-e)
    Hun=np.zeros((N,U,V));Hn=Hpose.swapaxes(-1,-2)@Wmeas
    Hnp=Hn@Hpose;gnp=np.einsum('tnij,tnj->tni',Hn,residual)
    cross=Hu@Hpose
    for k in range(W):
        ns=slice(6*k,6*k+6);us=slice(9*k,9*k+3)
        Hnn[:,ns,ns]+=Hnp[k];gn[:,ns]+=gnp[k];Hun[:,us,ns]=cross[k]
    diff=e-prior['mean']
    cost=float(x.ravel()@Hmotion@x.ravel()+np.einsum('tni,tnij,tnj->',residual,Wmeas,residual)+np.einsum('ni,nij,nj->',diff,prior['precision'],diff))
    return Huu,gu,Hnn,gn,Hun,cost


def eliminate(Huu,gu,Hnn,gn,Hun):
    if Hnn is None:
        factor=cho_factor(.5*(Huu+Huu.T),lower=True,check_finite=False)
        return cho_solve(factor,gu,check_finite=False),None,factor,None,None
    N,V=gn.shape;U=len(gu);Vsol=np.empty((N,V,U));nsol=np.empty_like(gn);factors=[]
    for i in range(N):
        factor=cho_factor(.5*(Hnn[i]+Hnn[i].T),lower=True,check_finite=False);factors.append(factor)
        both=cho_solve(factor,np.c_[Hun[i].T,gn[i]],check_finite=False)
        Vsol[i]=both[:,:U];nsol[i]=both[:,U]
    S=Huu-np.einsum('nuv,nvw->uw',Hun,Vsol,optimize=True)
    b=gu-np.einsum('nuv,nv->u',Hun,nsol,optimize=True)
    S=.5*(S+S.T);factor=cho_factor(S,lower=True,check_finite=False)
    dx=cho_solve(factor,b,check_finite=False)
    de=nsol-np.einsum('nvu,u->nv',Vsol,dx)
    return dx,de,factor,Vsol,factors


def fit(cfg,t,radar,nav,par=Parameters(),joint=True):
    T,N=radar.shape[:2];std=np.array([cfg.range_std,*([np.deg2rad(cfg.angle_std_deg)]*2)])
    states=[];covariances=[];nav_states=[];posterior=[];seconds=[];costs=[]
    xprev=None;eprev=None
    for k in range(T):
        tic=time.perf_counter();start=max(0,k-par.window+1);W=k-start+1
        prior=window_prior(nav,start,k) if joint else None
        if xprev is None:
            world=nav['p'][0]+rotate(nav['R'][0],cartesian(radar[0]));x=np.zeros((1,9));x[0,:3]=world.mean(0)
        else:
            F,_=gp_matrices(cfg.dt,par.jerk_psd)
            x=np.vstack([xprev,(F@xprev[-1])[None]])[-W:]
        e=prior['mean'].copy() if joint else None
        Hmotion=motion_precision(W,cfg.dt,par)
        args=(radar[start:k+1],nav['p'][start:k+1],nav['R'][start:k+1],nav['P'][start:k+1],prior,Hmotion,std**2,joint)
        first_cost=None
        for _ in range(par.iterations):
            system=build_system(x,e,*args)
            if first_cost is None:first_cost=system[-1]
            dx,de,*_=eliminate(*system[:5]);x+=dx.reshape(W,9)
            if joint:e+=de
        system=build_system(x,e,*args);_,_,factor,Vsol,factors=eliminate(*system[:5])
        Puu=cho_solve(factor,np.eye(9*W),check_finite=False)
        states.append(x[-1].copy());covariances.append(Puu[-9:,-9:]);costs.append([first_cost,system[-1]])
        if joint:
            delta=np.einsum('niv,nv->ni',prior['restore'],e-prior['mean'])
            p=nav['p'][k]+delta[:,:3];v=nav['v'][k]+delta[:,3:6]
            R=nav['R'][k]@exp(delta[:,6:9]);ba=nav['ba'][k]+delta[:,9:12]
            Cconditional=[]
            Lfull=-prior['restore']@Vsol
            for i in range(N):
                Cconditional.append(prior['conditional'][i]+prior['restore'][i]@cho_solve(factors[i],prior['restore'][i].T,check_finite=False))
            posterior.append(dict(Puu=Puu,nav_shared=Lfull,nav_conditional=np.asarray(Cconditional)))
        else:
            p=nav['p'][k];v=nav['v'][k];R=nav['R'][k];ba=nav['ba'][k]
            posterior.append(dict(Puu=Puu,nav_shared=np.zeros((N,16,9*W)),nav_conditional=nav['P'][k]))
        nav_states.append(dict(p=p.copy(),v=v.copy(),R=R.copy(),ba=ba.copy(),accel=nav['accel'][k],theta=delta[:,6:9].copy() if joint else np.zeros_like(p)))
        seconds.append(time.perf_counter()-tic);xprev=x.copy()
    return dict(state=np.asarray(states),covariance=np.asarray(covariances),nav=nav_states,posterior=posterior,seconds=np.asarray(seconds),cost=np.asarray(costs),parameters=asdict(par))


def predict(history,h,cfg,acceleration=True,covariance=True):
    """Joint future geometry and distance covariance. No future sensor data."""
    state=history['state'];hs=np.broadcast_to(np.asarray(h),(len(state),))
    p_u=state[:,:3]+hs[:,None]*state[:,3:6]+.5*hs[:,None]**2*state[:,6:9]
    buoy=[];distance_cov=[]
    for k,(n,post) in enumerate(zip(history['nav'],history['posterior'])):
        h=float(hs[k])
        a=rotate(n['R'],n['accel']-n['ba'])+np.array([0.,0.,-9.81]) if acceleration else np.zeros_like(n['p'])
        p_b=n['p']+h*n['v']+.5*h*h*a;buoy.append(p_b)
        if not covariance:continue
        F,Q=gp_matrices(h,history['parameters']['jerk_psd'])
        q=p_u[k]-p_b;u=q/np.linalg.norm(q,axis=-1,keepdims=True)
        U=post['Puu'].shape[-1];J_u=np.zeros((3,U));J_u[:,-9:]=F[:3]
        J_n=np.zeros((cfg.nodes,3,16));J_n[:,:,:3]=np.eye(3);J_n[:,:,3:6]=h*np.eye(3)
        if acceleration:
            J_n[:,:,6:9]=-.5*h*h*n['R']@skew(n['accel']-n['ba'])@right_jacobian(n['theta']);J_n[:,:,9:12]=-.5*h*h*n['R']
        shared=J_n@post['nav_shared'];D=np.einsum('ni,nij->nj',u,J_u[None]-shared)
        independent=J_n@post['nav_conditional']@J_n.swapaxes(-1,-2)
        independent+=cfg.accelerometer_density**2*h**3/3*np.eye(3)
        if acceleration:independent+=(.25*h**4*cfg.accelerometer_density**2/cfg.imu_dt)*np.eye(3)
        C=D@post['Puu']@D.T+u@Q[:3,:3]@u.T
        C+=np.diag(np.einsum('ni,nij,nj->n',u,independent,u))
        distance_cov.append(.5*(C+C.T))
    buoy=np.asarray(buoy);cmd=np.linalg.norm(p_u[:,None]-buoy,axis=-1)
    return p_u,buoy,cmd,np.asarray(distance_cov)
