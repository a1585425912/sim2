"""Raw range/azimuth/elevation IMM update with estimated navigation uncertainty."""
import numpy as np
from scipy.linalg import expm
from scipy.special import logsumexp
from imm_filter import Parameters,mix,transition,jacobian,process_covariance
from physics import rotate,skew,spherical


def wrap(a):return (a+np.pi)%(2*np.pi)-np.pi


def cartesian(z):
    r,a,e=np.moveaxis(z,-1,0)
    return np.stack([r*np.cos(e)*np.cos(a),r*np.cos(e)*np.sin(a),r*np.sin(e)],-1)


def cart_jac(z):
    r,a,e=np.moveaxis(z,-1,0);ca,sa,ce,se=np.cos(a),np.sin(a),np.cos(e),np.sin(e)
    return np.stack([np.stack([ce*ca,-r*ce*sa,-r*se*ca],-1),np.stack([ce*sa,r*ce*ca,-r*se*sa],-1),np.stack([se,np.zeros_like(r),r*ce],-1)],-2)


def spherical_jac(q):
    x,y,z=np.moveaxis(q,-1,0);r=np.linalg.norm(q,axis=-1);h=np.hypot(x,y)
    return np.stack([np.stack([x/r,y/r,z/r],-1),np.stack([-y/h**2,x/h**2,np.zeros_like(r)],-1),np.stack([-x*z/(h*r*r),-y*z/(h*r*r),h/r**2],-1)],-2)


def fit_radar(t,measurements,nav,radar_cov,enabled=(0,1,2),par=Parameters()):
    T,N=measurements.shape[:2];M=len(enabled)
    body=cartesian(measurements[0]);world=nav['p'][0]+rotate(nav['R'][0],body)
    J=np.zeros((N,3,16));J[:,:,:3]=np.eye(3);J[:,:,6:9]=-nav['R'][0]@skew(body)
    G=nav['R'][0]@cart_jac(measurements[0])
    initial_cov=G@radar_cov@G.swapaxes(-1,-2)+J@nav['P'][0]@J.swapaxes(-1,-2)
    x=np.zeros((N,M,7));x[...,:3]=world[:,None]
    P=np.zeros((N,M,7,7));P[...,:3,:3]=initial_cov[:,None]
    P[...,3:6,3:6]=par.initial_velocity_std**2*np.eye(3);P[...,6,6]=par.omega_initial_std**2
    mu=np.full((N,M),1/M);xs=np.empty((T,N,M,7));ps=np.empty((T,N,M,7,7));mus=np.empty((T,N,M))
    xs[0],ps[0],mus[0]=x,P,mu
    for k in range(1,T):
        dt=t[k]-t[k-1];Gmode=np.full((M,M),par.switching_rate);np.fill_diagonal(Gmode,-(M-1)*par.switching_rate)
        x0,P0,prior=mix(x,P,mu,expm(dt*Gmode));ll=np.empty_like(mu)
        Rbw=nav['R'][k].swapaxes(-1,-2)
        for j,mode in enumerate(enabled):
            F=jacobian(x0[:,j],dt,mode,par);xp=transition(x0[:,j],dt,mode,par)
            Pp=F@P0[:,j]@F.swapaxes(-1,-2)+process_covariance(dt,mode,par)
            local=rotate(Rbw,xp[:,:3]-nav['p'][k]);pred=spherical(local);B=spherical_jac(local)
            H=np.zeros((N,3,7));H[:,:,:3]=B@Rbw
            Jnav=np.zeros((N,3,16));Jnav[:,:,:3]=-B@Rbw;Jnav[:,:,6:9]=B@skew(local)
            R=radar_cov+Jnav@nav['P'][k]@Jnav.swapaxes(-1,-2)
            innovation=measurements[k]-pred;innovation[:,1:]=wrap(innovation[:,1:])
            S=H@Pp@H.swapaxes(-1,-2)+R;L=np.linalg.cholesky(S)
            white=np.linalg.solve(L,innovation[...,None])[...,0]
            ll[:,j]=-.5*(np.sum(white**2,-1)+2*np.log(np.diagonal(L,axis1=-2,axis2=-1)).sum(-1)+3*np.log(2*np.pi))
            K=np.linalg.solve(S,H@Pp).swapaxes(-1,-2);x[:,j]=xp+np.einsum('nij,nj->ni',K,innovation)
            J=np.eye(7)-K@H;P[:,j]=J@Pp@J.swapaxes(-1,-2)+K@R@K.swapaxes(-1,-2)
        lw=np.log(np.maximum(prior,np.finfo(float).tiny))+ll;mu=np.exp(lw-logsumexp(lw,axis=1,keepdims=True))
        P=.5*(P+P.swapaxes(-1,-2));xs[k],ps[k],mus[k]=x,P,mu
    state=np.sum(mus[...,None]*xs,axis=2);delta=xs-state[:,:,None]
    cov=np.sum(mus[...,None,None]*(ps+delta[..., :,None]*delta[...,None,:]),axis=2)
    return dict(branches=xs,branch_covariances=ps,probabilities=mus,state=state,covariance=cov,enabled=enabled)
