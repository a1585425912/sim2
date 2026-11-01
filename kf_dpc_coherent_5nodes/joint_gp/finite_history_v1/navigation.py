"""Time-correlated Gaussian navigation factors from raw sensors, no truth input.

Filtered navigation is a linearization reference, not independent observations.
Backward Gaussian conditionals preserve pose/time correlations for each window.
"""
import numpy as np
from scipy.spatial.transform import Rotation
from physics import skew,exp
from estimation import InertialNavigation

POSE=np.array([0,1,2,6,7,8])


def coordinates(nav):
    x=np.zeros((len(nav.p),16));x[:,:3]=nav.p;x[:,3:6]=nav.v
    x[:,9:12]=nav.ba;x[:,12:15]=nav.bg;x[:,15]=nav.clock
    return x


def navigation_factors(cfg,packets):
    nav=InertialNavigation(cfg,packets);step=round(cfg.dt/cfg.imu_dt)
    result={k:[] for k in ('p','v','R','P','omega','ba','bg','clock','accel')}
    n=cfg.nodes;Ftotal=np.broadcast_to(np.eye(16),(n,16,16)).copy()
    backward_B=[];backward_a=[];backward_D=[]
    for k in range(len(packets['imu_time'])):
        if k:
            dt=cfg.imu_dt;w=packets['gyro'][k]-nav.bg;f=packets['accel'][k]-nav.ba
            Rmid=nav.R@exp(.5*dt*w)
            A=np.zeros((n,16,16));A[:,:3,3:6]=np.eye(3)
            A[:,3:6,6:9]=-Rmid@skew(f);A[:,3:6,9:12]=-Rmid
            A[:,6:9,6:9]=-skew(w);A[:,6:9,12:15]=-np.eye(3)
            F=np.eye(16)+dt*A+.5*dt**2*(A@A);Ftotal=F@Ftotal
            nav.propagate(packets['accel'][k],packets['gyro'][k],dt)
        if k%step:continue
        j=k//step
        if k:
            prior_coordinates=coordinates(nav);prior_R=nav.R.copy();prior_P=nav.P.copy()
            nav.gnss_update(packets['pseudoranges'][j],packets['satellites'])
            nav.mag_update(packets['mag'][j])
            delta=coordinates(nav)-prior_coordinates
            delta[:,6:9]=Rotation.from_matrix(prior_R.swapaxes(-1,-2)@nav.R).as_rotvec()
            T=np.broadcast_to(np.eye(16),(n,16,16)).copy()
            K=skew(delta[:,6:9]);T[:,6:9,6:9]+= .5*K+(K@K)/12
            Aback=np.linalg.solve(prior_P,Ftotal@result['P'][-1]).swapaxes(-1,-2)
            D=result['P'][-1]-Aback@prior_P@Aback.swapaxes(-1,-2)
            backward_B.append(Aback@T);backward_a.append(np.einsum('nij,nj->ni',Aback,delta))
            backward_D.append(.5*(D+D.swapaxes(-1,-2)))
        values=(nav.p,nav.v,nav.R,nav.P,packets['gyro'][k]-nav.bg,nav.ba,nav.bg,nav.clock,packets['accel'][k])
        for name,value in zip(result,values):result[name].append(value.copy())
        Ftotal=np.broadcast_to(np.eye(16),(n,16,16)).copy()
    output={k:np.asarray(v) for k,v in result.items()}
    output.update(B=np.asarray(backward_B),a_back=np.asarray(backward_a),D=np.asarray(backward_D))
    return output


def window_prior(nav,start,end):
    """p(pose[start:end+1] | raw navigation[0:end]), one compressed factor.

    Returns conditional restoration of full current navigation from pose errors.
    Coordinates remain at filtered references, so no covariance is silently
    transplanted between different orientation tangents.
    """
    N=nav['p'].shape[1];W=end-start+1
    mean=np.zeros((N,W,16));C=np.zeros((N,W,16,W,16))
    C[:,-1,:,-1,:]=nav['P'][end]
    for j in range(W-2,-1,-1):
        k=start+j;B=nav['B'][k]
        mean[:,j]=nav['a_back'][k]+np.einsum('nij,nj->ni',B,mean[:,j+1])
        C[:,j,:,j,:]=nav['D'][k]+B@C[:,j+1,:,j+1,:]@B.swapaxes(-1,-2)
        cross=np.einsum('nij,njtk->nitk',B,C[:,j+1,:,j+1:,:])
        C[:,j,:,j+1:,:]=cross
        C[:,j+1:,:,j,:]=cross.transpose(0,2,3,1)
    selected=np.take(np.take(C,POSE,axis=2),POSE,axis=4)
    covariance=selected.reshape(N,6*W,6*W)
    covariance=.5*(covariance+covariance.swapaxes(-1,-2))+1e-12*np.eye(6*W)
    precision=np.linalg.solve(covariance,np.broadcast_to(np.eye(6*W),covariance.shape))
    cross_current=np.take(C[:,-1,:,:,:],POSE,axis=3).reshape(N,16,6*W)
    restore=cross_current@precision
    conditional=nav['P'][end]-restore@cross_current.swapaxes(-1,-2)
    return dict(mean=mean[:,:,POSE].reshape(N,6*W),covariance=covariance,precision=precision,
                restore=restore,conditional=.5*(conditional+conditional.swapaxes(-1,-2)))
