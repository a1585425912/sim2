"""Estimators: accept sensor packets only, never SeaTruth or evaluation arrays."""
import numpy as np
from scipy.linalg import expm
from physics import skew,exp,rotate


def pseudorange_fix(ranges,satellites,sigma):
    n=len(ranges);x=np.zeros((n,4))
    for _ in range(6):
        q=x[:,None,:3]-satellites[None];d=np.linalg.norm(q,axis=-1)
        H=np.concatenate([q/d[...,None],np.ones(d.shape+(1,))],axis=-1)
        residual=ranges-d-x[:,None,3]
        x+=np.linalg.solve(H.swapaxes(-1,-2)@H,(H.swapaxes(-1,-2)@residual[...,None]))[...,0]
    P=np.linalg.inv(H.swapaxes(-1,-2)@H)*sigma**2
    return x,P


class InertialNavigation:
    """Right-error ESKF [dp,dv,dtheta,dba,dbg,dc_clock], 16 dimensions.

    Position is solved from pseudoranges. Initial tilt/heading use measured
    specific force and magnetic vector, with broad attitude uncertainty.
    No true pose, true velocity, or true biases enter initialization.
    """
    def __init__(self,cfg,packet):
        self.cfg=cfg;n=cfg.nodes
        fix,Pfix=pseudorange_fix(packet['pseudoranges'][0],packet['satellites'],cfg.gnss_range_std)
        self.p=fix[:,:3];self.clock=fix[:,3];self.v=np.zeros((n,3))
        up=packet['accel'][0].copy();up/=np.linalg.norm(up,axis=-1,keepdims=True)
        north=packet['mag'][0]-up*np.sum(packet['mag'][0]*up,axis=-1,keepdims=True)
        north/=np.linalg.norm(north,axis=-1,keepdims=True)
        east=np.cross(up,north)
        self.R=np.stack([north,east,up],axis=-1).swapaxes(-1,-2)
        self.ba=np.zeros((n,3));self.bg=np.zeros((n,3))
        self.P=np.zeros((n,16,16))
        self.P[:,:3,:3]=Pfix[:,:3,:3];self.P[:,:3,15]=Pfix[:,:3,3];self.P[:,15,:3]=Pfix[:,3,:3];self.P[:,15,15]=Pfix[:,3,3]
        for start,var in ((3,1.),(6,np.deg2rad(10)**2),(9,cfg.accelerometer_bias_std**2),(12,cfg.gyro_bias_std**2)):
            self.P[:,start:start+3,start:start+3]=var*np.eye(3)

    def propagate(self,accel,gyro,dt):
        c=self.cfg;f=accel-self.ba;w=gyro-self.bg
        Rmid=self.R@exp(.5*dt*w)
        acceleration=rotate(Rmid,f)+np.array([0.,0.,-9.81])
        self.p+=dt*self.v+.5*dt**2*acceleration;self.v+=dt*acceleration
        n=len(self.p);A=np.zeros((n,16,16));A[:,:3,3:6]=np.eye(3)
        A[:,3:6,6:9]=-Rmid@skew(f);A[:,3:6,9:12]=-Rmid
        A[:,6:9,6:9]=-skew(w);A[:,6:9,12:15]=-np.eye(3)
        F=np.eye(16)+dt*A+.5*dt**2*(A@A)
        Q=np.zeros_like(A);qa=c.accelerometer_density**2
        Q[:,:3,:3]=qa*dt**3/3*np.eye(3);Q[:,:3,3:6]=qa*dt**2/2*np.eye(3);Q[:,3:6,:3]=Q[:,:3,3:6]
        Q[:,3:6,3:6]=qa*dt*np.eye(3)
        Q[:,6:9,6:9]=c.gyro_density**2*dt*np.eye(3)
        Q[:,9:12,9:12]=c.accelerometer_bias_rw**2*dt*np.eye(3)
        Q[:,12:15,12:15]=c.gyro_bias_rw**2*dt*np.eye(3)
        Q[:,15,15]=c.clock_filter_rw**2*dt
        self.P=F@self.P@F.swapaxes(-1,-2)+Q
        self.R=self.R@exp(dt*w)

    def update(self,residual,H,measurement_covariance):
        S=H@self.P@H.swapaxes(-1,-2)+measurement_covariance
        K=np.linalg.solve(S,H@self.P).swapaxes(-1,-2)
        delta=np.einsum('nij,nj->ni',K,residual)
        self.p+=delta[:,:3];self.v+=delta[:,3:6];self.R=self.R@exp(delta[:,6:9])
        self.ba+=delta[:,9:12];self.bg+=delta[:,12:15];self.clock+=delta[:,15]
        J=np.eye(16)-K@H
        self.P=J@self.P@J.swapaxes(-1,-2)+K@measurement_covariance@K.swapaxes(-1,-2)
        reset=np.broadcast_to(np.eye(16),self.P.shape).copy();reset[:,6:9,6:9]-=.5*skew(delta[:,6:9])
        self.P=reset@self.P@reset.swapaxes(-1,-2);self.P=.5*(self.P+self.P.swapaxes(-1,-2))

    def gnss_update(self,ranges,satellites):
        q=self.p[:,None]-satellites;d=np.linalg.norm(q,axis=-1)
        H=np.zeros((len(self.p),len(satellites),16));H[:,:,:3]=q/d[...,None];H[:,:,15]=1.
        self.update(ranges-d-self.clock[:,None],H,self.cfg.gnss_range_std**2*np.eye(len(satellites)))

    def mag_update(self,mag):
        predicted=rotate(self.R.swapaxes(-1,-2),np.asarray(self.cfg.magnetic_world))
        H=np.zeros((len(self.p),3,16));H[:,:,6:9]=skew(predicted)
        self.update(mag-predicted,H,self.cfg.mag_std**2*np.eye(3))


def navigate(cfg,packets):
    nav=InertialNavigation(cfg,packets);step=round(cfg.dt/cfg.imu_dt)
    result={k:[] for k in ('p','v','R','P','omega','ba','bg','clock')}
    for k,t in enumerate(packets['imu_time']):
        if k:nav.propagate(packets['accel'][k],packets['gyro'][k],cfg.imu_dt)
        if k%step:continue
        j=k//step
        if k:
            nav.gnss_update(packets['pseudoranges'][j],packets['satellites'])
            nav.mag_update(packets['mag'][j])
        values=(nav.p,nav.v,nav.R,nav.P,packets['gyro'][k]-nav.bg,nav.ba,nav.bg,nav.clock)
        for name,value in zip(result,values):result[name].append(value.copy())
    return {k:np.asarray(v) for k,v in result.items()}


class WavePredictor:
    """Causal slow/fast decomposition fitted to estimated navigation, not truth.

    Per-axis [center,center_velocity,wave,wave_velocity]. This secondary
    filter consumes temporally correlated navigation posteriors; its nominal
    covariance is not an independently calibrated joint posterior. Only means
    are used for TX. Current UAV confidence is evaluated separately.
    """
    def __init__(self,cfg):
        self.cfg=cfg;self.w=2*np.pi/cfg.tp;self.z=.15
        self.A=np.array([[0.,1.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,1.],[0.,0.,-self.w**2,-2*self.z*self.w]])
        self.H=np.array([[1.,0.,1.,0.],[0.,1.,0.,1.]])

    def matrices(self,h):
        F=expm(self.A*h);spectral=np.zeros((3,4,4))
        sig=self.cfg.hs/4*np.array([.7,.7,1.])
        spectral[:,1,1]=1e-5;spectral[:,3,3]=4*self.z*self.w**3*sig**2
        Q=[]
        for S in spectral:
            M=np.zeros((8,8));M[:4,:4]=-self.A;M[:4,4:]=S;M[4:,4:]=self.A.T
            E=expm(M*h);q=F@E[:4,4:];Q.append(.5*(q+q.T))
        return F,np.asarray(Q)

    def fit(self,nav):
        T,N=nav['p'].shape[:2];x=np.zeros((N,3,4));x[...,0]=nav['p'][0];x[...,1]=nav['v'][0]
        P=np.zeros((N,3,4,4));wave_var=(self.cfg.hs/4*np.array([.7,.7,1.]))**2
        for axis in range(3):
            S=np.diag([wave_var[axis],self.w**2*wave_var[axis]])
            P[:,axis,:2,:2]=nav['P'][0][:,[axis,axis+3]][:,:,[axis,axis+3]]+S
            P[:,axis,2:,2:]=S;P[:,axis,:2,2:]=-S;P[:,axis,2:,:2]=-S
        xs=[x.copy()];F,Q=self.matrices(self.cfg.dt)
        for k in range(1,T):
            x=x@F.T;P=F@P@F.T+Q
            measurement=np.stack([nav['p'][k],nav['v'][k]],axis=-1)
            R=np.empty((N,3,2,2))
            for axis in range(3):R[:,axis]=nav['P'][k][:,[axis,axis+3]][:,:,[axis,axis+3]]
            S=self.H@P@self.H.T+R+1e-10*np.eye(2)
            K=np.linalg.solve(S,self.H@P).swapaxes(-1,-2)
            x+=np.einsum('naij,naj->nai',K,measurement-x@self.H.T)
            J=np.eye(4)-K@self.H;P=J@P@J.swapaxes(-1,-2)+K@R@K.swapaxes(-1,-2)
            xs.append(x.copy())
        return np.asarray(xs)

    def predict(self,history,h):
        predicted=history@expm(self.A*h).T
        return predicted[...,0]+predicted[...,2]
