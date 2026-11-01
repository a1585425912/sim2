"""Simulator ONLY: smooth shared sea, analytic kinematics, raw sensor packets.

Adapted spectral/response mathematics from archived dfpc_experiment/maritime.py.
No import of archived algorithm or navigation-output simulators.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.spatial.transform import Rotation


@dataclass(frozen=True)
class Config:
    nodes: int = 20
    duration: float = 90.
    frame_budget_s: float = .2
    dt: float = .2
    imu_dt: float = .02
    hs: float = 1.
    tp: float = 6.
    wave_components: int = 48
    wave_gamma: float = 3.3
    sea_depth: float = 50.
    mean_current: tuple = (.3, .1, 0.)
    slow_velocity_std: float = .05
    gnss_range_std: float = .30
    accelerometer_density: float = .015
    gyro_density: float = .0003
    accelerometer_bias_std: float = .03
    gyro_bias_std: float = .001
    accelerometer_bias_rw: float = .0003
    gyro_bias_rw: float = .00001
    clock_stationary_std: float = 30.
    clock_tau: float = 30.
    clock_filter_rw: float = 10.
    mag_std: float = .15
    magnetic_world: tuple = (20., 0., -45.)
    range_std: float = .8
    angle_std_deg: float = .35
    processing_s: float = .012
    rounds: int = 25
    round_s: float = .005
    slot_s: float = .04
    slot_samples: int = 9
    total_power: float = 1.
    frequencies_mhz: tuple = (300., 1000.)

    def horizons(self):
        return self.processing_s+self.rounds*self.round_s+np.linspace(0,self.slot_s,self.slot_samples)


def skew(v):
    v=np.asarray(v);S=np.zeros(v.shape[:-1]+(3,3))
    S[...,0,1]=-v[...,2];S[...,0,2]=v[...,1]
    S[...,1,0]=v[...,2];S[...,1,2]=-v[...,0]
    S[...,2,0]=-v[...,1];S[...,2,1]=v[...,0]
    return S


def exp(v):
    return Rotation.from_rotvec(np.asarray(v).reshape(-1,3)).as_matrix().reshape(np.asarray(v).shape[:-1]+(3,3))


def rotate(R,v):
    return np.einsum('...ij,...j->...i',R,v)


def spherical(v):
    x,y,z=np.moveaxis(v,-1,0)
    return np.stack([np.linalg.norm(v,axis=-1),np.arctan2(y,x),np.arctan2(z,np.hypot(x,y))],-1)


class SeaTruth:
    """Continuous finite-spectral slow current + shared directional JONSWAP sea.

    Slow random coefficients are sampled once, not IID position jumps. All
    queries share one differentiable realization, including TX slot queries.
    Response coefficients are a stated surrogate, not a fitted buoy RAO.
    """
    def __init__(self,cfg,seed):
        self.cfg=cfg;rng=np.random.default_rng(seed)
        az=2*np.pi*np.arange(cfg.nodes)/cfg.nodes
        self.base=np.c_[100*np.cos(az),100*np.sin(az),np.full(cfg.nodes,1.5)]
        self.slow_w=2*np.pi/np.array([180.,120.,80.])
        scale=cfg.slow_velocity_std/np.sqrt(3)
        self.slow_a=rng.normal(0,scale,(cfg.nodes,3,2))
        self.slow_b=rng.normal(0,scale,(cfg.nodes,3,2))
        self.yaw=rng.uniform(-np.pi,np.pi,cfg.nodes)
        edges=np.linspace(.6*2*np.pi/cfg.tp,3*2*np.pi/cfg.tp,cfg.wave_components+1)
        self.w=edges[:-1]+np.diff(edges)*rng.uniform(.2,.8,cfg.wave_components)
        wp=2*np.pi/cfg.tp;sg=np.where(self.w<=wp,.07,.09)
        shape=self.w**-5*np.exp(-1.25*(wp/self.w)**4)*cfg.wave_gamma**np.exp(-(self.w-wp)**2/(2*sg**2*wp**2))
        spectral=shape*np.diff(edges)
        amp=np.sqrt(2*spectral/spectral.sum()*cfg.hs**2/16)
        direction=rng.normal(np.deg2rad(30),np.deg2rad(30),len(self.w))
        self.direction=np.c_[np.cos(direction),np.sin(direction)]
        k=self.w**2/9.81
        for _ in range(16):
            kh=k*cfg.sea_depth;th=np.tanh(kh)
            k-= (9.81*k*th-self.w**2)/(9.81*(th+kh*(1-th**2)))
        self.k=k[:,None]*self.direction
        self.phase=rng.uniform(-np.pi,np.pi,len(self.w))
        ratio=self.w/(2*np.pi/8.)
        heave=1/(1-ratio**2+2j*.2*ratio)
        self.Hp=np.c_[.7*amp*self.direction[:,0]*(-1j),.7*amp*self.direction[:,1]*(-1j),amp*heave]
        self.Ha=np.c_[-k*amp*self.direction[:,1],k*amp*self.direction[:,0],.03*k*amp].astype(complex)

    def at(self,times):
        t=np.asarray(times,float);sw=self.slow_w
        cs=np.cos(t[...,None]*sw);sn=np.sin(t[...,None]*sw)
        c=self.base+t[...,None,None]*np.asarray(self.cfg.mean_current)
        c=c.copy();c[...,:2]+=np.einsum('...m,nmj->...nj',sn/sw,self.slow_a)+np.einsum('...m,nmj->...nj',(1-cs)/sw,self.slow_b)
        v=np.broadcast_to(np.asarray(self.cfg.mean_current),c.shape).copy()
        v[...,:2]+=np.einsum('...m,nmj->...nj',cs,self.slow_a)+np.einsum('...m,nmj->...nj',sn,self.slow_b)
        a=np.zeros_like(c)
        a[...,:2]=np.einsum('...m,nmj->...nj',-sw*sn,self.slow_a)+np.einsum('...m,nmj->...nj',sw*cs,self.slow_b)
        phase=t[...,None,None]*self.w-np.einsum('...nj,mj->...nm',c[...,:2],self.k)+self.phase
        encounter=self.w-np.einsum('...nj,mj->...nm',v[...,:2],self.k)
        encounter_dot=-np.einsum('...nj,mj->...nm',a[...,:2],self.k)
        z=np.exp(1j*phase)
        wave=np.real(np.einsum('...nm,mj->...nj',z,self.Hp))
        velocity=v+np.real(np.einsum('...nm,mj->...nj',1j*encounter*z,self.Hp))
        acceleration=a+np.real(np.einsum('...nm,mj->...nj',(-encounter**2+1j*encounter_dot)*z,self.Hp))
        euler=np.real(np.einsum('...nm,mj->...nj',z,self.Ha));euler[...,2]+=self.yaw
        edot=np.real(np.einsum('...nm,mj->...nj',1j*encounter*z,self.Ha))
        R=Rotation.from_euler('xyz',euler.reshape(-1,3)).as_matrix().reshape(euler.shape[:-1]+(3,3))
        roll,pitch=euler[...,0],euler[...,1];rd,pd,yd=np.moveaxis(edot,-1,0)
        omega=np.stack([rd-yd*np.sin(pitch),pd*np.cos(roll)+yd*np.sin(roll)*np.cos(pitch),-pd*np.sin(roll)+yd*np.cos(roll)*np.cos(pitch)],-1)
        return dict(p=c+wave,v=velocity,a=acceleration,R=R,omega=omega,c=c,wave=wave)


def satellite_positions():
    az=np.deg2rad(np.arange(8)*45+12)
    el=np.deg2rad([25,50,35,70,30,60,40,55])
    return 2.1e7*np.c_[np.cos(el)*np.cos(az),np.cos(el)*np.sin(az),np.sin(el)]


def sensor_record(cfg,scene,seed):
    """Only this simulator receives truth. IMU packet describes prior interval.

    GNSS output is eight pseudoranges, with a latent drifting receiver clock.
    Range noise is a corrected-code idealization, not a carrier-ambiguity model.
    """
    ss=np.random.SeedSequence(seed).spawn(5)
    imu_rng,gnss_rng,mag_rng,bias_rng,clock_rng=[np.random.default_rng(s) for s in ss]
    t=np.arange(round(cfg.duration/cfg.imu_dt)+1)*cfg.imu_dt
    mid=np.r_[0.,.5*(t[1:]+t[:-1])]
    imu_truth=scene.at(mid);frame_step=round(cfg.dt/cfg.imu_dt)
    frame_t=t[::frame_step];frame_truth=scene.at(frame_t)
    nt=len(t);n=cfg.nodes;dt=cfg.imu_dt
    ba=cfg.accelerometer_bias_std*bias_rng.normal(size=(n,3))
    bg=cfg.gyro_bias_std*bias_rng.normal(size=(n,3))
    ba=ba[None]+np.cumsum(cfg.accelerometer_bias_rw*np.sqrt(dt)*bias_rng.normal(size=(nt,n,3)),axis=0)
    bg=bg[None]+np.cumsum(cfg.gyro_bias_rw*np.sqrt(dt)*bias_rng.normal(size=(nt,n,3)),axis=0)
    f=rotate(imu_truth['R'].swapaxes(-1,-2),imu_truth['a']-np.array([0.,0.,-9.81]))
    accel=f+ba+cfg.accelerometer_density/np.sqrt(dt)*imu_rng.normal(size=f.shape)
    gyro=imu_truth['omega']+bg+cfg.gyro_density/np.sqrt(dt)*imu_rng.normal(size=f.shape)
    nf=len(frame_t);rho=np.exp(-cfg.dt/cfg.clock_tau)
    clock=np.empty((nf,n));clock[0]=cfg.clock_stationary_std*clock_rng.normal(size=n)
    for k in range(1,nf):clock[k]=rho*clock[k-1]+cfg.clock_stationary_std*np.sqrt(1-rho*rho)*clock_rng.normal(size=n)
    sat=satellite_positions()
    ranges=np.linalg.norm(sat[None,None]-frame_truth['p'][:,:,None],axis=-1)+clock[...,None]
    ranges+=cfg.gnss_range_std*gnss_rng.normal(size=ranges.shape)
    mag=rotate(frame_truth['R'].swapaxes(-1,-2),np.asarray(cfg.magnetic_world))
    mag+=cfg.mag_std*mag_rng.normal(size=mag.shape)
    packets=dict(imu_time=t,frame_time=frame_t,accel=accel,gyro=gyro,pseudoranges=ranges,mag=mag,satellites=sat)
    evaluation=dict(frame=frame_truth,bias_accel=ba[::frame_step],bias_gyro=bg[::frame_step],clock=clock)
    return packets,evaluation
