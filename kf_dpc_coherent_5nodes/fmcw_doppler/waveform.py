"""Single-body FMCW dechirped complex ADC simulation, not truth + velocity noise.

Independent supplementary FMCW channel: existing radar range/angle data stay
unchanged. No truth used to choose a range or Doppler cell. Known one target,
global 2-D spectral peak, all frames retained, no CFAR or measurement gating.
"""
from dataclasses import dataclass,asdict
import time
import numpy as np
from scipy.fft import fft,fftshift

C=299792458.

@dataclass(frozen=True)
class Waveform:
    carrier_hz: float=77e9
    bandwidth_hz: float=5e6
    chirp_s: float=40e-6
    repetition_s: float=50e-6
    chirps: int=256
    samples: int=64
    doppler_fft: int=512
    reference_range_m: float=250.
    reference_adc_snr_db: float=5.
    # Conservative fixed likelihood scale, not calibrated against truth.
    velocity_std: float=.15

    @property
    def cpi_s(self):return self.chirps*self.repetition_s
    @property
    def half_span_s(self):return .5*((self.chirps-1)*self.repetition_s+self.chirp_s)
    @property
    def wavelength(self):return C/self.carrier_hz
    def metadata(self):
        return dict(**asdict(self),cpi_s=self.cpi_s,half_span_s=self.half_span_s,
                    velocity_fft_resolution=self.wavelength/(2*self.cpi_s),
                    velocity_unambiguous=self.wavelength/(4*self.repetition_s),
                    range_unambiguous_m=C*self.samples/(2*self.bandwidth_hz))


def buoy_bursts(scene,t,w):
    """Simulator-only continuous sea samples at chirp centers, once per sea."""
    slow=(np.arange(w.chirps)-(w.chirps-1)/2)*w.repetition_s
    p=[];v=[]
    for ti in t:
        truth=scene.at(ti+slow)
        p.append(truth['p']);v.append(truth['v'])
    return np.asarray(p),np.asarray(v)


def observe(cfg,t,kind,trajectory,buoy_p,buoy_v,seed,w=Waveform(),keep=False):
    """Truth appears only in forward propagation; output estimator gets ADC peaks.

    Stop-and-linearize within each 40us chirp, exact trajectory at each chirp
    center. All slow-time acceleration/wave motion is retained. Complex IQ,
    thermal white noise, point-body free-space amplitude proportional to r^-2.
    Zero rotor, sea clutter, phase noise and mutual interference: declared scope.
    """
    rng=np.random.default_rng(seed);lam=w.wavelength;S=w.bandwidth_hz/w.chirp_s
    slow=(np.arange(w.chirps)-(w.chirps-1)/2)*w.repetition_s
    fast=(np.arange(w.samples)+.5)*w.chirp_s/w.samples-w.chirp_s/2
    # Dechirping phase uses fast time measured from chirp start.
    fast_chirp=fast+w.chirp_s/2
    win=np.hanning(w.chirps)[:,None]*np.hanning(w.samples)[None,:]
    noise_std=np.sqrt(10**(-w.reference_adc_snr_db/10)/2)
    velocity=[];truth_center=[];truth_average=[];seconds=[];snr=[];sample={}
    for k,ti in enumerate(t):
        # Continuous physical states only generate the received signal.
        pu,vu=trajectory(ti+slow,kind,cfg)
        q=pu[:,None]-buoy_p[k];r=np.linalg.norm(q,axis=-1)
        vr=np.sum(q*(vu[:,None]-buoy_v[k]),axis=-1)/r
        rr=r[...,None]+vr[...,None]*fast
        delay=2*rr/C
        phase=2*np.pi*((w.carrier_hz+S*fast_chirp)*delay-.5*S*delay**2)
        amplitude=(w.reference_range_m/rr)**2
        adc=amplitude*np.exp(1j*phase)
        adc+=noise_std*(rng.normal(size=adc.shape)+1j*rng.normal(size=adc.shape))
        tic=time.perf_counter()
        # Axes slow-time, buoy, fast-time -> buoy, slow, fast for FFT.
        spectrum=fftshift(fft(fft(adc.transpose(1,0,2)*win,axis=-1),n=w.doppler_fft,axis=1),axes=1)
        power=np.abs(spectrum)**2
        peak=np.argmax(power.reshape(cfg.nodes,-1),axis=-1)
        kd,kr=np.divmod(peak,w.samples);nodes=np.arange(cfg.nodes)
        left=np.log(np.maximum(power[nodes,(kd-1)%w.doppler_fft,kr],1e-300))
        mid=np.log(np.maximum(power[nodes,kd,kr],1e-300))
        right=np.log(np.maximum(power[nodes,(kd+1)%w.doppler_fft,kr],1e-300))
        correction=.5*(left-right)/(left-2*mid+right)
        frequency=(kd-w.doppler_fft/2+correction)/(w.doppler_fft*w.repetition_s)
        velocity.append(.5*lam*frequency);seconds.append(time.perf_counter()-tic)
        pc,vc=trajectory(np.asarray(ti),kind,cfg)
        # Simulator-only center reference interpolates sea burst states at t.
        # Direct central query is provided by the caller through continuous sea
        # samples: two central chirps symmetrically bracket the nominal epoch.
        bp=.5*(buoy_p[k,w.chirps//2-1]+buoy_p[k,w.chirps//2])
        bv=.5*(buoy_v[k,w.chirps//2-1]+buoy_v[k,w.chirps//2])
        qc=pc-bp;truth_center.append(np.sum(qc*(vc-bv),axis=-1)/np.linalg.norm(qc,axis=-1))
        truth_average.append(np.mean(vr,axis=0))
        snr.append(w.reference_adc_snr_db+40*np.log10(w.reference_range_m/np.mean(r,axis=0)))
        if keep and k in (0,100,325):
            sample[f'frame{k}_adc']=adc[:,0].astype(np.complex64)
            sample[f'frame{k}_range_doppler']=power[0].astype(np.float32)
    velocity=np.asarray(velocity);truth_center=np.asarray(truth_center)
    error=velocity-truth_center
    return dict(velocity=velocity,std=np.full_like(velocity,w.velocity_std),
                seconds=np.asarray(seconds),truth_velocity=truth_center,
                truth_burst_mean=np.asarray(truth_average),error=error,
                snr_db=np.asarray(snr),sample=sample,metadata=w.metadata())
