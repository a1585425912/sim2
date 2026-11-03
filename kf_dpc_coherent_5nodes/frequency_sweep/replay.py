"""Evaluate arbitrary positive MHz values on ALL completed geometry caches."""
from pathlib import Path
import argparse,csv,json
import numpy as np

HERE=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mhz',type=float,nargs='+',required=True)
    ap.add_argument('--include-legacy',action='store_true',help='Include all eight completed legacy strategies')
    ap.add_argument('--output',type=Path,default=HERE/'outputs/custom_frequencies.csv');args=ap.parse_args()
    if any(f<=0 for f in args.mhz):raise ValueError('Frequency must be positive; f=0 is only an analytic limit')
    chunks=[HERE/'outputs/chunks'/str(x) for x in (20261028,20261033,20261038,20261043)]
    for p in chunks:
        if json.loads((p/'manifest.json').read_text())['status']!='completed':raise ValueError('Wait for full completion')
    files=sorted(p for chunk in chunks for p in (chunk/'geometry').glob('*.npz'))
    if len(files)!=400:raise ValueError('Incomplete full-trial geometry')
    if args.include_legacy:
        legacy=[HERE/'outputs/legacy_chunks'/str(x) for x in (20261028,20261033,20261038,20261043)]
        for p in legacy:
            if json.loads((p/'manifest.json').read_text())['status']!='completed':raise ValueError('Wait for complete legacy models')
        more=sorted(p for chunk in legacy for p in (chunk/'geometry').glob('*.npz'))
        if len(more)!=1600:raise ValueError('Incomplete legacy geometry')
        files+=more
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fields=['seed','hs','trajectory','method','timing','phase','frequency_mhz','eta_energy','eta_frame_mean','eta_frame_p05','gain_db']
    with args.output.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.writer(f);writer.writerow(fields)
        for file in files:
            seed,hs,kind,*method_parts=file.stem.split('_');method='_'.join(method_parts)
            with np.load(file) as d:
                t=d['time_s'];masks={'all':np.ones(len(t),bool),'late75_90':(t>=75)&(t<90)}
                for timing in ('fixed_delay','recorded_mean_cpu'):
                    if timing+'__error' not in d.files:continue
                    e=d[timing+'__error'];a=1/d[timing+'__distance'];ideal=a.sum(-1)**2;inc=(a*a).sum(-1)
                    for frequency in args.mhz:
                        rx=np.abs((a*np.exp(1j*2*np.pi*frequency*1e6/299792458.*e)).sum(-1))**2
                        frame=rx.sum(-1)/ideal.sum(-1)
                        for phase,mask in masks.items():
                            writer.writerow([int(seed),float(hs[2:]),kind,method,timing,phase,frequency,
                                float(rx[mask].sum()/ideal[mask].sum()),float(frame[mask].mean()),
                                float(np.quantile(frame[mask],.05)),float(10*np.log10(rx[mask].sum()/inc[mask].sum()))])
    print(str(args.output))

if __name__=='__main__':main()
