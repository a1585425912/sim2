"""Fixed-weight distributed covariance intersection, not independent summation.

The finite-round kernel is algebraically W**L. No central all-node mean
replaces a finite round. Marginalize omega BEFORE constructing information.
"""
import sys
from pathlib import Path
import numpy as np
IMM_DIR=Path(__file__).resolve().parent.parent/'imm_dpc'
sys.path.insert(0,str(IMM_DIR))
from imm_filter import transition, jacobian, process_covariance, Parameters


def forecast(history, h, par=Parameters()):
    """Frozen-mode EKF forecast moments, including inter-mode spread and Q."""
    means=[]; covs=[]
    for j,model in enumerate(history['enabled']):
        x=history['branches'][...,j,:]
        P=history['branch_covariances'][...,j,:,:]
        F=jacobian(x,h,model,par)
        means.append(transition(x,h,model,par)[...,:6])
        covs.append((F@P@F.swapaxes(-1,-2)+process_covariance(h,model,par))[...,:6,:6])
    means=np.stack(means,axis=-2); covs=np.stack(covs,axis=-3)
    mu=history['probabilities']
    mean=np.sum(mu[...,None]*means,axis=-2)
    delta=means-mean[...,None,:]
    P=np.sum(mu[...,None,None]*(covs+delta[..., :,None]*delta[...,None,:]),axis=-3)
    return mean,.5*(P+P.swapaxes(-1,-2))


def kernel(W, rounds):
    if int(rounds)!=rounds or rounds<0: raise ValueError('rounds must be a nonnegative integer')
    W=np.asarray(W)
    if W.ndim!=2 or W.shape[0]!=W.shape[1] or np.any(W<0) or not np.allclose(W.sum(1),1):
        raise ValueError('Need nonnegative row-stochastic network weights')
    return np.linalg.matrix_power(W,int(rounds))


def average(A, C, matrix=False):
    return np.einsum('ij,...jkl->...ikl',C,A,optimize=True) if matrix else np.einsum('ij,...jk->...ik',C,A,optimize=True)


def fuse(x, P, C, method):
    """x(...,N,6), P(...,N,6,6), C(N,N) convex finite-round weights."""
    if method=='local': return x.copy(),P.copy()
    if method=='dpc':
        # Unknown cross-correlation conservative bound; NOT sum(P)/N**2.
        return average(x,C),average(P,C,matrix=True)
    if method!='ci': raise ValueError(method)
    eye=np.broadcast_to(np.eye(P.shape[-1]),P.shape)
    Y=np.linalg.solve(P,eye)
    y=np.einsum('...ij,...j->...i',Y,x)
    Yf=average(Y,C,matrix=True); yf=average(y,C)
    Pf=np.linalg.solve(Yf,eye)
    xf=np.linalg.solve(Yf,yf[...,None])[...,0]
    return xf,.5*(Pf+Pf.swapaxes(-1,-2))
