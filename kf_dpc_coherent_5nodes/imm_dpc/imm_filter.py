"""Causal S/CV/CT IMM-EKF. Inputs contain observations only; SI units.

Common state: [px,py,pz,vx,vy,vz,omega]. S is a soft near-hover
prior, not an externally supplied zero-velocity observation.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.linalg import expm
from scipy.special import logsumexp

NAMES = ('S', 'CV', 'CT')


@dataclass(frozen=True)
class Parameters:
    beta: float = 1.0
    q_hover: float = .01
    q_cv: float = .5
    q_ct: float = .1
    q_omega: float = .0004
    omega_initial_std: float = .25
    switching_rate: float = .05  # per second to EACH other model
    initial_velocity_std: float = 20.


def transition(x, h, model, par=Parameters()):
    """Deterministic branch means, including stable CT zero-rate limit."""
    y = x.copy()
    if model == 0:
        rho = np.exp(-par.beta*h)
        y[..., :3] += (-np.expm1(-par.beta*h)/par.beta)*x[..., 3:6]
        y[..., 3:6] *= rho
        y[..., 6] *= np.exp(-h)
    elif model == 1:
        y[..., :3] += h*x[..., 3:6]
        y[..., 6] *= np.exp(-h)
    elif model == 2:
        a = x[..., 6]*h
        # np.sinc is analytic at zero, including for complex-step derivatives.
        A = h*np.sinc(a/np.pi)
        B = .5*h*a*np.sinc(a/(2*np.pi))**2
        c, s = np.cos(a), np.sin(a)
        vx, vy = x[..., 3], x[..., 4]
        y[..., 0] += A*vx-B*vy
        y[..., 1] += B*vx+A*vy
        y[..., 2] += h*x[..., 5]
        y[..., 3] = c*vx-s*vy
        y[..., 4] = s*vx+c*vy
    else:
        raise ValueError(model)
    return y


def jacobian(x, h, model, par=Parameters()):
    F = np.empty(x.shape[:-1]+(7, 7))
    for j in range(7):
        z = x.astype(complex)
        z[..., j] += 1e-20j
        F[..., :, j] = transition(z, h, model, par).imag/1e-20
    return F


def process_covariance(h, model, par):
    Q = np.zeros((7, 7))
    if model == 0:
        # Exact integrated OU velocity covariance.
        r = np.exp(-par.beta*h); b = par.beta; q = par.q_hover
        pp = q/b**2*(h-2*(1-r)/b+(1-r*r)/(2*b))
        pv = q*(1-r)**2/(2*b*b)
        vv = q*(1-r*r)/(2*b)
    else:
        # Explicit discrete model: deterministic CT then independent integrated
        # white-acceleration increment in world axes. Not exact continuous CT Q.
        q = par.q_cv if model == 1 else par.q_ct
        pp, pv, vv = q*h**3/3, q*h*h/2, q*h
    Q[:3, :3] = pp*np.eye(3)
    Q[:3, 3:6] = Q[3:6, :3] = pv*np.eye(3)
    Q[3:6, 3:6] = vv*np.eye(3)
    Q[6, 6] = par.q_omega*(h if model == 2 else -np.expm1(-2*h)/2)
    return Q


def mix(x, P, mu, Pi):
    prior = mu@Pi
    weight = mu[:, :, None]*Pi[None]/prior[:, None, :]
    mixed = np.einsum('nab,nad->nbd', weight, x)
    delta = x[:, :, None, :]-mixed[:, None, :, :]
    cov = np.einsum('nab,nabij->nbij', weight,
                    P[:, :, None]+delta[..., :, None]*delta[..., None, :])
    return mixed, cov, prior


def fit(t, z, R, par=Parameters(), enabled=(0, 1, 2), retain_branch_covariances=False):
    """Batched independent nodes; no consensus feedback. Returns branch history."""
    t, z, R = np.asarray(t), np.asarray(z), np.asarray(R)
    if len(t) < 2 or not np.all(np.diff(t) > 0):
        raise ValueError('Need increasing timestamps')
    if z.ndim != 3 or z.shape[-1] != 3 or R.shape != z.shape[:-1]+(3, 3):
        raise ValueError('Expected z(T,N,3), R(T,N,3,3)')
    nt, nn, _ = z.shape; nm = len(enabled)
    x = np.zeros((nn, nm, 7)); x[..., :3] = z[0, :, None]
    P = np.zeros((nn, nm, 7, 7)); P[..., :3, :3] = R[0, :, None]
    P[..., 3:6, 3:6] = par.initial_velocity_std**2*np.eye(3)
    P[..., 6, 6] = par.omega_initial_std**2
    mu = np.full((nn, nm), 1/nm)
    xs = np.empty((nt, nn, nm, 7)); ps = np.empty((nt, nn, nm, 7, 7))
    mus = np.empty((nt, nn, nm)); nis = np.zeros((nt, nn, nm))
    xs[0], ps[0], mus[0] = x, P, mu
    for k in range(1, nt):
        h = t[k]-t[k-1]
        G = np.full((nm, nm), par.switching_rate)
        np.fill_diagonal(G, -(nm-1)*par.switching_rate)
        Pi = expm(h*G)
        x0, P0, prior = mix(x, P, mu, Pi)
        ll = np.empty_like(mu)
        for j, model in enumerate(enabled):
            F = jacobian(x0[:, j], h, model, par)
            xp = transition(x0[:, j], h, model, par)
            Pp = F@P0[:, j]@F.swapaxes(-1, -2)+process_covariance(h, model, par)
            v = z[k]-xp[:, :3]; S = Pp[:, :3, :3]+R[k]
            chol = np.linalg.cholesky(S)
            whitened = np.linalg.solve(chol, v[..., None])[..., 0]
            ni = np.sum(whitened**2, -1)
            ll[:, j] = -.5*(ni+2*np.log(np.diagonal(chol, axis1=-2, axis2=-1)).sum(-1)+3*np.log(2*np.pi))
            K = np.linalg.solve(S, Pp[:, :3, :]).swapaxes(-1, -2)
            x[:, j] = xp+np.einsum('nij,nj->ni', K, v)
            A = np.broadcast_to(np.eye(7), Pp.shape).copy(); A[..., :3] -= K
            P[:, j] = A@Pp@A.swapaxes(-1, -2)+K@R[k]@K.swapaxes(-1, -2)
            nis[k, :, j] = ni
        logw = np.log(np.maximum(prior, np.finfo(float).tiny))+ll
        mu = np.exp(logw-logsumexp(logw, axis=1, keepdims=True))
        P = .5*(P+P.swapaxes(-1, -2))
        xs[k], ps[k], mus[k] = x, P, mu
    state = np.sum(mus[..., None]*xs, axis=2)
    delta = xs-state[:, :, None]
    covariance = np.sum(mus[..., None, None]*(ps+delta[..., :, None]*delta[..., None, :]), axis=2)
    result = dict(branches=xs, probabilities=mus, state=state, covariance=covariance,
                  nis=nis, enabled=enabled, parameters=asdict(par))
    if retain_branch_covariances:
        result['branch_covariances'] = ps
    return result


def predict(history, h, par=Parameters()):
    """Frozen-mode short-horizon mean, branch first then mixture; no future data."""
    return sum(history['probabilities'][..., j, None]*
               transition(history['branches'][..., j, :], h, model, par)[..., :6]
               for j, model in enumerate(history['enabled']))
