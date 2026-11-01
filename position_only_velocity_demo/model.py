"""Sensor geometry and causal estimators. No truth trajectory is available here.

World frame: right-handed ENU (east, north, up), metres and seconds.
Radar frame: right-handed; R maps radar coordinates into world coordinates.
Measurement: [range_m, azimuth_rad, elevation_rad]. Azimuth is atan2(y,x),
elevation is atan2(z,hypot(x,y)), positive above the radar xy plane.
State order: [px, py, pz, vx, vy, vz]. No Doppler/radial-velocity input.
"""
from __future__ import annotations

import numpy as np


def wrap_angle(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def spherical(q):
    q = np.asarray(q, dtype=float)
    x, y, z = np.moveaxis(q, -1, 0)
    return np.stack((np.linalg.norm(q, axis=-1), np.arctan2(y, x),
                     np.arctan2(z, np.hypot(x, y))), axis=-1)


def cartesian(z):
    r, a, e = np.moveaxis(np.asarray(z, dtype=float), -1, 0)
    return np.stack((r*np.cos(e)*np.cos(a), r*np.cos(e)*np.sin(a),
                     r*np.sin(e)), axis=-1)


def cartesian_jacobian(z):
    r, a, e = z
    return np.array([
        [np.cos(e)*np.cos(a), -r*np.cos(e)*np.sin(a), -r*np.sin(e)*np.cos(a)],
        [np.cos(e)*np.sin(a), r*np.cos(e)*np.cos(a), -r*np.sin(e)*np.sin(a)],
        [np.sin(e), 0., r*np.cos(e)]])


def spherical_jacobian(q):
    x, y, z = q
    rho2 = x*x + y*y
    r2 = rho2 + z*z
    if rho2 < 1e-12 or r2 < 1e-12:
        raise ValueError('Spherical measurement is singular on the radar z axis.')
    rho, r = np.sqrt(rho2), np.sqrt(r2)
    return np.array([[x/r, y/r, z/r], [-y/rho2, x/rho2, 0.],
                     [-x*z/(rho*r2), -y*z/(rho*r2), rho/r2]])


def measurement_prediction(state, radar_position, radar_rotation):
    """Actual nonlinear EKF h(x) and dh/dx; velocity columns are zero."""
    relative = radar_rotation.T @ (state[:3] - radar_position)
    h = np.zeros((3, 6))
    h[:, :3] = spherical_jacobian(relative) @ radar_rotation.T
    return spherical(relative), h


def convert_measurements(z, radar_positions, radar_rotations, sensor_cov):
    """World xyz and first-order covariance propagation at each noisy z.

    Navigation is assumed exact. The converted errors are only approximately
    zero-mean Gaussian; the linear KF does not remove this conversion bias.
    """
    p = radar_positions + np.einsum('nij,nj->ni', radar_rotations, cartesian(z))
    cov = np.empty((len(z), 3, 3))
    for k in range(len(z)):
        j = radar_rotations[k] @ cartesian_jacobian(z[k])
        cov[k] = j @ sensor_cov @ j.T
    return p, cov


def cv_matrices(dt, acceleration_psd):
    """Continuous white acceleration PSD q [m^2/s^3], not accel std.

    Q = q * [[dt^3/3 I, dt^2/2 I], [dt^2/2 I, dt I]].
    """
    if not np.isfinite(dt) or dt <= 0 or not np.isfinite(acceleration_psd) or acceleration_psd < 0:
        raise ValueError('dt must be positive and acceleration_psd nonnegative.')
    i = np.eye(3)
    f = np.eye(6)
    f[:3, 3:] = dt*i
    q = acceleration_psd * np.block([[dt**3/3*i, dt**2/2*i],
                                      [dt**2/2*i, dt*i]])
    return f, q


def validate_positions(t, p):
    t, p = np.asarray(t, dtype=float), np.asarray(p, dtype=float)
    if t.ndim != 1 or len(t) < 2 or p.shape != (len(t), 3):
        raise ValueError('Expected time (N,) and position (N,3), N>=2.')
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(p)) or np.any(np.diff(t) <= 0):
        raise ValueError('Inputs must be finite, with strictly increasing times.')
    return t, p


def difference_velocity(t, positions):
    t, positions = validate_positions(t, positions)
    v = np.full_like(positions, np.nan)
    v[1:] = np.diff(positions, axis=0) / np.diff(t)[:, None]
    # First-frame velocity is undefined. It is excluded explicitly in metrics.
    # At k>=1 this is an interval-average velocity assigned to the current frame.
    return v


def position_filter(t, positions, covariances, acceleration_psd=0.5,
                    initial_velocity_std=20.):
    """Linear CV Kalman filter using only xyz position observations."""
    return _filter(t, positions, covariances, acceleration_psd, initial_velocity_std)


def radar_ekf(t, measurements, radar_positions, radar_rotations, sensor_cov,
              acceleration_psd=0.5, initial_velocity_std=20.):
    """True EKF: nonlinear range/azimuth/elevation innovation + analytic H(x).

    xyz conversion is used once for initialization. Subsequent updates use raw
    spherical observations. Neither truth position nor truth velocity is input.
    """
    p, cov = convert_measurements(measurements, radar_positions, radar_rotations, sensor_cov)
    return _filter(t, p, cov, acceleration_psd, initial_velocity_std,
                   (measurements, radar_positions, radar_rotations, sensor_cov))


def _filter(t, positions, covariances, acceleration_psd, initial_velocity_std, radar=None):
    t, positions = validate_positions(t, positions)
    if covariances.shape != (len(t), 3, 3) or not np.all(np.isfinite(covariances)):
        raise ValueError('Expected finite position covariance (N,3,3).')
    if not np.isfinite(initial_velocity_std) or initial_velocity_std <= 0:
        raise ValueError('initial_velocity_std must be positive.')
    if np.min(np.linalg.eigvalsh(covariances)) <= 0:
        raise ValueError('Position covariances must be positive definite.')
    x = np.r_[positions[0], np.zeros(3)]
    p = np.zeros((6, 6))
    p[:3, :3] = covariances[0]
    p[3:, 3:] = initial_velocity_std**2 * np.eye(3)
    states, posterior_cov = np.empty((len(t), 6)), np.empty((len(t), 6, 6))
    states[0], posterior_cov[0] = x, p
    h_linear = np.c_[np.eye(3), np.zeros((3, 3))]
    # Frame 0 has already initialized position. Do not update twice with it.
    for k in range(1, len(t)):
        f, q = cv_matrices(t[k] - t[k-1], acceleration_psd)
        x, p = f @ x, f @ p @ f.T + q
        if radar is None:
            h, r = h_linear, covariances[k]
            innovation = positions[k] - h @ x
        else:
            z, radar_p, radar_r, r = radar
            predicted, h = measurement_prediction(x, radar_p[k], radar_r[k])
            innovation = z[k] - predicted
            innovation[1:] = wrap_angle(innovation[1:])
        s = h @ p @ h.T + r
        gain = np.linalg.solve(s, (p @ h.T).T).T
        x = x + gain @ innovation
        ikh = np.eye(6) - gain @ h
        p = ikh @ p @ ikh.T + gain @ r @ gain.T  # Joseph form
        p = (p + p.T) / 2
        states[k], posterior_cov[k] = x, p
    if not np.all(np.isfinite(states)) or not np.all(np.isfinite(posterior_cov)):
        raise FloatingPointError('Non-finite filter output.')
    return states, posterior_cov
