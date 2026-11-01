"""Truth-free estimators, arithmetic DPC, and transmit commands.

World coordinates are right-handed X/Y/Z (metres); Z is up.
Radar angles are radians: azimuth atan2(Y,X), elevation above the XY plane.
No function in this module receives UAV truth or a Doppler observation.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from position_only_velocity_demo.model import convert_measurements, position_filter
from dfpc_experiment.sim_core import weights_from_adjacency, consensus_linear


def ring_weights(n=5):
    if n < 3:
        raise ValueError('A ring needs at least three nodes.')
    adjacency = np.zeros((n, n), dtype=bool)
    idx = np.arange(n)
    adjacency[idx, (idx + 1) % n] = True
    adjacency[(idx + 1) % n, idx] = True
    return weights_from_adjacency(adjacency)


def local_estimates(time_s, measurements, nodes, rotations, sensor_cov,
                    acceleration_psd=.5, initial_velocity_std=20.):
    """Five independent causal KFs, no consensus feedback into their memories.

    measurements: (T,N,3), nodes: (N,3), rotations: (N,3,3).
    sensor_cov: (3,3), units m^2 / rad^2. All frames refer to the same time.
    """
    nt, nn, dim = measurements.shape
    if dim != 3 or len(time_s) != nt or nodes.shape != (nn, 3):
        raise ValueError('Bad measurement/time/node shape.')
    observed = np.empty((nt, nn, 3))
    obs_cov = np.empty((nt, nn, 3, 3))
    states = np.empty((nt, nn, 6))
    post_cov = np.empty((nt, nn, 6, 6))
    for i in range(nn):
        observed[:, i], obs_cov[:, i] = convert_measurements(
            measurements[:, i], np.broadcast_to(nodes[i], (nt, 3)),
            np.broadcast_to(rotations[i], (nt, 3, 3)), sensor_cov)
        states[:, i], post_cov[:, i] = position_filter(
            time_s, observed[:, i], obs_cov[:, i], acceleration_psd,
            initial_velocity_std)
    return observed, obs_cov, states, post_cov


def dpc(states, weights, rounds):
    """Freeze local posteriors at each frame, then perform L neighbour rounds.

    Every node retains its own finite-round result. No hidden central average.
    The original Metropolis/DPC functions are reused from dfpc_experiment.
    """
    states = np.asarray(states, dtype=float)
    if states.ndim != 3 or states.shape[-1] != 6 or rounds < 0 or int(rounds) != rounds:
        raise ValueError('Expected states (T,N,6) and nonnegative integer rounds.')
    nt, nn, _ = states.shape
    mixed = consensus_linear(weights, states.transpose(1, 0, 2).reshape(nn, -1), rounds)
    return mixed.reshape(nn, nt, 6).transpose(1, 0, 2).copy()


def distance_commands(states, nodes, delay_s, slot_times_s, use_velocity=True):
    """Compute geometry phase commands from ESTIMATES and known navigation.

    Predict position to transmission start by p_tx=p_hat+delay*v_hat.
    Within the short slot use d_hat(s)=d_hat(0)+s*(u_hat dot v_hat).
    The latter is an estimated geometric distance rate, never radar Doppler.
    Setting use_velocity=False holds current estimated position throughout,
    an ablation for the value of using the estimated velocity.
    Returns (T,S,N) predicted distances, in metres.
    """
    if delay_s < 0 or np.any(np.asarray(slot_times_s) < 0):
        raise ValueError('Prediction times must be nonnegative.')
    velocity = states[..., 3:] if use_velocity else np.zeros_like(states[..., 3:])
    q = states[..., :3] + delay_s * velocity - nodes[None]
    distance = np.linalg.norm(q, axis=-1)
    if np.any(distance <= 0):
        raise ValueError('Undefined line of sight.')
    rate = np.sum(q / distance[..., None] * velocity, axis=-1)
    return distance[:, None, :] + np.asarray(slot_times_s)[None, :, None] * rate[:, None, :]
