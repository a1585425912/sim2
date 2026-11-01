"""DPC 实验需要的最小仿真公共函数。

这个文件是从旧项目中迁移出来的精简核心，不包含旧实验脚本、
旧画图逻辑或其它无关方法。新项目内部只依赖这个文件，不再依赖
旧实验目录里的代码。
"""

from __future__ import annotations

import numpy as np


def wrap_to_pi(x: np.ndarray) -> np.ndarray:
    """把相位包裹到 [-pi, pi)。"""
    return (x + np.pi) % (2 * np.pi) - np.pi


def wrap_to_2pi(x: np.ndarray) -> np.ndarray:
    """把相位包裹到 [0, 2pi)。"""
    return x % (2 * np.pi)


def apply_polarity_correction(x: np.ndarray) -> np.ndarray:
    """Retired oracle API: truth-dependent polarity control is not executable."""
    raise ValueError("oracle polarity correction removed; supply an estimated phase command")


def effective_phase_errors(theta: np.ndarray, phi_true: np.ndarray, eps: np.ndarray) -> np.ndarray:
    """Actual residual, wrapped modulo 2*pi for reporting only."""
    return wrap_to_pi(np.asarray(theta) - np.asarray(phi_true) + np.asarray(eps))


def sample_phase_std(phases: np.ndarray) -> float:
    """计算相位残差的样本标准差。

    这里沿用旧实验中的 sigma_phi 口径：先减去样本均值，再计算
    普通样本标准差。它不是 circular std，目的是便于和之前图表比较。
    """
    vals = np.asarray(phases, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size < 2:
        return np.nan
    centered = vals - np.mean(vals)
    return float(np.sqrt(np.sum(np.abs(centered) ** 2) / (vals.size - 1)))


def phase_rmse(phases: np.ndarray) -> float:
    """计算残余相位相对于理想零相位的均方根误差。"""
    vals = np.asarray(phases, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan
    return float(np.sqrt(np.mean(vals**2)))


def disk_uniform(n: int, radius: float, rng: np.random.Generator) -> np.ndarray:
    """在半径为 radius 的圆盘中均匀部署二维点。"""
    rr = radius * np.sqrt(rng.random(n))
    aa = 2 * np.pi * rng.random(n)
    return np.stack([rr * np.cos(aa), rr * np.sin(aa)], axis=1)


def deploy_buoys_in_disk(n: int, area_radius: float, rng: np.random.Generator) -> np.ndarray:
    """部署 N 个节点/浮标的初始二维位置。"""
    return disk_uniform(n, area_radius, rng)


def metropolis_hastings_W(n: int, connectivity: float, rng: np.random.Generator) -> np.ndarray:
    """生成 Metropolis-Hastings 双随机共识矩阵。

    先根据 connectivity 生成随机无向图，再额外加入一个 ring backbone，
    避免出现孤立节点。W 用于节点间线性共识：X <- W @ X。
    """
    A = (rng.random((n, n)) < connectivity).astype(np.float64)
    A = np.triu(A, 1)
    A = A + A.T

    idx = np.arange(n)
    A[idx, (idx + 1) % n] = 1.0
    A[(idx + 1) % n, idx] = 1.0

    deg = A.sum(axis=1)
    W = np.zeros((n, n), dtype=np.float64)
    rows, cols = np.nonzero(np.triu(A, 1))
    for i, j in zip(rows, cols):
        wij = 1.0 / (max(deg[i], deg[j]) + 1.0)
        W[i, j] = wij
        W[j, i] = wij
    W[np.diag_indices(n)] = 1.0 - W.sum(axis=1)
    return W


def consensus_linear(W: np.ndarray, X: np.ndarray, steps: int) -> np.ndarray:
    """执行若干步线性共识。"""
    Y = X
    for _ in range(max(steps, 0)):
        Y = W @ Y
    return Y


def power_from_theta(theta: np.ndarray, phi_true: np.ndarray, amp: np.ndarray, eps: np.ndarray) -> float:
    """Actual received field power; never fold residual phases modulo pi."""
    residual = np.asarray(theta) - np.asarray(phi_true) + np.asarray(eps)
    return float(np.abs(np.sum(amp * np.exp(1j * residual))) ** 2)


def weights_from_adjacency(adjacency: np.ndarray) -> np.ndarray:
    """Metropolis weights for observed symmetric links, including isolates."""
    A = np.asarray(adjacency, dtype=bool).copy()
    if A.ndim != 2 or A.shape[0] != A.shape[1] or not np.array_equal(A, A.T):
        raise ValueError("adjacency must be square and symmetric")
    np.fill_diagonal(A, False)
    degree = A.sum(axis=1)
    W = A / (1 + np.maximum(degree[:, None], degree[None, :]))
    W[np.diag_indices(len(A))] = 1 - W.sum(axis=1)
    return W
