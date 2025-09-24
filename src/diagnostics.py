"""Legacy diagnostics module (shim).

Primary implementation lives in package module `redfield_education.diagnostics`.
Kept to allow `import diagnostics` in notebooks for brevity.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable
import numpy as np
from qutip import Qobj


@dataclass
class PositivityReport:
    min_eig: float
    negative: bool
    trace_dev: float

    def to_row(self) -> dict:
        return {
            "min_eig": self.min_eig,
            "negative": self.negative,
            "trace_dev": self.trace_dev,
        }


def density_matrix_positivity(rho: Qobj) -> PositivityReport:
    """Return positivity metrics for a density matrix.

    Parameters
    ----------
    rho : Qobj
        Density operator (should be Hermitian, trace ~ 1).
    """
    if not rho.isherm:
        raise ValueError("rho must be Hermitian")
    evals = rho.eigenenergies()  # eigenvalues since Hermitian
    min_eig = float(np.min(evals))
    trace_dev = float(abs(rho.tr() - 1.0))
    return PositivityReport(
        min_eig=min_eig, negative=min_eig < -1e-10, trace_dev=trace_dev
    )


def trajectory_min_eigs(rhos: Iterable[Qobj]) -> np.ndarray:
    """Vector of minimal eigenvalues over a trajectory."""
    return np.array([rho.eigenenergies().min() for rho in rhos], dtype=float)


def secular_filter(redfield_tensor, threshold: float) -> np.ndarray:
    """Apply a simple secular (|omega_i - omega_j| > threshold) mask to a Redfield tensor.

    This simplistic implementation expects a 4-index object shaped (N,N,N,N) as ndarray.
    Elements violating the secular condition are zeroed.
    """
    R = np.array(redfield_tensor, copy=True)
    if R.ndim != 4:
        raise ValueError("Expected rank-4 tensor")
    N = R.shape[0]
    # crude frequency ladder assumption: omega_i = i * delta (set delta=1) -> difference = |i-j|
    idx = np.arange(N)
    omega = idx.astype(float)  # pedagogical placeholder
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    dw = abs((omega[a] - omega[b]) - (omega[c] - omega[d]))
                    if dw < threshold:
                        continue
                    # keep; else zero? Actually secular keeps resonant terms -> set others zero
                    # so if dw != 0 and > threshold -> remove
                    if dw > threshold:
                        R[a, b, c, d] = 0.0
    return R
