"""Matplotlib plotting helpers (minimal) for consistent styling.

Utilities are intentionally lightweight to keep notebooks readable.
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from qutip import basis

PALETTE = [f"C{i}" for i in range(10)]


def style():
    plt.rcParams.update(
        {
            "figure.dpi": 130,
            "axes.grid": False,
            "font.size": 12,
            "lines.linewidth": 2,
        }
    )


def populations_from_states(states):
    """Return populations array shape (dim, n_times)."""
    dim = states[0].shape[0]
    projectors = [basis(dim, i) * basis(dim, i).dag() for i in range(dim)]
    pops = np.empty((dim, len(states)))
    for k, Pk in enumerate(projectors):
        pops[k] = [np.real((rho * Pk).tr()) for rho in states]
    return pops


def plot_populations(times, states, labels=None, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    pops = populations_from_states(states)
    for i, p in enumerate(pops):
        if np.allclose(p, 0):
            continue
        lab = labels[i] if labels else rf"$p_{{{i}}}$"
        ax.plot(times, p, color=PALETTE[i % len(PALETTE)], label=lab)
    ax.set_xlabel(r"Time $t$")
    ax.set_ylabel(r"Population")
    ax.legend(frameon=False, ncol=2)
    return ax


def plot_min_eig(times, min_eigs, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 3))
    ax.plot(times, min_eigs, color="C3", label=r"$\lambda_{\min}$")
    ax.axhline(0, color="k", lw=1, ls=":")
    ax.set_ylabel(r"$\lambda_{\min}$")
    ax.set_xlabel(r"Time $t$")
    ax.legend(frameon=False)
    return ax
