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


def coherences_from_states(states):
    """Return coherences array shape (n_coherences, n_times)."""
    dim = states[0].shape[0]
    n_coherences = dim * (dim - 1) // 2
    cohs = np.empty((n_coherences, len(states)))
    idx = 0
    for i in range(dim):
        for j in range(i + 1, dim):
            cohs[idx] = [np.abs(rho[i, j]) for rho in states]
            idx += 1
    return cohs


def plot_populations(times, states, labels=None, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    pops = populations_from_states(states)
    for i, p in enumerate(pops):
        if np.allclose(p, 0):
            continue
        lab = labels[i] if labels and i < len(labels) else rf"$p_{{{i}}}$"
        color_index = (len(ax.lines) + i) % len(PALETTE)
        ax.plot(times, p, color=PALETTE[color_index], label=lab)
    ax.set_xlabel(r"Time $t$")
    ax.set_ylabel(r"Population")
    ax.legend(frameon=False, ncol=2)
    return ax


def plot_coherences(times, states, labels=None, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    cohs = coherences_from_states(states)
    dim = states[0].shape[0]
    if labels is None:
        labels = [
            rf"$|\rho_{{{i}{j}}}|$" for i in range(dim) for j in range(i + 1, dim)
        ]
    for i, c in enumerate(cohs):
        if np.allclose(c, 0):
            continue
        lab = labels[i] if i < len(labels) else rf"$c_{{{i}}}$"
        color_index = (len(ax.lines) + i) % len(PALETTE)
        ax.plot(times, c, color=PALETTE[color_index], label=lab)
    ax.set_xlabel(r"Time $t$")
    ax.set_ylabel(r"Coherence $|\rho_{ij}|$")
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
