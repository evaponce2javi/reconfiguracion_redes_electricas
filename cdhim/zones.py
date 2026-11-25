"""Zones for aggregated uncertainty sources: a deterministic, nested partition of the buses.

The partition is a stand-in for real control areas or forecast zones. It uses only the topology:
the buses are embedded with the first eigenvectors of the graph Laplacian and split recursively at
the median of the coordinate with the largest spread, so zones are compact and have equal bus counts
(to within one bus). The zones for 2k are a refinement of the zones for k.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

from .network import Grid

N_VECTORS = 8          # Laplacian eigenvectors used for the embedding


def laplacian_embedding(grid: Grid, n_vec: int = N_VECTORS) -> np.ndarray:
    """First ``n_vec`` non-constant eigenvectors of the (unweighted) graph Laplacian, one row per bus."""
    n = grid.n_bus
    A = (grid.Ybus != 0).astype(float).tolil()
    A.setdiag(0.0)
    A = A.tocsr()
    A.eliminate_zeros()
    A = ((A + A.T) > 0).astype(float)
    L = sp.diags(np.asarray(A.sum(axis=1)).ravel()) - A
    n_vec = min(n_vec, n - 2)
    if n <= 400:
        _, vec = np.linalg.eigh(L.toarray())
    else:                                   # shift-invert just below zero: L - sigma I is positive definite
        v0 = np.cos(np.arange(n) + 1.0)     # fixed start vector, so the result does not depend on a random draw
        val, vec = eigsh(L.tocsc(), k=n_vec + 1, sigma=-1e-2, which="LM", v0=v0, tol=1e-10)
        vec = vec[:, np.argsort(val)]
    vec = vec[:, 1:n_vec + 1]
    lead = np.abs(vec).argmax(axis=0)       # fix the sign of every eigenvector
    return vec * np.sign(vec[lead, np.arange(vec.shape[1])])


def spectral_zones(grid: Grid, k: int, embedding: np.ndarray | None = None) -> np.ndarray:
    """Zone label (0..k-1) of every bus; ``k`` must be a power of two."""
    levels = int(round(np.log2(k)))
    if 2 ** levels != k or k < 1:
        raise ValueError("k must be a power of two")
    X = laplacian_embedding(grid) if embedding is None else embedding
    label = np.zeros(grid.n_bus, dtype=int)
    for _ in range(levels):
        new = np.empty_like(label)
        for z in np.unique(label):
            idx = np.flatnonzero(label == z)
            axis = int(np.argmax(X[idx].var(axis=0)))
            order = idx[np.lexsort((idx, X[idx, axis]))]       # ties broken by bus index
            half = (len(order) + 1) // 2
            new[order[:half]] = 2 * z
            new[order[half:]] = 2 * z + 1
        label = new
    return label
