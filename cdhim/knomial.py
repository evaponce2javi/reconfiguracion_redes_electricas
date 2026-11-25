"""RAA-K: k-nomial rewriting of a quadratic form by diagonal budgeting, and the bound it gives.

Port of the group's implementation, ``knomial_rewrite`` and ``knomial_quad_bound`` in
``ibex-lib/raak/affine/ibex_KNomial.h`` (itself a bit-for-bit port of ``raa_lib/knomial.py`` with
mode='sdsos', sign_mode='fixed'). The search is reproduced operation by operation;
``tests/test_knomial.py`` checks it against the golden cases shipped with that repository.

    x' Q x  =  x' Q_res x  +  sum_b c_b (u_b' x)^2 ,      c_b = +1 or -1

Each block takes a clique of variables whose diagonal entries have the same sign, with signs s_v such
that s_u s_v = sign(c Q_uv) on every edge, and removes from the diagonals (the "budget" C_i = |Q_ii|)
and from the edges of the clique the largest multiple lambda that all its edges can pay for.
Blocks are taken greedily by the gain in AF2 error, from a heap of pairs extended to cliques of at
most ``max_k`` variables by depth-first search.

The bound follows the k-nomial theorem, valid for any blocks: for eps in [-1, 1]^n,

    eps' Q eps  in  sum_b c_b [0, |u_b|_1^2] + sum_i R_ii [0, 1] + sum_{i<j} 2 R_ij [-1, 1],
    R = Q - sum_b c_b u_b u_b'.

Unlike the C++ version this one does not round outward: it is not a verified bound.

Scale. The search compares entries and gains with an ABSOLUTE tolerance (``TINY``, as in the original), so it
is not scale-invariant: a form with entries of 1e-6 loses most of its blocks and the bound falls back to the
term-by-term one, although x'(sQ)x = s x'Qx. ``knomial_bound(..., normalise=True)`` rewrites Q / max|Q| and
scales the bound back, which restores the invariance. The default stays False so that the experiments that
were run without it (X2 to X7) still regenerate their tables.
"""
from __future__ import annotations

import heapq

import numpy as np
from numba import njit

TINY = 1e-12


@njit(cache=True)
def _sign1(v):
    return np.copysign(1.0, v)


@njit(cache=True)
def _benefit(clique, m, Q, C, w):
    """(lambda, gain) of a clique: the scale every edge can pay and the AF2 error it removes."""
    lam = 1.0
    for x in range(m):
        for y in range(x + 1, m):
            u, v = clique[x], clique[y]
            cp = C[u] * C[v]
            if cp <= 1e-12:
                return 0.0, 0.0
            lam = min(lam, abs(Q[u, v]) / np.sqrt(cp))
    de = 0.0
    for x in range(m):
        for y in range(x + 1, m):
            u, v = clique[x], clique[y]
            de += lam * np.sqrt(C[u] * C[v]) * w[u] * w[v] / 2.0
    return lam, de


@njit(cache=True)
def _bfs_signs(S, m, Q, r, sign, seen, dq, tol):
    """Signs of a block r (sum s_v x_v)^2 with s_u s_v = sign(r Q_uv) on every edge; False if impossible."""
    for a in range(m):
        seen[S[a]] = False
    for a in range(m):
        root = S[a]
        if seen[root]:
            continue
        seen[root] = True
        sign[root] = 1.0
        dq[0] = root
        head, tail = 0, 1
        while head < tail:
            u = dq[head]
            head += 1
            for b in range(m):
                v = S[b]
                if seen[v] or abs(Q[u, v]) <= tol:
                    continue
                seen[v] = True
                sign[v] = sign[u] * _sign1(r * Q[u, v])
                dq[tail] = v
                tail += 1
    ok = True
    for a in range(m):
        for b in range(a + 1, m):
            u, v = S[a], S[b]
            if abs(Q[u, v]) > tol and sign[u] * sign[v] != _sign1(r * Q[u, v]):
                ok = False
    return ok


@njit(cache=True)
def _rewrite(Q, w, max_k, tol, budget):
    """In place: Q becomes Q_res. Returns (c, U, number of unbalanced cliques downgraded to their pair)."""
    n = Q.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            v = 0.5 * (Q[i, j] + Q[j, i])
            Q[i, j] = v
            Q[j, i] = v
    cap_b = 16
    blk_c = np.zeros(cap_b)
    blk_u = np.zeros((cap_b, n))
    nb = 0
    if n == 0 or max_k < 2:
        return blk_c[:0], blk_u[:0], 0

    C, rvec = np.zeros(n), np.zeros(n)
    P, N = np.zeros(n, np.bool_), np.zeros(n, np.bool_)
    for i in range(n):
        d = Q[i, i]
        C[i] = abs(d)
        P[i] = d > 0
        N[i] = d < 0
        rvec[i] = 1.0 if P[i] else (-1.0 if N[i] else 0.0)
    A = np.zeros((n, n), np.bool_)
    for i in range(n):
        for j in range(i + 1, n):
            if Q[i, j] != 0 and ((P[i] and P[j]) or (N[i] and N[j])):
                A[i, j] = True
                A[j, i] = True

    pair = np.zeros(2, np.int64)
    heap = [(0.0, 0, 0)]
    heap.pop()
    for i in range(n):
        for j in range(i + 1, n):
            if A[i, j]:
                pair[0], pair[1] = i, j
                de = _benefit(pair, 2, Q, C, w)[1]
                if de > 1e-12:
                    heapq.heappush(heap, (-de, i, j))

    kk = min(max_k, n)
    cap_s = kk * n + 4                                   # depth-first stack: at most n children per level
    st_clique = np.zeros((cap_s, kk), np.int64)
    st_sg = np.zeros((cap_s, kk))
    st_len = np.zeros(cap_s, np.int64)
    st_cands = np.zeros((cap_s, n), np.int32)
    st_nc = np.zeros(cap_s, np.int64)
    cur_clique, cur_sg, cur_cands = np.zeros(kk, np.int64), np.zeros(kk), np.zeros(n, np.int32)
    best_clique = np.zeros(kk, np.int64)
    clique, c_loc = np.zeros(kk, np.int64), np.zeros(kk)
    signs, seen, dq = np.zeros(n), np.zeros(n, np.bool_), np.zeros(n, np.int64)

    n_infeasible, guard, dfs_nodes = 0, 0, 0
    max_iter = 200 * max(1, n * n)
    while len(heap) > 0 and guard < max_iter:
        guard += 1
        neg_prev, bi, bj = heapq.heappop(heap)
        prev = -neg_prev
        if not A[bi, bj]:
            continue
        pair[0], pair[1] = bi, bj
        p_ij, de_ij = _benefit(pair, 2, Q, C, w)
        if de_ij < prev - 1e-9:                          # stale entry: its gain dropped, put it back
            if de_ij > 1e-12:
                heapq.heappush(heap, (-de_ij, bi, bj))
            continue
        if de_ij <= 1e-12:
            continue

        r = rvec[bi]
        best_clique[0], best_clique[1] = bi, bj
        best_len, best_de, best_p = 2, de_ij, p_ij

        # depth-first search of balanced cliques with at most max_k variables
        sp = 1
        st_clique[0, 0], st_clique[0, 1] = bi, bj
        st_len[0] = 2
        st_sg[0, 0], st_sg[0, 1] = 1.0, _sign1(r * Q[bi, bj])
        nc = 0
        for v in range(n):
            if A[bi, v] and A[bj, v]:
                st_cands[0, nc] = v
                nc += 1
        st_nc[0] = nc
        is_root = True
        while sp > 0:
            sp -= 1
            m = st_len[sp]
            nc = st_nc[sp]
            for t in range(m):
                cur_clique[t] = st_clique[sp, t]
                cur_sg[t] = st_sg[sp, t]
            for t in range(nc):
                cur_cands[t] = st_cands[sp, t]
            dfs_nodes += 1
            if not is_root:
                lam_c, de_c = _benefit(cur_clique, m, Q, C, w)
                if de_c > best_de and lam_c > 1e-12:
                    best_de, best_p, best_len = de_c, lam_c, m
                    for t in range(m):
                        best_clique[t] = cur_clique[t]
            is_root = False
            # each clique is enumerated once (descending order of the added variables); with the node
            # budget spent, pending nodes are still evaluated but not extended
            if m < max_k and (budget == 0 or dfs_nodes < budget):
                last = cur_clique[m - 1] if m > 2 else n
                for ci in range(nc):
                    v = cur_cands[ci]
                    if v >= last:
                        continue
                    s_v = cur_sg[0] * _sign1(r * Q[bi, v])
                    ok = True
                    for t in range(1, m):
                        if _sign1(r * Q[cur_clique[t], v]) != cur_sg[t] * s_v:
                            ok = False
                            break
                    if ok:
                        for t in range(m):
                            st_clique[sp, t] = cur_clique[t]
                            st_sg[sp, t] = cur_sg[t]
                        st_clique[sp, m] = v
                        st_sg[sp, m] = s_v
                        st_len[sp] = m + 1
                        k2 = 0
                        for cj in range(nc):
                            if A[v, cur_cands[cj]]:
                                st_cands[sp, k2] = cur_cands[cj]
                                k2 += 1
                        st_nc[sp] = k2
                        sp += 1

        m = best_len
        p = best_p
        for t in range(m):
            clique[t] = best_clique[t]
            c_loc[t] = C[clique[t]]
        if not _bfs_signs(clique, m, Q, r, signs, seen, dq, tol):
            n_infeasible += 1                            # unbalanced clique: fall back to the base pair
            m, p = 2, p_ij
            clique[0], clique[1] = bi, bj
            c_loc[0], c_loc[1] = C[bi], C[bj]
            _bfs_signs(clique, m, Q, r, signs, seen, dq, tol)
        lam = p
        if lam <= 1e-12:
            continue

        if nb == cap_b:
            cap_b *= 2
            blk_c = np.concatenate((blk_c, np.zeros(nb)))
            blk_u = np.concatenate((blk_u, np.zeros((nb, n))))
        blk_c[nb] = r
        for t in range(m):
            v = clique[t]
            blk_u[nb, v] = signs[v] * np.sqrt(lam * c_loc[t])
        nb += 1
        for t in range(m):
            v = clique[t]
            C[v] = max(0.0, C[v] * (1.0 - lam))
        for x in range(m):
            for y in range(x + 1, m):
                uu, vv = clique[x], clique[y]
                s = _sign1(Q[uu, vv])
                nv = max(0.0, abs(Q[uu, vv]) - p * np.sqrt(c_loc[x] * c_loc[y]))
                Q[uu, vv] = s * nv
                Q[vv, uu] = Q[uu, vv]
                if nv < 1e-12:
                    A[uu, vv] = False
                    A[vv, uu] = False
        if abs(Q[bi, bj]) > 1e-12:
            heapq.heappush(heap, (-de_ij, bi, bj))

    for k in range(n):
        Q[k, k] = C[k] if P[k] else (-C[k] if N[k] else 0.0)
    return blk_c[:nb], blk_u[:nb], n_infeasible


def knomial_rewrite(Q, w=None, max_k: int = 3, tol: float = TINY, node_budget: int = 0):
    """Rewrite x'Qx as x'Q_res x + sum_b c_b (u_b'x)^2. Returns (c, U, Q_res, n_downgraded).

    ``w``: widths of the variables (2 for noise symbols in [-1, 1], the default). ``node_budget``: cap on
    the nodes of the clique search over the whole call (0 = no cap); once spent, cliques stop growing.
    """
    Q = np.array(Q, dtype=float)
    w = np.full(len(Q), 2.0) if w is None else np.asarray(w, dtype=float)
    c, U, bad = _rewrite(Q, w, int(max_k), float(tol), int(node_budget))
    return c, U, Q, bad


def block_bound(H, c, U):
    """[lo, hi] enclosing eps'H eps on [-1, 1]^n from any blocks (c_b, u_b): the k-nomial theorem."""
    H = 0.5 * (H + H.T)
    R = H - (U.T * c) @ U
    s2 = np.abs(U).sum(axis=1) ** 2
    d = np.diag(R)
    off = np.abs(R).sum() - np.abs(d).sum()                   # = 2 sum_{i<j} |R_ij|
    lo = -s2[c < 0].sum() + np.minimum(d, 0.0).sum() - off
    hi = s2[c > 0].sum() + np.maximum(d, 0.0).sum() + off
    return float(lo), float(hi)


def knomial_bound(H, max_k: int = 3, node_budget: int = 0, normalise: bool = False):
    """Bound of eps'H eps on [-1, 1]^n with the RAA-K rewriting. Returns (lo, hi, number of blocks).

    ``normalise``: bound H / max|H| and scale the result back (see "Scale" in the module docstring)."""
    H = np.asarray(H, dtype=float)
    scale = float(np.abs(H).max()) if normalise and H.size else 1.0
    if scale == 0.0:
        return 0.0, 0.0, 0
    Q = H / scale if normalise else H
    c, U, _, _ = knomial_rewrite(Q, None, max_k, node_budget=node_budget)
    lo, hi = block_bound(Q, c, U)
    return lo * scale, hi * scale, len(c)
