"""The port of RAA-K against the group's golden cases and against brute force."""
import itertools
from pathlib import Path

import numpy as np
import pytest

from cdhim.second_order import bound_monomial
from cdhim.knomial import block_bound, knomial_bound, knomial_rewrite

GOLDEN = Path(__file__).parent / "data" / "knomial_golden.txt"      # copied from ibex-lib/raak/tests


def golden_cases():
    tok = iter(GOLDEN.read_text().split())
    take = lambda k: np.array([float(next(tok)) for _ in range(k)])
    for _ in range(int(next(tok))):
        n, k = int(next(tok)), int(next(tok))
        w, Q = take(n), take(n * n).reshape(n, n)
        nb = int(next(tok))
        c, U = np.zeros(nb), np.zeros((nb, n))
        for b in range(nb):
            c[b] = float(next(tok))
            U[b] = take(n)
        yield n, k, w, Q, c, U, take(n * n).reshape(n, n)


def test_port_reproduces_the_reference_implementation():
    """Same check as ibex-lib/raak/tests/knomial_golden.cpp: same blocks and same Q_res in every case."""
    n_cases = n_blocks = 0
    worst = 0.0
    for n, k, w, Q, c_ref, U_ref, Qres_ref in golden_cases():
        c, U, Qres, downgraded = knomial_rewrite(Q, w, k)
        assert len(c) == len(c_ref), f"case {n_cases}: {len(c)} blocks, expected {len(c_ref)}"
        assert downgraded == 0
        diff = max(np.abs(Qres - Qres_ref).max(), np.abs(c - c_ref).max() if len(c) else 0.0, np.abs(U - U_ref).max() if len(c) else 0.0)
        worst = max(worst, diff)
        n_cases += 1
        n_blocks += len(c)
    assert n_cases == 80 and n_blocks > 0
    assert worst <= 1e-12


@pytest.mark.parametrize("max_k", [0, 2, 3, 6])
def test_rewriting_is_an_identity_and_the_bound_encloses(max_k):
    rng = np.random.default_rng(11)
    for n in (3, 5, 8):
        for trial in range(20):
            if trial % 2:                                  # near rank one and negative, like the power-flow forms
                v = rng.uniform(0.2, 1.0, n) * rng.choice([-1.0, 1.0], n)
                H = -np.outer(v, v) + 0.05 * rng.standard_normal((n, n))
            else:
                H = rng.standard_normal((n, n))
            H = 0.5 * (H + H.T)
            c, U, Qres, _ = knomial_rewrite(H, None, max_k)
            np.testing.assert_allclose(Qres + (U.T * c) @ U, H, atol=1e-12)          # x'Qx = x'Q_res x + sum c (u'x)^2
            lo, hi, nb = knomial_bound(H, max_k)
            X = np.array(list(itertools.product([-1.0, 0.0, 1.0], repeat=n))) if n <= 5 else rng.uniform(-1, 1, (4000, n))
            val = np.einsum("ij,jk,ik->i", X, H, X)
            assert lo <= val.min() + 1e-12 and hi >= val.max() - 1e-12
            mono = block_bound(H, np.zeros(0), np.zeros((0, n)))                       # no blocks: term by term (AF2)
            assert hi - lo <= mono[1] - mono[0] + 1e-12                                # never wider than AF2
            if max_k < 2:
                assert nb == 0


def test_full_clique_is_exact_on_a_rank_one_form():
    v = np.array([0.5, -1.0, 2.0, 0.25])
    H = -np.outer(v, v)
    lo2, hi2, _ = knomial_bound(H, 2)
    lo4, hi4, nb = knomial_bound(H, 4)
    exact = -np.abs(v).sum() ** 2
    assert nb == 1 and lo4 == pytest.approx(exact) and hi4 == pytest.approx(0.0, abs=1e-12)
    assert lo2 <= lo4 and hi2 >= hi4 and hi2 - lo2 > hi4 - lo4                         # pairs alone leave cross terms


def test_node_budget_keeps_the_bound_valid():
    rng = np.random.default_rng(3)
    v = rng.uniform(0.1, 1.0, 12)
    H = -np.outer(v, v) + 0.02 * rng.standard_normal((12, 12))
    H = 0.5 * (H + H.T)
    X = rng.uniform(-1, 1, (5000, 12))
    val = np.einsum("ij,jk,ik->i", X, H, X)
    for budget in (0, 5, 50):
        lo, hi, _ = knomial_bound(H, 12, node_budget=budget)
        assert lo <= val.min() and hi >= val.max()


def test_normalised_bound_is_scale_invariant_and_valid():
    """x'(sQ)x = s x'Qx: the bound must scale with s. Without normalising it does not (absolute tolerance)."""
    rng = np.random.default_rng(11)
    A = rng.normal(size=(8, 3))
    H = -(A @ A.T) + 0.3 * np.diag(rng.normal(size=8))                 # mostly negative definite, like the voltage forms
    H = 0.5 * (H + H.T)
    lo1, hi1, _ = knomial_bound(H, 8, normalise=True)
    eps = np.array([[1.0 if (v >> j) & 1 else -1.0 for v in range(256)] for j in range(8)])
    val = np.einsum("iv,ij,jv->v", eps, H, eps)
    mono = bound_monomial(H)
    assert lo1 <= val.min() + 1e-12 and hi1 >= val.max() - 1e-12
    assert hi1 - lo1 < 0.9 * (mono[1] - mono[0])                         # the rewriting does something here
    for s in (1e-3, 1e-6, 1e-9):
        lo, hi, _ = knomial_bound(s * H, 8, normalise=True)
        np.testing.assert_allclose([lo / s, hi / s], [lo1, hi1], rtol=1e-9)
        raw = knomial_bound(s * H, 8)
        assert raw[0] <= lo + 1e-12 * abs(lo) and raw[1] >= hi - 1e-12 * abs(hi)     # the raw bound is valid, only wider
    raw = knomial_bound(1e-9 * H, 8)
    assert (raw[1] - raw[0]) / 1e-9 > 1.2 * (hi1 - lo1)                  # documents the loss at small scale
    assert knomial_bound(np.zeros((4, 4)), 4, normalise=True) == (0.0, 0.0, 0)
