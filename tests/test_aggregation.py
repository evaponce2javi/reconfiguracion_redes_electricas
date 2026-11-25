"""Aggregated uncertainty sources: several injections driven by one noise symbol (IEEE-33)."""
import warnings

import numpy as np
import pytest
import scipy.sparse as sp

from cdhim import load_grid, make_study, participation_factors, sample_eps, spectral_zones


@pytest.fixture(scope="module")
def grid():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_grid("case33bw")


def test_one_symbol_per_injection_reproduces_the_default(grid):
    """An identity aggregation must give exactly what the code gives without aggregation."""
    ref = make_study(grid)
    st = make_study(grid)
    st.src.Z = sp.identity(ref.src.n, format="csr")
    assert st.src.n == ref.src.n == 64
    eps = sample_eps(np.random.default_rng(1), 64, 5)
    np.testing.assert_array_equal(st.injections(eps, 0.15), ref.injections(eps, 0.15))
    one_per_bus = make_study(grid, zones=np.arange(grid.n_bus), couple_pq=False)      # every bus its own zone
    assert one_per_bus.src.n == 64
    np.testing.assert_allclose(np.sort(np.abs(one_per_bus.C1), axis=1), np.sort(np.abs(ref.C1), axis=1), rtol=1e-12, atol=1e-16)


@pytest.mark.parametrize("k", [4, 8, 16])
def test_zone_symbols_are_sums_of_elementary_symbols(grid, k):
    zones = spectral_zones(grid, k)
    ref = make_study(grid)
    st = make_study(grid, zones=zones)
    Z = st.src.Z
    assert st.src.n == Z.shape[1] <= k and st.src.n_injections == 64
    np.testing.assert_array_equal(np.asarray(Z.sum(axis=1)).ravel(), 1.0)             # each injection in exactly one zone
    np.testing.assert_allclose(st.C1, ref.C1 @ Z.toarray(), rtol=1e-10, atol=1e-15)
    eps = sample_eps(np.random.default_rng(2), st.src.n, 7)
    np.testing.assert_allclose(st.injections(eps, 0.15), ref.injections(Z @ eps, 0.15), rtol=0, atol=1e-15)
    # P and Q of a bus share the symbol: the power factor of every load is unchanged
    S = st.injections(eps, 0.15) - st.sbus[:, None]
    load = np.abs(grid.s_load) > 1e-6
    ratio = S[load] / (-grid.s_load[load, None])
    np.testing.assert_allclose(ratio.imag, 0.0, atol=1e-12)


def test_spectral_zones_are_balanced_nested_and_deterministic(grid):
    z4, z8, z16 = (spectral_zones(grid, k) for k in (4, 8, 16))
    for k, z in ((4, z4), (8, z8), (16, z16)):
        sizes = np.bincount(z, minlength=k)
        assert len(sizes) == k and sizes.max() - sizes.min() <= 1
        np.testing.assert_array_equal(z, spectral_zones(grid, k))
    np.testing.assert_array_equal(z8 // 2, z4)                                         # refinement
    np.testing.assert_array_equal(z16 // 2, z8)
    with pytest.raises(ValueError):
        spectral_zones(grid, 6)


def test_aggregated_linear_prediction_error_is_second_order(grid):
    st = make_study(grid, zones=spectral_zones(grid, 8))
    eps = sample_eps(np.random.default_rng(5), st.src.n, 200)
    errs = []
    for delta in (0.20, 0.10, 0.05):
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        assert conv.all()
        errs.append(np.abs(np.abs(V[grid.pq]) - st.vm_star[:, None] - st.C1 @ eps * delta).max())
    assert errs[0] / errs[1] == pytest.approx(4.0, rel=0.15)
    assert errs[1] / errs[2] == pytest.approx(4.0, rel=0.15)


@pytest.fixture(scope="module")
def grid118():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_grid("case118")


def test_generator_participation_balances_every_symbol(grid118):
    g = grid118
    alpha = participation_factors(g)
    assert alpha.sum() == pytest.approx(1.0) and (alpha >= 0).all() and alpha[g.pq].sum() == 0.0
    zones = spectral_zones(g, 8)
    single = make_study(g, zones=zones)
    shared = make_study(g, zones=zones, participation=alpha)
    eps = np.ones((shared.src.n, 1))                                     # every zone at +delta
    dS = shared.injections(eps, 0.05) - shared.sbus[:, None]
    assert abs(dS.real.sum()) < 1e-10                                    # specified injections stay balanced

    def slack_power(st):
        V = st.solver.solve(st.injections(eps, 0.05))[0][:, 0]
        return (V * np.conj(g.Ybus @ V)).real[g.ref].sum()

    base = (g.V0 * np.conj(g.Ybus @ g.V0)).real[g.ref].sum()
    d_single, d_shared = slack_power(single) - base, slack_power(shared) - base
    total = -(single.injections(eps, 0.05) - single.sbus[:, None]).real.sum()          # demand increase
    assert total < d_single < 1.15 * total                               # slack alone: the whole change plus extra losses
    # with participation the slack supplies its own share plus the change in losses (about 6 % of the increase here)
    assert alpha[g.ref].sum() * total < d_shared < alpha[g.ref].sum() * total + 0.10 * total
    assert d_shared < 0.25 * d_single


def test_participation_keeps_the_linear_model_second_order_accurate(grid118):
    g = grid118
    st = make_study(g, zones=spectral_zones(g, 8), participation=participation_factors(g))
    eps = sample_eps(np.random.default_rng(6), st.src.n, 200)
    errs = []
    for delta in (0.20, 0.10, 0.05):
        V, conv, _ = st.solver.solve(st.injections(eps, delta))
        assert conv.all()
        errs.append(np.abs(np.abs(V[g.pq]) - st.vm_star[:, None] - st.C1 @ eps * delta).max())
    assert errs[0] / errs[1] == pytest.approx(4.0, rel=0.2)
    assert errs[1] / errs[2] == pytest.approx(4.0, rel=0.2)
