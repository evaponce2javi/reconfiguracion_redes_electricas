"""Networks of experiments/x13_scale.py: the evaluation of the expressions of case533mt_hi and the source designs."""
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "experiments"))

from x13_scale import close_ties, evaluate_expressions, fixed_ratings, load  # noqa: E402


def test_expressions_are_evaluated_in_code_only():
    text = "mpc.baseMVA = 50/3; % base 50/3 MVA\n  1\t3\t0\t135/sqrt(3)\t1.05\t0.95;\n  7 -50/3 12/sqrt(3) 0.5 2;\n"
    out = evaluate_expressions(text).splitlines()
    assert out[0] == f"mpc.baseMVA = {50 / 3!r}; % base 50/3 MVA"
    assert out[1].split() == ["1", "3", "0", repr(135 / math.sqrt(3)), "1.05", "0.95;"]
    assert out[2].split() == ["7", repr(-50 / 3), repr(12 / math.sqrt(3)), "0.5", "2;"]


@pytest.fixture(scope="module")
def net533():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return load("numeric:case533mt_hi", ""), load("numeric:case533mt_hi", "ties")
        except OSError as e:                                  # no network access to download the case file
            pytest.skip(f"case533mt_hi could not be downloaded: {e}")


def test_case533_is_read_as_its_file_describes(net533):
    radial, meshed = net533
    assert radial.n_bus == meshed.n_bus == 533 and len(radial.pv) == 0 and len(radial.ref) == 1
    loops = lambda g: len({(min(a, b), max(a, b)) for a, b in zip(g.br_f, g.br_t) if a != b}) - g.n_bus + 1
    assert loops(radial) == 0 and loops(meshed) == 45
    # the file states a total net load of 45 MW (three-phase); its matrices hold single-phase MW on a 50/3 MVA base
    assert radial.base_mva == pytest.approx(50 / 3)
    assert 3 * radial.base_mva * radial.s_load.real.sum() == pytest.approx(45.0, rel=0.02)
    assert radial.base_kv[0] == pytest.approx(135 / math.sqrt(3)) and np.abs(radial.V0).min() > 0.9


def test_fixed_design_spreads_equal_sources(net533):
    g = net533[0]
    for n in (8, 32):
        r = fixed_ratings(g, n)
        dg = np.flatnonzero(r > 0)
        assert len(dg) == n and not np.isin(dg, g.ref).any()
        assert r.sum() == pytest.approx(abs(g.s_load.sum())) and np.allclose(r[dg], r[dg][0])
        lab = np.sort(np.asarray(g.bus_label)[dg])
        assert np.diff(lab).min() > 0.5 * g.n_bus / n         # evenly spread in label order
    assert close_ties is not None
