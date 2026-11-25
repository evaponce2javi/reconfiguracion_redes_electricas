"""Shared plumbing for the experiment scripts: configuration, cached grids, seeds, process pool."""
from __future__ import annotations

import os
import sys
import warnings
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402

from cdhim import load_config, load_grid, make_study, method_config  # noqa: E402
from cdhim.core import MethodConfig  # noqa: E402

CFG = load_config()
_GRIDS: dict = {}

# What the unmodified finalFinalCode.py does: flat-start Jacobian with the H sign error and the
# absolute sx*sy filter. Used only to reproduce the figures of the submitted manuscript.
LEGACY_LIN = dict(point="flat", jacobian="legacy_signbug")
LEGACY_METHOD = dict(min_sigma_product=1e-12, min_radius_pu=0.0)


def grid(case: str):
    if case not in _GRIDS:
        _GRIDS[case] = load_grid(case, tol=CFG["powerflow"]["tol_pu"])
    return _GRIDS[case]


def study(case: str, load_factor: float = 1.0, legacy: bool = False, threads=None, **kw):
    lin = dict(LEGACY_LIN) if legacy else dict(point=CFG["linearization"]["point"], jacobian=CFG["linearization"]["jacobian"])
    lin.update(kw)
    return make_study(grid(case), load_factor, permc_spec=CFG["linearization"]["permc_spec"],
                      tol=CFG["powerflow"]["tol_pu"], max_chord=CFG["powerflow"]["max_chord"], threads=threads,
                      min_injection_pu=CFG["uncertainty"]["min_injection_pu"], **lin)


def method(legacy: bool = False, **overrides) -> MethodConfig:
    return method_config(CFG, **({**LEGACY_METHOD, **overrides} if legacy else overrides))


def seed_for(*tags) -> int:
    """Deterministic per-experiment seed derived from the base seed in config.yaml."""
    key = zlib.crc32("|".join(str(t) for t in tags).encode())
    return int(np.random.SeedSequence([CFG["experiments"]["seed"], key]).generate_state(1)[0])


def pct(x) -> float:
    return 100.0 * float(x)


def parallel(func, tasks, workers: int | None = None):
    """Run ``func`` over ``tasks`` in separate processes (one BLAS/solver thread each)."""
    tasks = list(tasks)
    workers = min(len(tasks), workers or max(1, (os.cpu_count() or 2) // 2 - 1))
    if workers <= 1:
        return [func(t) for t in tasks]
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    os.environ["CDHIM_THREADS"] = "1"
    with ProcessPoolExecutor(workers) as pool:
        return list(pool.map(func, tasks))


def worker_threads():
    """None in the main process (use the default), 1 inside a pool worker."""
    return int(os.environ["CDHIM_THREADS"]) if "CDHIM_THREADS" in os.environ else None
