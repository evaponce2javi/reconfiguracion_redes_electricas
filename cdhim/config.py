"""Configuration file access and run provenance (versions, seed, commit)."""
from __future__ import annotations

import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from .core import MethodConfig

ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = ROOT / "config.yaml"
RESULTS_DIR = ROOT / "results"


def load_config(path: Path = CONFIG_PATH) -> dict:
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def method_config(cfg: dict, **overrides) -> MethodConfig:
    return MethodConfig(**{**cfg["method"], **overrides})


def versions() -> dict:
    import numpy
    import pandapower
    import scipy

    return {"python": platform.python_version(), "pandapower": pandapower.__version__,
            "numpy": numpy.__version__, "scipy": scipy.__version__, "platform": platform.platform()}


def _commit() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True, timeout=10)
        return out.stdout.strip()
    except Exception:
        return ""


def save_table(name: str, rows: list[dict], meta: dict | None = None, results_dir: Path = RESULTS_DIR) -> Path:
    """Write ``results/<name>.csv`` plus ``results/<name>.meta.json`` (seed, sample size, versions)."""
    import pandas as pd

    results_dir.mkdir(parents=True, exist_ok=True)
    path = results_dir / f"{name}.csv"
    pd.DataFrame(rows).to_csv(path, index=False, float_format="%.10g")
    info = {"table": name, "script": Path(sys.argv[0]).name, "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "git_commit": _commit(), "versions": versions()}
    info.update(meta or {})
    (results_dir / f"{name}.meta.json").write_text(json.dumps(info, indent=2, default=str), encoding="utf-8")
    return path
