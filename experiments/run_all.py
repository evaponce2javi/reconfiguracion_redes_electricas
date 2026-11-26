"""Regenerate every table in results/ from scratch.

    python experiments/run_all.py            # everything (roughly 2-2.5 h on 8 cores; C1 takes most of it)
    python experiments/run_all.py b1 c5      # only the named experiments (x1, x2 = affine-arithmetic studies)
    python experiments/build_report.py       # then rebuild RESULTS.md from the CSV files

The slowest single piece is the C1 cell at 95 % of the loadability limit with delta = 30 % on
PEGASE-2869 (about 30 min: most samples need full Newton-Raphson). The timing experiments (c7, c8)
run last and alone so that other jobs do not distort them.
"""
import subprocess
import sys
from pathlib import Path
from time import perf_counter

HERE = Path(__file__).resolve().parent
ORDER = ["a_reproduce_paper", "b1_ac_violation", "b2_conservatism", "b3_table1", "b4_second_order", "b5_affine_forms", "c1_ac_sweep", "c2_correlated_sources",
         "c3_nongaussian", "c4_rho_rule", "c5_system_level", "c6_gamma_transfer", "x1_aa_second_order", "x2_few_sources", "x3_worst_vertex", "x4_reference_checks",
         "x5_liao_baseline", "x6_fixed_point", "x7_block_choices", "x8_feeders",
         "x9_krawczyk", "x10_three_phase", "x12_verified", "x13_scale", "x14_pv",
         "c7_scalability",
         "c8_methods"]


def main():
    wanted = [a.lower() for a in sys.argv[1:]]
    scripts = [s for s in ORDER if not wanted or s.split("_")[0] in wanted or s in wanted]
    failed = []
    for name in scripts:
        t0 = perf_counter()
        print(f"=== {name}", flush=True)
        code = subprocess.run([sys.executable, str(HERE / f"{name}.py")], cwd=HERE.parent).returncode
        print(f"=== {name}: {'ok' if code == 0 else f'FAILED ({code})'} in {perf_counter() - t0:.0f} s", flush=True)
        if code:
            failed.append(name)
    if failed:
        sys.exit("failed: " + ", ".join(failed))


if __name__ == "__main__":
    main()
