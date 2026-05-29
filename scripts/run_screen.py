#!/usr/bin/env python3
"""CLI for the V-dominant high-entropy NVPF pre-screening pipeline.

Examples
--------
# fast coarse demo (tier-1 + tier-2 elements, equal split, step 0.1)
python scripts/run_screen.py --out candidates.csv

# full grid partition coarse scan (large pool; addendum 2 Stage-1)
python scripts/run_screen.py --full-grid --out coarse_pool.csv

# tier-1 only, entropy floor enforced
python scripts/run_screen.py --tiers 1 --entropy-floor --out t1.csv
"""
import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from henvpf import pipeline  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="V-dominant HE-NVPF pre-screening")
    ap.add_argument("--out", default="candidates.csv")
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("--tiers", nargs="+", type=int, default=[1, 2])
    ap.add_argument("--full-grid", action="store_true",
                    help="enumerate every grid partition (large Stage-1 pool)")
    ap.add_argument("--entropy-floor", action="store_true")
    ap.add_argument("--input", default=None,
                    help="score compositions from a CSV instead of auto-generating "
                         "(see examples/my_compositions.csv)")
    ap.add_argument("--max", type=int, default=None)
    ap.add_argument("--top", type=int, default=25, help="rows to print to stdout")
    args = ap.parse_args()

    rep = pipeline.run(
        composition_step=args.step,
        restrict_tiers=tuple(args.tiers),
        equal_split_only=not args.full_grid,
        max_candidates=args.max,
        enforce_entropy_floor=args.entropy_floor,
        input_csv=args.input,
    )
    pipeline.write_csv(rep, args.out)

    print(f"generated............. {rep.generated}")
    print(f"passed hard screen.... {rep.passed_hard_screen}")
    print(f"triage................ {rep.triage_counts}")
    print(f"top rejection reasons. "
          f"{dict(sorted(rep.rejected_reasons.items(), key=lambda kv: -kv[1])[:6])}")
    print(f"\nTop {args.top} by priority_score:")
    hdr = f"{'formula':<46}{'x':>4}{'nd':>3}{'PP':>6}{'FR':>6}{'VA':>6}{'cap':>7}  {'pri':>5}  triage"
    print(hdr)
    print("-" * len(hdr))
    for r in rep.ranked[:args.top]:
        b = r.bundle
        print(f"{r.candidate.formula():<46}{r.candidate.x:>4g}{r.candidate.n_dopants:>3}"
              f"{b.phase_purity_score:>6.2f}{b.f_retention_score:>6.2f}"
              f"{b.v_activity_retention_score:>6.2f}{b.theoretical_capacity:>7.1f}"
              f"  {r.priority_score:>5.3f}  {r.triage}")
    print(f"\nwrote -> {args.out}")


if __name__ == "__main__":
    main()
