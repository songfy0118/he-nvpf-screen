"""
End-to-end orchestrator for stages 1-3 + rule-calibrated ranking.

Stage 1  generate candidates (composition.py)
Stage 2  chemical hard screen  (screening.py)        <- deterministic rejection
Stage 3  charge solver + 3 proxy scores (scores.py)
Rank     priority score + triage (ranking.py)

Downstream (MLIP config search -> VASP labels -> Bayesian optimisation) is left
to the integration stubs; this orchestrator emits the ranked candidate pool that
feeds them, following the addendum 9 calculation order.
"""
from __future__ import annotations
import csv
from dataclasses import dataclass, field
from typing import Dict, List

from . import composition, screening, scores, ranking


@dataclass
class PipelineReport:
    generated: int = 0
    passed_hard_screen: int = 0
    ranked: List[ranking.Ranked] = field(default_factory=list)
    triage_counts: Dict[str, int] = field(default_factory=dict)
    rejected_reasons: Dict[str, int] = field(default_factory=dict)


def run(
    x_values=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    composition_step=0.1,
    restrict_tiers=(1, 2),
    equal_split_only=True,
    max_candidates=None,
    enforce_entropy_floor=False,
    input_csv=None,
) -> PipelineReport:
    if input_csv:
        from . import userinput
        cands = userinput.load_compositions_csv(input_csv)
    else:
        cands = composition.generate(
            x_values=x_values,
            composition_step=composition_step,
            restrict_tiers=restrict_tiers,
            equal_split_only=equal_split_only,
            max_candidates=max_candidates,
        )
    report = PipelineReport(generated=len(cands))

    for cand in cands:
        hs = screening.hard_screen(cand, enforce_entropy_floor=enforce_entropy_floor)
        if not hs.passed:
            for r in hs.reasons:
                key = r.split(" ")[0].split(":")[0]
                report.rejected_reasons[key] = report.rejected_reasons.get(key, 0) + 1
            continue
        report.passed_hard_screen += 1
        bundle = scores.score_candidate(cand.v_amount, cand.dopants, hs.charge)
        ps = ranking.priority_score(bundle, cand)
        tg = ranking.triage(bundle, cand)
        report.ranked.append(ranking.Ranked(cand, bundle, ps, tg))
        report.triage_counts[tg] = report.triage_counts.get(tg, 0) + 1

    report.ranked.sort(key=lambda r: r.priority_score, reverse=True)
    return report


def write_csv(report: PipelineReport, path: str) -> None:
    if not report.ranked:
        return
    rows = [r.to_row() for r in report.ranked]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
