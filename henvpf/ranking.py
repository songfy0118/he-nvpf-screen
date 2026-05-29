"""
Rule-calibrated ranking model + triage (addendum 8).

This is explicitly NOT a supervised predictor of formation energy, Ehull,
F-vacancy energy, cycle life, or real capacity retention. It is a weak-supervision
ranking layer that aggregates the three proxy scores into a single priority score
and triages each candidate into one of:
    priority | exploration | reject | active_learning_pool

The aggregate weighting encodes the addendum's priority order:
    phase purity (main-phase) and F retention are the gating concerns;
    V-activity retention is the practical floor;
    entropy is a tie-breaker, never a reason to beat a stronger NVPF candidate.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict

from .scores import ScoreBundle
from .composition import Candidate

# aggregate weights (priority order, addendum 1 & 8)
W_PHASE = 0.35
W_F = 0.30
W_V = 0.25
W_ENTROPY = 0.10   # tie-breaker only; small


def priority_score(bundle: ScoreBundle, cand: Candidate) -> float:
    # entropy contribution is normalised against the strict high-entropy line (1.5R)
    # and capped, so it can only ever act as a tie-breaker
    ent = min(cand.sconfig_over_R / 1.5, 1.0)
    s = (W_PHASE * bundle.phase_purity_score
         + W_F * bundle.f_retention_score
         + W_V * bundle.v_activity_retention_score
         + W_ENTROPY * ent)
    return round(s, 4)


def triage(bundle: ScoreBundle, cand: Candidate) -> str:
    gates = (bundle.phase_gate, bundle.f_gate, bundle.v_gate)
    if "reject" in gates:
        # rejected by a hard gate, but keep borderline / informative ones for AL
        if cand.low_dose_multi_dopant_risk or _near_miss(bundle):
            return "active_learning_pool"
        return "reject"
    if "exploration" in gates or cand.low_dose_multi_dopant_risk:
        return "exploration"
    return "priority"


def _near_miss(bundle: ScoreBundle) -> bool:
    return (bundle.phase_purity_score >= 0.50
            and bundle.f_retention_score >= 0.50
            and bundle.v_activity_retention_score >= 0.55)


@dataclass
class Ranked:
    candidate: Candidate
    bundle: ScoreBundle
    priority_score: float
    triage: str

    def to_row(self) -> Dict:
        row = {
            "formula": self.candidate.formula(),
            "x": self.candidate.x,
            "n_dopants": self.candidate.n_dopants,
            "v_amount": self.candidate.v_amount,
            "dopants": ";".join(f"{e}{a:g}" for e, a in sorted(self.candidate.dopants.items())),
            "sconfig_over_R": self.candidate.sconfig_over_R,
            "entropy_band": self.candidate.entropy_band,
            "low_dose_multi_dopant_risk": self.candidate.low_dose_multi_dopant_risk,
            "priority_score": self.priority_score,
            "triage": self.triage,
        }
        row.update(self.bundle.to_row())
        return row
