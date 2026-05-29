"""
Explicit charge solver (addendum 4).

Replaces the vague "charge compensation difficulty" judgement with an actual
search over allowed dopant oxidation states + mean-V valence + Na vacancy that
satisfies overall charge neutrality of Na_{3-d} (V,M..)_2 (PO4)2 F3.

Balance equation (per formula unit):
    Na_charge + M_site_charge + anion_framework_charge = 0
    (3 - na_delta) * (+1) + [V_amount*meanV + sum_i amt_i*ox_i] + (-9) = 0
=>  V_amount*meanV + sum_i amt_i*ox_i = 6 + na_delta
"""
from __future__ import annotations
import itertools
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .elements import ELEMENTS

# Default screening bounds (addendum 4)
NA_DELTA_MIN, NA_DELTA_MAX = 0.0, 0.3
V_VALENCE_MIN, V_VALENCE_MAX = 2.5, 4.5
V_NOMINAL = 3.0

# charge balance mode labels (addendum 4)
MODE_NEUTRAL = "neutral_without_compensation"
MODE_V_ADJUST = "v_valence_adjustment"
MODE_NA_VAC = "na_vacancy_compensation"
MODE_PAIRING = "low_high_valence_pairing"
MODE_MIXED = "mixed_compensation"
MODE_UNRESOLVED = "unresolved"


@dataclass
class ChargeResult:
    charge_balance_mode: str
    na_delta_required: float
    mean_v_valence: float
    v_valence_range: Tuple[float, float]
    valence_feasible_flag: bool
    charge_compensation_difficulty: float          # normalised 0..1
    dopant_valence_assignment: Dict[str, int] = field(default_factory=dict)
    difficulty_label: str = ""


def _difficulty_label(score: float) -> str:
    if score < 0.20:
        return "easy"
    if score < 0.50:
        return "moderate"
    if score < 0.80:
        return "difficult"
    return "reject_unless_exploration"


def solve(v_amount: float, dopants: Dict[str, float]) -> ChargeResult:
    """Find the lowest-difficulty feasible charge assignment for one candidate."""
    elems = list(dopants)
    state_lists = [ELEMENTS[e].ox_states for e in elems]

    best: Optional[ChargeResult] = None
    feasible_meanVs: List[float] = []

    # iterate over every combination of dopant oxidation states
    for states in itertools.product(*state_lists) if elems else [()]:
        dop_charge = sum(dopants[e] * s for e, s in zip(elems, states))
        # need: V_amount*meanV = (6 + na_delta) - dop_charge
        # scan na_delta grid in [0, 0.3]
        for k in range(0, 7):
            na_delta = round(k * 0.05, 4)
            target_v_charge = (6.0 + na_delta) - dop_charge
            if v_amount <= 0:
                continue
            mean_v = target_v_charge / v_amount
            if not (V_VALENCE_MIN <= mean_v <= V_VALENCE_MAX):
                continue
            feasible_meanVs.append(mean_v)
            assignment = dict(zip(elems, states))
            mode = _classify_mode(mean_v, na_delta, assignment, dopants)
            diff = _difficulty(mean_v, na_delta, assignment, dopants)
            res = ChargeResult(
                charge_balance_mode=mode,
                na_delta_required=na_delta,
                mean_v_valence=round(mean_v, 4),
                v_valence_range=(V_VALENCE_MIN, V_VALENCE_MAX),
                valence_feasible_flag=True,
                charge_compensation_difficulty=round(diff, 4),
                dopant_valence_assignment=assignment,
                difficulty_label=_difficulty_label(diff),
            )
            if best is None or res.charge_compensation_difficulty < best.charge_compensation_difficulty:
                best = res

    if best is None:
        return ChargeResult(
            charge_balance_mode=MODE_UNRESOLVED,
            na_delta_required=float("nan"),
            mean_v_valence=float("nan"),
            v_valence_range=(V_VALENCE_MIN, V_VALENCE_MAX),
            valence_feasible_flag=False,
            charge_compensation_difficulty=1.0,
            difficulty_label="reject_unless_exploration",
        )
    return best


def _classify_mode(mean_v, na_delta, assignment, dopants) -> str:
    has_low = any(s <= 2 for s in assignment.values())
    has_high = any(s >= 5 for s in assignment.values())
    v_shift = abs(mean_v - V_NOMINAL) > 0.05
    na_used = na_delta > 1e-6

    if not v_shift and not na_used and not has_low and not has_high:
        return MODE_NEUTRAL
    if has_low and has_high:
        return MODE_PAIRING
    if na_used and v_shift:
        return MODE_MIXED
    if na_used:
        return MODE_NA_VAC
    return MODE_V_ADJUST


def _difficulty(mean_v, na_delta, assignment, dopants) -> float:
    """0..1, rises with the compensation burden (addendum 4 'difficulty score')."""
    total = sum(dopants.values()) or 1.0
    high_val_frac = sum(dopants[e] for e, s in assignment.items() if s >= 5) / total
    has_low = any(s <= 2 for s in assignment.values())

    d = 0.0
    d += 0.9 * (na_delta / NA_DELTA_MAX)                      # larger Na vacancy
    d += 0.6 * (abs(mean_v - V_NOMINAL) / 1.5)               # more extreme mean V
    d += 0.6 * high_val_frac                                  # high-valence dopant load
    if high_val_frac > 0 and not has_low:
        d += 0.25                                            # no low-valence compensation
    # penalty for needing very disparate states simultaneously
    spread = (max(assignment.values()) - min(assignment.values())) if assignment else 0
    d += 0.05 * spread
    return max(0.0, min(1.0, d))
