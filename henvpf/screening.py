"""
Chemical hard screening (main design 2).

These are *deterministic* rejections applied before any scoring:
  2.1 charge balance feasibility       (uses the charge solver)
  2.2 M-site total occupancy == 2
  2.3 ionic radius deviation < 15%
  2.4 configurational entropy floor
  2.5 toxicity exclusion
  2.6 theoretical-capacity floor + V fraction floor
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from .elements import ELEMENTS, TOXIC_EXCLUDE, M_SITE_TOTAL
from .charge_solver import solve, ChargeResult
from .scores import theoretical_capacity, radius_mismatch_risk
from .composition import Candidate

# main design 2.x default thresholds
RADIUS_DEV_MAX = 0.15
SCONFIG_FLOOR_PRACTICAL = 1.0          # practical multi-doping: Sconfig > ~1.0-1.2 R
CAPACITY_FLOOR = 110.0                 # mAh/g (main design 2.6)
V_FRACTION_FLOOR = 0.60                # V occupies >= 60% of M site


@dataclass
class HardScreenResult:
    passed: bool
    reasons: List[str] = field(default_factory=list)
    charge: ChargeResult = None


def hard_screen(cand: Candidate, enforce_entropy_floor: bool = False) -> HardScreenResult:
    reasons: List[str] = []

    # 2.5 toxicity
    for e in cand.dopants:
        if e in TOXIC_EXCLUDE:
            reasons.append(f"toxic_element:{e}")

    # 2.2 M-site occupancy
    occ = cand.v_amount + sum(cand.dopants.values())
    if abs(occ - M_SITE_TOTAL) > 1e-6:
        reasons.append(f"m_site_occupancy!=2 ({occ:.4f})")

    # 2.6 V fraction floor
    if cand.v_amount / M_SITE_TOTAL < V_FRACTION_FLOOR:
        reasons.append("v_fraction<0.60")

    # 2.1 charge balance feasibility
    charge = solve(cand.v_amount, cand.dopants)
    if not charge.valence_feasible_flag:
        reasons.append("charge_infeasible")

    # 2.3 ionic radius deviation
    rmm = radius_mismatch_risk(cand.v_amount, charge, cand.dopants)
    if rmm >= 1.0:  # rmm normalised so that 15% deviation -> 1.0
        reasons.append("radius_deviation>=15%")

    # 2.6 theoretical capacity floor
    cap = theoretical_capacity(cand.v_amount, cand.dopants)
    if cap < CAPACITY_FLOOR:
        reasons.append(f"capacity<{CAPACITY_FLOOR} ({cap:.1f})")

    # 2.4 entropy floor (off by default: V-dominant design is medium-entropy)
    if enforce_entropy_floor and cand.sconfig_over_R < SCONFIG_FLOOR_PRACTICAL:
        reasons.append("sconfig<1.0R")

    return HardScreenResult(passed=(len(reasons) == 0), reasons=reasons, charge=charge)
