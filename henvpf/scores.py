"""
Risk sub-components and the three rule-calibrated proxy scores
(addendum 5 phase purity, 6 F-retention, 7 V-activity retention).

IMPORTANT: these are first-version rule-calibrated proxies, NOT physical
predictions. Weights and component definitions follow the addendum exactly.
They must be recalibrated against DFT / experimental labels (addendum 8).
"""
from __future__ import annotations
import statistics
from dataclasses import dataclass, field
from typing import Dict

from .elements import (
    ELEMENTS, atomic_mass, electronegativity, shannon_radius_oct,
    NVPF_BASELINE_CAPACITY,
)
from .charge_solver import ChargeResult

F_X = electronegativity("F")          # 3.98
O_X = electronegativity("O")          # 3.44
V3_RADIUS = shannon_radius_oct("V", 3)  # ~0.64 A


def clip01(v: float) -> float:
    return max(0.0, min(1.0, v))


# ---------------------------------------------------------------------------
# Risk components (addendum 5)
# ---------------------------------------------------------------------------
@dataclass
class RiskComponents:
    naf_risk: float
    mfx_risk: float
    mpo4_risk: float
    nasicon_competition_risk: float
    radius_mismatch_risk: float
    charge_compensation_difficulty: float


def _amount_weighted(dopants: Dict[str, float], attr: str) -> float:
    tot = sum(dopants.values()) or 1.0
    return sum(amt * getattr(ELEMENTS[e], attr) for e, amt in dopants.items()) / tot


def radius_mismatch_risk(v_amount: float, charge: ChargeResult, dopants: Dict[str, float]) -> float:
    """Concentration-weighted M-site radius mismatch (HEA delta parameter).

        delta = sqrt( sum_i f_i (1 - r_i / r_avg)^2 )

    where f_i is the M-site mole fraction (V included) and r_avg the
    concentration-weighted mean octahedral radius. Low-dose dopants therefore
    contribute little, matching how dilute substitution actually behaves.
    Normalised so that delta = 0.15 (the main-design 15% criterion) -> 1.0.
    """
    import math
    assign = charge.dopant_valence_assignment or {e: ELEMENTS[e].ox_states[0] for e in dopants}
    items = [(V3_RADIUS, v_amount)]
    for e, amt in dopants.items():
        items.append((shannon_radius_oct(e, assign.get(e, ELEMENTS[e].ox_states[0])), amt))
    total = sum(w for _, w in items) or 1.0
    r_avg = sum(r * w for r, w in items) / total
    if r_avg <= 0:
        return 1.0
    delta = math.sqrt(sum((w / total) * (1 - r / r_avg) ** 2 for r, w in items))
    return clip01(delta / 0.15)


def compute_risks(v_amount: float, dopants: Dict[str, float], charge: ChargeResult) -> RiskComponents:
    naf = clip01(_amount_weighted(dopants, "fluoride_former") * 0.6
                 + (1 - _amount_weighted(dopants, "f_stabilizer")) * 0.4)
    mfx = clip01(_amount_weighted(dopants, "fluoride_former"))
    mpo4 = clip01(_amount_weighted(dopants, "phosphate_former"))
    # NASICON (NVP-type, F-free) competition rises when F is poorly stabilised
    nasicon = clip01(0.5 * (1 - _amount_weighted(dopants, "f_stabilizer"))
                     + 0.5 * naf)
    rmm = radius_mismatch_risk(v_amount, charge, dopants)
    return RiskComponents(
        naf_risk=round(naf, 4),
        mfx_risk=round(mfx, 4),
        mpo4_risk=round(mpo4, 4),
        nasicon_competition_risk=round(nasicon, 4),
        radius_mismatch_risk=round(rmm, 4),
        charge_compensation_difficulty=charge.charge_compensation_difficulty,
    )


# ---------------------------------------------------------------------------
# Score 1: phase_purity_score (addendum 5)
# ---------------------------------------------------------------------------
def phase_purity_score(r: RiskComponents) -> float:
    s = (1.0
         - 0.20 * r.naf_risk
         - 0.20 * r.mfx_risk
         - 0.15 * r.mpo4_risk
         - 0.15 * r.nasicon_competition_risk
         - 0.15 * r.radius_mismatch_risk
         - 0.15 * r.charge_compensation_difficulty)
    return round(clip01(s), 4)


# ---------------------------------------------------------------------------
# Score 2: f_retention_score (addendum 6)
# ---------------------------------------------------------------------------
def _mf_bond_strength_score(dopants: Dict[str, float], charge: ChargeResult) -> float:
    """Proxy: ionic M-F interaction grows with cation charge/radius and with the
    M-F vs M-O electronegativity contrast."""
    assign = charge.dopant_valence_assignment or {e: ELEMENTS[e].ox_states[0] for e in dopants}
    tot = sum(dopants.values()) or 1.0
    acc = 0.0
    for e, amt in dopants.items():
        z = assign.get(e, ELEMENTS[e].ox_states[0])
        r = shannon_radius_oct(e, z)
        field_strength = z / (r * r)               # cation field strength
        ionicity = (F_X - electronegativity(e)) / F_X
        acc += amt * (0.5 * clip01(field_strength / 8.0) + 0.5 * clip01(ionicity))
    return clip01(acc / tot)


def _high_valence_stabilization_score(dopants: Dict[str, float], charge: ChargeResult) -> float:
    assign = charge.dopant_valence_assignment or {e: ELEMENTS[e].ox_states[0] for e in dopants}
    tot = sum(dopants.values()) or 1.0
    hv_frac = sum(amt for e, amt in dopants.items() if assign.get(e, 3) >= 5) / tot
    # capped: benefit saturates so excess high-valence content does not keep helping
    return clip01(min(hv_frac, 0.35) / 0.35 * 0.8)


def f_retention_score(dopants: Dict[str, float], charge: ChargeResult, r: RiskComponents) -> float:
    mf = _mf_bond_strength_score(dopants, charge)
    hv = _high_valence_stabilization_score(dopants, charge)
    low_naf = 1 - r.naf_risk
    low_mfx = 1 - r.mfx_risk
    low_dist = 1 - r.radius_mismatch_risk
    # desodiated_f_stability_proxy: WEAK before DFT (addendum 6 warning). Conservative.
    desod = 0.5 * mf + 0.5 * (1 - r.naf_risk)
    s = (0.25 * mf
         + 0.15 * hv
         + 0.15 * low_naf
         + 0.15 * low_mfx
         + 0.15 * low_dist
         + 0.15 * desod)
    return round(clip01(s), 4)


# ---------------------------------------------------------------------------
# Score 3: v_activity_retention_score (addendum 7, renamed from capacity_retention)
# ---------------------------------------------------------------------------
def formula_molar_mass(v_amount: float, dopants: Dict[str, float]) -> float:
    m = 3 * atomic_mass("Na") + 2 * atomic_mass("P") + 8 * atomic_mass("O") + 3 * atomic_mass("F")
    m += v_amount * atomic_mass("V")
    m += sum(amt * atomic_mass(e) for e, amt in dopants.items())
    return m


def theoretical_capacity(v_amount: float, dopants: Dict[str, float], reversible_na: float = 2.0) -> float:
    """mAh/g for `reversible_na` one-electron Na (de)insertion. 26801.5 = F/3.6."""
    m = formula_molar_mass(v_amount, dopants)
    return reversible_na * 26801.5 / m


def v_activity_retention_score(v_amount: float, dopants: Dict[str, float]) -> float:
    v_frac = v_amount / 2.0                                   # fraction of M site that is V
    # V >= 1.4 -> v_frac in [0.7, 1.0]; map to [0,1] with 0.7 -> 0.5, 1.0 -> 1.0
    v_fraction_score = clip01((v_frac - 0.6) / 0.4)
    cap = theoretical_capacity(v_amount, dopants)
    theoretical_capacity_score = clip01(cap / NVPF_BASELINE_CAPACITY)
    tot_m = 2.0
    inactive = sum(amt for e, amt in dopants.items() if ELEMENTS[e].is_inactive) / tot_m
    heavy = sum(amt for e, amt in dopants.items() if ELEMENTS[e].is_heavy) / tot_m
    low_inactive_fraction_score = clip01(1 - inactive / 0.3)
    low_heavy_element_fraction_score = clip01(1 - heavy / 0.1)   # main design: heavy cap ~5-10%
    reversible_na_proxy = clip01(0.9 - 0.5 * inactive)
    s = (0.35 * v_fraction_score
         + 0.25 * theoretical_capacity_score
         + 0.15 * low_inactive_fraction_score
         + 0.15 * low_heavy_element_fraction_score
         + 0.10 * reversible_na_proxy)
    return round(clip01(s), 4)


# ---------------------------------------------------------------------------
# Hard gates (addendum 5/6/7)
# ---------------------------------------------------------------------------
def gate(score: float, reject_below: float, exploration_below: float) -> str:
    if score < reject_below:
        return "reject"
    if score < exploration_below:
        return "exploration"
    return "priority"


@dataclass
class ScoreBundle:
    risks: RiskComponents
    phase_purity_score: float
    f_retention_score: float
    v_activity_retention_score: float
    theoretical_capacity: float
    molar_mass: float
    phase_gate: str = ""
    f_gate: str = ""
    v_gate: str = ""

    def to_row(self) -> Dict[str, float]:
        return {
            "phase_purity_score": self.phase_purity_score,
            "f_retention_score": self.f_retention_score,
            "v_activity_retention_score": self.v_activity_retention_score,
            "theoretical_capacity": round(self.theoretical_capacity, 2),
            "molar_mass": round(self.molar_mass, 2),
            "naf_risk": self.risks.naf_risk,
            "mfx_risk": self.risks.mfx_risk,
            "mpo4_risk": self.risks.mpo4_risk,
            "nasicon_competition_risk": self.risks.nasicon_competition_risk,
            "radius_mismatch_risk": self.risks.radius_mismatch_risk,
            "charge_compensation_difficulty": self.risks.charge_compensation_difficulty,
            "phase_gate": self.phase_gate,
            "f_gate": self.f_gate,
            "v_gate": self.v_gate,
        }


def score_candidate(v_amount: float, dopants: Dict[str, float], charge: ChargeResult) -> ScoreBundle:
    risks = compute_risks(v_amount, dopants, charge)
    pp = phase_purity_score(risks)
    fr = f_retention_score(dopants, charge, risks)
    va = v_activity_retention_score(v_amount, dopants)
    bundle = ScoreBundle(
        risks=risks,
        phase_purity_score=pp,
        f_retention_score=fr,
        v_activity_retention_score=va,
        theoretical_capacity=theoretical_capacity(v_amount, dopants),
        molar_mass=formula_molar_mass(v_amount, dopants),
        phase_gate=gate(pp, 0.55, 0.70),
        f_gate=gate(fr, 0.55, 0.70),
        v_gate=gate(va, 0.60, 0.75),
    )
    return bundle
