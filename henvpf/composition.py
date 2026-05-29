"""
Stage-1 / Stage-2 composition generation for V-dominant high-entropy NVPF.

Implements:
  * V >= 1.4 constraint  ->  0.1 <= x <= 0.6        (main design A)
  * coarse step 0.1, fine step 0.05 / 0.02          (addendum 2)
  * dopant-count coupled to total dopant amount x   (addendum 3)
  * configurational entropy Sconfig + entropy bands (addendum 1)

Formula:  Na3 (V_{2-x} D1_a D2_b ...) (PO4)2 F3,  sum of dopant amounts = x.
"""
from __future__ import annotations
import itertools
import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Tuple

from .elements import DOPANT_POOL, M_SITE_TOTAL, priority_tier

R = 8.314462618  # gas constant; Sconfig is reported in units of R below

# Entropy bands (addendum 1)
ENTROPY_MULTI_DOPING = "multi_element_co_doping"          # Sconfig < 1.0 R
ENTROPY_MEDIUM = "medium_entropy"                          # 1.0 R <= Sconfig < 1.5 R
ENTROPY_HIGH = "strict_high_entropy"                       # Sconfig >= 1.5 R


@dataclass
class Candidate:
    v_amount: float
    dopants: Dict[str, float]                  # symbol -> amount on the M site
    x: float                                   # total dopant amount
    composition_step: float
    sconfig_over_R: float = 0.0
    entropy_band: str = ""
    n_dopants: int = 0
    low_dose_multi_dopant_risk: bool = False
    flags: List[str] = field(default_factory=list)

    @property
    def site_occupancy(self) -> Dict[str, float]:
        d = {"V": self.v_amount}
        d.update(self.dopants)
        return d

    def formula(self) -> str:
        parts = [f"V{self.v_amount:g}"]
        parts += [f"{el}{amt:g}" for el, amt in sorted(self.dopants.items())]
        return "Na3(" + "".join(parts) + ")(PO4)2F3"


def sconfig_over_R(site_occupancy: Dict[str, float]) -> float:
    """Configurational entropy on the mixing (M) site, in units of R.

    Sconfig = -R * sum_i f_i ln f_i, where f_i is the mole fraction of element i
    on the M site (occupancy normalised so it sums to 1).
    """
    total = sum(site_occupancy.values())
    if total <= 0:
        return 0.0
    s = 0.0
    for amt in site_occupancy.values():
        f = amt / total
        if f > 0:
            s -= f * math.log(f)
    return s  # already in units of R (R cancels)


def entropy_band(sconfig_R: float) -> str:
    if sconfig_R >= 1.5:
        return ENTROPY_HIGH
    if sconfig_R >= 1.0:
        return ENTROPY_MEDIUM
    return ENTROPY_MULTI_DOPING


def allowed_dopant_counts(x: float) -> List[int]:
    """Dopant-count rule coupled to x (addendum 3)."""
    if x <= 0.2:
        return [3]
    if x >= 0.3:
        return [4, 5]
    return [3, 4]  # transition region 0.2 < x < 0.3


def _amount_partitions(x: float, n: int, step: float) -> Iterable[Tuple[float, ...]]:
    """All ways to split total amount x into n positive parts on a `step` grid."""
    units = round(x / step)
    n_units = n
    # compositions of `units` into `n` positive integers
    if units < n_units:
        return
    for combo in itertools.combinations(range(1, units), n_units - 1):
        cuts = (0,) + combo + (units,)
        yield tuple(round((cuts[i + 1] - cuts[i]) * step, 6) for i in range(n_units))


def generate(
    x_values: Iterable[float] = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    composition_step: float = 0.1,
    dopant_pool: List[str] = None,
    restrict_tiers: Tuple[int, ...] = (1, 2),
    max_candidates: int = None,
    equal_split_only: bool = True,
) -> List[Candidate]:
    """Generate the Stage-1 candidate pool.

    equal_split_only=True keeps the coarse scan tractable by giving each chosen
    dopant set an (near-)equal split of x; set False to enumerate every grid
    partition (this is what produces the ~10^5-10^6 coarse pool, addendum 2).
    """
    pool = dopant_pool or [d for d in DOPANT_POOL if priority_tier(d) in restrict_tiers]
    out: List[Candidate] = []
    for x in x_values:
        v_amount = round(M_SITE_TOTAL - x, 6)
        for n in allowed_dopant_counts(x):
            if n > len(pool):
                continue
            for dset in itertools.combinations(sorted(pool), n):
                if equal_split_only:
                    amt = round(x / n, 6)
                    partitions = [tuple([amt] * n)]
                    # fix rounding so the sum is exactly x
                    drift = round(x - amt * n, 6)
                    if abs(drift) > 1e-9:
                        partitions = [tuple([amt] * (n - 1) + [round(amt + drift, 6)])]
                else:
                    partitions = list(_amount_partitions(x, n, composition_step))
                for parts in partitions:
                    dopants = {el: a for el, a in zip(dset, parts)}
                    occ = {"V": v_amount, **dopants}
                    s = sconfig_over_R(occ)
                    min_amt = min(parts)
                    cand = Candidate(
                        v_amount=v_amount,
                        dopants=dopants,
                        x=round(x, 6),
                        composition_step=composition_step,
                        sconfig_over_R=round(s, 4),
                        entropy_band=entropy_band(s),
                        n_dopants=n,
                        low_dose_multi_dopant_risk=(n >= 5 and min_amt < 0.03),
                    )
                    out.append(cand)
                    if max_candidates and len(out) >= max_candidates:
                        return out
    return out
