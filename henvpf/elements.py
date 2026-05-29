"""
Element pool and curated chemistry table for V-dominant high-entropy NVPF screening.

Physical primitives (octahedral Shannon radius, atomic mass, electronegativity)
are pulled from pymatgen. Domain knowledge that pymatgen does NOT carry
(role classification, M-site screening oxidation states, fluoride / phosphate
formation tendency, cost / toxicity flags) is curated here as a first-version
table. Every heuristic here is a *screening prior*, not a physical truth, and is
meant to be recalibrated once DFT / experimental labels exist (see addendum 5-8).
"""
from __future__ import annotations
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, List, Tuple

from pymatgen.core import Element
from pymatgen.core.periodic_table import Species

# ---------------------------------------------------------------------------
# Roles (main design, section 1) -- the 24-element M-site dopant pool + V
# ---------------------------------------------------------------------------
ROLE_ACTIVE = "active"            # electrochemically active / quasi-active
ROLE_STRUCTURAL = "structural"    # high-valence structural stabilizer
ROLE_LIGHT = "light"              # light / mid-valence structural modulator
ROLE_ELECTRONIC = "electronic"    # electronic-structure / local-distortion modulator
ROLE_RARE_EARTH = "rare_earth"    # rare-earth modulator

# Anion framework (PO4)2 F3 fixed charge: 2*(-3) + 3*(-1) = -9
ANION_FRAMEWORK_CHARGE = -9.0
M_SITE_TOTAL = 2.0                # total V/M occupancy on the transition-metal site
NVPF_BASELINE_CAPACITY = 128.3    # mAh/g, 2 reversible Na in pristine NVPF


@dataclass(frozen=True)
class ElementInfo:
    symbol: str
    role: str
    # oxidation states allowed during *screening* (not a final physical assignment)
    ox_states: Tuple[int, ...]
    # first-version proxy tendencies in [0, 1]; higher = worse
    fluoride_former: float        # tendency to form a stable binary fluoride MFx
    phosphate_former: float       # tendency to form an isolated metal phosphate MPO4
    f_stabilizer: float           # tendency to STABILISE F inside the lattice (good)
    is_heavy: bool                # heavy element (mass penalty -> capacity dilution)
    is_inactive: bool             # does not contribute V-style redox capacity
    is_toxic: bool                # hard-exclude in chemical screening
    is_costly: bool               # allowed for computation, deprioritised for experiment


# Curated first-version table. Comments give the rationale used in the main design.
_RAW: List[ElementInfo] = [
    # --- active / quasi-active (main 1, group 1) ---
    ElementInfo("V",  ROLE_ACTIVE, (3, 4, 5), 0.20, 0.30, 0.70, False, False, False, False),
    ElementInfo("Cr", ROLE_ACTIVE, (3,),       0.25, 0.35, 0.55, False, False, False, False),
    ElementInfo("Mn", ROLE_ACTIVE, (2, 3, 4),  0.30, 0.35, 0.45, False, False, False, False),
    ElementInfo("Fe", ROLE_ACTIVE, (2, 3),     0.30, 0.55, 0.45, False, False, False, False),
    ElementInfo("Co", ROLE_ACTIVE, (2, 3),     0.30, 0.40, 0.45, False, False, False, False),
    ElementInfo("Ni", ROLE_ACTIVE, (2, 3),     0.35, 0.40, 0.40, False, False, False, False),
    ElementInfo("Ru", ROLE_ACTIVE, (3, 4, 5),  0.30, 0.30, 0.55, True,  False, False, True),
    # --- high-valence structural stabilizers (main 1, group 2) ---
    ElementInfo("Ti", ROLE_STRUCTURAL, (3, 4), 0.30, 0.35, 0.60, False, True,  False, False),
    ElementInfo("Zr", ROLE_STRUCTURAL, (4,),   0.35, 0.40, 0.55, False, True,  False, False),
    ElementInfo("Nb", ROLE_STRUCTURAL, (4, 5), 0.55, 0.40, 0.65, False, True,  False, False),
    ElementInfo("Mo", ROLE_STRUCTURAL, (4, 5, 6), 0.55, 0.45, 0.60, True, True, False, False),
    ElementInfo("W",  ROLE_STRUCTURAL, (4, 5, 6), 0.60, 0.45, 0.62, True, True, False, False),
    ElementInfo("Sn", ROLE_STRUCTURAL, (2, 4), 0.35, 0.45, 0.45, True,  True,  False, False),
    # --- light structural modulators (main 1, group 3) ---
    ElementInfo("Mg", ROLE_LIGHT, (2,),  0.45, 0.45, 0.45, False, True, False, False),
    ElementInfo("Al", ROLE_LIGHT, (3,),  0.55, 0.65, 0.50, False, True, False, False),
    ElementInfo("Sc", ROLE_LIGHT, (3,),  0.45, 0.55, 0.45, False, True, False, False),
    # --- electronic / local-distortion modulators (main 1, group 4) ---
    ElementInfo("Cu", ROLE_ELECTRONIC, (1, 2), 0.30, 0.40, 0.60, False, False, False, False),
    ElementInfo("Zn", ROLE_ELECTRONIC, (2,),   0.35, 0.45, 0.45, False, True,  False, False),
    ElementInfo("Ga", ROLE_ELECTRONIC, (3,),   0.50, 0.55, 0.45, False, True,  False, False),
    ElementInfo("In", ROLE_ELECTRONIC, (3,),   0.45, 0.55, 0.45, True,  True,  False, False),
    ElementInfo("Bi", ROLE_ELECTRONIC, (3,),   0.45, 0.60, 0.40, True,  True,  False, False),
    # --- rare-earth modulators (main 1, group 5) ---
    ElementInfo("Ce", ROLE_RARE_EARTH, (3, 4), 0.55, 0.60, 0.45, True, True, False, False),
    ElementInfo("La", ROLE_RARE_EARTH, (3,),   0.55, 0.65, 0.40, True, True, False, False),
    ElementInfo("Y",  ROLE_RARE_EARTH, (3,),   0.55, 0.60, 0.40, True, True, False, False),
    ElementInfo("Yb", ROLE_RARE_EARTH, (3,),   0.55, 0.60, 0.40, True, True, False, True),
]

ELEMENTS: Dict[str, ElementInfo] = {e.symbol: e for e in _RAW}

# Dopant pool = everything except V (the 24-element pool referenced in the addendum)
DOPANT_POOL: List[str] = [s for s in ELEMENTS if s != "V"]

# Hard-excluded toxic elements (main design 2.5)
TOXIC_EXCLUDE = {"Cd", "Pb", "Hg", "As"}

# Element priority tiers (main design "element priority")
PRIORITY_T1 = {"Cr", "Mn", "Fe", "Ti", "Zr", "Nb", "Mo", "W", "Al", "Mg"}
PRIORITY_T2 = {"Co", "Ni", "Cu", "Zn", "Sn", "Sc", "Ga", "Ce"}
PRIORITY_T3 = {"Ru", "In", "Bi", "La", "Y", "Yb"}


@lru_cache(maxsize=None)
def atomic_mass(symbol: str) -> float:
    return float(Element(symbol).atomic_mass)


@lru_cache(maxsize=None)
def electronegativity(symbol: str) -> float:
    x = Element(symbol).X
    return float(x) if x is not None else 1.8


@lru_cache(maxsize=None)
def shannon_radius_oct(symbol: str, valence: int) -> float:
    """Octahedral (CN VI) Shannon ionic radius in angstrom, with fallbacks."""
    try:
        r = Species(symbol, valence).get_shannon_radius("VI")
        if r:
            return float(r)
    except Exception:
        pass
    try:
        return float(Species(symbol, valence).ionic_radius)
    except Exception:
        # crude fallback from neutral atomic radius
        ar = Element(symbol).atomic_radius
        return float(ar) * 0.6 if ar else 0.75


def priority_tier(symbol: str) -> int:
    if symbol in PRIORITY_T1:
        return 1
    if symbol in PRIORITY_T2:
        return 2
    if symbol in PRIORITY_T3:
        return 3
    return 1  # V and any unlisted default to tier 1
