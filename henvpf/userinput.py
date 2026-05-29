"""
Load user-supplied compositions so they can be scored directly, instead of
(or in addition to) the auto-generated Stage-1 pool.

Accepted CSV formats (header required):

  A) compact  -- one column named `composition`, element:amount pairs:
        composition
        V:1.9;Cr:0.0333;Ti:0.0333;Zr:0.0334
        V:1.95;Mn:0.01;Ti:0.01;Mg:0.01;Cr:0.01;Zr:0.01

  B) wide     -- one column per element (V plus any dopants); blanks = 0:
        V,Cr,Ti,Zr
        1.9,0.0333,0.0333,0.0334

The M-site amounts should sum to 2 (a warning is printed otherwise; the
hard screen will flag it as m_site_occupancy!=2).
"""
from __future__ import annotations
import csv
from typing import Dict, List

from .composition import Candidate, sconfig_over_R, entropy_band


def _make_candidate(occ: Dict[str, float]) -> Candidate:
    v = occ.get("V", 0.0)
    dop = {e: a for e, a in occ.items() if e != "V" and a > 0}
    s = sconfig_over_R({"V": v, **dop})
    min_amt = min(dop.values()) if dop else 1.0
    return Candidate(
        v_amount=round(v, 6),
        dopants={e: round(a, 6) for e, a in dop.items()},
        x=round(sum(dop.values()), 6),
        composition_step=0.0,
        sconfig_over_R=round(s, 4),
        entropy_band=entropy_band(s),
        n_dopants=len(dop),
        low_dose_multi_dopant_risk=(len(dop) >= 5 and min_amt < 0.03),
    )


def load_compositions_csv(path: str) -> List[Candidate]:
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        fields = [c.strip() for c in (reader.fieldnames or [])]
        compact = "composition" in fields
        out: List[Candidate] = []
        for row in reader:
            occ: Dict[str, float] = {}
            if compact:
                spec = (row.get("composition") or "").strip()
                if not spec:
                    continue
                for pair in spec.replace(",", ";").split(";"):
                    pair = pair.strip()
                    if not pair:
                        continue
                    el, amt = pair.split(":")
                    occ[el.strip()] = float(amt)
            else:
                for el in fields:
                    val = (row.get(el) or "").strip()
                    if val:
                        occ[el.strip()] = float(val)
            if occ:
                out.append(_make_candidate(occ))
    return out
