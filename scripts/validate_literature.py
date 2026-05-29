#!/usr/bin/env python3
"""Sanity-check the proxy scores against real, published high-entropy NVPF.

These compositions were reported as successfully synthesised, single-main-phase,
and electrochemically good. A useful first-version proxy should NOT reject them.
(References are experimental papers, 2024-2026.)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import warnings; warnings.filterwarnings("ignore")

from henvpf import charge_solver, scores, ranking
from henvpf.composition import Candidate, sconfig_over_R, entropy_band

LITERATURE = [
    # (label, V_amount, {dopant: amount})
    ("ACS Nano 2025  V1.9(Fe,Ni,Co,Mg,Cr)0.02",
     1.9, {"Fe": 0.02, "Ni": 0.02, "Co": 0.02, "Mg": 0.02, "Cr": 0.02}),
    ("J Mater Chem A 2025  V1.9(Ca,Cr,Ti,Nb,Mo,Na)0.1/6 [Ca,Na->skip]",
     1.9, {"Cr": 0.1/6, "Ti": 0.1/6, "Nb": 0.1/6, "Mo": 0.1/6}),
    ("ScienceDirect 2025  V1.95(Ti,Mg,Mn,Cr,Zr)0.01",
     1.95, {"Ti": 0.01, "Mg": 0.01, "Mn": 0.01, "Cr": 0.01, "Zr": 0.01}),
    ("ZIB 2024  HE-NVPF (Al,Zn,Mn,Cr,Nb)",
     1.5, {"Al": 0.1, "Zn": 0.1, "Mn": 0.1, "Cr": 0.1, "Nb": 0.1}),
    ("pristine NVPF (control)", 2.0, {}),
]

print(f"{'composition':<52}{'Sconf':>6}{'PP':>6}{'FR':>6}{'VA':>6}{'cap':>7}{'pri':>6}  triage")
print("-" * 96)
for label, v, dop in LITERATURE:
    occ = {"V": v, **dop}
    s = sconfig_over_R(occ)
    cand = Candidate(v_amount=v, dopants=dop, x=round(sum(dop.values()), 4),
                     composition_step=0.0, sconfig_over_R=round(s, 4),
                     entropy_band=entropy_band(s), n_dopants=len(dop),
                     low_dose_multi_dopant_risk=(len(dop) >= 5 and (min(dop.values()) if dop else 1) < 0.03))
    ch = charge_solver.solve(v, dop)
    b = scores.score_candidate(v, dop, ch)
    ps = ranking.priority_score(b, cand)
    tg = ranking.triage(b, cand)
    print(f"{label:<52}{s:>6.2f}{b.phase_purity_score:>6.2f}{b.f_retention_score:>6.2f}"
          f"{b.v_activity_retention_score:>6.2f}{b.theoretical_capacity:>7.1f}{ps:>6.3f}  {tg}")
