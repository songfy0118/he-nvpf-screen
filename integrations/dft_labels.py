"""
VASP refinement + key-label computation  (addendum 9, steps 4-6).

Runs only on the low-energy configurations selected by mlip_relax.py. Computes
the labels that the rule-calibrated proxies should eventually be replaced by:

    formation energy            Eform
    energy above hull           Ehull          (needs MP phase diagram)
    competing-phase decomp E    (vs NaF, MFx, VF3, V2O3/5, MPO4, NVP, NaVPO4F, ...)
    F vacancy formation energy  (pristine AND de-sodiated -> the key F label)
    de-sodiated F stability
    average voltage             (from Na-full vs Na-empty energies)
    volume change               (de-sodiation)
    Na migration barrier        (NEB)

Recommended open-source stack:
    pymatgen + custodian + atomate2 / FireWorks   (VASP workflow management)
    pymatgen.analysis.phase_diagram               (Ehull vs Materials Project)
    pymatgen.analysis.transition_state            (NEB parsing)

These labels feed back into ranking.py to turn the rule-calibrated ranker into a
genuinely supervised model (addendum 8).
"""

def competing_phases():
    return ["NaF", "MFx", "VF3", "V2O3", "V2O5", "MPO4",
            "Na3V2(PO4)3", "NaVPO4F", "Na3M2(PO4)3", "MOx"]

def compute_labels(structure, mp_api_key=None):
    raise NotImplementedError(
        "Wire up atomate2 VASP workflows; parse with pymatgen; "
        "build phase diagram from Materials Project for Ehull.")
