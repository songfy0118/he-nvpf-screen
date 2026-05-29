"""
MLIP structure relaxation + configuration search  (addendum 9, steps 1-3).

This is the FIRST compute-heavy stage and runs BEFORE any VASP, per the addendum's
revised calculation order. It is a documented stub: it shows exactly where the
open-source universal MLIPs plug in. Install one of:

    pip install chgnet            # CHGNet
    pip install matgl             # M3GNet (and others)
    pip install mace-torch        # MACE

Workflow per surviving candidate composition:
    1. build a few symmetry-distinct supercell decorations of the NVPF
       Na3(V,M..)2(PO4)2F3 framework (V/M arrangements, Na-vacancy orderings,
       de-sodiated states, local F environments);
    2. relax each with the MLIP;
    3. keep the low-energy, chemically distinct configurations for VASP.

CAVEAT (surfaced during screening): universal MLIPs trained mostly on
near-equilibrium, neutral Materials Project relaxations can be unreliable for
F environments, charged/de-sodiated states, and F-vacancy energetics. Spot-check
MLIP ranking against DFT, or fine-tune on a small NVPF DFT set, before trusting
the lowest-energy configuration selection.
"""
from __future__ import annotations
from typing import Dict, List


def build_configurations(v_amount: float, dopants: Dict[str, float],
                         n_configs: int = 30) -> List["Structure"]:  # noqa: F821
    """Generate representative orderings on the NVPF framework.

    Implement with pymatgen + an ordering enumerator, e.g.:
        from pymatgen.core import Structure
        from pymatgen.transformations.advanced_transformations import (
            EnumerateStructureTransformation, OrderDisorderedStructureTransformation)
    Start from a pristine NVPF CIF (Materials Project mp-id for Na3V2(PO4)2F3),
    substitute partial occupancies on the V site, then enumerate.
    """
    raise NotImplementedError(
        "Provide a pristine NVPF Structure and an ordering enumerator; "
        "see docstring. Requires pymatgen + a CIF/mp-id.")


def relax_with_chgnet(structures: List["Structure"]):  # noqa: F821
    """Relax structures and return (relaxed_structure, energy_per_atom) pairs.

        from chgnet.model import StructOptimizer
        opt = StructOptimizer()
        return [(opt.relax(s)["final_structure"], opt.relax(s)["trajectory"].energies[-1])
                for s in structures]
    """
    raise NotImplementedError("pip install chgnet, then implement as in docstring.")


def select_low_energy(relaxed, per_composition: int = 5):
    """Keep the lowest-energy + chemically distinct configs per composition."""
    raise NotImplementedError
