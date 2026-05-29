"""
Three-tier training-data assembly + pretraining  (main design 5; addendum 8).

Tier 1  general electrode data (Materials Project)      -> base stability/voltage model
Tier 2  NASICON / polyanion / fluorophosphate data      -> transfer / pretrain
Tier 3  self-built NVPF-specific labels (from dft_labels)-> fine-tune

    from mp_api.client import MPRester
    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(
            elements=["Na","V","P","O","F"], fields=["material_id","formula_pretty",
            "formation_energy_per_atom","energy_above_hull","volume"])

NOTE: the label that matters most for NVPF -- F-vacancy formation energy in the
de-sodiated state -- has the LEAST transferable public training data. Expect weak
transfer on exactly that target; prioritise Tier-3 DFT there.

Models to compare (main design 6): CatBoost / XGBoost / LightGBM / RandomForest
for tabular features; GNN (matgl/MEGNet) once structures are available.
"""
def fetch_polyanion_training_set(api_key):
    raise NotImplementedError("pip install mp-api; implement query as in docstring.")
