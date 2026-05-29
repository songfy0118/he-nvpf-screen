"""
Multi-objective Bayesian optimisation of element ratios  (addendum 9, step 7;
main design 11). Runs INSIDE element combinations already validated by DFT.

Variables : xV, xM1..xM5, xNa     (subject to M-site sum = 2, charge balance)
Maximise  : main-phase prob, F-vacancy E, capacity, voltage, energy density,
            Sconfig, Na diffusivity, structural reversibility
Minimise  : Ehull, decomposition tendency, impurity risk, F-escape risk,
            volume change, Na migration barrier, polyhedral distortion,
            inactive-element fraction

Open-source options:
    BoTorch + Ax (qNEHVI / qEHVI for multi-objective)   <- recommended
    pymoo (NSGA-II) if you prefer evolutionary Pareto search
    Optuna (NSGAIISampler) for a lightweight start

Each round proposes a few points -> MLIP/DFT validate -> feed back (active loop).
Then locally refine composition step 0.1 -> 0.05 / 0.02 (addendum 2 Stage-2,
addendum 9 step 8) before secondary DFT on 5-10 candidates.
"""
def make_ax_experiment(validated_region):
    raise NotImplementedError("pip install ax-platform botorch; build qNEHVI loop.")
