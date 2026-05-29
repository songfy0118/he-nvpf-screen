# he-nvpf-screen

**V-dominant high-entropy NVPF cathode — ML pre-screening pipeline**
**V 主导高熵 NVPF 正极 — 机器学习预筛管线**

A runnable Python implementation of the *generate → chemical hard-screen →
proxy scoring → rule-calibrated ranking* front-end of a high-entropy
Na₃V₂(PO₄)₂F₃ (NVPF) cathode screening workflow. The compute-heavy back-end
(MLIP relaxation, VASP labels, Bayesian optimisation) is provided as documented
integration stubs.

本仓库实现了高熵 Na₃V₂(PO₄)₂F₃（NVPF）正极筛选流程中**现在就能跑、不需要
GPU/DFT** 的前端：组成生成 → 化学硬筛 → proxy 打分 → rule-calibrated 排序。
需要算力的后端（MLIP 弛豫、VASP 标签、贝叶斯优化）以接口桩形式给出。

> **What this is / 定位.** A **rule-calibrated ranking model**, *not* a supervised
> predictor. Outputs are **proxy / priority scores** — it does **not** claim to
> predict formation energy, energy-above-hull, F-vacancy formation energy, cycle
> life, or real capacity retention. All scores must be re-calibrated once DFT /
> experimental labels exist.
>
> 这是一个 **rule-calibrated ranking model**，不是监督预测器。输出是
> **proxy / priority 分**，不声称预测 formation energy、Ehull、F 空位形成能、
> 循环寿命或真实容量保持率。所有打分都应在拿到 DFT / 实验标签后重新校准。

---

## Quick start / 快速开始

```bash
pip install pymatgen

# 1) auto-generate and score the candidate pool  /  自动生成并打分候选池
python scripts/run_screen.py --out candidates.csv --top 25

# 2) score YOUR OWN list of compositions          /  给你自己的配方列表打分
python scripts/run_screen.py --input examples/my_compositions.csv --out my_scored.csv

# 3) sanity-check against real published HE-NVPF   /  用真实已发表配方做验证
python scripts/validate_literature.py
```

---

## Where do the thousands of compositions come from? / 那成千上万种配方哪来的？

**You do not type them.** `henvpf/composition.py` *enumerates* them
combinatorially from three inputs:

1. the element pool (24 M-site dopants + V),
2. the total-dopant grid `x ∈ {0.1 … 0.6}` (because **V ≥ 1.4**),
3. the dopant-count rule coupled to `x`.

The default run generates ~48,000 candidates in a few seconds. You only choose
*parameters* (which element tiers, the grid step) — never the formulas themselves.

**你不需要手动输入配方。** `henvpf/composition.py` 会按规则用排列组合**自动枚举**，
输入只有三样：元素池（24 个 M 位掺杂元素 + V）、总掺杂量网格 `x ∈ {0.1…0.6}`
（因为 **V ≥ 1.4**）、以及掺杂元素数随 `x` 的耦合规则。默认跑几秒生成约 48,000 个候选。
你只选**参数**（用哪几档元素、网格步长），从不手写配方。

If you *do* have specific compositions, feed a CSV with `--input` (see
`examples/my_compositions.csv`).
如果你确实有特定配方，用 `--input` 喂一个 CSV 即可（见 `examples/my_compositions.csv`）。

---

## What is implemented (runs now) / 已实现（可运行）

| Stage / 阶段 | Module / 模块 | Spec / 对应文档 |
|---|---|---|
| 1. Composition generation / 组成生成 | `henvpf/composition.py` | main §1A / addendum §1–3 |
| 2. Chemical hard screen / 化学硬筛 | `henvpf/screening.py` | main §2 |
| — Explicit charge solver / 显式电荷求解器 | `henvpf/charge_solver.py` | addendum §4 |
| 3. Proxy scores / 三个 proxy 打分 | `henvpf/scores.py` | addendum §5/§6/§7 |
| Ranking + triage / 排序 + 分流 | `henvpf/ranking.py` | addendum §8 |
| Orchestrator / 编排 | `henvpf/pipeline.py` | main §1–8 |
| User input / 用户输入 | `henvpf/userinput.py` | — |

Key rules / 关键规则:

- **V ≥ 1.4** → `0.1 ≤ x ≤ 0.6`; V occupies ≥ 70 % of the M site.
  V 占 M 位 ≥ 70%。
- **Dopant count coupled to x** / 掺杂数耦合 x (addendum §3):
  `x ≤ 0.2 → 3`, `x ≥ 0.3 → 4–5`; low-dose 5-dopant sets are flagged
  `low_dose_multi_dopant_risk`.
- **Entropy bands** / 构型熵分级 (addendum §1):
  `<1.0R co-doping / 1.0–1.5R medium / ≥1.5R strict high-entropy`. Entropy is a
  ranking tie-breaker only. 熵只作排序 tie-breaker。
- **Explicit charge solver** / 显式电荷求解器 (addendum §4): searches allowed
  oxidation states × mean-V valence (2.5–4.5) × Na vacancy (0–0.3); outputs
  `charge_balance_mode`, `na_delta_required`, `mean_v_valence`,
  `valence_feasible_flag`, `charge_compensation_difficulty`.
- **Three proxy scores** follow the addendum formulas/weights exactly:
  三个 proxy 分完全按 addendum 公式权重：
  `phase_purity_score`, `f_retention_score`,
  `v_activity_retention_score` (renamed from `capacity_retention_score`).
- **Radius mismatch** uses the concentration-weighted HEA δ parameter, so dilute
  dopants contribute little; δ = 0.15 maps to the 15 % criterion.
  离子半径失配用浓度加权 δ 参数。
- **Hard gates** / 硬门槛: `reject / exploration / priority`; gate-rejected but
  borderline candidates go to `active_learning_pool`.

## Not yet implemented (integration stubs) / 未实现（接口桩）

In `integrations/`, ordered per **addendum §9 (MLIP before VASP)**:
`integrations/` 下，严格按 **addendum §9「先 MLIP 后 VASP」** 顺序：

- `mlip_relax.py` — CHGNet / M3GNet (matgl) / MACE configuration search (steps 1–3)
- `dft_labels.py` — VASP refinement + key labels incl. F-vacancy energy (steps 4–6)
- `mp_pretrain.py` — three-tier training data (Materials Project / NASICON / NVPF)
- `bayes_opt.py` — multi-objective Bayesian optimisation (step 7; BoTorch+Ax / pymoo)

---

## Output columns / 输出列

`candidates.csv`: `formula, x, n_dopants, v_amount, dopants, sconfig_over_R,
entropy_band, low_dose_multi_dopant_risk, priority_score, triage,
phase_purity_score, f_retention_score, v_activity_retention_score,
theoretical_capacity, molar_mass, naf_risk, mfx_risk, mpo4_risk,
nasicon_competition_risk, radius_mismatch_risk, charge_compensation_difficulty,
phase_gate, f_gate, v_gate`.

---

## Known calibration limitations / 已知校准缺陷

Honest notes, to be fixed once DFT labels exist. 待 DFT 标签出来后修正。

1. **The F-retention proxy rejects pristine NVPF** as a control, because it
   rewards dopant-driven M–F stabilisation and lacks a baseline term for the
   intrinsic V–F bond. De-sodiated F stability is the weakest proxy before DFT
   (addendum §6 warns of this).
   **F-retention proxy 会把原始 NVPF 对照判成 reject**：它缺一个本征 V–F 键基线项；
   脱钠态 F 稳定性 DFT 前本就最弱。
2. **The dopant-count rule deprioritises a known-best experimental material**:
   ACS Nano 2025's V₁.₉(Fe,Ni,Co,Mg,Cr)₀.₀₂ is a 5-dopant x = 0.1 formula, which
   addendum §3 (x ≤ 0.2 → 3 dopants) flags and demotes to `exploration`. The rule
   is conservative; real winners cluster at very low dose, so the x ceiling and
   dopant-count rule may need loosening toward the low-dose end.
   **掺杂数规则会让已知最佳实验材料降级**：那个 5 元 0.02 配方被打 risk 标记降到
   exploration。规则偏保守，实测 work 的都在极低剂量端。
3. **Main-spec §7 vs addendum §8 tension**: §7 applies physical-unit thresholds
   (Eform < 0, Ehull ≤ 0.025 …) to model outputs, but round 1 has no calibrated
   predictor. This implementation follows the addendum — round 1 ranks on proxy
   scores only; physical-unit thresholds apply after `dft_labels.py`.
   **主文档 §7 与 addendum §8 矛盾**：本实现按 addendum，第一轮只用 proxy 排序。

---

## Repository structure / 仓库结构

```
he-nvpf-screen/
├── henvpf/              # core package (runs now) / 核心包（可运行）
│   ├── elements.py          # element pool + curated chemistry table
│   ├── composition.py       # Stage-1 generation + Sconfig
│   ├── charge_solver.py     # explicit charge solver (addendum §4)
│   ├── screening.py         # chemical hard screen (main §2)
│   ├── scores.py            # 3 proxy scores (addendum §5/6/7)
│   ├── ranking.py           # rule-calibrated ranking + triage (§8)
│   ├── userinput.py         # score user-supplied CSV
│   └── pipeline.py          # orchestrator
├── integrations/        # documented stubs (need compute/data) / 接口桩
├── scripts/
│   ├── run_screen.py        # CLI
│   └── validate_literature.py
├── examples/my_compositions.csv
├── requirements.txt
└── README.md
```

## Open-source stack / 开源栈

pymatgen · CHGNet / matgl / mace-torch · mp-api · atomate2 / custodian ·
CatBoost / XGBoost / LightGBM / scikit-learn · BoTorch + Ax / pymoo / Optuna.

## License

MIT (see `LICENSE`).
