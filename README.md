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

## 1. Get the code / 获取代码

You only need **git** (or just a browser) and **Python ≥ 3.9**.
只需要 **git**（或一个浏览器）和 **Python ≥ 3.9**。

### Option A — clone with git (recommended) / 用 git 克隆（推荐）

```bash
git clone https://github.com/songfy0118/he-nvpf-screen.git
cd he-nvpf-screen
```

### Option B — download the ZIP (no git needed) / 下载 ZIP（不用装 git）

1. Open <https://github.com/songfy0118/he-nvpf-screen>
2. Click the green **Code** button → **Download ZIP**.
3. Unzip it, then `cd` into the unzipped `he-nvpf-screen` folder.

打开上面的链接 → 点绿色 **Code** 按钮 → **Download ZIP** → 解压 → 进入解压出的
`he-nvpf-screen` 文件夹。

---

## 2. Install & run / 安装并运行

The runnable core (stages 1–3 + ranking) needs **exactly one** package: `pymatgen`.
可运行的核心（1–3 阶段 + 排序）**只需要一个**依赖：`pymatgen`。

```bash
pip install pymatgen        # or: pip install -r requirements.txt

# 1) auto-generate and score the candidate pool  /  自动生成并打分候选池
python scripts/run_screen.py --out candidates.csv --top 25

# 2) score YOUR OWN list of compositions          /  给你自己的配方列表打分
python scripts/run_screen.py --input examples/my_compositions.csv --out my_scored.csv

# 3) sanity-check against real published HE-NVPF   /  用真实已发表配方做验证
python scripts/validate_literature.py
```

> A harmless `UserWarning: No ionic radius for Sn2+` may print — pymatgen lacks
> that one Shannon radius and the code falls back automatically. It does not
> affect results. 会打印一句无害的 `Sn2+` 半径警告，代码有 fallback，不影响结果。

---

## 3. How to verify it works / 如何确认跑成功

This was run end-to-end on a clean machine (Python 3.14, pymatgen 2026.5).
Your numbers should match closely. 下面是在干净环境实测的输出，你的结果应基本一致。

**`python scripts/run_screen.py` should print something like:**

```
generated............. 48144
passed hard screen.... 47347
triage................ {'reject': 25956, 'active_learning_pool': 17340, 'exploration': 4051}

Top 25 by priority_score:
formula                                          x nd    PP    FR    VA    cap    pri  triage
Na3(V1.9Cr0.033333Ti0.033333Zr0.033334)(PO4)2F3 0.1  3  0.73  0.63  0.93  127.9  0.694  exploration
...
wrote -> candidates.csv
```

**`python scripts/validate_literature.py` should print a table where the four
real published materials all get scored (not crashed) and the pristine NVPF
control is correctly demoted to `reject`:**

```
composition                                  Sconf    PP    FR    VA    cap   pri  triage
ACS Nano 2025  V1.9(Fe,Ni,Co,Mg,Cr)0.02       0.28  0.67  0.55  0.94  128.3 0.655  exploration
...
pristine NVPF (control)                       0.00  0.81  0.43  0.99  128.3 0.663  reject
```

✅ **Success checklist / 成功标志:**
- both scripts exit without a traceback (the `Sn2+` warning is fine);
- `candidates.csv` is created in the current folder;
- `run_screen.py` reports ~48k generated and ~47k passing the hard screen.

两个脚本不报错（`Sn2+` 警告除外）、当前目录生成了 `candidates.csv`、生成数约 4.8 万、
过筛约 4.7 万 —— 就说明跑通了。

---

## 4. What each file does / 每个文件干什么

### `henvpf/` — core package (runs now) / 核心包（现在就能跑）

| File / 文件 | What it does / 作用 |
|---|---|
| `__init__.py` | Package init; exposes submodules and `__version__`. 包入口。 |
| `elements.py` | **The chemistry knowledge base.** The 24-element M-site dopant pool + V, with curated per-element data: allowed oxidation states, fluoride/phosphate-forming tendency, F-stabilising tendency, heavy/toxic/costly flags, priority tier. Physical radii/masses come from pymatgen. 元素化学知识表：24 个掺杂元素的氧化态、成氟化物/磷酸盐倾向、F 稳定性、毒性/成本/优先级。 |
| `composition.py` | **Stage 1 — generate candidates.** Combinatorially enumerates `Na₃(V₂₋ₓ D…)(PO₄)₂F₃` formulas under the rules (`V ≥ 1.4` ⇒ `0.1 ≤ x ≤ 0.6`, dopant-count coupled to `x`), and computes configurational entropy `Sconfig` + entropy band. ~48k candidates in seconds. 组成生成：按规则枚举配方并算构型熵。 |
| `charge_solver.py` | **Explicit charge-balance solver.** For each candidate, searches dopant oxidation states × mean-V valence (2.5–4.5) × Na vacancy (0–0.3) for the lowest-difficulty neutral assignment; outputs the compensation mode + a 0–1 difficulty score. 显式电荷求解器：搜索可行的电荷平衡方案并打难度分。 |
| `screening.py` | **Stage 2 — chemical hard screen.** Deterministic pass/reject *before* any scoring: toxicity, M-site occupancy = 2, V-fraction floor, charge feasibility, ionic-radius deviation < 15%, theoretical-capacity floor. 化学硬筛：打分前的确定性淘汰。 |
| `scores.py` | **Stage 3 — the three proxy scores.** Computes risk sub-components (NaF / MFₓ / MPO₄ / NASICON-competition / radius-mismatch) and combines them into `phase_purity_score`, `f_retention_score`, `v_activity_retention_score`, plus theoretical capacity & molar mass, with hard gates. 三个 proxy 打分 + 风险子项 + 闸门。 |
| `ranking.py` | **Rule-calibrated ranking + triage.** Aggregates the three scores (weights 0.35/0.30/0.25 + 0.10 entropy tie-breaker) into one `priority_score` and sorts each candidate into `priority` / `exploration` / `reject` / `active_learning_pool`. 加权排序 + 分流。 |
| `userinput.py` | **Score your own compositions.** Parses a user CSV (compact `V:1.9;Cr:0.03;…` form *or* one-column-per-element form) into candidates so they run through the same scoring. 读用户自定义配方 CSV。 |
| `pipeline.py` | **Orchestrator.** Ties stages 1→2→3→rank together and writes the output CSV. Called by the CLI. 端到端编排器，CLI 调它。 |

### `scripts/` — command-line entry points / 命令行入口

| File / 文件 | What it does / 作用 |
|---|---|
| `run_screen.py` | The main CLI. Auto-generates+scores the pool, or scores a `--input` CSV; prints a summary + top-N table and writes the full CSV. Flags: `--out --step --tiers --full-grid --entropy-floor --input --max --top`. 主命令行。 |
| `validate_literature.py` | Runs the scorer on real published HE-NVPF compositions as a sanity check — a good proxy should not reject known-good materials. 用已发表配方自检。 |

### `integrations/` — documented stubs (need compute/data) / 接口桩（需算力/数据，未实现）

Ordered per the spec's "MLIP-before-VASP" calculation order. 按「先 MLIP 后 VASP」排序。

| File / 文件 | What it would do / 计划做什么 |
|---|---|
| `mlip_relax.py` | CHGNet / M3GNet (matgl) / MACE structure relaxation & config search. MLIP 弛豫。 |
| `dft_labels.py` | VASP refinement + key labels incl. F-vacancy formation energy. VASP 标签。 |
| `mp_pretrain.py` | Three-tier training data (Materials Project / NASICON / NVPF). 训练数据。 |
| `bayes_opt.py` | Multi-objective Bayesian optimisation (BoTorch+Ax / pymoo). 多目标贝叶斯优化。 |

### Other files / 其它文件

| File / 文件 | What it is / 是什么 |
|---|---|
| `examples/my_compositions.csv` | Example input for `--input`: a handful of compositions in the accepted CSV format. `--input` 的示例输入。 |
| `candidates.csv` | A pre-generated sample output (the ranked pool) so you can inspect results without running anything. 预生成的示例输出，不跑也能看结果。 |
| `requirements.txt` | Dependencies. Only `pymatgen` is needed for the runnable core; the rest are commented optional deps for the integration stubs. 依赖清单（核心只需 pymatgen）。 |
| `.gitignore` | Ignores caches, virtualenvs, and regenerable outputs (`my_scored.csv`, `coarse_pool.csv`, `*.zip`). |
| `LICENSE` | MIT. |

---

## 5. Where do the thousands of compositions come from? / 那成千上万种配方哪来的？

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

## 6. Key rules / 关键规则

- **V ≥ 1.4** → `0.1 ≤ x ≤ 0.6`; V occupies ≥ 70 % of the M site. V 占 M 位 ≥ 70%。
- **Dopant count coupled to x**: `x ≤ 0.2 → 3`, `x ≥ 0.3 → 4–5`; low-dose 5-dopant
  sets are flagged `low_dose_multi_dopant_risk`.
- **Entropy bands**: `<1.0R co-doping / 1.0–1.5R medium / ≥1.5R strict high-entropy`.
  Entropy is a ranking tie-breaker only. 熵只作排序 tie-breaker。
- **Explicit charge solver**: searches allowed oxidation states × mean-V valence
  (2.5–4.5) × Na vacancy (0–0.3).
- **Three proxy scores**: `phase_purity_score`, `f_retention_score`,
  `v_activity_retention_score`.
- **Hard gates**: `reject / exploration / priority`; gate-rejected but borderline
  candidates go to `active_learning_pool`.

---

## 7. Output columns / 输出列

`candidates.csv`: `formula, x, n_dopants, v_amount, dopants, sconfig_over_R,
entropy_band, low_dose_multi_dopant_risk, priority_score, triage,
phase_purity_score, f_retention_score, v_activity_retention_score,
theoretical_capacity, molar_mass, naf_risk, mfx_risk, mpo4_risk,
nasicon_competition_risk, radius_mismatch_risk, charge_compensation_difficulty,
phase_gate, f_gate, v_gate`.

---

## 8. Known calibration limitations / 已知校准缺陷

Honest notes, to be fixed once DFT labels exist. 待 DFT 标签出来后修正。

1. **The F-retention proxy rejects pristine NVPF** as a control, because it
   rewards dopant-driven M–F stabilisation and lacks a baseline term for the
   intrinsic V–F bond. De-sodiated F stability is the weakest proxy before DFT.
   **F-retention proxy 会把原始 NVPF 对照判成 reject**：缺一个本征 V–F 键基线项。
2. **The dopant-count rule deprioritises a known-best experimental material**:
   ACS Nano 2025's V₁.₉(Fe,Ni,Co,Mg,Cr)₀.₀₂ is a 5-dopant x = 0.1 formula, which
   the `x ≤ 0.2 → 3 dopants` rule flags and demotes to `exploration`. Real winners
   cluster at very low dose, so the x ceiling / dopant-count rule may need loosening.
   **掺杂数规则会让已知最佳实验材料降级**：规则偏保守，实测 work 的都在极低剂量端。
3. **Round 1 ranks on proxy scores only**; physical-unit thresholds (Eform < 0,
   Ehull ≤ 0.025 …) apply only after `dft_labels.py` produces calibrated labels.
   **第一轮只用 proxy 排序**，物理单位阈值要等 DFT 标签出来才用。

---

## 9. Open-source stack / 开源栈

pymatgen · CHGNet / matgl / mace-torch · mp-api · atomate2 / custodian ·
CatBoost / XGBoost / LightGBM / scikit-learn · BoTorch + Ax / pymoo / Optuna.

## License

MIT (see `LICENSE`).
