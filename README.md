# Sorghum SDI Pipeline (Hyperspectral 400–1000 nm + GWAS + Robustness)

Reproducible analysis pipeline for separating visible pigmentation from an NIR-defined structural axis in sorghum seeds and testing associations with head smut outcomes.  
Core outputs include the **Structural Defense Index (SDI)**, GWAS results, robustness checks (permutation; covariate sensitivity), effect-space quadrant enrichment (χ² / Fisher OR), and optional LD clumping summaries (PLINK).

> Manuscript nickname: “The Invisible Fortress” (code/analysis focus is SDI + validation).

---
## Important note on terminology & variable names

To preserve exact computational reproducibility (and to avoid introducing accidental changes during refactoring), some internal variable names and filenames retain early working labels. These legacy terms map to the manuscript’s finalized effect-sign quadrant terminology as follows:

| Legacy term (code/output) | Manuscript terminology | Effect-sign definition |
|---|---|---|
| `Trap` | Q++ | β_SDI > 0 and β_HS > 0 |
| `FreeLunch` / `Free_lunch` | Q+- | β_SDI > 0 and β_HS < 0 |
| `Sacrifice` | Q-+ | β_SDI < 0 and β_HS > 0 |
| `Fortress` | Q-- | β_SDI < 0 and β_HS < 0 |
| `Weak` / `LOWSDI` / `LowEffect` (if present) | low-effect group | variants not assigned to quadrants under the chosen filtering rule |

> **Note:** Scripts or filenames containing **`Maslow`** utilize a **legacy metaphorical label** to describe the **Effect-size Landscape Analysis**. The mathematical logic is unchanged regardless of this naming convention.


## What this repo does

1. **Spectral decomposition**
   - Performs PCA separately on VIS (400–700 nm) and NIR (700–1000 nm) reflectance domains.
   - Defines **SDI** as the standardized residual of NIR_PC1 after regressing on VIS_PC1.

2. **Association analyses**
   - GWAS for VIS_PC1 and SDI, with genotype-PC covariates and optional race covariates.
   - Sensitivity checks via genotype-PC sweeps and robustness reporting.

3. **Validation / robustness modules**
   - Permutation-based checks for SDI–disease association robustness.
   - Effect-space quadrant enrichment in the top-hit tail (χ² tests; Fisher’s exact tests / odds ratios).
   - Optional LD clumping workflow (PLINK v1.9) for independent lead-signal summaries (typically reported as a Supplementary Table).

---

## Repository contents

- `kimura_pipeline_public.py`  
  Main pipeline: prepares phenotype table with spectral PCs and SDI, and writes GWAS-ready outputs.

- `kimura_addon_robustness_public.py`  
  Robustness add-on: SDI–disease models, within-group checks, permutation tests, threshold scans, and plots.

- `kimura_addon_maslow_enrichment_public.py`  
  Effect-space quadrant enrichment: summarizes quadrant counts and enrichment statistics across top-hit thresholds (χ²; Fisher OR).  
  *Note: the manuscript uses effect-sign quadrants (Q++/Q+-/Q-+/Q--) rather than metaphor labels.*

- `kimura_addon_gwas_pc_sweep_public.py`  
  Genotype-PC sensitivity sweep: evaluates stability/inflation across covariate settings and reports summary metrics.

- `kimura_go_keyword_addon.py`  
  Exploratory enrichment utilities for candidate gene annotations (keyword-based / hypothesis-generating).

- `kimura_in_silico_ridge_public.py`  
  Exploratory genotype-to-phenotype projection (illustrative sensitivity analysis; **not** an evolutionary simulation).

---

## Installation

Recommended: Python 3.9+ (conda environment)

```bash
conda create -n gwas_env python=3.10 -y
conda activate gwas_env
pip install numpy pandas scipy scikit-learn matplotlib openpyxl
Quick start (typical workflow)
1) Run main pipeline
python kimura_pipeline_public.py --out outputs
Expected outputs (examples):

outputs/Phenotypes_with_Spectral_PCs.csv

outputs/GWAS_*.csv

2) Run robustness add-on
python kimura_addon_robustness_public.py \
  --out outputs \
  --sdi_col SDI \
  --traits headsmut_greenhouse_avg,headsmut_highest_score \
  --n_pcs 0 \
  --n_perm 2000 \
  --make_plots
3) Run quadrant enrichment add-on
python kimura_addon_maslow_enrichment_public.py \
  --out outputs \
  --disease_pattern "GWAS_headsmut*.csv" \
  --top_pcts "0.5,1,2,5" \
  --use_fdr \
  --make_plots
4) Run genotype-PC sweep add-on (optional)
python kimura_addon_gwas_pc_sweep_public.py \
  --out outputs \
  --trait SDI \
  --pc_list "0,1,2,3,4,6,8,10" \
  --include_race "0,1"
Terminology (effect-sign quadrants)
Quadrants are defined by the signs of GWAS beta coefficients:

Q++: β_SDI > 0 and β_HS > 0

Q+-: β_SDI > 0 and β_HS < 0

Q-+: β_SDI < 0 and β_HS > 0

Q--: β_SDI < 0 and β_HS < 0

If any legacy labels exist in intermediate outputs, they should be interpreted as sign-quadrant categories rather than biological mechanism claims.

LD clumping (optional; PLINK)
If you have PLINK v1.9 available, you can generate LD clumping summaries (e.g., r²=0.2 within 250 kb) to report independent lead signals.
This step is not required to run the Python pipeline but can be used to produce a Supplementary Table summarizing independent loci.

Notes / limitations
This repository provides reproducible computation and statistical validation; biological mechanism claims require targeted experimental validation.

In silico ridge projection is provided as an illustrative sensitivity analysis and should not be interpreted as evolution through time.
