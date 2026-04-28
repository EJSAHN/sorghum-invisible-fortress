# Sorghum SDI Analysis Pipeline

Reproducible analysis code for separating visible pigmentation from an NIR-defined structural spectral axis in sorghum seeds and testing associations with head smut outcomes.

The core outputs are:

- Structural Defense Index (SDI), defined as the standardized residual of NIR PC1 after regressing on VIS PC1.
- GWAS outputs for VIS PC1, NIR-derived SDI, and disease-response traits.
- Robustness checks for SDI–disease association, including covariate sensitivity and permutation tests.
- Effect-sign quadrant summaries using GWAS beta coefficients for SDI and disease-response traits.
- Optional candidate-gene annotation and keyword/GO enrichment summaries.

## Terminology

This repository uses formal effect-sign quadrant labels rather than informal labels.

| Code label | Effect-sign definition | Interpretation |
|---|---|---|
| `Qpp` | beta_SDI > 0 and beta_disease > 0 | positive SDI effect, positive disease-score effect |
| `Qpn` | beta_SDI > 0 and beta_disease < 0 | positive SDI effect, negative disease-score effect |
| `Qnp` | beta_SDI < 0 and beta_disease > 0 | negative SDI effect, positive disease-score effect |
| `Qnn` | beta_SDI < 0 and beta_disease < 0 | negative SDI effect, negative disease-score effect |

These labels are descriptive statistical categories. They should not be interpreted as causal mechanisms without additional validation.

## Repository contents

- `sdi_pipeline.py`  
  Main pipeline. Builds the phenotype table with spectral PCs and SDI, runs GWAS, annotates candidate SNPs, and writes supplementary-data tables.

- `sdi_robustness_checks.py`  
  Robustness module for SDI–disease models, within-group checks, permutation tests, threshold scans, and correlations with physical-proxy traits when available.

- `effect_sign_enrichment.py`  
  Effect-sign quadrant enrichment module. Summarizes quadrant counts and enrichment statistics across top-hit thresholds using chi-square and Fisher exact tests.

- `gwas_pc_sensitivity.py`  
  Genotype-PC sensitivity sweep. Evaluates GWAS stability and genomic inflation across covariate settings.

- `candidate_gene_keyword_enrichment.py`  
  Exploratory candidate-gene annotation and keyword/GO enrichment utilities. These outputs are intended as hypothesis-generating summaries.

- `genotype_to_spectral_projection.py`  
  Exploratory genotype-to-spectral projection using ridge regression. This is a sensitivity/projection exercise, not a time-series or causal model.

## Installation

Recommended: Python 3.9+ in a conda environment.

```bash
conda create -n sdi_env python=3.10 -y
conda activate sdi_env
pip install numpy pandas scipy scikit-learn matplotlib openpyxl
```

## Quick start

### 1. Run the main pipeline

Run from the directory containing the raw input files.

```bash
python sdi_pipeline.py --out outputs
```

Expected outputs include:

```text
outputs/Phenotypes_with_Spectral_PCs.csv
outputs/GWAS_Vis_PC1.csv
outputs/GWAS_SDI.csv
outputs/GWAS_<disease_trait>.csv
outputs/EffectSign_<disease_trait_short>.csv
outputs/Candidate_SNPs_Annotated.csv
outputs/Supplementary_Data_S1.xlsx
```

### 2. Run robustness checks

```bash
python sdi_robustness_checks.py \
  --out outputs \
  --sdi_col SDI \
  --traits headsmut_greenhouse_avg,headsmut_highest_score \
  --n_pcs 3 \
  --n_perm 2000 \
  --make_plots
```

### 3. Run effect-sign quadrant enrichment

```bash
python effect_sign_enrichment.py \
  --out outputs \
  --disease_pattern "GWAS_headsmut*.csv" \
  --top_pcts "0.5,1,2,5" \
  --use_fdr \
  --make_plots
```

### 4. Run genotype-PC sensitivity sweep

```bash
python gwas_pc_sensitivity.py \
  --out outputs \
  --trait SDI \
  --pc_list "0,1,2,3,4,6,8,10" \
  --include_race "0,1"
```

### 5. Run exploratory genotype-to-spectral projection

```bash
python genotype_to_spectral_projection.py \
  --input outputs/Supplementary_Data_S1.xlsx \
  --sheet S1_Phenotypes \
  --out outputs/addons_projection
```

## Optional LD clumping

If PLINK v1.9 is available, LD clumping can be used to summarize independent lead signals, for example with r2 = 0.2 within 250 kb windows. This step is optional and is not required for the Python pipeline.

## Notes and limitations

This repository provides reproducible computation and statistical validation for the SDI framework. Biological mechanism claims require independent experimental validation, such as tissue-level structural measurements, biochemical profiling, or functional genetics.

The genotype-to-spectral projection module is included only as an exploratory sensitivity analysis. It does not model change through time and should not be interpreted as a validated breeding trajectory.
