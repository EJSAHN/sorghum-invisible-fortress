# Sorghum SDI Analysis Pipeline

Analysis code for separating visible pigmentation from an NIR-defined structural spectral axis in sorghum seeds and testing associations with head smut outcomes.

This repository is intended as a clean, reproducible code package for manuscript review and reuse. It contains analysis scripts, environment files, documentation, and lightweight sanity tests. Raw input data and generated outputs are not tracked by git by default.

## Main outputs

- Structural Defense Index (`SDI`), defined as the residual of `NIR_PC1` after regressing on `Vis_PC1`.
- GWAS tables for `Vis_PC1`, `SDI`, and available head smut traits.
- Robustness checks for SDI–disease association, including covariate sensitivity and permutation tests.
- Effect-sign enrichment summaries using GWAS beta coefficients for SDI and disease-response traits.
- Optional candidate-gene annotation and keyword/GO enrichment summaries.
- Exploratory genotype-to-spectral projection outputs.

## Repository layout

```text
.
├── README.md
├── LICENSE
├── requirements.txt
├── requirements-dev.txt
├── environment.yml
├── pyproject.toml
├── Makefile
├── sdi_pipeline.py
├── sdi_robustness_checks.py
├── effect_sign_enrichment.py
├── gwas_pc_sensitivity.py
├── candidate_gene_keyword_enrichment.py
├── genotype_to_spectral_projection.py
├── data/
│   └── README.md
├── docs/
│   ├── inputs.md
│   ├── outputs.md
│   └── statistical_notes.md
├── examples/
│   └── run_full_workflow.sh
└── tests/
    └── test_imports.py
```

## Installation

Python 3.10 or newer is recommended.

### Option A: conda

```bash
conda env create -f environment.yml
conda activate sorghum-sdi
```

### Option B: pip

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

For development checks:

```bash
python -m pip install -r requirements-dev.txt
```

The project can also be installed in editable mode:

```bash
python -m pip install -e .
```

This enables command-line entry points such as `sorghum-sdi-pipeline`, but direct script execution also works.

## Required inputs

The main pipeline auto-detects files in the current working directory. See `data/README.md` and `docs/inputs.md` for details.

Minimum required inputs:

1. phenotype CSV with an `accession` column;
2. hyperspectral accession-mean CSV with reflectance columns named like `R_650`;
3. HapMap genotype file.

Optional inputs include a GFF3 file and sorghum annotation/defline text files for candidate-gene annotation.

## Quick start

Run from the directory containing the raw input files:

```bash
python sdi_pipeline.py --out outputs
```

or, after editable installation:

```bash
sorghum-sdi-pipeline --out outputs
```

Expected core outputs include:

```text
outputs/Phenotypes_with_Spectral_PCs.csv
outputs/GWAS_Vis_PC1.csv
outputs/GWAS_SDI.csv
outputs/GWAS_<disease_trait>.csv
outputs/EffectSign_<disease_trait_short>.csv
outputs/Candidate_SNPs_Annotated.csv
outputs/Supplementary_Data_S1.xlsx
outputs/run_manifest.json
```

## Add-on analyses

### Robustness checks

```bash
python sdi_robustness_checks.py \
  --out outputs \
  --sdi_col SDI \
  --traits headsmut_greenhouse_avg,headsmut_highest_score \
  --n_pcs 3 \
  --n_perm 2000
```

### Effect-sign enrichment

```bash
python effect_sign_enrichment.py \
  --out outputs \
  --disease_pattern "GWAS_headsmut*.csv" \
  --top_pcts "0.5,1,2,5" \
  --use_fdr
```

### Genotype-PC sensitivity sweep

```bash
python gwas_pc_sensitivity.py \
  --out outputs \
  --trait SDI \
  --pc_list "0,1,2,3,4,6,8,10" \
  --include_race "0,1"
```

### Candidate-gene keyword enrichment

```bash
python candidate_gene_keyword_enrichment.py \
  --out outputs \
  --resources .
```

### Exploratory genotype-to-spectral projection

```bash
python genotype_to_spectral_projection.py \
  --input outputs/Supplementary_Data_S1.xlsx \
  --sheet S1_Phenotypes \
  --out outputs/projection
```

A full example workflow is provided in `examples/run_full_workflow.sh`.

## Effect-sign terminology

This repository uses formal quadrant labels rather than informal labels.

| Code label | Effect-sign definition | Interpretation |
|---|---|---|
| `Qpp` | beta_SDI > 0 and beta_disease > 0 | positive SDI effect, positive disease-score effect |
| `Qpn` | beta_SDI > 0 and beta_disease < 0 | positive SDI effect, negative disease-score effect |
| `Qnp` | beta_SDI < 0 and beta_disease > 0 | negative SDI effect, positive disease-score effect |
| `Qnn` | beta_SDI < 0 and beta_disease < 0 | negative SDI effect, negative disease-score effect |

These labels are descriptive statistical categories. They should not be interpreted as causal mechanisms without additional validation.

## Quality checks

```bash
make test
```

This compiles the scripts and runs import tests. It does not require raw input data.

Optional formatting/linting commands are included:

```bash
make lint
make format
```

## Optional LD clumping

If PLINK v1.9 is available, LD clumping can be used to summarize independent lead signals, for example with r² = 0.2 within 250 kb windows. This step is optional and is not required for the Python pipeline.

## Notes and limitations

SDI is an NIR-defined spectral index. It is not a direct measurement of seed coat density, porosity, thickness, or mechanical strength without independent validation.

Candidate-gene and enrichment outputs are intended as hypothesis-generating summaries. Gene-level mechanism claims require follow-up validation such as tissue-level structural measurements, biochemical profiling, or functional genetics.

The genotype-to-spectral projection module is included only as an exploratory sensitivity analysis. It does not model change through time and should not be interpreted as a validated breeding trajectory.
