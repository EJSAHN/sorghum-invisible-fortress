#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-outputs}"

python sdi_pipeline.py --out "${OUT_DIR}"

python sdi_robustness_checks.py \
  --out "${OUT_DIR}" \
  --sdi_col SDI \
  --traits headsmut_greenhouse_avg,headsmut_highest_score \
  --n_pcs 3 \
  --n_perm 2000

python effect_sign_enrichment.py \
  --out "${OUT_DIR}" \
  --disease_pattern "GWAS_headsmut*.csv" \
  --top_pcts "0.5,1,2,5" \
  --use_fdr

python gwas_pc_sensitivity.py \
  --out "${OUT_DIR}" \
  --trait SDI \
  --pc_list "0,1,2,3,4,6,8,10" \
  --include_race "0,1"

python candidate_gene_keyword_enrichment.py \
  --out "${OUT_DIR}" \
  --resources .

python genotype_to_spectral_projection.py \
  --input "${OUT_DIR}/Supplementary_Data_S1.xlsx" \
  --sheet S1_Phenotypes \
  --out "${OUT_DIR}/projection"
