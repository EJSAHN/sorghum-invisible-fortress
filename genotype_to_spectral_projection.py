#!/usr/bin/env python3
"""
genotype_to_spectral_projection.py

Exploratory genotype-to-spectral projection using ridge regression.

This script maps genotype PCs to selected spectral traits and applies a user-defined
shift along one genotype PC. The output is intended as a sensitivity/projection exercise,
not as a time-series model or a validated breeding prediction.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


def parse_cols(value: str) -> List[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="outputs/Supplementary_Data_S1.xlsx", help="Input Excel workbook or CSV file")
    ap.add_argument("--sheet", default="S1_Phenotypes", help="Excel sheet name if --input is an .xlsx file")
    ap.add_argument("--targets", default="Vis_PC1,NIR_PC1,SDI", help="Comma-separated spectral target columns")
    ap.add_argument("--geno_pcs", default="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10", help="Comma-separated genotype PC columns")
    ap.add_argument("--pc_index", type=int, default=1, help="One-based genotype PC index to perturb")
    ap.add_argument("--shift", type=float, default=2.0, help="Perturbation size in standardized genotype-PC units")
    ap.add_argument("--alpha", type=float, default=1.0, help="Ridge regression penalty")
    ap.add_argument("--out", default="outputs/addons_projection", help="Output directory")
    args = ap.parse_args()

    input_path = Path(args.input)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if input_path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(input_path, sheet_name=args.sheet)
    else:
        df = pd.read_csv(input_path)

    target_cols = parse_cols(args.targets)
    geno_pc_cols = parse_cols(args.geno_pcs)

    missing_targets = [c for c in target_cols if c not in df.columns]
    missing_pcs = [c for c in geno_pc_cols if c not in df.columns]
    if missing_targets:
        raise ValueError(f"Missing target columns: {missing_targets}")
    if missing_pcs:
        raise ValueError(f"Missing genotype PC columns: {missing_pcs}")
    if not 1 <= args.pc_index <= len(geno_pc_cols):
        raise ValueError(f"--pc_index must be between 1 and {len(geno_pc_cols)}")

    S_data = df[target_cols].astype(float).fillna(0).to_numpy()
    G_data = df[geno_pc_cols].astype(float).fillna(0).to_numpy()

    scaler_G = StandardScaler()
    scaler_S = StandardScaler()
    G_scaled = scaler_G.fit_transform(G_data)
    S_scaled = scaler_S.fit_transform(S_data)

    model = Ridge(alpha=args.alpha)
    model.fit(G_scaled, S_scaled)
    r2 = float(model.score(G_scaled, S_scaled))

    G_shifted = G_scaled.copy()
    G_shifted[:, args.pc_index - 1] += float(args.shift)
    S_pred_scaled = model.predict(G_shifted)
    S_pred = scaler_S.inverse_transform(S_pred_scaled)

    pred_df = pd.DataFrame(S_pred, columns=[f"projected_{c}" for c in target_cols])
    base_cols = [c for c in ["accession", "Accession", "sample", "Sample"] if c in df.columns]
    if base_cols:
        pred_df.insert(0, base_cols[0], df[base_cols[0]].values)
    pred_df.to_csv(out_dir / "genotype_pc_projection.csv", index=False)

    summary = pd.DataFrame([{
        "input": str(input_path),
        "sheet": args.sheet if input_path.suffix.lower() in {".xlsx", ".xls"} else "",
        "n_samples": int(G_data.shape[0]),
        "n_genotype_pcs": int(G_data.shape[1]),
        "target_columns": ",".join(target_cols),
        "perturbed_pc": geno_pc_cols[args.pc_index - 1],
        "shift_sd_units": float(args.shift),
        "ridge_alpha": float(args.alpha),
        "in_sample_r2": r2,
    }])
    summary.to_csv(out_dir / "genotype_pc_projection_summary.csv", index=False)

    print(f"[OK] Wrote projection outputs to: {out_dir}")


if __name__ == "__main__":
    main()
