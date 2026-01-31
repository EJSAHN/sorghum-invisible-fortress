# =============================================================================
# GO & Keyword Enrichment Analysis
#Analyzes the functional enrichment of candidate genes identified by GWAS.
#Matches genes to Gene Ontology (GO) terms and specific keywords (Defense/Cell Wall).
#Usage: Ensure 'GWAS_Results.csv', 'Annotation_Info.txt', and 'GFF.gff3' are in the directory.
# - Run inside a directory that contains the public pipeline outputs (default: ./outputs)
#- Optional resources (place in working dir or provide --resources):
#    - Sbicolor_454_v3.1.1.gene_exons.gff3
#    - Sbicolor_454_v3.1.1.P14.annotation_info.txt(.gz)
#    - Sbicolor_454_v3.1.1.P14.defline.txt(.gz)
# =============================================================================
# Outputs (in OUT_DIR):
#  CandidateGenes_Annotated.csv
#  GO_Enrichment_CandidateGenes.csv
#  Keyword_Enrichment_Defense_CellWall.csv
#  Supplementary_Data_S1_with_GO.xlsx  (if Supplementary_Data_S1.xlsx exists)
# =============================================================================

from __future__ import annotations
import os, re, glob, gzip, shutil, argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

GO_ID_RE = re.compile(r"GO:\d{7}")

DEFENSE_KWS = [
    "defense", "defence", "immune", "immunity", "pathogen", "disease", "resistance",
    "receptor-like kinase", "nbs-lrr", "nlr", "lrr", "pr protein",
    "chitinase", "glucanase", "wrky", "npr", "salicylic", "jasmon", "ethylene",
    "mapk", "kinase", "peroxidase", "oxidase", "ros"
]
CELLWALL_KWS = [
    "cell wall", "cellwall", "cellulose", "hemicellulose", "xylan", "pectin",
    "lignin", "phenylpropanoid", "laccase", "expansin",
    "callose", "cuticle", "cutin", "wax", "suberin"
]

def _mkdir(p: str) -> str:
    Path(p).mkdir(parents=True, exist_ok=True)
    return p

def _read_text_maybe_gz(path: str, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="ignore")
    return open(path, mode, encoding="utf-8", errors="ignore")

def _find_first(base: str, patterns: list[str]) -> str | None:
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(base, pat)))
        if hits:
            return hits[0]
    return None

def _norm_sobic(g: str) -> str | None:
    if g is None:
        return None
    m = re.search(r"(Sobic\.\d{3}G\d{6})", str(g))
    return m.group(1) if m else None

def _bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n+1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out

def _kw_match(text: str, kws: list[str]) -> bool:
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return False
    s = str(text).lower()
    return any(k.lower() in s for k in kws)

def load_defline_map(defline_path: str) -> dict[str, str]:
    m = {}
    with _read_text_maybe_gz(defline_path, "rt") as f:
        for ln in f:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            parts = ln.split("\t")
            if len(parts) < 2:
                continue
            gid = _norm_sobic(parts[0])
            if gid:
                m[gid] = parts[1].strip()
    return m

def load_p14_annotation_table(path: str) -> pd.DataFrame:
    hdr = pd.read_csv(path, sep="\t", nrows=0, compression="infer")
    cols = list(hdr.columns)

    gene_col = None
    for c in cols:
        if c.lower() in {"gene_id", "gene", "id", "name"}:
            gene_col = c
            break
    if gene_col is None:
        gene_col = cols[0]

    def interesting(c: str) -> bool:
        lc = c.lower()
        return (
            "uniprot" in lc or "interpro" in lc or "pfam" in lc or "go" in lc or
            "domain" in lc or "function" in lc or "description" in lc or "product" in lc
        )

    usecols = [gene_col] + [c for c in cols if c != gene_col and interesting(c)]
    usecols = usecols[:25]
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, compression="infer", usecols=usecols)
    df = df.rename(columns={gene_col: "gene_id"}).copy()
    df["gene_id"] = df["gene_id"].map(_norm_sobic)
    df = df.dropna(subset=["gene_id"]).drop_duplicates("gene_id")
    return df

def extract_go_map(df_annot: pd.DataFrame) -> dict[str, set[str]]:
    go_cols = [c for c in df_annot.columns if "go" in c.lower()]
    go_map = {}
    for _, r in df_annot.iterrows():
        gid = r["gene_id"]
        terms = set()
        for c in go_cols:
            val = r.get(c, None)
            if val is None:
                continue
            for m in GO_ID_RE.findall(str(val)):
                terms.add(m)
        go_map[gid] = terms
    return go_map

def go_enrichment(candidate_genes: set[str], go_map: dict[str, set[str]]) -> pd.DataFrame:
    universe = set(go_map.keys())
    cand = candidate_genes & universe
    if len(cand) == 0:
        return pd.DataFrame()

    term_genes = {}
    for g, tset in go_map.items():
        for t in tset:
            term_genes.setdefault(t, set()).add(g)

    U = len(universe); C = len(cand)
    rows = []
    for term, genes_with in term_genes.items():
        a = len(cand & genes_with)
        if a < 1:
            continue
        b = C - a
        c = len(genes_with - cand)
        d = (U - len(genes_with)) - b
        if (a + c) < 5:
            continue
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "GO": term,
            "cand_with": a, "cand_total": C,
            "bg_with": a + c, "bg_total": U,
            "oddsratio_approx": ((a+0.5)*(d+0.5))/((b+0.5)*(c+0.5)),
            "p": p
        })
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    out["fdr_bh"] = _bh_fdr(out["p"].to_numpy(float))
    out["neglog10_fdr"] = -np.log10(np.clip(out["fdr_bh"].to_numpy(float), 1e-300, 1.0))
    return out.sort_values(["fdr_bh", "p"]).reset_index(drop=True)

def keyword_enrichment(candidate_genes: set[str], universe_genes: set[str], cat_genes: dict[str, set[str]]) -> pd.DataFrame:
    cand = candidate_genes & universe_genes
    U = len(universe_genes); C = len(cand)
    rows = []
    for cat, gset in cat_genes.items():
        a = len(cand & gset)
        b = C - a
        c = len((gset & universe_genes) - cand)
        d = (U - len(gset & universe_genes)) - b
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "category": cat,
            "cand_with": a, "cand_total": C,
            "bg_with": a + c, "bg_total": U,
            "oddsratio_approx": ((a+0.5)*(d+0.5))/((b+0.5)*(c+0.5)),
            "p": p
        })
    out = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    out["fdr_bh"] = _bh_fdr(out["p"].to_numpy(float))
    out["neglog10_fdr"] = -np.log10(np.clip(out["fdr_bh"].to_numpy(float), 1e-300, 1.0))
    return out

def append_sheets_copy(src_xlsx: str, dst_xlsx: str,
                       genes_tbl: pd.DataFrame, go_tbl: pd.DataFrame, kw_tbl: pd.DataFrame):
    shutil.copy2(src_xlsx, dst_xlsx)
    with pd.ExcelWriter(dst_xlsx, engine="openpyxl", mode="a", if_sheet_exists="replace") as xl:
        genes_tbl.to_excel(xl, sheet_name="S1_CandGenes_Annot", index=False)
        (go_tbl if go_tbl is not None else pd.DataFrame()).to_excel(xl, sheet_name="S1_GO_Enrich", index=False)
        (kw_tbl if kw_tbl is not None else pd.DataFrame()).to_excel(xl, sheet_name="S1_KW_Enrich", index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="outputs", help="Directory containing public pipeline outputs")
    ap.add_argument("--resources", default=".", help="Directory containing GFF/annotation files")
    ap.add_argument("--cand_logp_min", type=float, default=4.5)
    args = ap.parse_args()

    out_dir = _mkdir(args.out)
    res_dir = args.resources
    cand_logp = float(args.cand_logp_min)

    # Find candidate gene list source:
    # Option A) CandidateGenes_Annotated.csv already exists (from local); use it.
    # Option B) Otherwise, try to derive genes from Candidate SNP annotation CSV if present.
    cand_csv = os.path.join(out_dir, "CandidateGenes_Annotated.csv")
    if os.path.exists(cand_csv):
        cand_tbl = pd.read_csv(cand_csv)
        cand_genes = set(cand_tbl["gene_id"].dropna().map(_norm_sobic))
        cand_genes = {g for g in cand_genes if g}
    else:
        # fallback: derive from disease GWAS using nearest-gene mapping only if GFF is available
        gff_path = _find_first(res_dir, ["Sbicolor_454_v3.1.1.gene_exons.gff3", "*.gff3"])
        if not gff_path:
            raise FileNotFoundError("Need CandidateGenes_Annotated.csv OR a GFF3 file in --resources.")
        # (Public mode: we stop here to avoid re-implementing full SNP->gene mapping without your run context.)
        raise FileNotFoundError("For public add-on, please provide CandidateGenes_Annotated.csv in outputs/.")

    # Load annotations
    annot_info_path = _find_first(res_dir, [
        "Sbicolor_454_v3.1.1.P14.annotation_info.txt",
        "Sbicolor_454_v3.1.1.P14.annotation_info.txt.gz",
        "*annotation_info*.txt", "*annotation_info*.txt.gz",
    ])
    defline_path = _find_first(res_dir, [
        "Sbicolor_454_v3.1.1.P14.defline.txt",
        "Sbicolor_454_v3.1.1.P14.defline.txt.gz",
        "*defline*.txt", "*defline*.txt.gz",
    ])
    if not annot_info_path and not defline_path:
        raise FileNotFoundError("No annotation_info/defline found in --resources.")

    df_annot = load_p14_annotation_table(annot_info_path) if annot_info_path else pd.DataFrame({"gene_id": []})
    dmap = load_defline_map(defline_path) if defline_path else {}

    # Build/refresh CandidateGenes_Annotated.csv
    base_tbl = pd.DataFrame({"gene_id": sorted(cand_genes)})
    base_tbl["defline"] = base_tbl["gene_id"].map(dmap)
    if not df_annot.empty:
        base_tbl = base_tbl.merge(df_annot, on="gene_id", how="left")

    # keyword flags
    text_cols = [c for c in base_tbl.columns if any(k in c.lower() for k in ["defline","uniprot","description","function","product","interpro","pfam","domain","go"])]
    base_tbl["defense_hit"] = [ _kw_match(" | ".join([str(r.get(c,"")) for c in text_cols]), DEFENSE_KWS) for _, r in base_tbl.iterrows() ]
    base_tbl["cellwall_hit"] = [ _kw_match(" | ".join([str(r.get(c,"")) for c in text_cols]), CELLWALL_KWS) for _, r in base_tbl.iterrows() ]

    # GO enrichment
    go_tbl = pd.DataFrame()
    kw_tbl = pd.DataFrame()
    if not df_annot.empty and annot_info_path:
        go_map = extract_go_map(df_annot)
        base_tbl["go_n"] = base_tbl["gene_id"].map(lambda g: len(go_map.get(g, set())))
        go_tbl = go_enrichment(set(base_tbl["gene_id"]), go_map)

        universe = set(go_map.keys())
        # keyword gene sets in universe
        uni_tbl = df_annot[["gene_id"]].copy()
        uni_tbl["defline"] = uni_tbl["gene_id"].map(dmap)
        desc_cols = [c for c in df_annot.columns if any(k in c.lower() for k in ["description","function","product","uniprot","interpro","pfam","domain","go"])]
        for c in desc_cols[:8]:
            uni_tbl[c] = df_annot.set_index("gene_id").loc[uni_tbl["gene_id"], c].values

        def_set = set()
        cw_set = set()
        for _, r in uni_tbl.iterrows():
            blob = " | ".join([str(r.get(c,"")) for c in uni_tbl.columns if c != "gene_id"])
            if _kw_match(blob, DEFENSE_KWS): def_set.add(r["gene_id"])
            if _kw_match(blob, CELLWALL_KWS): cw_set.add(r["gene_id"])

        kw_tbl = keyword_enrichment(set(base_tbl["gene_id"]), universe, {
            "defense/immune keywords": def_set,
            "cell wall keywords": cw_set,
        })

    # write outputs
    out_genes = os.path.join(out_dir, "CandidateGenes_Annotated.csv")
    out_go = os.path.join(out_dir, "GO_Enrichment_CandidateGenes.csv")
    out_kw = os.path.join(out_dir, "Keyword_Enrichment_Defense_CellWall.csv")

    base_tbl.to_csv(out_genes, index=False)
    (go_tbl if go_tbl is not None else pd.DataFrame()).to_csv(out_go, index=False)
    (kw_tbl if kw_tbl is not None else pd.DataFrame()).to_csv(out_kw, index=False)

    # merge into S1 if exists
    s1 = os.path.join(out_dir, "Supplementary_Data_S1.xlsx")
    if os.path.exists(s1):
        dst = os.path.join(out_dir, "Supplementary_Data_S1_with_GO.xlsx")
        append_sheets_copy(s1, dst, base_tbl, go_tbl, kw_tbl)
        print("Wrote:", dst)

    print("Done.")
    print("Outputs:")
    print(" -", out_genes)
    print(" -", out_go)
    print(" -", out_kw)

if __name__ == "__main__":
    main()
