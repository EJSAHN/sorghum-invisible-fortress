# ============================================================
# Assumes raw inputs are in current working directory.
#1) Run `kimura_pipeline_public.py` in the directory containing the raw inputs to reproduce the main GWAS / SDI / Maslow results and Supplementary_Data_S1.xlsx.
#2) (Optional) Run `kimura_go_keyword_addon.py --out outputs --resources .` to generate GO_Enrichment_CandidateGenes.csv and keyword enrichment tables, and to append them as new sheets to Supplementary_Data_S1_with_GO.xlsx.
# Same core logic as local:
#  - Vis/NIR PC1 + SDI
#  - HapMap -> dosage int8 -> QC
#  - OLS GWAS with PCs (+race if present)
#  - Dual disease traits (highest + avg)
#  - Candidate SNP annotation via GFF + P14 annotation_info/defline if present
# ============================================================

from __future__ import annotations
import os, re, json, gzip, glob, tempfile, warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.simplefilter("ignore", category=pd.errors.PerformanceWarning)

MAF_MIN  = 0.05
MISS_MAX = 0.10

N_GENO_PCS_MAIN = 3
PC_SWEEP_LIST   = [0, 1, 2, 3, 4, 6, 8, 10]
TOP_PCT_FOR_MASLOW = 0.01

DISEASE_TRAITS_ORDER = ["headsmut_highest_score", "headsmut_greenhouse_avg", "headsmut_score"]
CANDIDATE_TOPN_PER_TRAIT = 200
CANDIDATE_GENE_WINDOW_BP = 100_000

def _mkdir(p: str | Path) -> str:
    Path(p).mkdir(parents=True, exist_ok=True)
    return str(p)

def _find_first(patterns: list[str]) -> str:
    base = os.getcwd()
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(base, pat)))
        if hits:
            return hits[0]
    raise FileNotFoundError(f"Could not find any of: {patterns} in {os.getcwd()}")

def _read_maybe_gz(path: str, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="ignore")
    return open(path, mode, encoding="utf-8", errors="ignore")

def _parse_chr_int(x) -> int | None:
    s = str(x).strip()
    s = re.sub(r"(?i)^chr", "", s).strip()
    m = re.search(r"(\d+)", s)
    return int(m.group(1)) if m else None

IUPAC_PAIR_TO_CODE = {
    frozenset(("A","G")): "R",
    frozenset(("C","T")): "Y",
    frozenset(("G","C")): "S",
    frozenset(("A","T")): "W",
    frozenset(("G","T")): "K",
    frozenset(("A","C")): "M",
}
def _split_alleles(a: str):
    s = str(a).strip().upper().replace(",", "/")
    parts = s.split("/")
    if len(parts) != 2:
        return None, None
    a1, a2 = parts[0].strip(), parts[1].strip()
    if len(a1) != 1 or len(a2) != 1:
        return None, None
    return a1, a2

def _lambda_gc_from_p(pvals: np.ndarray) -> float:
    p = np.clip(np.asarray(pvals, float), 1e-300, 1.0)
    chi2 = stats.chi2.isf(p, df=1)
    return float(np.median(chi2) / 0.4549364)

def _sheet(name: str) -> str:
    name = re.sub(r"[\[\]\*\?/\\:]", "_", name)
    return name[:31]

def trait_short(tr: str) -> str:
    if tr == "headsmut_highest_score":
        return "HS_max"
    if tr == "headsmut_greenhouse_avg":
        return "GH_avg"
    tr2 = tr.replace("headsmut_", "")
    tr2 = re.sub(r"[^A-Za-z0-9]+", "_", tr2)
    return tr2[:12]

def count_hapmap_snps(hmp_path: str) -> int:
    with open(hmp_path, "rb") as f:
        n = sum(1 for _ in f)
    return max(0, n - 1)

def hapmap_to_dosage(hmp_path: str, keep_samples: list[str], cache_dir: str):
    _mkdir(cache_dir)
    base = os.path.basename(hmp_path).strip()
    safe = re.sub(r"[^A-Za-z0-9]+", "_", base)[:40].strip("_")
    tag = f"{safe}_n{len(keep_samples)}"
    cache_npy  = os.path.join(cache_dir, f"geno_dosage_{tag}.npy")
    cache_meta = os.path.join(cache_dir, f"snp_meta_{tag}.csv")
    cache_samp = os.path.join(cache_dir, f"geno_samples_{tag}.txt")

    if os.path.exists(cache_npy) and os.path.exists(cache_meta) and os.path.exists(cache_samp):
        G = np.load(cache_npy)
        meta = pd.read_csv(cache_meta)
        with open(cache_samp, "r", encoding="utf-8", errors="ignore") as f:
            samples = [ln.strip() for ln in f if ln.strip()]
        return G, samples, meta, {"cache_npy": cache_npy, "cached": True}

    with open(hmp_path, "r", encoding="utf-8", errors="ignore") as f:
        header = f.readline().rstrip("\n").split("\t")
    rs_col, alleles_col, chr_col, pos_col = header[0], header[1], header[2], header[3]
    all_samples = [str(s).strip() for s in header[11:]]
    all_set = set(all_samples)
    keep_in_file = [s for s in keep_samples if s in all_set]
    n_snps = count_hapmap_snps(hmp_path)
    n = len(keep_in_file)

    usecols = [rs_col, alleles_col, chr_col, pos_col] + keep_in_file
    G = np.empty((n, n_snps), dtype=np.int8)
    meta_rows = []
    off = 0

    for chunk in pd.read_csv(hmp_path, sep="\t", usecols=usecols, dtype=str, low_memory=False, chunksize=2000):
        k = len(chunk)
        rs = chunk[rs_col].astype(str).to_numpy()
        al = chunk[alleles_col].astype(str).to_numpy()
        chr_raw = chunk[chr_col].astype(str).to_numpy()
        pos_raw = chunk[pos_col].to_numpy()
        chr_int = np.array([_parse_chr_int(x) for x in chr_raw], dtype=float)
        pos_int = pd.to_numeric(pd.Series(pos_raw), errors="coerce").to_numpy(dtype=float)
        geno = chunk[keep_in_file].to_numpy(dtype=str)
        Gc = np.full((n, k), -1, dtype=np.int8)

        for i in range(k):
            a1, a2 = _split_alleles(al[i])
            meta_rows.append({"rs": rs[i], "alleles": str(al[i]).replace(",", "/"),
                              "a1": a1, "a2": a2, "chr": chr_int[i], "pos": pos_int[i]})
            if a1 is None or a2 is None:
                continue
            het = IUPAC_PAIR_TO_CODE.get(frozenset((a1, a2)))
            s12 = f"{a1}/{a2}"; s21 = f"{a2}/{a1}"
            s12b = f"{a1}{a2}"; s21b = f"{a2}{a1}"
            row = np.array([str(x).strip().upper() for x in geno[i,:]], dtype=object)
            d = np.full(n, -1, dtype=np.int8)
            d[row == a1] = 0
            d[row == a2] = 2
            if het is not None:
                d[row == het] = 1
            d[(row==s12)|(row==s21)|(row==s12b)|(row==s21b)] = 1
            Gc[:, i] = d

        G[:, off:off+k] = Gc
        off += k

    if off != n_snps:
        G = G[:, :off]
        meta_rows = meta_rows[:off]

    meta = pd.DataFrame(meta_rows)

    # best-effort cache
    try:
        meta.to_csv(cache_meta, index=False)
        with open(cache_samp, "w", encoding="utf-8") as f:
            f.write("\n".join(keep_in_file) + "\n")
    except:
        pass

    def _atomic_save(path: str) -> bool:
        tmp = None
        try:
            fd, tmp = tempfile.mkstemp(prefix="__tmp_geno_", suffix=".npy", dir=cache_dir)
            os.close(fd)
            with open(tmp, "wb") as fh:
                np.save(fh, G, allow_pickle=False)
            os.replace(tmp, path)
            return True
        except:
            if tmp and os.path.exists(tmp):
                try: os.remove(tmp)
                except: pass
            return False

    ok = _atomic_save(cache_npy)
    if not ok:
        fallback = os.path.join(cache_dir, f"geno_dosage_n{len(keep_in_file)}.npy")
        _atomic_save(fallback)

    return G, keep_in_file, meta, {"cache_npy": cache_npy, "cached": False}

def geno_qc(G_int8: np.ndarray, meta: pd.DataFrame, maf_min=0.05, miss_max=0.1):
    X = G_int8.astype(np.float32)
    miss = np.mean(X < 0, axis=0)
    X[X < 0] = np.nan
    p = np.nanmean(X, axis=0) / 2.0
    maf = np.minimum(p, 1 - p)
    chr_ok = pd.to_numeric(meta["chr"], errors="coerce").notna().to_numpy()
    pos_ok = pd.to_numeric(meta["pos"], errors="coerce").notna().to_numpy()
    valid = np.isfinite(maf) & (maf >= maf_min) & (miss <= miss_max) & chr_ok & pos_ok
    meta_q = meta.loc[valid].copy().reset_index(drop=True)
    meta_q["maf"] = maf[valid]
    meta_q["miss"] = miss[valid]
    G_q = G_int8[:, valid]
    return G_q, meta_q

def compute_geno_pcs(G_int8: np.ndarray, n_pcs=3, seed=0, max_snps=30000):
    X = G_int8.astype(np.float32)
    X[X < 0] = np.nan
    mu = np.nanmean(X, axis=0)
    nan_idx = np.where(~np.isfinite(X))
    if nan_idx[0].size:
        X[nan_idx] = np.take(mu, nan_idx[1])
    m = X.shape[1]
    rng = np.random.default_rng(seed)
    use = rng.choice(m, size=min(max_snps, m), replace=False)
    Xs = X[:, use]
    sd = Xs.std(axis=0, ddof=0); sd[sd==0]=1.0
    Xs = (Xs - Xs.mean(axis=0)) / sd
    pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=seed)
    pcs = pca.fit_transform(Xs)
    return pcs

def spectral_pcs(spec: pd.DataFrame, accessions: list[str]):
    r_cols = [c for c in spec.columns if str(c).startswith("R_")]
    def wl(c):
        try: return float(str(c).split("_", 1)[1])
        except: return np.nan
    wls = np.array([wl(c) for c in r_cols], float)
    ok = np.isfinite(wls)
    r_cols = [c for c,o in zip(r_cols, ok) if o]
    wls = wls[ok]
    vis_cols = [c for c,w in zip(r_cols, wls) if 400 <= w < 700]
    nir_cols = [c for c,w in zip(r_cols, wls) if 700 <= w <= 1000]
    spec_sub = spec.loc[accessions, :]

    def pc1(cols):
        X = spec_sub[cols].to_numpy(dtype=float)
        Xs = StandardScaler().fit_transform(X)
        pca = PCA(n_components=1, random_state=0)
        scores = pca.fit_transform(Xs).ravel()
        m = np.nanmean(X, axis=1)
        r = np.corrcoef(scores, m)[0,1]
        if np.isfinite(r) and r < 0:
            scores = -scores
        return scores

    vis = pc1(vis_cols)
    nir = pc1(nir_cols)
    return vis, nir

def compute_sdi(vis_pc1: np.ndarray, nir_pc1: np.ndarray):
    m = np.isfinite(vis_pc1) & np.isfinite(nir_pc1)
    slope, intercept, r, p, _ = stats.linregress(vis_pc1[m], nir_pc1[m])
    sdi = nir_pc1 - (slope*vis_pc1 + intercept)
    return sdi, float(slope), float(intercept), float(r)

def gwas_ols(y: np.ndarray, G_int8: np.ndarray, meta: pd.DataFrame, cov: np.ndarray | None, block=20000):
    y = np.asarray(y, float)
    n0 = len(y)
    if cov is not None:
        cov = np.asarray(cov, float)
        if cov.shape[0] != n0:
            raise ValueError("cov rows do not match y length")

    msk = np.isfinite(y)
    if cov is not None:
        msk = msk & np.all(np.isfinite(cov), axis=1)
    y = y[msk]
    G = G_int8[msk, :]
    n = len(y)

    X = np.ones((n,1), float) if cov is None else np.column_stack([np.ones((n,1), float), cov[msk,:]])
    y = (y - np.mean(y)) / (np.std(y, ddof=0) + 1e-12)
    beta_y = np.linalg.lstsq(X, y, rcond=None)[0]
    y_res = y - X @ beta_y

    XtX_inv = np.linalg.inv(X.T @ X); Xt = X.T
    df = max(1, n - X.shape[1] - 1)
    ssy = float(np.sum(y_res*y_res))

    m = G.shape[1]
    betas = np.empty(m,float); pvals = np.empty(m,float)

    for s in range(0,m,block):
        e = min(m, s+block)
        Gb = G[:, s:e].astype(np.float32)
        Gb[Gb<0]=np.nan
        mu = np.nanmean(Gb, axis=0)
        nan_idx = np.where(~np.isfinite(Gb))
        if nan_idx[0].size:
            Gb[nan_idx] = np.take(mu, nan_idx[1])
        Gb = Gb - mu
        sd = Gb.std(axis=0, ddof=0); sd[sd==0]=1.0
        Gb = Gb / sd
        B = XtX_inv @ (Xt @ Gb)
        Gr = Gb - X @ B
        num = Gr.T @ y_res
        den = np.sum(Gr*Gr, axis=0); den[den==0]=np.nan
        b = num/den
        rss = np.maximum(ssy - b*num, 1e-30)
        se = np.sqrt((rss/df)/den)
        t = b/se
        p = 2*stats.t.sf(np.abs(t), df=df)
        betas[s:e]=b
        pvals[s:e]=p

    out = meta.copy()
    out["beta"]=betas
    out["p"]=np.clip(pvals,1e-300,1.0)
    out["-log10p"]=-np.log10(out["p"].values)
    lam = _lambda_gc_from_p(out["p"].values)
    lead_i = int(out["-log10p"].values.argmax())
    lead = out.iloc[lead_i][["chr","pos","rs","-log10p"]].to_dict()
    return out, lam, lead

def load_genes_from_gff(gff_path: str):
    genes = []
    with _read_maybe_gz(gff_path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, src, ftype, start, end, score, strand, phase, attr = parts
            if ftype != "gene":
                continue
            chr_i = _parse_chr_int(seqid)
            if chr_i is None:
                continue
            try:
                st = int(start); en = int(end)
            except:
                continue
            m = re.search(r"ID=([^;]+)", attr)
            gene_id = m.group(1) if m else None
            if not gene_id:
                continue
            genes.append({"chr": chr_i, "gene_id": gene_id, "start": st, "end": en})
    return pd.DataFrame(genes)

def load_defline(defline_path: str):
    d = {}
    with _read_maybe_gz(defline_path, "rt") as f:
        for ln in f:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            parts = ln.split("\t")
            if len(parts) < 2:
                parts = re.split(r"\s+", ln, maxsplit=1)
            if len(parts) >= 2:
                d[parts[0].strip()] = parts[1].strip()
    return d

def load_annotation_info(annot_path: str):
    df = pd.read_csv(annot_path, sep="\t", dtype=str, low_memory=False)
    keycol = df.columns[0]
    for c in df.columns:
        if c.lower() in {"gene_id", "gene", "id", "name"}:
            keycol = c
            break
    df = df.rename(columns={keycol: "gene_id"})
    return df

def annotate_candidates(candidate_df: pd.DataFrame, genes: pd.DataFrame,
                        annot_info_path: str|None, defline_path: str|None,
                        window_bp: int):
    if genes is None or genes.empty:
        return candidate_df.assign(nearest_gene=np.nan, gene_dist_bp=np.nan, gene_dist_kb=np.nan)

    genes_by_chr = {}
    for c, sub in genes.groupby("chr"):
        genes_by_chr[int(c)] = sub.sort_values("start").reset_index(drop=True)

    nearest_gene = []
    dist_bp = []
    for _, r in candidate_df.iterrows():
        chr_i = int(r["chr"]); pos = int(r["pos"])
        sub = genes_by_chr.get(chr_i, None)
        if sub is None or len(sub) == 0:
            nearest_gene.append(np.nan); dist_bp.append(np.nan); continue
        st = sub["start"].to_numpy(); en = sub["end"].to_numpy()
        d = np.where(pos < st, st-pos, np.where(pos > en, pos-en, 0))
        j = int(np.argmin(d))
        if d[j] > window_bp:
            nearest_gene.append(np.nan); dist_bp.append(np.nan)
        else:
            nearest_gene.append(sub.loc[j, "gene_id"]); dist_bp.append(float(d[j]))

    out = candidate_df.copy()
    out["nearest_gene"] = nearest_gene
    out["gene_dist_bp"] = dist_bp
    out["gene_dist_kb"] = out["gene_dist_bp"] / 1000.0

    if defline_path and os.path.exists(defline_path):
        dmap = load_defline(defline_path)
        out["defline"] = out["nearest_gene"].map(dmap)
    else:
        out["defline"] = np.nan

    if annot_info_path and os.path.exists(annot_info_path):
        ai = load_annotation_info(annot_info_path)
        keep_cols = ["gene_id"]
        for c in ai.columns:
            lc = c.lower()
            if any(k in lc for k in ["uniprot", "interpro", "go", "pfam", "domain", "function", "description"]):
                keep_cols.append(c)
        keep_cols = list(dict.fromkeys(keep_cols))
        ai2 = ai[keep_cols].copy().rename(columns={"gene_id":"nearest_gene"})
        out = out.merge(ai2, on="nearest_gene", how="left")

    return out

def main():
    OUT_DIR = _mkdir("outputs")
    CACHE_DIR = _mkdir(os.path.join(OUT_DIR, "cache"))

    pheno_path = _find_first(["senegal_master_phenotypes.csv", "*phenotype*.csv"])
    spec_path  = _find_first(["hyperspec_accession_mean.csv", "*hyperspec*mean*.csv", "*hyperspec*.csv"])
    hapmap_path = _find_first(["kimura_geno_aligned.hmp.txt", "*.hmp.txt", "*.hmp*.txt"])

    gff_path = None
    gff_hits = sorted(glob.glob("*.gff3"))
    if gff_hits:
        prefer = [p for p in gff_hits if "gene_exons" in os.path.basename(p)]
        gff_path = prefer[0] if prefer else gff_hits[0]

    annot_info = None
    defline = None
    for cand in ["Sbicolor_454_v3.1.1.P14.annotation_info.txt", "Sbicolor_454_v3.1.1.P14.annotation_info.txt.gz"]:
        if os.path.exists(cand):
            annot_info = cand; break
    for cand in ["Sbicolor_454_v3.1.1.P14.defline.txt", "Sbicolor_454_v3.1.1.P14.defline.txt.gz"]:
        if os.path.exists(cand):
            defline = cand; break

    ph = pd.read_csv(pheno_path)
    ph["accession"] = ph["accession"].astype(str).str.strip()
    spec = pd.read_csv(spec_path, index_col=0)
    spec.index = spec.index.astype(str).str.strip()

    with open(hapmap_path, "r", encoding="utf-8", errors="ignore") as f:
        header = f.readline().rstrip("\n").split("\t")
    hap_samples = [str(s).strip() for s in header[11:]]
    common = sorted(list(set(ph["accession"]) & set(spec.index) & set(hap_samples)))
    if len(common) < 50:
        raise RuntimeError(f"Too few common accessions: {len(common)}")

    df = ph.set_index("accession").loc[common].copy()
    vis_pc1, nir_pc1 = spectral_pcs(spec, common)
    sdi, slope, intercept, r = compute_sdi(vis_pc1, nir_pc1)
    df["Vis_PC1"] = vis_pc1
    df["NIR_PC1"] = nir_pc1
    df["SDI"] = sdi
    df_ix = df.copy()

    disease_traits = [t for t in DISEASE_TRAITS_ORDER if t in df_ix.columns]
    if not disease_traits:
        raise KeyError("No disease trait found among: " + ", ".join(DISEASE_TRAITS_ORDER))

    # genotype
    G_raw, geno_samples, meta_raw, cache_info = hapmap_to_dosage(hapmap_path, keep_samples=common, cache_dir=CACHE_DIR)
    G, meta = geno_qc(G_raw, meta_raw, maf_min=MAF_MIN, miss_max=MISS_MAX)

    pcs_g = compute_geno_pcs(G, n_pcs=max(N_GENO_PCS_MAIN, max(PC_SWEEP_LIST)), seed=0)
    for i in range(pcs_g.shape[1]):
        df_ix[f"PC{i+1}"] = pcs_g[:, i]

    if "race" in df_ix.columns:
        race = df_ix["race"].astype(str)
        dummies = pd.get_dummies(race, drop_first=True)
        for c in dummies.columns:
            df_ix[f"race_{c}"] = dummies[c].values

    def cov_matrix(npcs: int, use_race: bool=True):
        cols = [f"PC{i+1}" for i in range(npcs)]
        if use_race:
            cols += [c for c in df_ix.columns if c.startswith("race_")]
        return df_ix[cols].to_numpy(float) if cols else None

    C_main = cov_matrix(N_GENO_PCS_MAIN, use_race=True)

    # core GWAS
    gwas_vis, lam_vis, lead_vis = gwas_ols(df_ix["Vis_PC1"].to_numpy(float), G, meta, C_main)
    gwas_sdi, lam_sdi, lead_sdi = gwas_ols(df_ix["SDI"].to_numpy(float), G, meta, C_main)

    # disease loops
    gwas_dis_map = {}
    maslow_hits_map = {}
    maslow_sum_map = {}

    for tr in disease_traits:
        gwas_dis, lam_dis, lead_dis = gwas_ols(df_ix[tr].to_numpy(float), G, meta, C_main)
        gwas_dis_map[tr] = gwas_dis

        merged = pd.merge(
            gwas_sdi[["chr","pos","rs","beta","-log10p","maf"]].rename(columns={"beta":"beta_sdi","-log10p":"logp_sdi"}),
            gwas_dis[["chr","pos","rs","beta","-log10p"]].rename(columns={"beta":"beta_dis","-log10p":"logp_dis"}),
            on=["chr","pos"], how="inner"
        )
        merged["H"] = 2*merged["maf"]*(1-merged["maf"])
        top_n = max(1, int(len(merged) * TOP_PCT_FOR_MASLOW))
        hits = merged.nlargest(top_n, "logp_dis").copy()
        cut = float(np.nanmedian(np.abs(hits["beta_dis"].to_numpy(float))))

        def cls(row):
            bd = float(row["beta_dis"]) if pd.notna(row["beta_dis"]) else np.nan
            bs = float(row["beta_sdi"]) if pd.notna(row["beta_sdi"]) else np.nan
            if (not np.isfinite(bd)) or (abs(bd) < cut):
                return "Weak"
            if bs > 0 and bd > 0: return "Trap"
            if bs > 0 and bd < 0: return "Free_lunch"
            if bs < 0 and bd < 0: return "Fortress"
            if bs < 0 and bd > 0: return "Sacrifice"
            return "Mixed"
        hits["class"] = hits.apply(cls, axis=1)

        maslow_hits_map[tr] = hits
        maslow_sum_map[tr] = hits.groupby("class").size().reset_index(name="n")

    # candidate SNPs
    cand_rows = []
    def add_top(df_gwas: pd.DataFrame, trait_name: str):
        top = df_gwas.sort_values("-log10p", ascending=False).head(CANDIDATE_TOPN_PER_TRAIT)
        for _, rrr in top.iterrows():
            cand_rows.append({"set":"top","trait":trait_name,"chr":int(rrr["chr"]),"pos":int(rrr["pos"]),
                              "rs":rrr.get("rs",""),"beta":float(rrr["beta"]),"-log10p":float(rrr["-log10p"])})
    add_top(gwas_sdi, "SDI")
    for t in disease_traits:
        add_top(gwas_dis_map[t], t)

    candidates = pd.DataFrame(cand_rows).drop_duplicates(subset=["chr","pos","trait","set"]).reset_index(drop=True)
    genes = load_genes_from_gff(gff_path) if gff_path else pd.DataFrame()
    candidates_annot = annotate_candidates(candidates, genes, annot_info_path=annot_info, defline_path=defline,
                                          window_bp=CANDIDATE_GENE_WINDOW_BP)

    # write outputs
    df_ix.reset_index().to_csv(os.path.join(OUT_DIR, "Phenotypes_with_Spectral_PCs.csv"), index=False)
    gwas_vis.to_csv(os.path.join(OUT_DIR, "GWAS_Vis_PC1.csv"), index=False)
    gwas_sdi.to_csv(os.path.join(OUT_DIR, "GWAS_SDI.csv"), index=False)
    for tr in disease_traits:
        gwas_dis_map[tr].to_csv(os.path.join(OUT_DIR, f"GWAS_{tr}.csv"), index=False)
        maslow_hits_map[tr].to_csv(os.path.join(OUT_DIR, f"Maslow_{trait_short(tr)}.csv"), index=False)
    candidates_annot.to_csv(os.path.join(OUT_DIR, "Candidate_SNPs_Annotated.csv"), index=False)

    # S1
    s1_path = os.path.join(OUT_DIR, "Supplementary_Data_S1.xlsx")
    with pd.ExcelWriter(s1_path, engine="openpyxl") as xl:
        df_ix.reset_index().to_excel(xl, sheet_name="S1_Phenotypes", index=False)
        gwas_vis.to_excel(xl, sheet_name="S1_GWAS_Vis", index=False)
        gwas_sdi.to_excel(xl, sheet_name="S1_GWAS_SDI", index=False)

        for tr in disease_traits:
            sh = trait_short(tr)
            gwas_dis_map[tr].to_excel(xl, sheet_name=_sheet(f"S1_GWAS_{sh}"), index=False)
            maslow_hits_map[tr].to_excel(xl, sheet_name=_sheet(f"S1_Maslow_{sh}"), index=False)

        candidates_annot.to_excel(xl, sheet_name="S1_Candidates", index=False)

    manifest = {
        "cwd": os.getcwd(),
        "inputs": {"pheno": pheno_path, "spec": spec_path, "hapmap": hapmap_path, "gff": gff_path,
                   "annot_info": annot_info, "defline": defline},
        "n_common": len(common),
        "qc": {"maf_min": MAF_MIN, "miss_max": MISS_MAX, "n_snps_qc": int(G.shape[1])},
        "disease_traits": disease_traits,
        "outputs_dir": OUT_DIR,
        "S1": s1_path,
    }
    with open(os.path.join(OUT_DIR, "run_manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print("[DONE] public pipeline complete. Outputs ->", OUT_DIR)

if __name__ == "__main__":
    main()
