# Data directory

This repository is configured as a code-only package. Raw input files are not tracked by git by default.
Place the required input files either in the repository root or in a working directory from which the scripts are run.

The main pipeline auto-detects the following files:

| Input type | Preferred filename or pattern | Notes |
|---|---|---|
| Phenotype table | `senegal_master_phenotypes.csv` or `*phenotype*.csv` | Must contain `accession`. Optional columns include `race`, `headsmut_greenhouse_avg`, `headsmut_highest_score`, and `headsmut_score`. |
| Hyperspectral table | `hyperspec_accession_mean.csv`, `*hyperspec*mean*.csv`, or `*hyperspec*.csv` | Row index should be accession IDs. Reflectance columns should use the `R_<wavelength>` convention, for example `R_650`. |
| Genotype table | `*_geno_aligned.hmp.txt`, `*.hmp.txt`, or `*.hmp*.txt` | HapMap format. Sample IDs are expected from column 12 onward. |
| Gene annotation | `*.gff3` | Optional. Used for candidate-gene annotation. |
| Annotation text | `Sbicolor_454_v3.1.1.P14.annotation_info.txt(.gz)` and/or `Sbicolor_454_v3.1.1.P14.defline.txt(.gz)` | Optional. Used to add functional annotation to candidate-gene outputs. |
