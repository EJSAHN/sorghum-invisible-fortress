# Output guide

The main pipeline writes outputs to the directory provided by `--out`.

## Core outputs

| File | Description |
|---|---|
| `Phenotypes_with_Spectral_PCs.csv` | Phenotype table after accession matching, with `Vis_PC1`, `NIR_PC1`, `SDI`, and genotype PCs. |
| `GWAS_Vis_PC1.csv` | GWAS table for the visible spectral component. |
| `GWAS_SDI.csv` | GWAS table for SDI. |
| `GWAS_<disease_trait>.csv` | GWAS table for each detected disease-response trait. |
| `EffectSign_<trait_short>.csv` | Effect-sign quadrant table for top disease-associated SNPs merged with SDI effects. |
| `Candidate_SNPs_Annotated.csv` | Top SNP/candidate-gene table, with annotation when available. |
| `Supplementary_Data_S1.xlsx` | Combined workbook for manuscript supplementary data. |
| `run_manifest.json` | Input paths, QC settings, detected disease traits, and output summary. |

## Add-on outputs

Additional scripts write robustness checks, effect-sign enrichment summaries, PC-sensitivity tables, candidate-gene keyword summaries, and genotype-to-spectral projection outputs.
