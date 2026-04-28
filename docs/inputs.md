# Input requirements

## Phenotype table

The phenotype table must include an `accession` column. The main pipeline intersects accessions across the phenotype table, hyperspectral table, and HapMap genotype file.

Commonly used columns:

- `race`: optional population/race covariate.
- `headsmut_greenhouse_avg`: primary disease-response trait when available.
- `headsmut_highest_score`: alternative disease-response trait for sensitivity checks.
- `headsmut_score`: fallback disease-response trait name.

## Hyperspectral table

The hyperspectral table should be accession-indexed. Reflectance columns are detected by the prefix `R_`. The numeric suffix is interpreted as wavelength in nm.

Examples:

```text
R_400, R_406, R_412, ..., R_1000
```

Visible PC1 is computed from 400 <= wavelength < 700 nm. NIR PC1 is computed from 700 <= wavelength <= 1000 nm.

## Genotype table

The genotype file is expected in HapMap format. The pipeline converts allele calls to dosage, filters SNPs by MAF and missingness, computes genotype PCs, and then runs OLS GWAS with covariates.

Default filters:

- MAF >= 0.05
- missing rate <= 0.10

## Optional annotation files

A GFF3 file and sorghum annotation/defline files can be provided for candidate-gene annotation. These are not required for SDI calculation or GWAS.
