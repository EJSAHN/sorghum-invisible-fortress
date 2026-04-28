# Statistical notes

## SDI definition

SDI is defined as the standardized residual of NIR PC1 after regression on visible PC1. This makes SDI orthogonal to the leading visible pigmentation component by construction.

SDI should be interpreted as an NIR-defined structural spectral index. It is not a direct measurement of seed coat density, porosity, thickness, or mechanical strength without independent validation.

## Effect-sign quadrants

The effect-sign analysis classifies SNPs by the signs of their SDI and disease-score beta coefficients:

- `Qpp`: beta_SDI > 0 and beta_disease > 0
- `Qpn`: beta_SDI > 0 and beta_disease < 0
- `Qnp`: beta_SDI < 0 and beta_disease > 0
- `Qnn`: beta_SDI < 0 and beta_disease < 0

These are descriptive categories for testing whether top disease-associated variants have a non-random distribution of effect directions. They are not causal labels.

## Genotype-to-spectral projection

The projection module fits a ridge regression from genotype PCs to spectral PCs and visualizes predicted shifts under perturbations of a genotype PC. This is a sensitivity/projection analysis, not a model of temporal change or a validated breeding simulation.
