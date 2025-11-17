# scDC: Single-cell Differential Composition Analysis

## 📋 Overview
scDC performs differential abundance analysis with bootstrap resampling to quantify uncertainty in cell type proportions, using both GLM and GLMM frameworks.

## 🎯 Core Principle
- Bootstrap resampling for uncertainty estimation
- Generalized linear (mixed) models for testing
- Handles complex experimental designs
- Provides confidence intervals for abundance changes

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: GLM/GLMM
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: scDC, lme4

## ⚙️ Key Parameters
NULL

## 💡 Biological Applications
- Robust abundance estimation
- Complex experimental designs
- Clinical sample analysis with covariates

## 🔗 Official Resources
- **GitHub**: https://github.com/SydneyBioX/scDC
- **Publication**: [BMC Bioinformatics](https://link.springer.com/article/10.1186/s12859-019-3211-9)

## ⚠️ Important Note
**Current Implementation Status**: This package has known issues in the current version. Only the `scDC_noClustering` function is working reliably. Other functions may produce errors or unexpected results.

## ⚙️ Installation

```r
## Some CRAN packages required by scDC
install.packages(c("parallel", "DescTools", "lme4", "reshape2", "ggridges", 
"lme4", "mice", "broom.mixed"))

## Some BioConductor packages required by scDC
BiocManager::install(c("scran"))

## Some Github packages required by scDC
devtools::install_github("SydneyBioX/scClustBench")

## Installing scDC 
devtools::install_github("SydneyBioX/scDC")