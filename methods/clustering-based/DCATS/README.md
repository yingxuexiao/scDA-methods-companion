
# DCATS: Differential Composition Analysis for single-cell data with uncertainty

## 📋 Overview
DCATS performs differential abundance analysis while accounting for classification uncertainty in cell type assignments using beta-binomial regression models.This R package contains methods to detect the differential composition abundances between multiple conditions in singel-cell experiments.

## 🎯 Core Principle
- Models cell type proportions with beta-binomial distribution
- Accounts for classification uncertainty in cell typing
- Uses similarity matrices to correct misclassification bias
- Handles over-dispersion in count data

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Beta-binomial Regression
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: DCATS, lme4, Matrix


## 💡 Biological Applications
- Cell type abundance analysis with classification uncertainty
- Single-cell data with ambiguous cell type assignments
- Studies requiring robust confidence intervals
- Multi-sample comparison with technical variation

## 🔗 Official Resources
- **GitHub**: https://github.com/holab-hku/DCATS
- **Publication**: [Genome Biology](https://link.springer.com/article/10.1186/s13059-023-02980-3)


## ⚙️ Installation

### GitHub Installation
```r
## install dependencies
install.packages(c("MCMCpack", "matrixStats", "robustbase", "aod", "e1071"))
## dependencies for vignette
install.packages(c("SeuratObject", "Seurat", "robustbase", "aod", "e1071"))
devtools::install_github('satijalab/seurat-data')

# install.packages("devtools")
devtools::install_github("holab-hku/DCATS", build_vignettes = TRUE)

#You can also install DCATS without building the vignette:

devtools::install_github("holab-hku/DCATS")

#From Biocounductor (required R >= 4.3.0)
if (!requireNamespace("BiocManager"))
# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("DCATS")