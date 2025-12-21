
# propeller: Proportion-based differential abundance testing with arcsin transformation

## 📋 Overview
propeller performs differential abundance analysis by testing for differences in cell type proportions using linear models on arcsin square-root transformed proportions, providing robust statistical testing for compositional data.The propeller, propeller.ttest and propeller.anova functions perform statistical tests for differences in cell type composition in single cell data. In order to test for differences in cell type proportions between multiple experimental conditions at least one of the groups must have some form of biological replication (i.e. at least two samples). For a two group scenario, the absolute minimum sample size is thus three. Since there are many technical aspects which can affect cell type proportion estimates, having biological replication is essential for a meaningful analysis.

## 🎯 Core Principle
- Transforms cell type proportions using arcsin square-root
- Uses linear models for hypothesis testing
- Handles compositional nature of proportion data
- Provides FDR-corrected p-values and confidence intervals

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Linear Model
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: speckle, limma, stats


## 💡 Biological Applications
- Simple to complex experimental designs
- Large-scale screening studies
- Clinical sample comparisons

## 🔗 Official Resources
- **GitHub**: https://github.com/phipsonlab/speckle
- **Publication**: [Bioinformatics](https://academic.oup.com/bioinformatics/article/38/20/4720/6675456)

## ⚙️ Installation

### Bioconductor Installation
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("speckle")