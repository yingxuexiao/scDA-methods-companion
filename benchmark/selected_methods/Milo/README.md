# Milo: Differential abundance testing on KNN graphs

## 📋 Overview
Milo identifies differentially abundant populations of cells by defining overlapping neighborhoods on a k-nearest neighbor (k-NN) graph and testing for abundance changes using negative binomial generalized linear models.

## 🎯 Core Principle
- Constructs k-NN graph in low-dimensional space
- Samples representative, overlapping neighborhoods
- Counts cells per neighborhood per sample
- Tests differences with NB-GLM while controlling spatial FDR

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Negative Binomial GLM
- **Dependencies**: miloR, SingleCellExperiment

## ⚙️ Key Parameters
- `k`: Number of neighbors for graph construction (default: 30)
- `d`: PCA dimensions for embedding (default: 30)
- `prop`: Proportion of cells to sample as neighborhoods (default: 0.1)

## 💡 Biological Applications
- Fine-grained cell state transitions
- Local abundance changes within cell types
- Developmental trajectory analysis

## 🔗 Official Resources
- **GitHub**: https://github.com/MarioniLab/miloR
- **Publication**: [Nature Biotechnology](https://www.nature.com/articles/s41587-021-01033-z)


## 🔗 Installation
### Stable Release (Bioconductor)
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("miloR")