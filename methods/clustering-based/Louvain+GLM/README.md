
# Louvain+GLM: Graph-based clustering with generalized linear models

## 📋 Overview
Louvain+GLM is a two-step differential abundance analysis pipeline that first identifies cell communities using Louvain graph clustering, then tests for abundance differences using negative binomial generalized linear models.

## ⚠️ Important Note
**Custom Pipeline**: This is not a standalone package but a commonly used analysis pipeline that combines Louvain clustering from Seurat/Scanpy with generalized linear models for differential abundance testing.

## 🎯 Core Principle
- Uses Louvain algorithm for graph-based community detection
- Applies negative binomial GLM to cluster counts
- Leverages cell-cell similarity graphs for clustering
- Tests for condition effects on cluster abundances

## 🔧 Technical Implementation
- **Language**: R/Python
- **Core Model**: Louvain Clustering + Negative Binomial GLM
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: Seurat/Scanpy, edgeR/DESeq2, igraph

## ⚙️ Key Parameters
- `resolution`:  Louvain clustering resolution
- `k`: Number of neighbors for graph construction


## 💡 Biological Applications
- Studies requiring community detection
- Single-cell data with clear graph structure


## 🔗 Official Resources
- **GitHub**: https://github.com/sunlab/Scissor
- **Publication**: [Nature Biotechnology](https://www.nature.com/articles/s41587-021-01187-w)

## ⚠️ Important Considerations
- Results depend on graph construction parameters
- Clustering resolution affects DA results
- Requires careful parameter tuning
- Not a unified method but a custom pipeline
- Interpretation should consider graph topology

## ⚙️ Installation

### R Implementation
```r
# Install required packages
install.packages(c("Seurat", "igraph", "MASS"))
BiocManager::install("edgeR")