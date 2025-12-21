
# TreeCorTreat: Tree-based correlation screening for differential abundance

## 📋 Overview
TreeCorTreat performs differential abundance analysis by screening correlations in tree-structured cell type hierarchies to identify significant abundance changes across multiple resolutions.TreeCorTreat takes a gene expression matrix (raw count), cell-level metadata and sample-level metadata as input. It provides a whole pipeline to integrate data across samples, identify cell clusters and their hierarchical structure, evaluate the association between sample phenotype and cell type at different resolution levels in terms of both cell type proportion and gene expression, and summarize and visualize the results in a tree structured TreeCorTreat plot.

## 🎯 Core Principle
- Uses tree-based correlation screening algorithms
- Analyzes hierarchical cell type structures
- Screens for significant correlations with conditions
- Identifies abundance changes at multiple tree levels

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Tree-based Correlation Screening
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: TreeCorTreat, tree, stats


## 💡 Biological Applications
- Hierarchical cell type analysis
- Multi-resolution abundance screening
- Developmental lineage studies
- Cell differentiation hierarchy analysis

## 🔗 Official Resources
- **GitHub**: https://github.com/byzhang23/TreeCorTreat
- **Publication**: [biorxiv](https://www.biorxiv.org/content/10.1101/2021.10.27.466024v1.full)


## ⚠️ Important Note
**Current Implementation Status**: This package has known issues in the current version and may produce errors during execution. Users should verify functionality with their specific data.

## ⚙️ Installation

### GitHub Installation
```r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("byzhang23/TreeCorTreat")