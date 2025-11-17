# ELVAR: Attribute-aware community detection for differential abundance

## 📋 Overview
ELVAR is an R-package for differential abundance (DA) testing of cell-types in single-cell RNA-Seq data. It implements an Extended Louvain clustering Algorithm (EVA) that takes cell attribute information into acccount when inferring cellular communities from the cell-cell nearest neighbour graph. By taking cell attribute information into account, improved community detection is possible, which improves the power to detect DA-shifts of cell-types in relation to disease risk factors or disease itself.

## 🎯 Core Principle
- Attribute-aware Louvain clustering algorithm
- Negative binomial generalized linear models
- Integrates cell metadata into clustering
- Tests for differential abundance in refined communities

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Negative Binomial GLM
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: ELVAR, igraph, MASS

## ⚙️ Key Parameters
- `alpha`: the most important parameter as it controls the relative importance of purity (defined by how homogeneous the inferred communities are in relation to the cell attribute of interest) and modularity (the Louvain modularity compares the edge density of the communities to that of a null distribution obtained by randomising edges keeping the degrees fixed)
- `threshold`: the tolerance threshold for the optimization function

## 💡 Biological Applications
- Attribute-informed cell clustering
- Multi-modal single-cell data integration
- Cell community detection with metadata
- Complex experimental designs

## 🔗 Official Resources
- **GitHub**: https://github.com/aet21/ELVAR
- **Publication**: [nature communications ](https://www.nature.com/articles/s41467-023-39017-z)

## ⚙️ Installation

### GitHub Installation
```r
library(devtools)
devtools::install_github("aet21/ELVAR")