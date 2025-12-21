
# DA-seq: Multi-scale differential abundance analysis

## 📋 Overview
DA-seq computes multi-scale differential abundance scores at single-cell resolution by analyzing local neighborhood compositions across multiple scales.DA-seq first computes a score vector for each cell to represent the DA behavior in the neighborhood to select cells in the most DA neighborhoods; then groups these cells into distinct DA cell subpopulations.

## 🎯 Core Principle
- Calculates DA scores at multiple k-NN resolutions
- Identifies cells with consistent abundance changes
- Uses logistic regression for statistical testing
- Provides cell-level DA probability scores

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Logistic Regression
- **Multi-group**: ❌ Binary comparisons only
- **Dependencies**: DAseq, RANN, glmnet, caret, Seurat, e1071, reticulate, ggplot2, cowplot, scales, ggrepel

## ⚙️ Key Parameters
- `k.vector`: values of k to use for the calculation of score vector with kNN
- `cell.labels`: sample label for every cell in the data

## 💡 Biological Applications
- Multi-resolution abundance analysis
- Cell state transition identification
- Rare cell population detection

## 🔗 Official Resources
- **GitHub**: https://github.com/KlugerLab/Daseq
- **Publication**: [PNAS](https://www.pnas.org/doi/abs/10.1073/pnas.2100293118)

## ⚙️ Installation

```r
devtools::install_github("KlugerLab/DAseq")