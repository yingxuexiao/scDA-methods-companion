# Cydar: Using mass cytometry for differential abundance analysis

## 📋 Overview
Cydar performs differential abundance analysis by counting cells in fixed-radius hyperspheres in high-dimensional marker space and testing for abundance changes using negative binomial models.

## 🎯 Core Principle
- Defines M-dimensional hyperspheres with fixed radius
- Counts cells within overlapping hyperspheres
- Tests differential abundance with NB-GLM
- Controls spatial FDR across multiple testing

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Negative Binomial GLM
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: cydar, limma

## ⚙️ Key Parameters
- `tol`: A numeric scalar to be used as the scaling factor for the hypersphere radius.
- `filter`: An integer scalar specifying the minimum count sum required to report a hypersphere.
- `downsample`: An integer scalar specifying the frequency with which cells are sampled to form
hyperspheres.

## 💡 Biological Applications
- Mass cytometry data analysis
- High-dimensional abundance changes
- Immune cell population dynamics

## 🔗 Official Resources
- **GitHub**: https://github.com/MarioniLab/cydar
- **Publication**: [Nature Methods](https://www.nature.com/articles/nmeth.4295)

## ⚙️ Installation

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cydar")
