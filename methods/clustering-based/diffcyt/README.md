
# diffcyt: Differential discovery in high-dimensional cytometry data

## 📋 Overview
The diffcyt package implements statistical methods for differential discovery analyses in high-dimensional cytometry data (including flow cytometry, mass cytometry or CyTOF, and oligonucleotide-tagged cytometry), based on a combination of high-resolution clustering and empirical Bayes moderated tests adapted from transcriptomics.

## 🎯 Core Principle
- Performs high-resolution clustering of cytometry data
- Tests for differential abundance using generalized linear models
- Handles complex experimental designs with multiple factors
- Provides flexible workflow for cytometry data analysis

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Generalized Linear Models
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: diffcyt, flowCore, SummarizedExperiment


## 💡 Biological Applications
- High-dimensional cytometry data analysis
- Immune cell subset discovery
- Clinical biomarker identification


## 🔗 Official Resources
- **GitHub**: https://github.com/lmweber/diffcyt
- **Publication**: [ communications biology](https://www.nature.com/articles/s42003-019-0415-5)


## ⚙️ Installation

### Bioconductor Installation
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("diffcyt")