
# CNA: Covarying neighborhood analysis is a method for finding structure in- and conducting association analysis with multi-sample single-cell datasets.

## 📋 Overview
cna does not require a pre-specified transcriptional structure such as a clustering of the cells in the dataset. It aims instead to flexibly identify differences of all kinds between samples. cna is fast, does not require parameter tuning, produces measures of statistical significance for its association analyses, and allows for covariate correction.

## 🎯 Core Principle
- Constructs phenotypic networks from single-cell data
- Analyzes network connectivity patterns
- Identifies regions with reorganized cellular relationships
- Uses linear models for statistical testing

## 🔧 Technical Implementation
- **Language**: Python/R
- **Core Model**: Linear Model
- **Dependencies**: cna-tools, networkx

## ⚙️ Key Parameters
- `k`: the number of NAM PCs used for the association test (automatically selected in a data-dependent way)
- `ncorrs`: the vector of neighborhood coefficients

## 💡 Biological Applications
- Network-level cellular reorganization
- Cell-cell interaction changes
- Systems-level abundance analysis

## 🔗 Official Resources
- **GitHub**: https://github.com/immunogenomics/cna
- **Publication**: [Nature Biotechnology](https://www.nature.com/articles/s41587-021-01066-4)

## ⚙️ Installation

### PyPI 
```bash
pip install cna
```
### R
```r
library(devtools)
install_github('korsunskylab/rcna')