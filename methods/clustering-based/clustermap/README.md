
# ClusterMap: Multi-sample cluster alignment for differential abundance

## 📋 Overview
ClusterMap is an R package designed to analyze and compare two or more single cell expression datasets. ClusterMap performs differential abundance analysis by aligning and matching cell clusters across multiple samples using tree-based clustering and pruning algorithms.

## 🎯 Core Principle
- Constructs hierarchical clustering trees for each sample
- Aligns clusters across samples using tree matching
- Identifies conserved and sample-specific cell populations
- Tests for abundance differences in matched clusters

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Tree Pruning Algorithm
- **Multi-group**: ✅ Supports multiple conditions


## 💡 Biological Applications
- Multi-sample cluster alignment studies
- Cross-dataset cell type comparison
- Conserved cell population identification
- Batch effect correction in clustering

## 🔗 Official Resources
- **GitHub**: https://github.com/xgaoo/ClusterMap
- **Publication**: [Bioinformatics](https://academic.oup.com/bioinformatics/article/35/17/3038/5289328?login=false)


## ⚙️ Installation

### PyPI Installation
```r
install_github('devtools')  
library('devtools')  
install_github("xgaoo/ClusterMap")
library('ClusterMap')  