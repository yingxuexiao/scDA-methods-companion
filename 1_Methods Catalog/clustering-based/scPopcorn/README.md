
# scPopCorn: Single-cell Population Correspondence analysis

## 📋 Overview
A python tool to do comparative analysis of mulitple single cell RNA-seq datasets. scPopCorn performs analysis by matching corresponding cell populations across different samples or conditions using Google's personalized PageRank algorithm.

## 🎯 Core Principle
- Uses personalized PageRank for cell-to-cell matching
- Identifies corresponding populations across samples
- Computes population correspondence scores
- Tests for abundance differences in matched populations

## 🔧 Technical Implementation
- **Language**: Python
- **Core Model**: Personalized PageRank
- **Dependencies**: scpopcorn, networkx, numpy


## 💡 Biological Applications
- Clinical outcome association
- Biomarker discovery
- Treatment response prediction

## 🔗 Official Resources
- **GitHub**: https://github.com/ncbi/scPopCorn
- **Publication**: [Cell Systems](https://www.sciencedirect.com/science/article/pii/S2405471219301887)

## ⚠️ Important Note
**Current Implementation Status**: This package has known issues in the current version and may produce errors during execution. Users should verify functionality with their specific data.

## ⚙️ Installation

### PyPI Installation
```bash
pip install scpopcorn