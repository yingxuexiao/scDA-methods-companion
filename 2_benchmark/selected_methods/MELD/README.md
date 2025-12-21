
# MELD: Quantifying the effect of experimental perturbations at single-cell resolution

## 📋 Overview
MELD uses kernel density estimation to compute relative likelihoods of cell states under different experimental conditions, enabling differential abundance analysis on continuous manifolds.The goal of MELD is to identify populations of cells that are most affected by an experimental perturbation. Rather than clustering the data first and calculating differential abundance of samples within clusters, MELD provides a density estimate for each scRNA-seq sample for every cell in each dataset. Comparing the ratio between the density of each sample provides a quantitative estimate the effect of a perturbation at the single-cell level. We can then identify the cells most or least affected by the perturbation.

## 🎯 Core Principle
- Estimates probability density functions for each condition
- Computes relative likelihoods using kernel density estimation
- Identifies regions with significant abundance differences
- Works on continuous cell state manifolds

## 🔧 Technical Implementation
- **Language**: Python
- **Core Model**: Kernel Density Estimation

## ⚙️ Key Parameters
- `knn` : Number of nearest neighbors (default: 5)
- `decay` : Alpha decay (default: 15)


## 💡 Biological Applications
- Continuous cell state analysis
- Developmental trajectory comparison
- Drug response characterization

## 🔗 Official Resources
- **GitHub**: https://github.com/KrishnaswamyLab/MELD
- **Publication**: [Nature Biotechnology](https://www.nature.com/articles/s41587-020-00803-5)


## ⚙️ Installation

### PyPI 
```bash
pip install meld