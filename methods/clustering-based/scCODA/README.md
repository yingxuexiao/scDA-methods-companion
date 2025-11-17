
# scCODA: Bayesian model for compositional differential abundance

## 📋 Overview
scCODA is a Bayesian hierarchical model for compositional differential abundance analysis that accounts for the compositional nature of single-cell data and provides credible intervals for effect sizes. It also provides a framework for integration of cell-type annotated data directly from scanpy and other sources. 

## 🎯 Core Principle
- Uses Dirichlet-Multinomial model for compositional data
- Bayesian inference with Markov Chain Monte Carlo sampling
- Automatically adjusts for compositional effects
- Provides credible intervals and Bayes factors

## 🔧 Technical Implementation
- **Language**: Python
- **Core Model**: Bayesian Hierarchical Model
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: scCODA, arviz, pymc3

## ⚙️ Key Parameters
- `reference_cell_type`: Reference cell type for compositionality

## 💡 Biological Applications
- Compositional differential abundance analysis
- Small sample size studies
- Bayesian uncertainty quantification
- Cell type proportion changes

## 🔗 Official Resources
- **GitHub**: https://github.com/theislab/scCODA
- **Publication**: [Nature Communications](https://www.nature.com/articles/s41467-021-27150-6)


## ⚙️ Installation

### PyPI Installation
```bash
pip install sccoda