# Clustering-Free Differential Abundance Methods

This directory contains implementations and examples for **clustering-free differential abundance (DA) analysis methods**.

## 📋 Overview

Clustering-free methods identify abundance changes without relying on pre-defined cell clusters. Instead, they operate directly on the continuous cellular manifold, detecting local neighborhoods or density changes that differ between conditions.

### 🔍 **Core Paradigm**:
- **No Discrete Clusters** → Bypasses hard clustering boundaries
- **Local Neighborhoods** → Analyzes small, overlapping cell communities  
- **Continuous Resolution** → Captures subtle state transitions and gradients
- **Data-Driven Discovery** → Identifies novel, unannotated cell states

## 🗂️ Available Methods

| Method | Core Model | Language | Multi-group | Key Feature |
|--------|---------------|-----------------|----------|
| [Milo](https://www.nature.com/articles/s41587-021-01033-z) | k-NN Graph + NB-GLM | R/Python | ✅ | Overlapping neighborhoods on k-NN graph for high-resolution local changes |
| [Cydar](https://www.nature.com/articles/nmeth.4295) | Hyperspheres + NB-GLM | R | ✅ | Fixed-radius hyperspheres for density-based abundance analysis |
| [DA-seq](https://www.pnas.org/doi/abs/10.1073/pnas.2100293118) | Multi-scale DA Scoring + Logistic Regression | R | ❌ | Multi-scale differential abundance scoring at single-cell level |
| [MELD](https://www.nature.com/articles/s41587-020-00803-5) | Kernel Density Estimation | Python | ✅ | Relative likelihood estimation for continuous probability surfaces |
| [CNA](https://www.nature.com/articles/s41587-021-01066-4) | Network Connectivity + Linear Model | Python/R | ✅ | Phenotypic network connectivity analysis for relationship changes |
| [scDC](https://link.springer.com/article/10.1186/s12859-019-3211-9) | Bootstrap + GLM/GLMM | R | ✅ | Bootstrap resampling for uncertainty quantification in proportions |
| [Dawnn](https://www.biorxiv.org/content/10.1101/2023.05.05.539427.abstract) | Deep Neural Networks | R | ❌ | Deep learning for end-to-end condition probability learning |
| [CellCnn](https://www.nature.com/articles/ncomms14825) | Convolutional Neural Networks | Python | ✅ | Filter-based CNN for rare cell population detection |
| [MrVI](https://www.nature.com/articles/s41592-025-02808-x) | Variational Autoencoder + Linear Model| Python | ✅ | Variational autoencoder with multi-head attention for batch effects |
| [Kompot](https://www.biorxiv.org/content/10.1101/2025.06.03.657769v3.full) | Gaussian Process Regression | Python | ✅ | Gaussian process regression for continuous manifold modeling |
| [Scissor](https://www.nature.com/articles/s41587-021-01091-3) | Regularized Sparse Regression | R | ✅ | Network-regularized regression linking single-cell to bulk phenotypes |

### ✅ **Advantages:**
- **High Resolution**: Detects changes within annotated cell types
- **Discovery Power**: Identifies novel, unannotated cell states
- **Continuous Transitions**: Captures gradient and intermediate states
- **Cluster-Free**: Avoids biases from clustering parameter choices
- **Subtle Effects**: Sensitive to small but consistent population shifts

### ⚠️ **Limitations:**
- **Parameter Sensitivity**: Results depend on neighborhood definitions
- **Computational Cost**: Generally more intensive than clustering-based
- **Interpretation Challenge**: Results may not map to traditional cell types
- **Complex Workflow**: Requires careful validation and interpretation

**Note**: Clustering-free methods represent the cutting edge of DA analysis, offering unprecedented resolution but requiring careful application and interpretation.