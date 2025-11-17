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

| Method | Core Approach | Key Innovation | Best For |
|--------|---------------|-----------------|----------|
| [Milo](./Milo/) | k-NN Graph + NB-GLM | Overlapping neighborhoods on graph | High-resolution local changes |
| [Cydar](./Cydar/) | Hyperspheres + NB-GLM | Fixed-radius density neighborhoods | Density-based abundance shifts |
| [DA-seq](./DA-seq/) | Multi-scale DA Scoring | Cell-level differential abundance scores | Multi-resolution analysis |
| [MELD](./MELD/) | Kernel Density Estimation | Relative likelihood estimation | Continuous probability surfaces |
| [CNA](./CNA/) | Network Connectivity | Phenotypic network reorganization | Cell-cell relationship changes |
| [scDC](./scDC/) | Bootstrap + GLM/GLMM | Uncertainty quantification with resampling | Robust confidence intervals |
| [Dawnn](./Dawnn/) | Deep Neural Networks | End-to-end condition probability learning | Complex non-linear patterns |
| [CellCnn](./CellCnn/) | Convolutional Neural Networks | Rare cell population detection | Phenotype-associated rare subsets |
| [MrVI](./MrVI/) | Variational Autoencoder | Nonlinear latent space modeling | Complex batch integration |
| [Kompot](./Kompot/) | Gaussian Process Regression | Continuous manifold regression | Smooth abundance fields |

## 🎯 When to Use Clustering-Free Methods?

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