# Clustering-Based Differential Abundance Methods

This directory contains implementations and examples for **clustering-based differential abundance (DA) analysis methods**.

## 📋 Overview

Clustering-based methods operate by first grouping cells into discrete clusters (cell types/states) using unsupervised learning algorithms, then testing for abundance differences of these pre-defined clusters across experimental conditions.

### 🔍 **Core Workflow**:
1. **Cell Clustering** → Identify cell types (Louvain, k-means, etc.)
2. **Count Aggregation** → Summarize cells per cluster per sample  
3. **Statistical Testing** → Test for abundance differences (GLM, Bayesian, etc.)

## 🗂️ Available Methods

| Method | Core Model | Language | Multi-group | Key Feature |
|--------|------------|----------|-------------|-------------|
| [scCODA](./scCODA/) | Bayesian | Python | ✅ | Compositional analysis with credible intervals |
| [propeller](./propeller/) | Linear Model | R | ✅ | Arcsin transformation for proportions |
| [diffcyt](./diffcyt/) | GLM | R | ✅ | Designed for cytometry data, generalizable to scRNA-seq |
| [MASC](./MASC/) | GLMM | R | ✅ | Mixed models for repeated measures |
| [DCATS](./DCATS/) | Beta-binomial | R | ✅ | Handles classification uncertainty |
| [sccomp](./sccomp/) | Beta-binomial | R | ✅ | Compositional constraints |
| [scanpro](./scanpro/) | Linear Model | R | ✅ | Proportions-based testing |
| [scellpam](./scellpam/) | Negative Binomial | R | ❌ | PAM clustering integration |
| [DirichletReg](./DirichletReg/) | Dirichlet | R | ✅ | Dirichlet regression framework |
| [Citrus](./Citrus/) | Regularized Regression | R | ✅ | High-dimensional feature selection |

## 🎯 When to Use Clustering-Based Methods?

### ✅ **Advantages:**
- **Interpretability**: Results map directly to biologically meaningful cell types
- **Stability**: Less sensitive to parameter tuning than clustering-free methods
- **Integration**: Easily combines with existing cell type annotation pipelines
- **Efficiency**: Generally faster for datasets with clear cell type structure

### ⚠️ **Limitations:**
- **Resolution dependent**: Results influenced by clustering granularity
- **Discrete boundaries**: May miss subtle state transitions within clusters
- **Annotation quality**: Dependent on accurate cell type identification


