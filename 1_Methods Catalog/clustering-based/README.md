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
| [scCODA](scCODA/README.md) | Bayesian | Python | ✅ | Compositional analysis with credible intervals |
| [tascCODA](tascCODA/README.md) | Bayesian| Python |✅ |Tree-structured hierarchical modeling|
| [propeller](propeller/README.md) | Linear Model | R | ✅ | Arcsin transformation for proportions |
| [diffcyt](diffcyt/README.md) | GLM | R | ✅ | Designed for cytometry data, generalizable to scRNA-seq |
| [MASC](MASC/README.md) | GLMM | R | ✅ | Mixed models for repeated measures |
| [DCATS](DCATS/README.md) | Beta-binomial | R | ✅ | Handles classification uncertainty |
| [sccomp](sccomp/README.md) | Beta-binomial | R | ✅ | Compositional constraints |
| [scanpro](scanpro/README.md) | Linear Model | R | ✅ | Proportions-based testing |
| [scellpam](scellpam/README.md) | Negative Binomial | R | ✅ | PAM clustering integration |
| [DirichletReg](dirichletReg/README.md) | Dirichlet | R | ✅ | Dirichlet regression framework |
| [Citrus](CITRUS/README.md) | Regularized Regression | R | ✅ | High-dimensional feature selection |
| [TreecorTreat](TreeCorTreat/README.md)| Tree-based Correlation | R | ✅ | Hierarchical correlation screening |
| [ELVAR](RLVAR/README.md)| NB-GLM | R | ✅ | Attribute-aware community detection |
| [Clustermap](Clustermap/README.md)| Tree Pruning | Python | ✅ | Multi-sample cluster alignment |
| [scPopCorn](scPopCorn/README.md)|Personalized PageRank | Python | ✅ | Cross-sample cell-to-cell matching |
| [Louvain+GLM](Louvain+GLM/README.md)|NB-GLM | R | ✅ | Graph-based clustering with GLM testing |
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

## 📝 Implementation Notes
### Code Availability
- **Methods with Custom Code**: We provide executable code examples and reproducible workflows for selected methods
- **Methods without Custom Code**: For certain methods, we do not provide custom implementations due to:
  - Technical complexity requiring specialized expertise
  - Platform-specific dependencies (e.g., flow cytometry data formats)
  - Known stability issues in current versions
 
### Alternative Resource Access
For methods without our custom code, each README file includes:
- **Direct GitHub Links**: Access to official implementations and source code
- **Complete Documentation**: Installation guides and usage instructions
- **Reference Publications**: Links to original methodological papers


