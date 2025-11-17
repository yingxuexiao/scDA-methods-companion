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
| [scCODA](https://www.nature.com/articles/s41467-021-27150-6) | Bayesian | Python | ✅ | Compositional analysis with credible intervals |
| [tascCODA](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2021.766405/full) | Bayesian| Python |✅ |Tree-structured hierarchical modeling|
| [propeller](https://academic.oup.com/bioinformatics/article/38/20/4720/6675456) | Linear Model | R | ✅ | Arcsin transformation for proportions |
| [diffcyt](https://www.nature.com/articles/s42003-019-0415-5) | GLM | R | ✅ | Designed for cytometry data, generalizable to scRNA-seq |
| [MASC](https://www.science.org/doi/abs/10.1126/scitranslmed.aaq0305) | GLMM | R | ✅ | Mixed models for repeated measures |
| [DCATS](https://link.springer.com/article/10.1186/s13059-023-02980-3) | Beta-binomial | R | ✅ | Handles classification uncertainty |
| [sccomp](https://www.pnas.org/doi/abs/10.1073/pnas.2203828120) | Beta-binomial | R | ✅ | Compositional constraints |
| [scanpro](https://www.nature.com/articles/s41598-024-66381-7) | Linear Model | R | ✅ | Proportions-based testing |
| [scellpam](https://link.springer.com/article/10.1186/s12859-023-05569-6) | Negative Binomial | R | ✅ | PAM clustering integration |
| [DirichletReg](https://www.researchgate.net/profile/Marco-Maier-4/publication/260191096_DirichletReg_Dirichlet_Regression_for_Compositional_Data_in_R/links/02e7e5300fd1412f55000000/DirichletReg-Dirichlet-Regression-for-Compositional-Data-in-R.pdf) | Dirichlet | R | ✅ | Dirichlet regression framework |
| [Citrus](https://www.pnas.org/doi/abs/10.1073/pnas.1408792111) | Regularized Regression | R | ✅ | High-dimensional feature selection |
| [TreecorTreat](https://www.biorxiv.org/content/10.1101/2021.10.27.466024v1.full)| 
| [CTDS](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.855076/full)|Tree-based Correlation $ regression model | R | ✅ | Hierarchical correlation screening | 
| [ELVAR](https://www.nature.com/articles/s41467-023-39017-z)| NB-GLM | R | ✅ | Attribute-aware community detection |
| [clustermap](https://academic.oup.com/bioinformatics/article/35/17/3038/5289328?login=false)| Tree Pruning | Python | ✅ | Multi-sample cluster alignment |
| [scPopCorn](https://www.sciencedirect.com/science/article/pii/S2405471219301887)|Personalized PageRank | Python | ✅ | Cross-sample cell-to-cell matching |
| [Louvain+GLM]|NB-GLM | R | ✅ | Graph-based clustering with GLM testing |
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


