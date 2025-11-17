# scellpam: Single-cell analysis using PAM clustering

## 📋 Overview
scellpam performs differential abundance analysis using Partitioning Around Medoids (PAM) clustering followed by negative binomial regression to identify differentially abundant cell populations.

## ⚠️ Important Note
**Primary Purpose**: This package is primarily designed for single-cell clustering using PAM algorithm, not specifically for differential abundance analysis. DA functionality is achieved by combining PAM clustering with post-hoc negative binomial regression.


## 🎯 Core Principle
- Uses PAM clustering for robust cell type identification
- Applies negative binomial regression for count data
- Handles over-dispersed single-cell count data
- Tests for abundance differences in PAM-derived clusters

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Negative Binomial Regression
- **Dependencies**: scellpam, cluster, MASS


## 💡 Biological Applications
- Robust cell clustering with PAM algorithm
- Studies requiring stable cluster assignments
- Single-cell data with clear cluster structure

## ⚠️ Important Considerations
- High memory requirements for large datasets
- Long execution times for large-scale problems
- Requires prior determination of cluster number k
- May not scale well to very large single-cell datasets

## 🔗 Official Resources
- **GitHub**: https://cran.r-project.org/web/packages/scellpam/index.html
- **Publication**: [BMC Bioinformatics](https://link.springer.com/article/10.1186/s12859-023-05569-6#availability-of-data-and-materials)


## ⚙️ Installation

```r
install.packages("scellpam")