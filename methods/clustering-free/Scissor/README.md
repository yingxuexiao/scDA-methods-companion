
# Scissor: Single-cell identification of subpopulations with bulk sample phenotypes

## 📋 Overview
Scissor identifies phenotype-associated cell subpopulations by linking single-cell data to bulk sample phenotypes using network-regularized sparse regression. Scissor is a novel approach that utilizes the phenotypes, such as disease stage, tumor metastasis, treatment response, and survival outcomes, collected from bulk assays to identify the most highly phenotype-associated cell subpopulations from single-cell data. 

## 🎯 Core Principle
- Network-regularized sparse regression
- Links single-cell to bulk phenotypes
- Identifies phenotype-associated cells
- Uses correlation-based feature selection

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Regularized Sparse Regression
- **Multi-group**: ❌ Binary comparisons only
- **Dependencies**: Scissor, glmnet

## ⚙️ Key Parameters
- `alpha`: Parameter α balances the effect of the L1-norm and the network-based penalties
- `family`: Model family ("cox", "binomial")


## 💡 Biological Applications
- Clinical outcome association
- Biomarker discovery
- Treatment response prediction

## 🔗 Official Resources
- **GitHub**: https://github.com/sunduanchen/Scissor
- **Publication**: [Nature Biotechnology](https://www.nature.com/articles/s41587-021-01091-3)

## ⚙️ Installation

```r
devtools::install_github('sunduanchen/Scissor')