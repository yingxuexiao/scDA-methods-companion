
# MASC: Mixed-effects association testing for single cells

## 📋 Overview
MASC is a novel reverse single cell association strategy for testing whether a specified covariate influences the membership of single cells in any of multiple cellular subsets while accounting for technical confounds and biological variation.

## 🎯 Core Principle
- Uses mixed-effects logistic regression models
- Accounts for donor-level variability and repeated measures
- Handles hierarchical experimental designs
- Tests for condition effects while controlling for random effects

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Mixed-effects Logistic Regression
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: MASC, lme4, stats


## 💡 Biological Applications
- Multi-donor single-cell studies
- Clinical cohort studies with multiple samples
- Hierarchical experimental designs

## 🔗 Official Resources
- **GitHub**: https://github.com/immunogenomics/MASC
- **Publication**: [Science Translational Medicine](https://www.science.org/doi/abs/10.1126/scitranslmed.aaq0305)


## ⚙️ Installation

### GitHub Installation
```r
if (!require(devtools)) install.packages("devtools")
library(devtools)
install_github("immunogenomics/masc")