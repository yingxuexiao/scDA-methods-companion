
# DirichletReg: Dirichlet regression for compositional data

## 📋 Overview
DirichletReg implements Dirichlet regression for modeling compositional data, allowing differential abundance analysis of cell type proportions while respecting the compositional nature of the data.

## 🎯 Core Principle
- Uses Dirichlet distribution for compositional data modeling
- Fits regression models to composition data
- Handles multivariate proportional outcomes
- Provides interpretable coefficients for proportion changes

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Dirichlet Regression
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: DirichletReg, stats

## ⚠️ Important Considerations
- Requires data where proportions sum to 1
- Choice of reference category affects interpretation
- May require extension for overdispersion/zero inflation
- Model interpretation can be complex

## 💡 Biological Applications
- Multivariate proportional outcome modeling
- Cell type interdependence studies
- Proportional data with multiple components

## 🔗 Official Resources
- **GitHub**: https://cran.r-project.org/web/packages/DirichletReg/index.html
- **Publication**: (https://www.academia.edu/28801683/DirichletReg_Dirichlet_Regression_for_Compositional_Data_in_R)


## ⚙️ Installation

### CRAN Installation
```r
install.packages("DirichletReg")