
# sccomp: Compositional differential abundance with sum constraints

## 📋 Overview
sccomp performs differential abundance analysis using a beta-binomial model with sum-to-one constraints to explicitly account for the compositional nature of single-cell data across cell types. sccomp, a generalised method for differential composition and variability analyses capable of jointly modelling data count distribution, compositionality, group-specific variability, and proportion mean-variability association, while being robust to outliers.

## 🎯 Core Principle
- Uses beta-binomial regression for overdispersed count data
- Applies sum-to-one constraints for compositional data
- Models cell type proportions with hierarchical shrinkage
- Provides Bayesian-style credible intervals

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Beta-binomial Regression
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: sccomp, rstan, brms


## 💡 Biological Applications
- Studies requiring explicit sum constraints
- Cell type proportion changes with dependencies
- Multi-condition experimental designs

## 🔗 Official Resources
- **GitHub**: https://github.com/MangiolaLaboratory/sccomp
- **Publication**: [pnas](https://www.pnas.org/doi/full/10.1073/pnas.2203828120)


## ⚠️ Important Considerations
- Results interpretation can be complex due to compositionality

## ⚙️ Installation

```r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")

# Step 1
BiocManager::install("sccomp")

# Step 2
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))

# Step 3
cmdstanr::check_cmdstan_toolchain(fix = TRUE) # Just checking system setting
cmdstanr::install_cmdstan()