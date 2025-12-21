# CTDS: An Entropy-Based Metric to Compare Overall Cell Type Composition Across Samples

## 📋 Overview
CTDS quantifies cellular diversity and dynamic changes within cell populations using entropy-based metrics to identify differentially abundant trajectories and states.

## 🎯 Core Principle
- Uses entropy-based diversity metrics
- Quantifies cellular state heterogeneity
- Analyzes trajectory dynamics and transitions
- Identifies diversity changes across conditions

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Entropy-based Metric
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: Custom R scripts, entropy


## 🔗 Official Resources
- **GitHub**: https://github.com/tanya-karagiannis/Cell-Type-Diversity-Statistic
- **Publication**: [frontiers](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.855076/full)


## ⚙️ Installation

### Manual Installation
```r
# CTDS is implemented as custom R scripts
# Download from publication supplementary materials
source("CTDS_functions.R")
#Before running the CTDS.score function, make sure that R is installed.
install.packages("tidyverse")
install.packages("Seurat")
install.packages("SingleCellExperiment")
install.packages("hablar")
install.packages("broom")
install.packages("patchwork")