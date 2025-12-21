# scanpro: Proportion-based differential abundance testing

## 📋 Overview
scanpro performs differential abundance analysis by testing for differences in cell type proportions across conditions using linear models on transformed proportion data.

## 🎯 Core Principle
- Analyzes cell type proportions directly
- Uses linear models for statistical testing
- Handles proportion data with appropriate transformations
- Supports both real and pseudo-replicates

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Linear Model
- **Multi-group**: ✅ Supports multiple conditions

## Parameters
- clusters: name of clusters/celltypes column in obs table
- sample: name of sample column in obs table
- cond: name of condition/group column in obs table
- transform: type of transformation; logit or arcsin, default is logit

## 💡 Biological Applications
- Cell type proportion analysis
- Simple experimental designs
- Large-scale screening studies
- Proportion-based differential testing

## 🔗 Official Resources
- **GitHub**: https://github.com/loosolab/scanpro
- **Publication**: [scientific reports](https://www.nature.com/articles/s41598-024-66381-7)


## ⚙️ Installation

### 
```bash
pip install scanpro