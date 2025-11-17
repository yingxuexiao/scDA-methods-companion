
# MrVI: Multi-resolution Variational Inference

## 📋 Overview
MrVI uses hierarchical variational autoencoders with multi-head attention to model single-cell data, enabling differential abundance analysis in a continuous latent space with sophisticated batch effect correction.

## 🎯 Core Principle
- Variational autoencoder architecture
- Multi-head attention for batch effects
- Continuous latent space modeling
- Sample-level and cell-level effect separation

## 🔧 Technical Implementation
- **Language**: Python
- **Core Model**: Variational Autoencoder + Linear Model
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: scvi-tools, pytorch

## ⚙️ two latent representations
- `u`: designed to capture broad cell states invariant to sample and nuisance covariates
- `z`: “sample-aware” representation of a cell, invariant to nuisance covariates. z augments u with sample-specific effects but remains corrected for nuisance covariate effects.


## 💡 Biological Applications
- Complex batch effect integration
- Multi-condition experimental designs
- High-dimensional data analysis

## 🔗 Official Resources
- **GitHub**: https://scvi-tools.readthedocs.io/en/latest/user_guide/models/mrvi.html
- **Publication**: [Nature Methods](https://www.nature.com/articles/s41592-025-02808-x)

##🚨 This standalone package will receive no further updates or support.Use instead the implementation in the scvi-tools package:

```
!pip install --quiet scvi-colab
from scvi_colab import install

install()

import os
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
from scvi.external import MRVI

scvi.settings.seed = 0  # optional: ensures reproducibility
print("Last run with scvi-tools version:", scvi.__version__)
save_dir = tempfile.TemporaryDirectory()

