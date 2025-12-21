
# Kompot: Gaussian process regression for differential abundance

## 📋 Overview
Kompot uses Gaussian process regression to model cellular density changes across conditions directly on phenotypic manifolds, providing smooth abundance field estimates.Kompot implements methodologies from the Mellon package for computing differential abundance and gene expression, with a focus on using Mahalanobis distance as a measure of differential expression significance. It leverages JAX for efficient computations and provides a scikit-learn like API with .fit() and .predict() methods.

## 🎯 Core Principle
- Gaussian process regression framework
- Continuous manifold modeling
- Smooth abundance field estimation
- Uncertainty quantification

## 🔧 Technical Implementation
- **Language**: Python
- **Core Model**: Gaussian Process Regression
- **Multi-group**: ✅ Supports multiple conditions


## 💡 Biological Applications
- Smooth abundance field estimation
- Continuous differentiation analysis
- Developmental trajectory comparison

## 🔗 Official Resources
- **GitHub**: https://github.com/settylab/Kompot
- **Publication**: [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.06.03.657769v3)

## ⚙️ installation:
```bash
pip install kompot

#For using the default diffusion map cell state representation:
pip install palantir

#For additional plotting functionality with scanpy integration:
pip install kompot[plot]

#For disk-backed storage with Dask support (recommended for large datasets):
pip install kompot[dask]

#To install all optional dependencies:
pip install kompot[all]