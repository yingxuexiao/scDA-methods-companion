# Pipeline for Differential Abundance (DA) Method Benchmarking

This directory contains the pipeline scripts and Docker configuration files for running differential abundance (DA) methods on simulated and real datasets. The pipeline is designed to be executed within a Docker container to ensure a reproducible environment for evaluating DA methods.

## 🛠️ **Pipeline Overview**

The pipeline consists of several stages:

1. **Data Preprocessing**: Raw data is processed to ensure compatibility with DA methods (normalization, feature selection, etc.).
2. **DA Method Execution**: Different DA methods are applied to the processed data.
3. **Performance Evaluation**: The output of DA methods is evaluated using various metrics (e.g., precision, recall, F1-score).

The pipeline is designed to be run inside a Docker container, ensuring a consistent environment for the analyses.

