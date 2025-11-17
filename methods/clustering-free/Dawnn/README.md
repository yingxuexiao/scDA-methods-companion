
# Dawnn: Deep learning for differential abundance analysis

## 📋 Overview
Dawnn uses deep neural networks to learn complex, non-linear mappings from single-cell data to experimental conditions, enabling end-to-end differential abundance analysis.

## 🎯 Core Principle
- Deep neural network architecture
- Learns condition probability for each cell
- Handles complex, non-linear patterns
- End-to-end training without manual feature engineering

## 🔧 Technical Implementation
- **Language**: R
- **Core Model**: Neural Networks
- **Multi-group**: ❌ Binary comparisons only
- **Dependencies**: tensorflow

## ⚙️ Key Parameters
- `nn_model`: String containing the path to the model's .hdf5 file (default ~/.dawnn/dawnn_nn_model.h5).
- `recalculate_graph`: Boolean whether to recalculate the KNN graph. If FALSE, then the one stored in the ‘cells’ object will be used (default TRUE).
- `alpha`: Numeric target false discovery rate supplied to the Benjamini–Yekutieli procedure (default 0.1, i.e. 10%).
- `verbosity`: Integer how much output to print. 0: silent; 1: normal output; 2: display messages from predict() function. (default 2)

## 💡 Biological Applications
- Complex non-linear abundance patterns
- High-dimensional data integration
- Novel cell state discovery

## 🔗 Official Resources
- **GitHub**: https://github.com/george-hall-ucl/dawnn
- **Publication**: [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.05.05.539427v1)

## ⚠️ Important Note
**Configuration Requirement**: Dawnn requires a pre-trained model file at `~/.dawnn/dawnn_nn_model.h5`. If this file is missing, the method will fail to run. You may need to manually create this directory or download the model file.


## ⚙️ Installation

```r
# Step 1: Install Dawnn package (may need to install `remotes` package first)
remotes::install_github("george-hall-ucl/dawnn")

# Step 2: Download Dawnn's model
# By default, model stored at ~/.dawnn/dawnn_nn_model.h5
dawnn::download_model()

# Step 3: Install Tensorflow Python package in Reticulate environment
reticulate::py_install("tensorflow")