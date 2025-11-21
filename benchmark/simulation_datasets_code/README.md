# Simulated Data for DA Method Benchmarking

This directory contains simulated datasets created to benchmark various differential abundance (DA) methods for single-cell RNA-seq data. These datasets allow for controlled experimentation with different scenarios, such as varying cell proportions, noise levels, and gene expression gradients.

## 🛠️ **Simulated Scenarios**

### 1. **Scenario 1: Varying Cell Proportions**
- **Description**: In this scenario, we simulate datasets with different proportions of cell types.
- **Objective**: Evaluate DA methods' ability to detect changes in cell-type proportions.

### 2. **Scenario 2: Varying Gene Expression**
- **Description**: Simulate datasets with different levels of gene expression variability.
- **Objective**: Assess how well DA methods can detect differentially expressed genes.

### 3. **Scenario 3: Noise and Data Sparsity**
- **Description**: Simulate datasets with added noise or sparse gene expression to evaluate the robustness of DA methods under less-than-ideal conditions.
- **Objective**: Determine the impact of noisy data on DA method performance.

### 4. **Scenario 4: Mixed Cell Types**
- **Description**: Generate datasets where cell types overlap or are difficult to separate.
- **Objective**: Test the performance of DA methods on datasets with challenging cell-type separability.
