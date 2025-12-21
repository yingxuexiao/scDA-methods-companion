# Differential Abundance Analysis Methods

This repository hosts the complete implementation, documentation, and benchmarking framework for differential abundance (DA) analysis methods in single-cell RNA-seq data. It supports the review article *"Uncovering Cellular Composition Changes: Cutting-Edge Differential Abundance Methods for Single-Cell RNA-Sequencing"* and provides two core resources:



## 📁 Repository Structure

### 1. **Methods Catalog** - Comprehensive Collection of DA Methods
This section systematically compiles and documents all 28 differential abundance analysis methods evaluated in our study. Each method is accompanied by implementation details, algorithm descriptions, and example usage code, providing a valuable reference for researchers seeking to understand or apply these tools.

| Directory | Description |
|-----------|-------------|
| **[`methods/`](methods/)** | Root directory containing all method implementations and documentation |
| **[`clustering_free/`](methods/clustering-free/)** | Methods operating at single-cell or neighborhood resolution without relying on pre-clustered cell types (e.g., Milo, CNA, DA-seq) |
| **[`clustering_based/`](methods/clustering-based/)** | Methods that require predefined cell type annotations and analyze differential abundance at the cell-type level (e.g., scCODA, DCATS, propeller) |

### 2. **Benchmark Framework** - Systematic Evaluation Suite
This section contains the complete experimental framework for our comprehensive performance assessment. It enables the reproduction of all analyses presented in the manuscript, from data preprocessing to final performance metrics.

| Directory | Description |
|-----------|-------------|
| **[`benchmark/`](benchmark/)** | Root directory for the entire benchmarking workflow |
| **[`chunk1_analysis_pipeline/`](benchmark/chunk1_analysis_pipeline/)** | Core execution engine: Docker-based pipeline for R methods and automation scripts for Python methods, ensuring reproducible execution across all evaluated tools |
| **[`chunk2_real_datasets/`](benchmark/chunk2_real_datasets/)** | Processing pipelines for 29 curated real-world scRNA-seq datasets, including standardized preprocessing, quality control, and annotation procedures |
| **[`chunk3_simulation_datasets/`](benchmark/chunk3_simulation_datasets/)** | Code for generating synthetic benchmark data, including scalability test series and pattern-specific simulations with known ground truth |
| **[`chunk4_usability/`](benchmark/chunk4_usability/)** | Implementation of the systematic usability assessment framework based on 27 evaluation criteria covering code quality, documentation, and maintenance |
| **[`chunk5_accuracy/`](benchmark/chunk5_accuracy/)** | Performance metric calculators and result aggregation scripts for evaluating method accuracy across all benchmark datasets |

