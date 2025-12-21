# Differential Abundance Analysis Methods Benchmark

This repository hosts the complete code and data for benchmarking differential abundance (DA) analysis methods for single-cell RNA-seq data, supporting the review article *"Uncovering Cellular Composition Changes: Cutting-Edge Differential Abundance Methods for Single-Cell RNA-Sequencing"*.


## 📁 Repository Structure
| Directory | Description |
|-----------|-------------|
| **[`methods/`](methods/)** | Implementation and documentation for all 28 evaluated DA methods |
| **[`clustering_free/`](methods/clustering_free/)** | Methods analyzing cellular neighborhoods without pre-clustering |
| **[`clustering_based/`](methods/clustering_based/)** | Methods operating on predefined cell type clusters |
| **[`benchmark/`](benchmark/)** | Complete benchmarking framework and evaluation code |
| **[`chunk1_analysis_pipeline/`](benchmark/chunk1_analysis_pipeline/)** | Core execution framework (Docker for R, automation scripts for Python) |
| **[`chunk2_real_datasets/`](benchmark/chunk2_real_datasets/)** | Processing code for 29 curated real datasets |
| **[`chunk3_simulation_datasets/`](benchmark/chunk3_simulation_datasets/)** | Synthetic data generation and simulation code |
| **[`chunk4_usability/`](benchmark/chunk4_usability/)** | Systematic usability assessment implementation |
| **[`chunk5_accuracy/`](benchmark/chunk5_accuracy/)** | Performance metric calculation and result aggregation |
