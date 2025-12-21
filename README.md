# Differential Abundance Analysis Methods Benchmark

This repository hosts the complete code and data for benchmarking differential abundance (DA) analysis methods for single-cell RNA-seq data, supporting the review article *"Uncovering Cellular Composition Changes: Cutting-Edge Differential Abundance Methods for Single-Cell RNA-Sequencing"*.


## 📁 Repository Structure
| Directory | Description |
|-----------|-------------|
| **`methods/`** | Implementation and documentation for all 28 evaluated DA methods |
| `methods/clustering_free/` | Methods analyzing cellular neighborhoods without pre-clustering |
| `methods/clustering_based/` | Methods operating on predefined cell type clusters |
| **`benchmark/`** | Complete benchmarking framework and evaluation code |
| `benchmark/pipeline/` | Core execution framework (Docker for R, automation scripts for Python) |
| `benchmark/real_datasets/` | Processing code for 29 curated real datasets |
| `benchmark/simulation_datasets/` | Synthetic data generation and simulation code |
| `benchmark/usability/` | Systematic usability assessment implementation |
| `benchmark/accuracy/` | Performance metric calculation and result aggregation |
| `benchmark/efficiency/` | Computational resource profiling (integrated in pipeline) |