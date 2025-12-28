# Pipeline for Differential Abundance (DA) Method Benchmarking

## 📋 Overview

This repository provides a comprehensive benchmarking pipeline for evaluating Differential Abundance (DA) methods on single-cell RNA sequencing data. The pipeline systematically compares multiple DA methods across both simulated and real datasets to assess their performance in detecting cell type composition changes.

The pipeline employs a **dual-environment architecture** to handle the diverse computational requirements of R-based and Python-based methods:

### 🐳 **Docker Environment (R Methods)**
- **Purpose**: Execute R-based DA methods in an isolated, reproducible environment
- **Base Image**: `r-based:4.4.0` with comprehensive R package ecosystem
- **Included R Methods**:
  - `Milo`
  - `Dawnn`
  - `DCATS`
  - `CNA`
  - `Louvain+GLM`
  - `DirichletReg`
  - `ELVAR`
  - `sccomp`
  - `DA-seq`
  - `cydar`
  - `propeller`
  

### 🐍 **Conda Environments (Python Methods)**
- **Purpose**: Run Python-based DA methods in separate, optimized environments
- **Rationale**: Python methods often have conflicting dependencies and require different Python versions
- **Managed Environments**:
  - `scCODA` Python version: 3.10.19
  - `MELD` Python version: 3.10.19  
  - `CNA` Python version: 3.8.20
  - `scanpro` Python version: 3.10.19
  

