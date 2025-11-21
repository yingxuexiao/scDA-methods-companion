# Real Data for DA Method Benchmarking

This directory contains real-world datasets used to evaluate the performance of various differential abundance (DA) methods for single-cell RNA-seq data.

## 🧬 **Datasets**

### 1. **Dataset 1: Mouse Cortex Data**
- **Description**: A real single-cell RNA-seq dataset from mouse cortex tissue, used for analyzing neuronal cell types.
- **Source**: [DOI or publication link]
- **Contents**: Gene expression data with metadata for cell types and experimental conditions.

### 2. **Dataset 2: Human PBMC Data**
- **Description**: A dataset from human peripheral blood mononuclear cells (PBMCs), used for immune cell analysis.
- **Source**: [DOI or publication link]
- **Contents**: Gene expression data with metadata for cell types and treatment conditions.

## 🛠️ **Using the Real Data**

1. **Clone the repository**:
    ```bash
    git clone https://github.com/yourusername/scDA-methods-companion.git
    cd scDA-methods-companion/benchmark/real_data
    ```

2. **Load Data**:
    Use the `Seurat` package to load the real data:
    ```r
    library(Seurat)
    data <- readRDS("mouse_cortex_data.rds")
    ```

3. **Run the Benchmarking Pipeline**:
    Real data can be directly used in the benchmarking pipeline by following the steps outlined in the `pipeline/README.md` file.

## 📝 **Data Format**

Each dataset is provided in `.rds` format, with the following key components:
- **Counts matrix**: Gene expression data.
- **Metadata**: Information about the samples, cell types, and conditions.

