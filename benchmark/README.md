Single-Cell Differential Abundance (DA) Analysis Benchmark

📊 Overview
This repository contains the comprehensive benchmark implementation for evaluating 16 differential abundance (DA) analysis methods on single-cell RNA sequencing data, as described in our publication. The benchmark includes evaluation code, workflows, simulated datasets, real datasets, and metrics computation for thorough method assessment.


🎯 Evaluated Methods
We systematically evaluated 16 state-of-the-art DA methods for single-cell data analysis. The complete list of methods and their implementations are detailed in our publication.


🔍 Method Selection Criteria
Methods were selected based on the following rigorous criteria:


Criteria
✅ Method must be currently functional and maintained	
✅ Specifically designed for single-cell data analysis	
✅ Primary intent must be differential abundance analysis	

🚫 Exclusion Criteria
Several methods were excluded from our benchmark:

❌ scellpam: Primarily a clustering method, not designed for DA analysis

❌ TascCODA: Requires construction of cell-type hierarchy trees, dependent on well-defined cell-type relationships

❌ Other exclusions: Methods that were deprecated, not maintained, or primarily designed for other analytical purposes

📁 Repository Structure
benchmark/
├── 📂 pipeline/
│   ├── evaluation_pipeline/      # Main evaluation workflows
│   
├── 📂 data/
│   ├── simulated_datasets/       # Synthetic data for controlled evaluation
│   └── real_datasets/           # Biological datasets for validation
├── 📂 results/                  # Pre-computed benchmark results
└── 📂 config/                   # Configuration files

⚠️ Important Notes & Limitations
📝 Disclaimer: This benchmark provides an initial comprehensive evaluation of DA methods, but several limitations should be considered:

🔬 Current Limitations
📏 Evaluation Depth: The current analysis provides a broad overview rather than deep, method-specific optimization

🔧 Parameter Sensitivity: Limited exploration of parameter spaces for each method

🌐 Dataset Scope: Evaluation on a curated but not exhaustive set of biological scenarios

⏱️ Computational Resources: Performance metrics may vary with different computational environments



🎯 Primary Goals
This benchmark aims to:

Help researchers select appropriate DA methods for their specific biological questions

Provide transparent and reproducible evaluation workflows

Highlight method strengths and weaknesses across different scenarios

Serve as a foundation for future, more specialized benchmarks


⭐ If you find this benchmark useful, please consider starring this repository!



