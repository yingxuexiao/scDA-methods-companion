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

scellpam: Primarily a clustering method, not designed for DA analysis

TascCODA: Requires construction of cell-type hierarchy trees, dependent on well-defined cell-type relationships

Other exclusions: Methods that were deprecated, not maintained, or primarily designed for other analytical purposes



⚠️ Important Notes & Limitations

This benchmark provides an initial comprehensive evaluation of DA methods, but several limitations should be considered:

🔬 Current Limitations

1. Evaluation Depth: The current analysis provides a broad overview rather than deep, method-specific optimization

2. Parameter Sensitivity: Limited exploration of parameter spaces for each method

3. Dataset Scope: Evaluation on a curated but not exhaustive set of biological scenarios

4. Computational Resources: Performance metrics may vary with different computational environments



🎯 This benchmark aims to:

1. Help researchers select appropriate DA methods for their specific biological questions

2. Provide transparent and reproducible evaluation workflows

3. Highlight method strengths and weaknesses across different scenarios

4. Serve as a foundation for future, more specialized benchmarks


⭐ If you find this benchmark useful, please consider starring this repository!



