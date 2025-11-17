
# tascCODA: Tree-aggregated compositional differential abundance analysis

## 📋 Overview
tascCODA extends scCODA by incorporating tree-structured hierarchies of cell types, enabling differential abundance analysis across multiple resolutions of cellular organization. Hereby, tascCODA can infer credible changes on the features (i.e. cell types/ASVs/taxa) of a high-throughput sequencing dataset, as well as effects on subgroups of the set of features, which are defined by a tree structure, for example a cell lineage hierarchy or a taxonomic tree.

## 🎯 Core Principle
- Incorporates cell type hierarchies (lineage trees)
- Bayesian model with tree-based aggregation
- Multi-resolution abundance testing
- Handles structured cell type relationships

## 🔧 Technical Implementation
- **Language**: Python
- **Core Model**: Bayesian Hierarchical Model
- **Multi-group**: ✅ Supports multiple conditions
- **Dependencies**: tasccoda, ete3, scCODA


## 💡 Biological Applications
- Hierarchical cell type analysis
- Developmental lineage studies
- Multi-resolution abundance testing
- Cell differentiation trajectory analysis


## 🔗 Official Resources
- **GitHub**: https://github.com/bio-datascience/tascCODA
- **Publication**: [Frontiers in Genetics](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2021.766405/full)


## ⚙️ Installation

### PyPI Installation
```bash
pip install tasccoda