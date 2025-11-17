
# CellCnn: Convolutional neural networks for rare cell population detection

## 📋 Overview
CellCnn uses convolutional neural networks in a multi-instance learning framework to detect rare cell populations associated with specific phenotypes or conditions.

## 🎯 Core Principle
- Convolutional neural network architecture
- Multi-instance learning framework
- Max-pooling for rare cell detection
- Phenotype-supervised feature learning

## 🔧 Technical Implementation
- **Language**: Python
- **Core Model**: Convolutional Neural Networks
- **Multi-group**: ❌ Binary comparisons only
- **Dependencies**: cellcnn, tensorflow


## 💡 Biological Applications
- Rare cell population detection
- Biomarker discovery
- Phenotype-associated cell identification

## 🔗 Official Resources
- **GitHub**: https://github.com/eiriniar/CellCnn
- **Publication**: [Nature Communications](https://www.nature.com/articles/ncomms14825)

## ⚙️ Installation

### PyPI
```bash
Clone the CellCnn repository and checkout the python3 branch:
git clone -b python3 https://github.com/eiriniar/CellCnn.git
Go to the CellCnn root directory:
cd CellCnn
Install CellCnn and its dependencies:
pipenv install
The above steps have to be performed only once. Then, each time you want to perform a CellCnn analysis, go to the CellCnn root directory and activate the pipenv virtual environment by running:

pipenv shell