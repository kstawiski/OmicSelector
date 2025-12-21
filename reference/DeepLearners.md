# Deep Learning Learners for Omics Data

Provides deep learning models optimized for high-dimensional omics data
using mlr3torch. Includes TabTransformer, GNNs, and other architectures.

## Details

These learners require the 'mlr3torch' and 'torch' packages to be
installed. They are optional components that provide state-of-the-art
deep learning capabilities for biomarker discovery.

Available architectures: - \*\*MLP\*\*: Multi-Layer Perceptron with
dropout regularization - \*\*TabTransformer\*\*: Attention-based model
for tabular data - \*\*GNN\*\*: Graph Neural Network for pathway-aware
classification

## Installation

“\`r install.packages("torch") torch::install_torch() \# Downloads
LibTorch install.packages("mlr3torch") “\`
