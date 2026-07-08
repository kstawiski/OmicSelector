# Deep Learning Learners for Omics Data

Provides deep learning models optimized for high-dimensional omics data
using mlr3torch and torch. Includes TabTransformer and autoencoder-based
representation learning.

## Details

These learners require the 'mlr3torch' and 'torch' packages to be
installed. They are optional components that provide state-of-the-art
deep learning capabilities for biomarker discovery.

Available architectures: - \*\*MLP\*\*: Multi-Layer Perceptron with
dropout regularization - \*\*TabTransformer\*\*: Attention-based model
for tabular data - \*\*Autoencoder\*\*: Unsupervised feature compression
(PipeOp)

## Installation

“\`r install.packages("torch") torch::install_torch() \# Downloads
LibTorch install.packages("mlr3torch") “\`
