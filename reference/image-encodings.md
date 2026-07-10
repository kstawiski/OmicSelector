# Image Encoding Methods for Biomarker Panels

Converts tabular biomarker expression data into 2D image representations
suitable for CNN classification. Four encoding methods are supported:

1.  `encode_simple_grid`: Fixed row-major layout in rows×cols grid

2.  `encode_corr_grid`: Features ordered by hierarchical correlation
    clustering

3.  `encode_deepinsight`: Feature positions learned via t-SNE on
    training data (Sharma et al., 2019, Scientific Reports)

4.  `encode_ratio_image`: Pairwise log-ratio matrix (see
    ratio-image-cnn.R)

These complement the within-sample normalization approaches in
`within-sample.R` and can be combined with `train_ratio_cnn` and related
CNN training utilities.

## References

Sharma A et al. (2019). DeepInsight: A methodology to transform a
non-image data to an image for convolution neural network architecture.
Scientific Reports 9:11399.

Aitchison J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

Quinn TP, Erb I, Richardson MF, Crowley TM. (2020). Understanding
sequencing data as compositions: an outlook and review.
*Bioinformatics*, 36(16), 4424-4432.
