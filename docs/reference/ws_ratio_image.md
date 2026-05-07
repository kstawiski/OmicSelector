# Within-Sample Pairwise Ratio Image

Creates a p x p pairwise ratio matrix for each sample, suitable as input
to a CNN. pixel(i,j) = feature_i - feature_j (on log scale = log-ratio).

This encoding is analogous to color normalization in histopathology
imaging: batch effects (brightness/contrast) cancel in the pairwise
ratios, just as stain normalization preserves relative channel
intensities.

## Usage

``` r
ws_ratio_image(x)
```

## Arguments

- x:

  A numeric matrix (samples x features) or vector (single sample).
  Assumed log-scale.

## Value

A 3D array (samples x features x features) for matrix input, or a 2D
matrix (features x features) for vector input.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single sample -> 14x14 ratio image
img <- ws_ratio_image(c(5.2, 3.1, 8.7, 2.4, 6.1, 4.3, 7.2,
                         1.9, 5.5, 9.0, 3.8, 4.7, 6.8, 2.1))
image(img, main = "miRNA ratio image")

# Matrix -> array of images
imgs <- ws_ratio_image(matrix(rnorm(140), nrow = 10, ncol = 14))
dim(imgs)  # 10 x 14 x 14
} # }
```
