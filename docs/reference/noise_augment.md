# Gaussian Noise Augmentation

Adds Gaussian noise to samples for data augmentation. Simpler than SMOTE
but maintains sample independence.

## Usage

``` r
noise_augment(task, n_augment = 1L, noise_sd = 0.1, class = NULL)
```

## Arguments

- task:

  An mlr3 classification task

- n_augment:

  Number of augmented samples per original (default: 1)

- noise_sd:

  Standard deviation of noise as fraction of feature SD (default: 0.1)

- class:

  Target class to augment (default: NULL = minority)

## Value

A new task with augmented data
