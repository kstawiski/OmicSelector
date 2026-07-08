# Create a TabPFN Learner

Create an \`mlr3\` classification learner backed by \`tabpfn\` through
\`reticulate\`. This is a zero-shot tabular foundation model suited to
small and medium tabular problems. TabPFN may require gated model access
through Hugging Face depending on the selected model version.

## Usage

``` r
make_tabpfn_learner(
  model_version = NULL,
  class_weights = "balanced",
  device = c("auto", "cpu", "cuda"),
  seed = 20260331L
)
```

## Arguments

- model_version:

  Optional TabPFN model version identifier exposed by
  \`tabpfn.constants.ModelVersion\`, for example \`"V2_5"\`. Default:
  \`NULL\`, which uses the package default.

- class_weights:

  Class weighting scheme. Included for interface consistency, but
  ignored by the current wrapper because TabPFN does not expose
  training-time class weighting here.

- device:

  \`"auto"\` (default), \`"cpu"\`, or \`"cuda"\`.

- seed:

  Random seed used inside the backend. Default: \`20260331\`.

## Value

An \`mlr3::LearnerClassif\` with \`predict_type = "prob"\`.

## References

Hollmann, N., et al. (2025). TabPFN-2.
