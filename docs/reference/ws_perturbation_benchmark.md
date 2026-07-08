# Within-Sample Perturbation Benchmark

Benchmarks a set of classification learners under within-sample
perturbation conditions: baseline, per-sample additive shift, per-sample
scaling, and the combined perturbation. Learners are trained on the
original training fold and evaluated on a perturbed copy of the held-out
fold.

## Usage

``` r
ws_perturbation_benchmark(
  task,
  learners,
  conditions = c("baseline", "shift", "scale", "both"),
  resampling = mlr3::rsmp("cv", folds = 5L),
  shift_sd = 2,
  scale_range = c(0.5, 2),
  seed = 42L
)
```

## Arguments

- task:

  An \`mlr3::TaskClassif\` with numeric features.

- learners:

  A list of \`mlr3\` learners or learner IDs understood by
  \[mlr3::lrn()\].

- conditions:

  Character vector of perturbation conditions. Defaults to
  \`c("baseline", "shift", "scale", "both")\`.

- resampling:

  An \`mlr3\` resampling object. Defaults to 5-fold CV.

- shift_sd:

  Standard deviation of the per-sample additive perturbation. Defaults
  to \`2\`.

- scale_range:

  Length-2 numeric vector defining the uniform scaling range. Defaults
  to \`c(0.5, 2.0)\`.

- seed:

  Integer random seed controlling perturbation draws.

## Value

A list with class \`"ws_perturbation_benchmark"\` containing per-fold
scores, aggregated AUC by learner and condition, and per-fold
predictions.

## References

Aitchison J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. (2017). Microbiome
datasets are compositional: and this is not optional. *Frontiers in
Microbiology*, 8, 2224.

Quinn TP, Erb I, Richardson MF, Crowley TM. (2020). Understanding
sequencing data as compositions: an outlook and review.
*Bioinformatics*, 36(16), 4424-4432.

## Examples

``` r
if (FALSE) { # \dontrun{
task <- mlr3::tsk("sonar")
learners <- list("classif.rpart", "classif.clr_mlp")
ws_perturbation_benchmark(task, learners, seed = 1)
} # }
```
