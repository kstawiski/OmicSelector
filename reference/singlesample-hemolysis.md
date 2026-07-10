# Paper 3 robust-regression hemolysis correction (Module B)

Frozen robust-regression hemolysis correction introduced in Paper 3
(Module B; Stawiski et al., in preparation). Performs per-feature robust
regression of log-expression on a hemolysis proxy score (and optionally
a platelet score), fitting only on training controls. At deployment, the
predicted nuisance contribution is subtracted from each feature.

Note: the existing package file `hemolysis-correction.R` provides
[`fit_hemolysis_prefilter()`](https://kstawiski.github.io/OmicSelector/reference/fit_hemolysis_prefilter.md)
(a different approach based on marker-ratio gating). This file provides
the Module-B fit-and-freeze approach and is named
`singlesample-hemolysis.R` to avoid collision.

Methods provided:

- [`fit_hemolysis_rr`](https://kstawiski.github.io/OmicSelector/reference/fit_hemolysis_rr.md):
  fit the per-feature robust-regression model.

- [`apply_hemolysis_rr`](https://kstawiski.github.io/OmicSelector/reference/apply_hemolysis_rr.md):
  apply the frozen model to new samples.

## References

Huber PJ. (1964) Robust Estimation of a Location Parameter. *The Annals
of Mathematical Statistics* 35(1): 73–101.

Stawiski K. (in preparation) Paper 3 of the OmicSelector programme.
