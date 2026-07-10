# Paper 3 additional within-sample methods (Module A, P2)

Two additional Module-A within-sample methods introduced in Paper 3
(Stawiski et al., in preparation), complementing the five core methods
in `singlesample-within-sample.R`.

Methods provided:

- [`fit_logistic_normal_eb`](https://kstawiski.github.io/OmicSelector/reference/fit_logistic_normal_eb.md)
  /
  [`apply_logistic_normal_eb`](https://kstawiski.github.io/OmicSelector/reference/apply_logistic_normal_eb.md):
  frozen-reference empirical-Bayes denoiser using logistic-normal
  posterior shrinkage (Aitchison & Shen 1980; Efron 2010). Per-feature
  depth-weighted CLR shrinkage toward a training-cohort prior.

- [`fit_frozen_quantile`](https://kstawiski.github.io/OmicSelector/reference/fit_frozen_quantile.md)
  /
  [`apply_frozen_quantile`](https://kstawiski.github.io/OmicSelector/reference/apply_frozen_quantile.md):
  monotone quantile calibrator for cross-platform mapping. Frozen
  empirical-CDF mapping from test-sample feature values to
  training-distribution quantile values (Bolstad et al. 2003; Hicks et
  al. 2018).

## References

Aitchison J, Shen SM. (1980) Logistic-normal distributions: some
properties and uses. *Biometrika* 67(2): 261–272.

Efron B. (2010) Large-Scale Inference. Cambridge University Press.

Bolstad BM, Irizarry RA, Astrand M, Speed TP. (2003) A comparison of
normalization methods for high density oligonucleotide array data based
on variance and bias. *Bioinformatics* 19(2): 185–193.

Hicks SC, Okrah K, Paulson JN, Quackenbush J, Irizarry RA, Bravo HC.
(2018) Smooth quantile normalization. *Biostatistics* 19(2): 185–198.

Stawiski K. (in preparation) Paper 3 of the OmicSelector programme.
