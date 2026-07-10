# Fit the DA-cVAE single-sample discriminator

Trains a domain-adversarial conditional VAE on the standardised Stage-1
structural-ILR balances (encoder \\\to\mu\\, FIT-only conditional
decoder \\p(\mathrm{balances}\mid z,\mathrm{platform})\\, KL prior,
FIT-only gradient-reversal platform adversary, supervised disease label
head on \\\mu\\; via reticulate-python torch, full-batch Adam), EXPORTS
the encoder (trunk up to and including the \\\mu\\ head) and the disease
LABEL head to R as plain numeric arrays, fits a frozen kNN-distance
novelty reference over the train control latents, and DISCARDS the
python module. The fitted model holds no external pointer; scoring is
pure base-R with no platform/cohort input.

## Usage

``` r
fit_dacvae(X_train, y_train, meta_train = NULL, annotation, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples x features) of non-negative abundances with
  unique, non-empty feature names (colnames).

- y_train:

  0/1 disease labels (1 = case), length `nrow(X_train)`, with at least
  one case and one control.

- meta_train:

  Optional per-sample metadata carrying a platform/assay column (the
  conditional-decoder + adversary domain). With \\\ge 2\\ platforms the
  conditional decoder is always used; the adversary is engaged only when
  `grl_lambda > 0` (so `grl_lambda = 0` here is the no-GRL ablation:
  conditional decoder preserved, adversary off). Must have one row per
  row of `X_train`.

- annotation:

  data.frame with columns `feature` (miRNA id), `cluster`
  (genomic-cluster id), `gc` (GC fraction) — the Stage-1 SBP priors.

- hp:

  Optional list of hyperparameters (all frozen into the model). Fields:
  `enc_hidden` (encoder trunk widths; default `c(64L,32L)`),
  `dec_hidden` (FIT-only decoder widths; default `c(32L,64L)`),
  `latent_dim` (\\L\\; default `16L`), `activation` (`"relu"`/`"tanh"`),
  `epochs` (default `200L`), `lr` (default `1e-3`), `weight_decay`
  (default `1e-4`), `kl_beta` (\\\beta\\; default `1.0`), `recon_weight`
  (\\\rho\\; default `1.0`), `grl_lambda` (adversary strength; default
  `1.0`; `0` = no-GRL ablation), `domain_hidden` (FIT-only adversary
  widths; default `c(32L)`), `platform_col` (`NULL` = auto-detect),
  `novelty_k` (kNN novelty neighbours; default `5L`), `input_type`
  (`"abundance"` (default) or `"ct"` for qPCR cycle-threshold input,
  converted to relative abundance via a frozen per-miRNA train-median
  reference), `ct_efficiency` (assay efficiency \\E\\ for \\E^{-\Delta
  Ct}\\; default `2`), `ct_floor` (censored-Ct MZR floor abundance;
  `NULL` = fit-computed), `pseudocount` / `aggregate` /
  `min_balance_coverage` / `min_part` (Stage-1 balance controls),
  `min_features` (score-time overlap floor; default `3L`), `device`
  (`"cpu"`/`"cuda"`/`"auto"`), `seed` (default `42L`).

## Value

Object of class `dacvae_model`: a list with `feature_universe`, `sbp`,
`balance_names`, balance controls (`pseudocount`, `aggregate`,
`min_balance_coverage`), `center`/`scale` (frozen balance
standardisation), `weights` (exported encoder layers), `activation`,
`head_w`/`head_b` (frozen disease head), `latent_dim`,
`platform_levels`, `n_platforms` (platforms used by the conditional
decoder; `1` = unconditional), `adv_engaged` (`FALSE` = no-GRL
ablation), `grl_lambda`, `novelty_ref` (control-latent reference
matrix), `novelty_k`, `device`, `seed`, `hp`.

## See also

[`score_dacvae`](https://kstawiski.github.io/OmicSelector/reference/score_dacvae.md),
[`score_dacvae_novelty`](https://kstawiski.github.io/OmicSelector/reference/score_dacvae_novelty.md),
[`fit_ss_struct_ilr`](https://kstawiski.github.io/OmicSelector/reference/fit_ss_struct_ilr.md)
