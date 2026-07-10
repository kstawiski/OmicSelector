# Advanced single-sample learned and transport scorers

Thin exported aliases for the learned kit-aware, Group-DRO, and
Sinkhorn-OT scorer implementations used by the OmicSelector
single-sample scoring bank. The aliases keep experimental method names
under the \`os\_\` namespace while preserving the train/apply separation
in the underlying fit and score functions.

## Usage

``` r
os_fit_dann_kit_extended(...)

os_score_dann_kit_extended(...)

os_transform_dann_kit_extended(...)

os_fit_kit_conditional_vae(...)

os_score_kit_conditional_vae(...)

os_transform_kit_conditional_vae(...)

os_fit_icp_per_kit(...)

os_score_icp_per_kit(...)

os_predict_icp_per_kit_details(...)

os_fit_group_dro_scorer(...)

os_score_group_dro_scorer(...)

os_fit_sinkhorn_ot_scorer(...)

os_apply_sinkhorn_ot_scorer(...)

os_score_sinkhorn_ot_scorer(...)
```

## Arguments

- ...:

  Arguments passed to the underlying method-specific fit, transform,
  apply, or score function.
