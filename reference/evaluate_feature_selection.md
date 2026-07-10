# Evaluate Feature Selection for Leakage

Assesses whether a feature selection method has data leakage by checking
if TRAP features were selected.

## Usage

``` r
evaluate_feature_selection(selected_features, gold_standard)
```

## Arguments

- selected_features:

  Character vector of selected feature names

- gold_standard:

  Output from \[generate_gold_standard()\]

## Value

A list of class \`OmicSelectorEvaluation\` containing:

- n_causal_selected:

  Number of true causal features found

- n_trap_selected:

  Number of trap features found (should be 0)

- sensitivity:

  Proportion of causal features recovered

- precision:

  Proportion of selected features that are causal

- leakage_score:

  Proportion of selected features that are traps

- verdict:

  Human-readable assessment

- pass:

  Logical indicating if selection passed (no traps)

## Examples

``` r
if (FALSE) { # \dontrun{
gs <- generate_gold_standard(n_features = 200, n_causal = 10, n_trap = 15)

# A good selection: only causal features
good <- evaluate_feature_selection(gs$ground_truth$causal_features[1:5], gs)
print(good)  # PASS

# A leaky selection: includes trap features
bad <- evaluate_feature_selection(c("CAUSAL_01", "TRAP_01"), gs)
print(bad)  # FAIL: Data leakage detected
} # }
```
