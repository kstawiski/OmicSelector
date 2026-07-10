# Provenance-aware within-sample scoring methods

These functions implement scoring methods that explicitly account for
cohort and specimen-provenance structure. All functions expect an
expression matrix with samples in rows and features in columns. \`meta\`
must have one row per sample in the same order as \`expr\`.

## Details

The mixed-effects scorer is split into a training step and a prediction
step: disease labels are used only by \`fit_mixed_effects_scorer()\` on
the training pool. \`predict_mixed_effects_scorer()\` and
\`score_mixed_effects_scorer()\` with a supplied \`fit\` object do not
read disease labels from test metadata.
