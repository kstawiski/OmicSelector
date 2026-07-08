# Modern Tabular Learners for OmicSelector

Optional mlr3 learner wrappers for modern tabular methods commonly used
as strong baselines on biomarker-style classification problems. The
wrappers are designed for small-\\P\\, medium-\\N\\, and
class-imbalanced omics tasks.

## Details

The learners in this file are intentionally lazy about backend
availability: they can be constructed without the Python or CatBoost
runtime being present, but training will fail with an actionable error
if the required dependency is missing.

Backends: - \`make_tabnet_learner()\`: \`pytorch_tabnet\` through
\`reticulate\` - \`make_tabm_learner()\`: \`tabm\` or
\`rtdl_revisiting_models\` through \`reticulate\` -
\`make_fttransformer_learner()\`: \`rtdl_revisiting_models\` through
\`reticulate\` - \`make_catboost_learner()\`: \`catboost\` R package -
\`make_tabpfn_learner()\`: \`tabpfn\` through \`reticulate\`

GPU selection follows \`torch.cuda.is_available()\` when \`device =
"auto"\`.
