# Technology-aware within-sample scoring methods

These scorers keep the platform technology axis separate from kit
labels. \`technology_col\` must identify one of \`NGS\`, \`Toray\`,
\`qPCR\`, \`NanoString\`, or \`unknown\` for every sample, including the
aliases below. Training helpers fit only on the supplied training pool;
prediction helpers do not read test labels.

Alias matching ignores case and surrounding whitespace. TaqMan, TLDA,
OpenArray, Exiqon/LNA, SYBR and relative-quantity or delta-delta-Ct qPCR
labels map to \`qPCR\`. Small-RNA/miRNA sequencing and HTG EdgeSeq
counts map to \`NGS\`; 3D-Gene arrays map to \`Toray\`. These are assay
classes, not claims that measurement scales are interchangeable. No
abundance conversion is performed here. Kit and biofluid labels are not
changed.

The package has no general microarray class: \`microarray\`,
\`Affymetrix miRNA microarray\` and \`five-platform repeated assay\` are
refused. The historical \`array-other\` alias for \`NanoString\` is
retained only for compatibility; it does not classify other array
labels. Other unrecognized labels are also refused (deployment code
\`UNSEEN_METADATA_LEVEL\`). Direct helpers retain the \`unknown\`
sentinel for missing labels; deployment requires known metadata.
