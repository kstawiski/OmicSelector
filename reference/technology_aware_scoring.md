# Technology-aware within-sample scoring methods

These scorers keep the platform technology axis separate from kit
labels. \`technology_col\` must identify one of \`NGS\`, \`Toray\`,
\`qPCR\`, \`NanoString\`, or \`unknown\` for every sample. Training
helpers fit only on the supplied training pool; prediction helpers do
not read test labels.
