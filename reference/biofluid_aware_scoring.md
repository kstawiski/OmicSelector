# Biofluid-aware within-sample scoring methods

These scorers treat serum, plasma anticoagulant class, exosome-enriched
fractions, and unknown biofluid status as first-class pre-analytical
strata. All functions expect expression matrices with samples in rows
and miRNAs in columns. Fitting functions use training data only;
prediction functions do not read disease labels from test metadata.
