# Cohort-Provenance Bias Audit

Diagnostic tools for detecting and quantifying non-biological signal in
biomarker classifier evaluation. These are generic miRNA/omics QC checks
that every biomarker study should pass before any AUC claim. They expose
confounding from:

- cohort identity (institution, batch, dataset),

- demographic covariates (age, sex, BMI),

- pre-analytical asymmetry (storage, collection site, tube type),

- specimen duplication across accessions.

## Details

The motivating problem is simple: on the four public 3D-Gene
ovarian-cancer cohorts (GSE106817, GSE211692, GSE113486, GSE113740), a
logistic regression using *only* dataset-identity dummy variables
reaches AUC 0.72 against cancer/healthy labels without touching a single
miRNA value. Any classifier trained on the same cohort structure
therefore inherits this "bias floor" even if its apparent performance is
higher.

The functions here make that floor visible before any biological claim
is made, and extend to:

- covariate-only AUCs (age, sex, etc.),

- clustered bootstrap CIs that respect specimen duplication,

- feature-level batch-signal ANOVA within healthy samples,

- case/control pre-analytical asymmetry tables.

## References

Leek JT, Scharpf RB, Bravo HC, et al. (2010). Tackling the widespread
and critical impact of batch effects in high-throughput data. *Nature
Reviews Genetics*, 11(10), 733-739.

Pritchard CC, Kroh E, Wood B, et al. (2012). Blood cell origin of
circulating microRNAs: a cautionary note for cancer biomarker studies.
*Cancer Prevention Research*, 5(3), 492-497.

Kirschner MB, Kao SC, Edelman JJ, et al. (2011). Haemolysis during
sample preparation alters microRNA content of plasma. *PLoS ONE*, 6(9),
e24145.
