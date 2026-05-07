# Paper 3 qPCR non-detect imputation (Module B add-on)

Two qPCR non-detect handlers used by the Paper-3 pipeline:

- [`qpcr_nondetect_impute`](https://kstawiski.github.io/OmicSelector/reference/qpcr_nondetect_impute.md):
  Bayesian hierarchical imputation via the `nondetects` Bioconductor
  package (McCall et al. 2014).

- [`qpcr_nondetect_lod_fallback`](https://kstawiski.github.io/OmicSelector/reference/qpcr_nondetect_lod_fallback.md):
  limit-of-detection fallback that fills undetermined Cts with a fixed
  value (default 40).

Both operate on a features-by-samples Ct matrix. NAs in the input
represent undetermined Cts. The Bayesian path requires both `nondetects`
and `HTqPCR` (Bioconductor); if either is unavailable, the function
falls back to LOD imputation with a logged message. This mirrors the
manuscript Methods section "Non-detect handling for quantitative PCR
data".

## References

McCall MN, McMurray HR, Land H, Almudevar A. (2014) On non-detects in
qPCR data. *Bioinformatics* 30(16): 2310-2316.
