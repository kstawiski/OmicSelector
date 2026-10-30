# Threshold selection and numerical score precision

The [scoring comparison table](../../scoring_agreement/out/deployed_vs_benchmark.tsv) reports score differences, rank correlation, area under the receiver operating characteristic curve (AUC) and changed calls for each method and fold.
The [factorial coverage table](../../scoring_agreement/out/factorial_coverage.tsv) identifies the evaluated scoring conditions.
Threshold selection and thread-dependent preprocessing have different numerical effects.

## Threshold selection

The two threshold discrepancies concern GSE210329 TabPFN fold 4 and TabDPT fold 3.
The Youden rule selects the smallest midpoint that maximizes `TP*n_negative + TN*n_positive`.
This integer objective preserves ties between thresholds with identical class counts.
Floating-point sums of sensitivity and specificity can differ by one unit in the last place for these tied counts.
Consequently, `pROC::coords(..., "best")` can exclude tied midpoints.
The integer-based R and Python checks agree on threshold selection for the same score vectors.

## Thread-dependent preprocessing

TabPFN 8.0.8 applies `AddSVDFeaturesStep` followed by `AddFingerprintFeaturesStep`.
The first step adds singular value decomposition (SVD) features using `StandardScaler(with_mean=False)` and `TruncatedSVD(algorithm="arpack")`.
Changing CPU threads from four to six slightly changes the appended SVD features.
The second step rounds each row to 12 decimal places and hashes its bytes into an auxiliary fingerprint feature.
A small SVD difference near a rounding boundary can therefore produce a substantially different fingerprint.

GSE259327_DIAGNOSTIC fold 1 contains 161 training specimens, 42 held-out specimens and 2,050 features.
Both scoring configurations use the same CPU server.
Deployment uses four threads and one `score_specimen()` call per specimen.
Benchmark scoring uses six threads and batches of specimens.
The following table summarizes the comparison.

| Quantity | Value |
|---|---|
| Deployment package commit | `e68e243b068d0584211c9f25b8993072e391e983` |
| Benchmark package commit | `a55f6c2b99c986cb4ac897cced9c0996f5014c82` |
| Maximum absolute score difference | 0.0004056096076969773 |
| Spearman correlation | 0.9998379385787214 |
| AUC and calls | Identical between configurations |

Controlled preprocessing uses identical training-input hashes.
Only the 17 appended SVD columns differ in this fold.
The maximum difference before hashing is 1.2083667400020204e-12.
The largest fingerprint difference is 0.9801999482151503.
Training fingerprints differ for 73, 50, 58 and 72 rows in the four ensemble members with SVD features.
The four members without SVD features agree.
Native and forced AVX2 preprocessing at six threads produce identical input, rounded-input and fingerprint hashes across all five GSE259327 folds.

## Model replay and interpretation

Model replay uses the benchmark package on the same CPU server.
Four threads reproduce deployment scores within 6.556510919741854e-7.
Six threads reproduce stored benchmark scores exactly.
The two package implementations differ in cached query workspace and chunk size but use the same SVD and fingerprint implementation.
The thread-controlled replay supports a preprocessing explanation for this fold.

Comparisons above 1e-5 include their scoring configurations, metrics and direct four-versus-six-thread preprocessing evidence.
Preprocessing probes localize differences in all training inputs and the first eight query inputs.
Hashes cover every query event.
These probes omit the transformer cache and do not provide full probability replays for the other folds.
Attribution in those folds rests on the shared mechanism and the representative full replay.
Some comparisons also use different computing hosts, whose separate contributions have not been quantified.
The comparison table includes all score differences, including differences in AUC and calls.
