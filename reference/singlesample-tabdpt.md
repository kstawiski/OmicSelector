# TabDPT in-context tabular-foundation discriminator (conditional single-sample)

Conditional-single-sample discriminator (family O) built on TabDPT, the
pretrained tabular foundation model of Ma et al. (2024). Like TabPFN and
TabICL, TabDPT classifies a query row by *in-context learning*: a
forward pass of a frozen pretrained transformer conditions on a provided
labelled TRAINING context (optionally retrieval-augmented) and emits the
query's class posterior. There is no per-dataset gradient training –
"fitting" only loads the training context into the model. This scorer is
the direct sibling of
[`fit_tabpfn`](https://kstawiski.github.io/OmicSelector/reference/fit_tabpfn.md)
and
[`fit_tabicl`](https://kstawiski.github.io/OmicSelector/reference/fit_tabicl.md),
swapping the in-context forward calls for TabDPT's.

**Why this is single-sample-deployable (conditional-in-N).** A specimen
is single-sample-classifiable iff the TRAINING context that conditions
the transformer is FROZEN into the model at fit time, so that at
inference the query's posterior depends only on that query row and the
frozen context – never on which other specimens are scored alongside it.
This implementation freezes the per-sample-rCLR-transformed training
context (rows + labels) and the feature universe at fit, and at score
conditions every query against that one frozen context. The frozen
context is the TRAINING rows only, so it is leakage-safe (no scored /
test data ever enters the context). The scorer pins the
context-reduction mode to `"subsample"` (a single context shared across
queries, not per-query retrieval) and passes the full frozen context (no
per-query subsampling), so the conditioning set is identical for every
query and never depends on the scored batch.

**Row-independent scoring modes.** With the frozen shared (subsample)
context and the full context passed to every query, TabDPT's
`predict_proba` is invariant to the row ORDER of the scored batch and to
the scored subset – a query's posterior does not depend on the other
queries (measured companion maxdiff \\0\\). The default deployment path
still evaluates `predict_proba` ONE QUERY ROW AT A TIME, so the \\n=1\\
forward path is used uniformly. The optional benchmark-only
`score_batch` path evaluates balanced chunks of query rows. That batched
path is leakage-free and AUC-faithful (\\\|dAUC\| = 0\\; specimen
ranks/verdicts unchanged), but it is not bit-identical to the \\n=1\\
row-by-row scorer: the installed build has a deterministic batch-size
numerical effect (observed `max|b-r| = 1.97e-3`), not a cross-row
coupling.

**Single-sample transform.** Each specimen is mapped to the
self-contained per-sample robust CLR (rCLR) over its OWN
strictly-positive support (geometric-mean centring on `v > 0`) in the
FROZEN feature universe; this is exactly invariant to per-specimen
scaling and uses no cross-row statistic. The SAME rCLR transform is
applied to the frozen training context at fit and to each query at
score. Universe features absent from a specimen carry the neutral rCLR
value `0` (centred log-ratio of an absent coordinate), matching the
missing-feature rule of the other roster scorers.

**Token-free pinned checkpoint.** The default TabDPT checkpoint
(`tabdpt1_2.safetensors`, the TabDPT v1.2 model of Ma et al. 2024) is
hosted in the public Hugging Face repository `Layer6/TabDPT` and is
fetched non-interactively by `hf_hub_download` – there is no license
acceptance or interactive token step (unlike TabPFN's gated V2.6). This
method pins that checkpoint version explicitly (via the installed
`tabdpt` package's pinned `_VERSION = "1_2"`) so the model is
reproducible and the download path is token-free.

**Orientation and output.** The score is the frozen-context posterior
\\P(\mathrm{case}) = P(y = 1)\\ for the query, in \\\[0, 1\]\\; larger =
more case-like. Degenerate queries (fewer than `hp$min_features`
universe features present, or an all-zero / empty specimen) return the
neutral probability `0.5`.

## References

Ma J, Thomas V, Hosseinzadeh R, Kamkari H, Labach A, Cresswell JC,
Golestan S, Yu G, Volkovs M, Caterini AL. (2024) TabDPT: Scaling Tabular
Foundation Models. arXiv:2410.18164.

Qu J, Holzmuller D, Varoquaux G, Le Morvan M. (2025) TabICL: A Tabular
Foundation Model for In-Context Learning on Large Data. *International
Conference on Machine Learning (ICML)*; arXiv:2502.05564.

Hollmann N, Muller S, Purucker L, et al. (2025) Accurate predictions on
small data with a tabular foundation model. *Nature* 637:319-326.
