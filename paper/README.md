# OmicSelector study reproduction

This deposit contains analysis scripts, public-data accession manifests and expected tables for OmicSelector 2.6.10. `result_map.tsv` connects the scientific quantities to their expected tables and computing scripts. The release object archives provide saved fold scores, model states, cohort annotations and pinned software dependencies. Expected tables are separate from regenerated results.

## Summary reconstruction

Summary reconstruction uses saved per-fold scores and the specified statistical producers. It reproduces benchmark AUC and matched-baseline lift, companion dependence, normalization, re-assay, unusable-input outcomes, transfer, calibration, computing costs and figure inputs. Use a Linux x86-64 CPU with four threads and at least 36 GiB RAM. Budget approximately 30–90 minutes, including software extraction. Elapsed time depends on storage and the selected analyses.

```sh
export OMICSELECTOR_PAPER_ROOT="$PWD/work"
(cd ../bundle && sha256sum -c SHA256SUMS)
bash scripts/run_all_summary.sh --bundle ../bundle
```

## Representative refits

Refits reconstruct raw matrices from public deposits, prepare each training fold and fit the named methods. They cover GSE117063 benchmark and normalization cells, GSE210329 calibration, GSE188627 case mix, GSE268524 re-assay, three transfer targets and the published pancreatic-cancer rule in GSE85589. The published rule applies fixed weights. `config/refit_routes.json` gives the complete route and representative specification.

```sh
bash scripts/run_all_refit.sh --bundle ../bundle
bash scripts/run_all_replay.sh --bundle ../bundle
```

The representative benchmark uses three classical methods and five seeds. Normalization uses nine arms and five seeds. Calibration uses k-TSP and seed 101. Case mix evaluates eleven arms. Re-assay trains on array 1 and evaluates array 2. Transfer uses GSE210329, GSE320031 and GSE41922. The classical representatives take minutes to tens of minutes after setup. Allow up to two hours for the complete set on four CPU threads. Frozen replay covers five methods and three folds per method. Foundation-model inference supports CUDA. The numerical tolerance is 1e-6. Other declared comparisons use 1e-10 unless a route states a stricter criterion.

## Full recomputation

`data/accessions.tsv`, `data/raw_downloads.tsv` and `data/preparation_dependencies.tsv` identify public cohorts and their SHA-256-pinned inputs. Acquire a selected accession with `python3 data/acquire.py --accession GSE117063`. Prepare cohorts with `Rscript --vanilla data/prepare_all.R . GSE117063`. `scripts/benchmark/refit_registered.py --help` lists registered-cell fitting options. The manifest contains 831 cohort-source files. GSE71008 and PRJEB30542 include 532 FASTQs. `data/quantify_public.py --help` provides the nf-core/smrnaseq 2.4.1 route and requires Java, Nextflow and a supported container engine. Full quantification and all registered model fits require substantially more compute than the representative tier. Their total wall time and peak storage have not been measured. Plan a multi-core server, project storage sized for the FASTQs and quantification outputs, and a CUDA GPU for foundation-model fits. These routes require an independent resource estimate before execution.

GSE198834 label assignment uses hash-verified count matching of sample files to GEO identities. 154 profiles match exactly and GSM5957628 is the unique remaining pair. Counts differ for that pair. Deposited expression values are retained. The assigned label set contains 37 cases among 155 samples. Toray preparation follows the deposited sample and feature selection annotations and the 79-feature intersection. It does not independently validate biological annotations.

## Software and interpretation

`container/README.md` gives Docker and Apptainer recipes for CPU and CUDA execution. Software versions, source hashes and analysis dependencies are pinned in `container/` and the bundle manifest. Summary producers use the deposited R 4.6.1 runtime. Model execution uses R 4.5.2 and Python 3.12.3. Runtime and hardware measurements are in the scientific cost tables. Computational agreement establishes the declared reproduction contracts. Clinical validity requires separate evidence.
