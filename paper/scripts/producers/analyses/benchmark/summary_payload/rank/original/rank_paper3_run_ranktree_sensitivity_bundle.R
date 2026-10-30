#!/usr/bin/env Rscript

.public_prose_columns <- c('Aggregation','Control release binding','Control source record','Native scoring output or omission','Outcome through deployment adapter','Parity comparison scope','Parity records','Parity tolerance basis','Plain-language method description','Producer summary record','Real-data reason','Real-data records','Real-data split','Supplemental parity records','Synthetic assay records','Synthetic history','Synthetic release binding','Training preflight outcome','baseline','cause','ci_method','contrast_scope','criterion','decision','definition','derivation','detail','execution','execution_environment','interpretation','interval','limitation','mechanism','model_policy','note','orientation','orientation_note','path_description','reason','repair','resolution','scale_handling','score','score_package','score_source','secondary_score_rule','source_host','source_population','summary')
.public_schema <- c('old_calibrator'='reference_calibrator','old_n_seeds'='reference_n_seeds','old_status'='reference_status','old_complete'='reference_complete','new_complete'='analysis_complete','old_estimate'='reference_estimate','old_ci_low'='reference_ci_low','old_ci_high'='reference_ci_high','estimate_change'='estimate_minus_reference','transition'='coverage_comparison','Synthetic history'='Synthetic assay records','canonical_repair_session_seconds'='canonical_state_session_seconds','canonical_repair_peak_rss_bytes'='canonical_state_peak_rss_bytes','old_tp'='reference_tp','old_fn'='reference_fn','old_tn'='reference_tn','old_fp'='reference_fp','old_auc'='reference_auc','historical_attempts'='reference_implementation_attempts','historical_rejections'='reference_implementation_rejections','historical_silent_scores'='reference_implementation_silent_scores','sensitivity_rerun_minus_shipped_auc'='sensitivity_replay_minus_deposited_auc','max_score_difference_sensitivity_rerun_minus_shipped'='max_score_difference_sensitivity_replay_minus_deposited','primary_rerun_minus_shipped_auc'='primary_replay_minus_deposited_auc','max_score_difference_primary_rerun_minus_shipped'='max_score_difference_primary_replay_minus_deposited','n_matched_sensitivity_reruns'='n_matched_sensitivity_replays','rerun_reason'='replay_outcome','new_status'='analysis_status','old_test_profiles'='reference_test_profiles','new_test_profiles'='analysis_test_profiles','old_rows'='reference_rows','new_rows'='analysis_rows','old_label'='reference_label','old_evaluable'='reference_evaluable','new_evaluable'='analysis_evaluable','new_value'='analysis_value','historical_value'='reference_value','previous_path'='reference_path','previous_sha256'='reference_sha256','rerun_difference'='replay_difference','rerun_mean_auc'='replay_mean_auc')
.public_notes <- c('2a4617576ea43f672b367ee4308117ea66b53278f904ffc173bf601828e92df9'='Only the matched baseline is included. Other method contrasts belong to the benchmark analysis.','9e50328b7c7db559753dd182761a37ff7cefe4f3bac739e87e90679286074edf'='Intervals use the normalization REML/Knapp-Hartung algorithm and stopping rules with 2000 joint person draws.','177d91f71eceeef6ea083b23ddaec4cb6b6e171ab70cc24d8fd047e52c8e3bc4'='Normalization arms and the rCLR reference use the same methods in both comparisons.','21b5de0a1c57dff42fd7efff3a0803b2289c434205e738adf7517c11d3e2f811'='The computing host for the local session is unavailable.','ac55f980ece6922b975d33727ad6b1ff2a880ad15b330285379535c472cde8eb'='The local computing host is identified in the device benchmark evidence.','3807984e5bf12da75c7d7a2f8942c054f196d53b384e2023a90e1b9cbdfb1334'='The analysis uses the local computing host.','198b2e0ff57a6edfeb2900e176b41c7baef0386f4a8026c69841446006d94903'='The executed functions match the recorded release.','5ff88db66ada08f625d3df773988bad16804ce26b5b1b92cc161ff3f7efee46b'='Execution uses the recorded release.','e91d7e4e94d90e8583348f6854e7ccc568297db53ea87d9d03af605475ada9bd'='The result uses identical executed functions.','ad6354c009ceadefbde3084d5ce95eaeb39c83482b059085013ef2d38009889b'='The grouped person split uses seed-101 and frozen real_fixture.rds.','26716332fde421d24fbf9fa75b4cf6256c28f9de1707695f19f7b59b97186892'='The grouped person split uses seed-101 and real_fixture.rds.','35fc825a9940a82cf04a739a79d70980c1eab4397dd9d82504b92fc64b8e3068'='The analysis uses the held-out scores.','5ccc67591ef21f45731a5435b622b0f7a5d4023e9c5e04c9b00cd112f53bc9f2'='The analysis uses saved case-mix predictions.','1e2868f69080b53089c8317a51852ba9d6e3507520d358911d4d92f71c6036c6'='Training-only preparation uses OmicSelector 2.6.7.','3aa22e35dfe1d5d7a9e6691f9db9e53646360ff5504dbc235d0270af16be6a8a'='Saved companion-analysis fits use retained-intercept reconstruction with OmicSelector 2.6.8.','72f1e642d961d1b52fe9d339aa582699d2edd0ba01bb262bf91536143ad35f83'='Re-assay scores combine frozen comparators and retained-intercept reconstruction.','838d67a9b4d1595511dbb487aead5ec39bc55dce90e5e34de71f61b567433ff0'='Panel scoring uses training-only preparation on CPU. Per-method sessionInfo gives the environment.','a7ef312a21b41540693d08fad47815d3cfc6956d9f8376d363ef81e3568692fa'='Clinical operating points use training-only preparation with OmicSelector 2.6.8 on the GPU workstation. Per-method sessionInfo gives the environment.','9379350f467f8b1d34d03bc7e4d8bad3367e224cd9fb4f0454a817d06ea512b4'='Unusable-input scoring uses training-only preparation with OmicSelector 2.6.8 on the CPU server with 4 threads. The environment image is 73cbfec2. Analysis provenance gives the environment.','445ea9c058ce4c4e2516dee9c02284579460c0fe7a9386d586ee3709b7694544'='Across-unit summaries give the median and range. The active-cancer synthesis excludes GSE145176.','3a708009ec8c57d88acfc284d793473793e613d2b147d7dbfd34a88ee8cbb829'='Across-unit summaries give the median and range. Absolute AUCs are not meta-analyzed.','d7cf4a38ee6bb14c10e7748bc68dd55a90ed3971573dc7b04176b0605e289c82'='Across-unit summaries give the median and range.','b7bff465477f746bdf3cf47dcce88c6f29e18e68fb5122236449c37c397e7820'='OmicSelector 2.6.9 refits use the saved train/test allocation. Label-sensitivity assignments apply to training and test specimens. Balanced batches use those assignments.','22182a3e6cee8f25881cf724472b8bb260d66d25addf20a6e7cc16be69bf1c44'='Fitted models, thresholds and memberships are fixed. Label-sensitivity assignments and exclusions apply only to evaluation rows.','9921ff96144695105b45992f03a663796f8625ce01f7911123eb92fae9c2b27e'='Transfer summaries use transfer_performance_r3.tsv. Outcome counts use validated_cells.tsv.','3775ce4dab990d8287fa751a3d3053e3b4f88ed633da66c22ea18efd63227144'='The display uses a score-difference threshold. Source classification requires identical AUC and a strictly increasing score mapping.','04f34ee27ea83383c50313297df0685fe1189cb3b25819f2b41fa8c66246617e'='The re-assay analysis has 1350 evaluations per policy. The 7970 count does not describe this analysis.','08866403a3feabc5e5861b4215717b3c8c78fbed67634c791abfd99d61145357'='The score is 1 minus the published three-of-four cancer score. No model is fitted.','3815e63efc80cb808c0b4ea002a7ed0b23d6b7288ccfd82ab04acf9b72ea4d23'='Intervals use 2000 joint person-bootstrap percentile draws. Repeated appearances and splits stay with each person. Models, thresholds and batches are fixed.','630f1a471587ae82f758d73a1016060fade7b01a72dee0db4afd3c81aa63e15d'='Intervals use 2000 joint person-bootstrap percentile draws and the equal mean of weighted batch rates. Scores and batches are fixed.','36245e2550b2b2c6cf55b14f3c726e72fdf936a4c9d0ee53e01ac53a67a32feb'='Intervals use 2000 person-bootstrap percentile draws and the appearance-weighted rate. Scores and batches are fixed.','f5c6dacd06266f32a2650bffe4216cad025e654868bd4b3b3a087637013d6fdc'='The analysis uses 5 grouped stratified folds. Folds 1+2 are held out. The seed is 101.','92e79a4c40ebd1f295a57fea601ba710167a10b28f3a4d861c9bb268d6f40c49'='All methods use specimen metadata without imputation.','798c6e404e94bc24a01499e45bb232279e06610dcc61549dffd9065673743caa'='At least one mature anchor is absent from the input schema. The index is NA.','7b85b4252a87cfd474e11c1b33b5cb4c1f7bc01c700200be3841e14c7c4888de'='Both anchors are present in the unimputed input. Missing anchor values remain NA.','0ee3600723bcae97ca70c0c5e87c7208bf54821c2f91a26db65ff89b2839fb98'='Discordance includes positive and negative calls across batch appearances of the same specimen within a split. Differences between fitted split models are excluded.','e1d2a32cf8db5f1a6b1c5a4aa6ece9a4789faaab198e245f3c6236c47f769721'='The source name and stage 2 support CRC. A missing disease field does not imply a healthy label.','c709419bc8968dc93b1c5cc2b99122ad1f31c699de9cb17976447b105618b1fe'='Primary annotation fields conflict. Source-label sensitivity and exclusion analyses address this uncertainty.','789eb71414f4d504f407bac1494025fb7e7c01d4afa7c8cf37df19324cd8c849'='Intervals use DeLong 95% Wald estimates. Score ties receive half credit. Limits are clipped to [0,1].','44557162e5a7d82688bf23b8a26981c4d2527c1ed00203bf102772c2c449d589'='Direction follows the observed case-control association. This orientation describes the same data and does not estimate externally selected diagnostic performance.','a4a8faa788d8b815424410a22c208f9ce7145241fde72a938c700d4c591b21b5'='Part A contrasts and single-sample versus cohort-normalized B1 are descriptive without a directional expectation.','da1f06769e0dee77a6bee3f66e670da64005fe431c45ce130e8f2fd3dbc4d86b'='The default registered fit path does not use split-based tuning. Training-only fitting and selection use the full outer-training data.','77791a264e89dc188739b282e7358718158042c1ee139383a985c88d22c79e5b'='Non-evaluable method-target-seed cells include structural exclusions and uncoded outcomes.','fd5787e000ad0c243a802285c18e5322a52324c54f56007c26fb7dae8b504b2c'='Target estimates are arithmetic means across all five eligible seeds. Medians give each target equal weight.','8af4831ccff749b3845d9ba6969ad738bd6821a679c479d939c855528802e33e'='Cell counts include all refusals. Structurally eligible targets have at least two source cohorts.','9db59e8968c4bea19da4e451b823022cb28f5b36983d17f6229a8d0386bdf889'='The model has a constant intercept and AUC 0.5.','1915b6381bdecd6c3d719fb6b4994f7cdea2f30e06a64f09509e7e2f679b275a'='The model is not fitted. No intercept-only model exists.','83333466dae6e9a7bd1e25dc0b425de8f79807f69e67ee80c95536e53a5d861b'='Held-out profiles belong to different people from the training profiles.','f3afa242cbf209440a87484c211e14f8c596e07ba064292fab41fbd3da66b712'='Both scoring paths use the same held-out split.','957bd9e2eaf90a60c57adec35e3c15689a39b0d94cfe24cc759142afba37b5dd'='Both scoring paths use the supplemental held-out split.','8022362ff2a6aa8888b4052030076202f78713f87ecb74a4f97db14b823fec88'='Commands execute sequentially without classifier refits.','8789beb9884807b59ff46002ceb7f0ec4c6566b3b833deadb005813b97e26ae4'='The transfer figure producer supplies the reference values.','a24dd159e0e423249b1b79ccdbe10fc0edacb4b858f0fbd1dcdd3f2366fa2635'='The score sums strict cancer inequalities multiplied by the published normalized weights.','d153bd31298c91caf840ee8fd85176c30c13a32a2b275005b683836a45963f8b'='Higher scores indicate cancer. The direction is fixed without optimization.','3be548e09feb3479b127a9e80f348db523ee56c0edc745136a3bd129d534b9c8'='Ten public accessions form one sample-deduplicated analysis unit. Pooled inference does not treat the accessions as independent replicates.','45605e7782b735b133139d52606380d92957fa83f59eec075ff73fa70b096040'='The transfer runner fits and scores through the public deployment interface. These results do not include an independent second-path numerical comparison.','e0fb063eea030760f281f0dd8f5c2137cf4262b1074a161157b8103dd821b001'='The analysis uses real profiles and a stand-in covariate. The searched metadata lack specimen-level RIN.','0a14ca1fd14116ba91b7155b16529196c3c07e8c03763e014c73b249bd7fd039'='The registration receipt has no separate commit field. The accompanying positive-control receipt identifies the release.','2559eafde521269666a41108ef50e80c2c781beefc878610499038ce1eee2d06'='Shared preparation uses the stated analysis specification.','4c282e3317cb9d498c230d4e856a1459b83e17991eca320fe66b3e11a00f9b3f'='Test-batch operations follow the registered designs. The preparation specification defines training and test transformations.','fb406d9e5b89f36223c297c9b05d11f4e7ccad52a924e42f32e3af3c4b29bd1b'='Transfer summaries use the specified grid.','04e38393a43c424f5c693b4445e54a1c72d375a9b390700be0a29edf71e4898e'='The method roster excludes sel-stablemate. The transfer grid contains 2070 cells and 414 method-target summaries.','11b4689338fe5e7ec0a35680783ee250de3f48854067b9cef366ead849b161f5'='Four decision-curve routes use training-only preparation.','dbd712c3cbf8eee6207c562e048014538171418969ba29ab1983e3ca8d0e0bbb'='Merged transfer predictions lack paired training scores for calibrators and exclude four decision-curve routes. The decision-curve partitions include GSE204951 JM08-JM10 sources. Transfer partitions exclude them. Test-batch designs use frozen fitted states.','51332f2ea76a8dcd4c74cf742910886c09180953625e3563c2243fcb100e89fb'='PLCO evaluation uses frozen scores, models and cutoffs.','5b2f0f923ac730e8802d9b6c433ba180dd14f0f3867302b9ff6baca49e9c174f'='The diagnostic training cohort and full preparation are identical across the compared calculations. The target contains no non-finite cells. Its preparation reference is CODE_REPLAY, without an independent stored abundance matrix.','97ccfd7601e78ea95b099b1cb6d76e5b79794c86a2ba5faa174853e3456a9b05'='Within-cohort decision-curve estimates are unavailable for the affected units.','579b7c65b88ae93dd4dfb46df44e1b4be20be445c98d1a6df046efdbc01e1a58'='Within-cohort decision-curve estimates are missing for affected units. Unaffected estimates are available.','8bc47d6081b7e574bff3b83a102a1b05c33e940c9e599724f116b4ef7a61ee0a'='The paired baseline contains non-finite values. The exclusion reason is recorded.','dbb836a5ec02d95e11f5f182edf70ad3c9364223b967c91ce1fcf77b1356ae72'='Selected binary chunk predictions include GSE204951.','d55b02462c643dcf2c982d04984acc24d3242db10814a4cc4cd472cecb9bd144'='Both caches and the producer entry contain 126 features and lack named alignment-counter columns. The 1884-feature reference is unavailable.','8d07bec773a4dfe10efc878c74690ab617f0c300118bf6f8cf5628aba0006ccc'='The exact reference path is unavailable. The 126-feature matrix does not substitute for the required DCA reference.','29b7eda566e311c0f2becce1c9b288c7469bfc62ff2b9c53137d30dcc0f922fb'='The analysis uses 100 pre-surgery profiles and 20 people in both roles.','287c8efcd16127bc7aeebd0461d366c7f45ec05604becc2a4c3c7be9fa015aac'='Decimal TSV serialization introduces rounding differences. Matrix labels require feature-identity alignment.','7208012e7915fe157a7dd278f9663ce49b071c8464c16ad9d687f3a7b21edc32'='Full deposited cohort rows permuted, benchmark preprocessing, held-out rows scored together in permuted order.','ffaff7000cd5f6dd8416052004196837552353215a1c0d2487fdc32489868e32'='Benchmark preprocessing applied to the held-out rows only (a new batch), scored together.','5344664127eabc0cfefbd84685f9b9022bd7a49109bbab6dcb6aca7a65ffed0e'='Each held-out specimen alone through the benchmark preprocessing code, scored alone.','f545119890ac3afbf29d6673409213b13d3095f9ea86ad4f2ca3340cec5ebe91'='Each held-out specimen alone through per-specimen preprocessing with the frozen feature universe and training-only fill constants, scored alone.','0dd9ecd200abe81fb79630a14d5adfcc7fa26890110986c2ee564a341bf31f93'='Control. As S2 but fill constants replayed from the full deposited cohort (the benchmark state).','faf28ce4b7c2341344e871025a86773af8e9c5d5d09f10829b6658f14c75c212'='Held-out specimen plus four extreme companions through the benchmark preprocessing, whole batch scored.','475da2ee274a6d2f3eb5487dea1d83f16ad840b3cedbb03842da3a6e339296c2'='Held-out specimen plus four extreme companions through the frozen per-specimen preprocessing, whole batch scored.','b2d769593704d494fde1269f3af64d6734ab94e1025d89495c8fef0083720976'='Exploratory analysis dated 2026-10-01. The analysis is result-informed and exploratory. It is outside the primary benchmark and FDR family. 48 PLCO_Pre-PDAC/48 PLCO_Control from GEO disease field (label_audit_v2.tsv). Unpaired. No lead-time strata. Model and direction frozen on GSE259327_DIAGNOSTIC n=203 (plco_locked_models.rds). DeLong 95% CI. Delta = frozen 5-seed grouped-CV diagnostic AUC 0.9245 minus PLCO AUC (approx. 95% CI 0.267 to 0.512). Diagnostic lift vs ws-rclr-trimmed +0.0000 (below +0.05).','2b81172106f1bba57305cb0be3be5fe267104daae4022fc5591882aaaba1f900'='Exploratory analysis dated 2026-10-01. The analysis is result-informed and exploratory. It is outside the primary benchmark and FDR family. 48 PLCO_Pre-PDAC/48 PLCO_Control from GEO disease field (label_audit_v2.tsv). Unpaired. No lead-time strata. Model and direction frozen on GSE259327_DIAGNOSTIC n=203 (plco_locked_models.rds). DeLong 95% CI. Delta = frozen 5-seed grouped-CV diagnostic AUC 0.9403 minus PLCO AUC (approx. 95% CI 0.282 to 0.482). Diagnostic lift vs ws-rclr-trimmed +0.0158 (below +0.05).','a669a0ef8818181aae624388324893b0880b93b395beccc794a05bbb6ccab4c6'='Exploratory analysis dated 2026-10-01. The analysis is result-informed and exploratory. It is outside the primary benchmark and FDR family. 48 PLCO_Pre-PDAC/48 PLCO_Control from GEO disease field (label_audit_v2.tsv). Unpaired. No lead-time strata. Model and direction frozen on GSE259327_DIAGNOSTIC n=203 (plco_locked_models.rds). DeLong 95% CI. Delta = frozen 5-seed grouped-CV diagnostic AUC 0.9605 minus PLCO AUC (approx. 95% CI 0.302 to 0.540). Diagnostic lift vs ws-rclr-trimmed +0.0359 (below +0.05).','110cbc2766223ee8ce2153adc731d8a5a6195cf059aa8a83d269b1ec1a4f9ce9'='Near-identical profile deposited under another accession (residual duplicate after identifier and exact-profile deduplication).','507462d3f5cfe4eb884dd0de5dc5099a5c58294a510da0c08ad0b0f83e51c632'='Exact duplicate profile under two sample accessions.','d7075de93e0bc6a01b2f75df27ab5ef180c0be74af49ad4174f0d028a426a9db'='Unimputed specimen-local linear conversion. No cohort statistic. Log(x+0.5) anchor difference.','da98ae58307038efe5a61b2b9c880e797e7b0e804d6ed4a5e5e419bdac8edd65'='Deposited log2 converted to 2^x. No filling. Log(x+0.5) anchor difference.','5a30b61e371a073f89c6f578dc2071ffa4e5192462a537a6b60508f92eabbdc3'='Maximum score difference <=1e-12. C1=0. C2=0.','47904dfbc6f155e82b727ba6c52d16c572ffffcd961f7762d24a82eb8f4486be'='Median C1 at size 80 >= 0.5 times median C1 at size 20. Available datasets at each size.','7d853487dfb27f0be399b986bbd96602a53319787f6326637091b4b7152dd279'='Median within-dataset change (90% minus 10% cases). Sensitivity <0 and specificity >0.','0da5aae5b83ac9bc06e61f25153ab98fddeb257646937a3b62c49b4b18ed2759'='Reference method expressing each miRNA relative to the specimen\'s trimmed typical miRNA level before classification.','4c28989c08368c6ccbf13ee2753eac4182abd2bb2d6416947e781325b826d099'='Expresses selected miRNAs relative to a training-chosen pivot miRNA and classifies the resulting log-ratios.','3b3512acd2a4aa14685c914200766e143135d0b743ebebab7e0d42c7edefeb90'='Classifies from within-specimen log-ratios scaled by training-derived median absolute deviation.','299f9cd40dc80e4678fe8815e4ff935d87fc1e86bd33f388a20eba1ef4a0a9da'='Scores whether selected miRNAs dominate others within the same specimen.','5c5c6a37b98c1ef121a8a1b34ae28b4b8323977685ff53ce38578be901d2a9db'='Removes unwanted variation estimated only from training controls, then applies a frozen classifier.','1fe33bc7cd24a17bf46ecf997f1b36accfd3939d9ffcc70eda2c00706cc86ba6'='Classifies after removing a robust principal-component representation of background variation.','a2fb48b1dff05452920fbb836428509b86598b43f705af12a78104bb1c3eb978'='Removes variation associated with a fitted hemolysis score before classification.','eda0bc9d72351d3f3fc93ab2502f25005f159411dc54ba8182a31c671af009a9'='Fits robust log-ratio normalization separately by assay kit and freezes the kit-specific parameters.','bc88d214a36b424799f9d72ebe015fe7fcd55a677f0cae86c063d7ef73456f9e'='Whole-source constrained stratification. Target 60/40. Seed 101. Exhaustive minimum imbalance.','9c146e7b1f2422310c32bb043fb0866a02c691314c5d1b01947a835eb6569fed'='Removes training-estimated assay-kit effects from additive log-ratios before classification.','1d87aaebe863bb33c1c61d69cb8427268d2f83de981cbd919c554499a8eebecf'='Constructs log-ratio coordinates designed to be orthogonal to assay-kit effects learned in training.','2cac8e6953e6028acedf676e0e8ae0733e63de387a86d7beea96298172ab5c91'='Removes assay-kit effects and classifies from robust residual log-ratio dispersion.','85f13dbd598cca62bb16a83d7bb05190560cddb5ac64e015a968d4b52b0f9ac0'='Normalizes to miRNA anchors selected for stability across assay kits before classification.','5f5461a928ce7694764e85c5d9c9db8b8116a012f6047a3d223d81b15b630ab2'='Anchor status=ok. N_anchors=8.','c175f6da942b24539bc5867c5a0d7d3e8f2d5bc6975104d304c0eeb78923c062'='Normalizes to stable anchors while accounting jointly for hemolysis and assay-kit effects learned in training.','905b1d378e032bc02fafd4b48655d98da96571ae9c53472a67805c3461d7ccb1'='Hemo_ok=FALSE. Anchor_only_hemo_proxy_missing_or_weak.','d24534b3aee1a419fd62bd0ed81e4df79733e3128c09995a2bae58dc26aa6019'='Combines training-selected miRNAs after normalization to rank-invariant anchor miRNAs.','71631d5597c0a04445d5bf40e5da90bf5e186fe10b30b0e03bb3e78eec5b9a3f'='Fits robust log-ratio normalization separately by measurement technology and freezes the stratum parameters.','ab68f2efd456f2670d38b7c82dcfc7efa61ef3fc5ac1d94b50f7a0d46e4b9f85'='Uses training estimates to align robust log-ratio features across measurement technologies before classification.','3485069d877f8459846edfb753c2dc362f27cf2c520b5f2e30cfa94c5932d26b'='Removes training-estimated measurement-technology effects from additive log-ratios before classification.','e13ccac2a91cd50c2e756655fde067ad10b216ed3c6ff7b34d6b0a846ac35cd3'='Fits the robust log-ratio reference separately by biofluid stratum, then freezes the stratum-specific parameters.','22d23dd4d204251af43b4b8afcdbab039e6d1275d1f987589a7c5cb6b8eff622'='Normalizes to training-selected stable miRNA anchors within each biofluid before classification.','6b1d96df1955e00fa1de7de7d312ad81c9865227eae63597a763246a7ba07539'='Removes training-estimated biofluid effects from additive log-ratios before classification.','e0f283aa4eb1fdda15b7b68a7b73908f832f4b6bc315364a8a4f9ea45a29152d'='Excludes training-specified blood-cell marker miRNAs before robust log-ratio normalization and classification.','eb9ca23fb13d45adf39a916b699abd3a156bc5c59b2109a38d2b5f98f5bf6ab9'='Learns features that retain cancer/control information while suppressing training-cohort information.','2a7deb36667827a0053950224d83f7d20a66444cfddf383625a2bd8231a109ee'='Fits separate class-conditional generative models and scores whether a specimen is better explained as cancer or control.','23fdcfda17a654b43a382f012c31e391341d334d5fb2c43514bbae1db4817054'='Selects predictors whose relationships remain stable across training environments.','45012cd37072ab02f9b00b4c50eba032f8a68e583e453afc60326d28e740cb15'='Voting classifier based on which miRNA is higher within each selected pair.','af9fab1d397aac03e30092baecc20daa646adbb9ad8abd48295c6a033613fe6b'='Penalized classifier using selected log-ratios between miRNA pairs.','b5fd6fb7882f8b46edb99c79d7a96d78b3dcbe954ffa87ea5e630c72da27fcb9'='N_selected_pairs=14. Lambda_rule=1se.','10da3e2ad8dc4a414425b50c75cb427999eba6820538880c521b22ebe94dca41'='Scores enrichment of training-selected up- and down-miRNA sets within each specimen\'s ranks.','96d02ec466dfcd58dba9ee6b37a20225f751304e80553e6cdb4f1f746bfcadba'='Uses a Mann--Whitney rank statistic to score enrichment of a training-selected miRNA set within one specimen.','682f0dc382f9b8fcd21d885f34372dd6545e910336766eaf6a50ac42c131758a'='Combines selected miRNA-pair ordering rules across training cohorts.','1e51e29c6c331dc222717fc2415029e91a6b07949fa525dc099c636279b594a2'='Classifies from a selected balance between geometric means of two miRNA sets.','365f3094aa38b842895004840abc9a097c0b29dfbca2da0afa9b48a2de08017e'='Classifies from a training-selected isometric log-ratio balance between two miRNA groups.','498a219cab73718183f5ff2740fef08cf52e2b3d7d33aca5f46fd867db0444ce'='Compares fitted cancer and control distributions with a Bures--Wasserstein likelihood-ratio score.','4c205609daea0c7194cdef21fb1fd355f4b2e4d86312db32ae09c3e7cc0008d1'='Multiplies feature-wise cancer-versus-control density evidence estimated by kernel smoothing.','66a3118da5b094796e5bd048107ac19b6e8db22295ce0b644782adc6e2a1bd76'='Compares fitted cancer and control Gaussian-copula distributions with a likelihood ratio.','bb0ce340376b929161288bd3ba6c9790a0034fffed199d6da75e278bd8fb93e1'='Compares fitted heavy-tailed cancer and control copula distributions with a likelihood ratio.','ccf183b668c185751dd48c30729e6add43f025a0c631baace90cb22bb8a44cff'='Compares class-specific conformity of a specimen and returns a case-versus-control evidence score.','105e0f8bdedd09ff77d946c2b190e109ef8fcf01eb6a9dc60adc26702df6f39d'='Learns a specimen representation by distinguishing genuine training profiles from corrupted copies, then applies a frozen classifier.','fbfbeaa8fe1df296619b033b8beb3ace9b2270079498fa56b7c422680909485b'='Uses a pretrained tabular model conditioned on frozen training examples to score each new specimen.','649ac5be41723382fb013665436177ed4a4b6450866f447df3491d74a770f1e7'='Fits a classifier while emphasizing the worst-performing recorded training group.','742e0dfc682ab5b7ba5fd88a79adb0e2c412ac434217b17d40c47f8da53a7156'='Fits a classifier while penalizing variation in prediction risk across training groups.','2a3f8cd2a1a371dacc4c6051de3b0f96bed535d56e4281f096bdd180e9603ff3'='Combines several frozen expert classifiers using weights determined only from the incoming specimen.','e6e4859c04fb59f20b5afa217219c8cbc5b2a4796991ef5244c24c2c68f44c15'='Uses multi-scale wavelet summaries of the within-specimen profile as inputs to a frozen classifier.','cb4886fdcbc837c591ad21f2e75adb75a5ceeef08a8a7aa9ddde4abc14431a26'='Summarizes local ordinal patterns in the sorted within-specimen abundance profile before classification.','3d3090659ff1987d04d8de5e36799b5dea0ef0ad0933d5ee49814f0e1fed5558'='Summarizes persistent peaks and valleys across thresholds of the within-specimen profile before classification.','11ba8911b0e911517bf8c978b3a5b7ce51d7bc62028daf1967fd64dd9a225b2e'='Summarizes ordered changes along the fixed within-specimen miRNA profile with path-signature features.','905cbfea9358cc321719e54d134d03e8c0eb2dfe29eced26e385d4c7f38e7898'='Summarizes multi-scale fluctuation patterns along the fixed-order within-specimen log-ratio profile before classification.','1f669613c18cbfd6342863cc16daa86341735d0f38c6fe3135ee22105f6a9266'='Encodes pairwise within-specimen relations as an image and applies a frozen convolutional classifier.','fcb3ca6da16f1b92afac2d82432abb54ac7b928f832f99ab54a6ab9c46b5fc1b'='Maps each specimen to transport displacements from a training-derived reference distribution before classification.','1aeea3a98bfb56f6eeac240e3b9b616dc70e9d72e23344527b0895d791570ee9'='Combines cancer-versus-control distribution evidence across one-dimensional transport projections.','3f8481191d47eb84707dc86c816e9fd8eeda8c3daf90ddb8998351bb44f4f4dd'='Scores whether a specimen is closer to the training cancer or control distribution in a kernel representation.','51d3abd3908390198a32e694e5222f3c1ba59031a82a349efc68a8750a8334d2'='Scores relative proximity to fitted cancer and control distributions using Fisher--Rao geometry.','bf8448cc57f49a0c0dfc5eec280abd37e46c68dc24d8277793198671dc880d3c'='Directly estimates the cancer-to-control density ratio and uses it as the specimen score.','2a6890c3cd2899a3afca028270fbefc8e69c0e566a022cfb6ff6a44d0d8e0250'='Projects each specimen onto a training-derived miRNA co-expression graph and classifies from graph-frequency features.','eb21252bf99b2d0ab268df76986a6de911aedf1017b3c3700f69ff6bc877bd5e'='Places a new specimen into a training-derived diffusion map and classifies from its fixed coordinates.','b4fe08ac68a257ed0fdf82a2fc20c444b2fce051075e0fad45b500b617d3dc21'='Uses a regularized neural classifier with a Gaussian-process output to produce a frozen specimen score.','ccd20c28972a6631b9abed16c6654c52e3485cc2a69127c58acddbd0ef065e16'='Learns a self-supervised specimen representation by preserving information while reducing redundancy, then freezes a classifier.','d131e19b2ff9a08abd308aa62cef704c71779bde3753c07da1049353445e4bec'='Learns a nonlinear classifier constrained to use log-ratios between miRNAs.','6aa9b0ab4aeca7a96a1a9e250446f5ddcbe794ce73dda9514728bf60e7da5864'='Uses a pretrained in-context tabular model with frozen training examples to score each specimen.','06c111576ffb8da0eb1266ba523b8974f74f2a5c4e41da92c0ed27d16eefb2c5'='Uses a pretrained retrieval model with a frozen bank of training examples to score each specimen.','c1abc88d153fc69d125770e4c0a47781b091f570f85f0aadcd8eef0a4384a3bc'='Summarizes texture-like co-occurrence patterns in the fixed-order within-specimen log-ratio profile.','4421c19a06073ddcf2a24367084f11281d55b5b48b021dd0b77fc5b2b296e4d0'='Represents each specimen by its abundance quantile function and classifies from training-derived principal components.','fec1336e95cad159503a1a58ebf1f3dc08189e4385df2838d19a460cbfea32fb'='Describes curvature across smoothing scales of the fixed-order within-specimen profile before classification.','1bfaee03e5c939680d97a1d310dc36ce9849b5dc73bba1528fe48adba065d37c'='Compares fitted cancer and control vine-copula distributions with a likelihood ratio.','f169b420f6889d016622bb960649bab1898e439a0ec585fc37f445493b976992'='Scores relative proximity to frozen cancer and control prototype profiles in a learned representation.','df5034b0901714dc2eb7cca1a8b9a2ea99f28940fcf3be2d82ada6061d1c02a2'='Learns a frozen representation and scores relative Mahalanobis distance to the cancer and control classes.','14b24520ebde72ac52bb52fcb12503d769b85ed5651b8b6042a4022254f9377a'='Learns a sparse ratio between two miRNA groups and scores the resulting fixed balance.','2f4fde319ad7b6bbf209c7dab47f2b48523e10c4e586e797dffb4d4c8c736948'='Fits a classifier while penalizing differences in gradient variability across training groups.','cfe803209c70a55896c0fcd70dba10f79b5e5e0866e87364fe8c4397e7d11c08'='Learns a compact representation whose prediction rule is constrained to remain stable across training groups.','90bfa0db960f32479574d1e156c8e4102c9605306f425d45478cfa6394f0f342'='Selects predictors that remain stable across training groups, then freezes the resulting classifier.','732f14b49c6738dbf6b432096190a4f4745c01f0091ebb223acea3f51f303d64'='Negative control based on robust distance from the centre of the training distribution.','81b59b099df980062fbc8428d120a4e4da282899e4731928f356c1ccce703d64'='Negative control that scores how atypical a specimen is relative to a frozen one-class training distribution.','d3989d0d04761050268e5d33179623e2197bb4697a2c5bd8fa67065f89979c0f'='Negative control that scores how readily a specimen is isolated from the training distribution.','8532fdc644ba0fcc6f95025f76a32bbfbdbfc5d659e275cf5e72035385f74953'='Negative control that combines feature-wise tail probabilities into a global outlier score.','cb166d78f12736415fe062f20da6a8a36e3564ccb54856660be0d34b953c99fa'='Negative control based on transport distance from a single specimen to the frozen training reference.','349affd56b64fed6c149fc858c1ec3bf91361d6f1f5983e2fc2ef87e92a4b9c6'='Deliberately uses companion specimens to centre scores.','a118b9acdc504d22333fe24cc56feb2a6c32a2d9db645a8971f41b2d7d2ea2e9'='Attempts audited registration of the deliberately batch-dependent scorer.','0de95dd0b31f4eb63a73d07e90365e0f7bcc77639622eb285728c3e0001d56b0'='Requires metadata, a frozen training panel and panel weights, which a matrix-only call omits.','6cf29e2857af2aa3445b634bc84b301c218ba090f9b63d2bdcbe6c6a07e72bf3'='Stops because metadata are missing.','f4e4f158bcc73059bce61c2db3a382519c2897aec02456f6f39e471376996783'='Returns one score per specimen and passes the synthetic and real-profile audits.','1b0f2bbe8018f96231fb5d4472808ecc685874ff1f691fdea2ef85aeae025712'='Returns a denoised matrix and takes specimen values before the fitted object.','a150560aaea1befd56c3d8e2dfbed8ad80a817446402f8c65799c045b8dc1fb4'='Stops at the scorer-registration guard before native scoring.','08091e76f7401013c353740d2db61a52ba270061477c5e0828570520c77e8852'='Returns a list containing scores, empirical P values, thresholds and flags.','cd4f64e3191e852747f87bfd265cd0359027d57b193b23022357c57ee011e8a7'='Returns a list containing squared distances, distances and chi-square P values.','661a3a5d9286fa472f972ffcf1a6dbd0bc129656b993961218f6e37dfb95a635'='Returns a matrix of residuals rather than one score per specimen.','366cd4c257e44b03cf9d4fe69f9da9f7b51f1dfc8b3f77576a00c47cb7fd05df'='Returns an ALR matrix and requires sufficient pivot features plus a frozen selected panel.','00319b7e2942b9b1a80904ab19025ae93e50ebc54a00db1a2ed92c323c5b6ef0'='Stops because too few default pivot features are present.','32a55cb0f32e2800c5e07b4143901cee5c4b02336a0a099ab7db9bb9f13ba446'='Omits the frozen fitted top-feature count and pseudocount when called without an adapter.','057f3bad901c4a914582bd2904ddc695a469cb1633a5df34a0eb2ee7228a0fb9'='Returns a within-specimen MAD-scaled log-ratio matrix and omits the frozen selected panel.','d4e36562b78b7df5c793a48910463733e47c136e917d6589ac5ad61e1f8a6e3a'='Shape validation rejects the matrix returned for the training probe.','837a6139e2d9089156a39e4490d7a5e3efc14238d125fea26bcbf4616a9d441c'='Fitted intercept retained.','2ae68edc5d86f948f1ee4697bcb8e9ea7fca6a4aea4375c6b43894b1c61236a4'='Singleton/batch score identity, zero discordance, and singleton AUC over identical appearances.','34d71b1bf93db120ba82c5f57c4e804563116189157813d6606c6c635aef24f0'='Positive discordance expectation applies only to cohort-normalized arms with selected features. Intercept-only arms remain valid fitted predictors.','2115e628bd0d0439c714a3ed1631978aa7c6a881bbcce1e979a2128383ffb8c0'='Equal-dataset mean deployment AUC minus full-set AUC.','54eb0cfeb556779e6c977853c097f6743162e6423bd50b1bd9de109606f40862'='Survivor cohort. Outside cancer-discrimination synthesis.','c23ddee30406b84a293eee0abda73ca990381734360afe841a8d217ec3e1e215'='Donor/pool membership unresolved. Independent-group threshold unverifiable.','9a96f936e2d621a1c13a6a9b612d0658ad1fba0771c2a03203fe009c04cd84cd'='Imaging screen status without pathological confirmation.','9172a48509d156a0be17268f786609d866126e9b2f7d4f95dfa91e81277a57b8'='Whole blood. Outside primary biofluid synthesis.','8ae5b6dc5b69abb9ef400a4f8aa8e8bf768dfc19ed2bb51854af2ce504bac628'='Cross-validation selects glmnet penalty lambda. .reo_pairratio_foldid -> .ss_grouped_inner_folds (minimum 3 folds).','ba7889f7895d17c14a59f6a11ed106e5553578b938537169ad32b754a695351e'='Internal cross-validation selects ridge lambda when hp$lambda is NULL (default). .inv_scatter_cv_lambda -> .inv_scatter_cv_folds -> .ss_grouped_inner_folds.','27b388145b1d44b3dd79029e04beabf30e8f57fcaafdc3abaa8a622881e8398f'='Internal cross-validation selects ridge lambda when hp$lambda is NULL (default). .sig_path_cv_lambda -> .sig_path_cv_folds -> .ss_grouped_inner_folds.','6d5aa630cf05aca5cc37684eff2a8f28e8fa5838c14cd5752110687f7e4bf39e'='Episodic representation training separates support and query sets. .proto_net_episode_plan -> .ss_grouped_episode_plan -> .ss_grouped_inner_folds. This is internal resampling, not a hyperparameter grid.','44adc6c8cd6713d77c2e57b3ad36f08f4c1b67e67dcb98a409e1c8960027dd2d'='For each specimen observed at both extremes. P10*(1-p90)+(1-p10)*p90. Equal specimen mean within split, then equal split mean.','fed072798de2178e9fd60bbcda3b1d92a325993a861a7825537f0931b6c6f59c'='95% DeLong.','0a6f89751620e0b2eaea41a1a2453ef90dbea71b186dbb93f5cac52f3925f5c1'='2000 joint person-bootstrap percentile draws. Saved scores/models/batches fixed.','95fd6982763f2bc3aa267cbe28c5319efc39226768192fd39e422e72064602bc'='Joint person bootstrap (2000 draws). Pooled REML/Knapp-Hartung. Scores and batches fixed.','d06dc36e25bc2dea75e45bc918f81cfa6abda034f157acfa761205233b8269a8'='Released code. Larger of cancer and non_cancer weighted scores. Cancer first on score ties.','8ea2d3f0bdc7b921ef62eb4bc227fc8e8eca136288e7cd067dd67bb5d65428e1'='N_scenario_comparisons=6090. Max_abs_change=0.','fb086b4300d3e2fcc76e44964d3f30284b6bc92b98fe4077081181abc43e6673'='N_examples=6. Max_abs_difference=6.0715321659188248e-18.')
.public_table_values <- c('controls_cpu_final'='cpu_four_thread_controls','repair_cpu_1_8060eb0ba3'='cpu_four_thread_controls_128gib','gate_gpu_s8_tm1c_matched_v3'='cuda_seed_101_replay','gate_gpu_s8_tm1c_matched_v2'='cuda_seed_101_replay','sensitivity_gpu_s8_tm1b'='cuda_cv_seed_sensitivity','tm1c_seed_101_s8'='cuda_matched_seed_101','cpd-konrad-worker'='gpu_server_01','zbimt-s8'='gpu_workstation_01','tm1c_seed_202_mgb_v3'='cuda_seed_202_device_7b1d7eaf','tm1c_seed_303_mgb_v3'='cuda_seed_303_device_f6af6362','tm1c_seed_404_mgb'='cuda_seed_404_device_eefe4b23','tm1c_seed_404_mgb_v2'='cuda_seed_404_device_f6af6362','tm1c_seed_505_mgb'='cuda_seed_505_device_f6af6362','tm1c_seed_505_mgb_v2'='cuda_seed_505_device_f6af6362','D2_corrected_labels'='count_matched_labels','D2_corrected_outer_training_refit'='count_matched_outer_training_refit','R15b_outer_training_refit_historical_bank_missing'='outer_training_refit_without_score_bank','cpu-release'='cpu_measurements','s8-release'='cuda_measurements','cpu-R-dependencies'='cpu_dependency_measurements','cpu-hemo'='cpu_hemolysis_measurements','transfer-2610'='transfer_measurements','cpu-canonical-2610'='cpu_foundation_measurements','mgb-gpu-2610'='cuda_foundation_measurements','cpu-refusals-pb1c'='cpu_refusal_measurements','s8-gpu-pb1c'='cuda_scoring_measurements','s8-host-cache-pb1c'='cuda_host_cache_measurements','gpu-retention'='gpu_model_storage','source-loaded-s8-full'='source_loaded_cuda_measurements','R15_outer_training_refit_historical_bank_missing'='reference_outer_training_refit_without_score_bank')
.public_path_pattern <- '(?i)(?:^|_)(?:path|file|filename|directory|root|script|entry|archive|source_bundle|predictions|output)(?:_|$)'
.public_join_pattern <- '(?i)(?:^|_)(?:sample|specimen|patient|person|donor|subject|row|group|cohort|dataset|unit|fold|seed|method|model|arm|contrast|analysis|family|table|metric|endpoint|setting|tool|scorer|feature|id|key|label|task)(?:_|$)'
.public_path_fields <- c('execution_source','evidence','artifact','Source path','Producer summary source','bundle','primary_bundle','sensitivity_bundle','control_bundle','analysis_primary_bundle','analysis_sensitivity_bundle','GPU_measurement_adapter','source','target','output_member','script_member','manuscript_file','matrix_path','checkpoint','reference_cache','source_bundle','model_path','prediction_path','expected_path','input_path','output_path','predictions')
.public_table_labels <- function(x) {
 if(!is.data.frame(x))return(x)
 x <- data.table::copy(x)
 paths <- unlist(jsonlite::fromJSON(Sys.getenv('OMICSELECTOR_PUBLIC_PATHS','{}')),use.names=TRUE)
 for(column in names(x)) {
  value <- x[[column]]
  if(!is.character(value) && !is.factor(value))next
  value <- as.character(value)
  index <- match(value,names(.public_table_values));keep <- !is.na(index)
  value[keep] <- unname(.public_table_values[index[keep]])
  if(column %in% .public_prose_columns) {
   index <- match(vapply(value,function(v)if(is.na(v))NA_character_ else digest::digest(enc2utf8(v),algo='sha256',serialize=FALSE),character(1)),names(.public_notes));keep <- !is.na(index)
   value[keep] <- unname(.public_notes[index[keep]])
  }
  path <- !column %in% c('source_host','baseline') && (column %in% .public_path_fields || grepl(.public_path_pattern,column,perl=TRUE) && !grepl(.public_join_pattern,column,perl=TRUE))
  if(path && length(paths))value <- vapply(value,function(v) {
   if(is.na(v) || grepl('^[A-Za-z][A-Za-z0-9+.-]*://',v))return(v)
   if(!grepl('/',v,fixed=TRUE) && !grepl('[.](tsv|csv|json|rds|R|py|gz|txt|pkl|tar|xlsx|xls)$',v,ignore.case=TRUE))return(v)
   v <- sub('^/(.*/)?OmicSelector_paper/','./',v)
   root <- Sys.getenv('OMICSELECTOR_PAPER_ROOT','')
   if(nzchar(root) && startsWith(v,paste0(root,'/')))v <- paste0('./',substring(v,nchar(root)+2))
   for(prefix in c('delegated_tasks/review_response_20261020/','delegated_tasks/complete_all_routes_20261010/'))v <- gsub(prefix,'analyses/',v,fixed=TRUE)
   v <- gsub('delegated_tasks/single_sample_reframe_20261017/','analyses/single_specimen/',v,fixed=TRUE)
   trailing <- endsWith(v,'/');parts <- strsplit(v,'/',fixed=TRUE)[[1]];i <- match(parts,names(paths));k <- !is.na(i)
   parts[k] <- unname(paths[i[k]]);paste0(paste(parts,collapse='/'),if(trailing)'/'else'')
  },character(1),USE.NAMES=FALSE)
  x[[column]] <- value
 }
 index <- match(names(x),names(.public_schema));keep <- !is.na(index)
 names(x)[keep] <- unname(.public_schema[index[keep]])
 x
}
.public_input_path <- function(file) {
 if(!is.character(file) || length(file)!=1L || is.na(file))return(file)
 root <- Sys.getenv('OMICSELECTOR_PAPER_ROOT','')
 rooted <- grepl('^(?:[.]/)?(?:delegated_tasks|analyses|data_bucket|processed|reference|external)/',file,perl=TRUE) || (nzchar(root) && startsWith(file,paste0(root,'/'))) || grepl('/OmicSelector_paper/',file,fixed=TRUE)
 value <- .public_table_labels(data.frame(source_bundle=file,stringsAsFactors=FALSE))$source_bundle[[1]]
 if(rooted && nzchar(root) && !startsWith(value,'/'))value <- file.path(root,sub('^[.]/','',value))
 value
}
.public_readRDS <- function(file,...)base::readRDS(.public_input_path(file),...)
.public_write_table <- function(x,...)utils::write.table(.public_table_labels(x),...)
.public_fwrite <- function(x,...)data.table::fwrite(.public_table_labels(x),...)

# Exact, fail-closed runner for the OmicSelector paper's separately registered
# external rank-tree competitor sensitivity. Scientific results are produced
# only through the pre-existing grouped-CV benchmark engine after every source,
# cache, environment, package, and provenance input has been hash-bound.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

.ranktree_methods <- function() {
  c("ws-rclr-trimmed", "reo-ktsp", "reo-pairratio", "reo-singscore",
    "reo-ucell", "reo-rankforest", "reo-rankboost")
}

.ranktree_seeds <- function() c(101L, 202L, 303L, 404L, 505L)

.ranktree_arg_value <- function(args, flag) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) NA_character_ else args[[hit + 1L]]
}

.ranktree_runner_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  flags <- c(
    package_root = "--package-root", paper_root = "--paper-root",
    analysis_plan = "--analysis-plan",
    analysis_amendment = "--analysis-amendment",
    analysis_amendment2 = "--analysis-amendment2",
    analysis_amendment3 = "--analysis-amendment3",
    cache = "--cache", output_dir = "--output-dir", methods = "--methods",
    provenance_manifest = "--provenance-manifest",
    provenance_union_inventory = "--provenance-union-inventory",
    runtime_image = "--runtime-image",
    runtime_attestation_manifest = "--runtime-attestation-manifest",
    engine_manifest = "--engine-manifest",
    environment_manifest = "--environment-manifest",
    expected_package_version = "--expected-package-version",
    expected_analysis_plan_sha256 = "--expected-analysis-plan-sha256",
    expected_analysis_amendment_sha256 =
      "--expected-analysis-amendment-sha256",
    expected_analysis_amendment2_sha256 =
      "--expected-analysis-amendment2-sha256",
    expected_analysis_amendment3_sha256 =
      "--expected-analysis-amendment3-sha256",
    expected_package_commit = "--expected-package-commit",
    expected_cache_commit = "--expected-cache-commit",
    expected_cache_sha256 = "--expected-cache-sha256",
    expected_cache_manifest_sha256 = "--expected-cache-manifest-sha256",
    expected_cache_digest_sha256 = "--expected-cache-digest-sha256",
    expected_runtime_image_sha256 = "--expected-runtime-image-sha256",
    expected_runtime_attestation_manifest_sha256 =
      "--expected-runtime-attestation-manifest-sha256",
    expected_container_launcher_sha256 =
      "--expected-container-launcher-sha256",
    expected_r_version = "--expected-r-version",
    expected_engine_manifest_sha256 = "--expected-engine-manifest-sha256",
    expected_environment_manifest_sha256 =
      "--expected-environment-manifest-sha256",
    expected_provenance_script_sha256 =
      "--expected-provenance-script-sha256",
    expected_provenance_manifest_sha256 =
      "--expected-provenance-manifest-sha256",
    expected_provenance_union_inventory_sha256 =
      "--expected-provenance-union-inventory-sha256",
    expected_analysis_code_id = "--expected-analysis-code-id"
  )
  out <- lapply(flags, function(flag) .ranktree_arg_value(args, flag))
  missing <- names(out)[vapply(out, function(x) is.na(x) || !nzchar(x), logical(1L))]
  if (length(missing)) {
    stop("Missing required argument(s): ",
         paste(unname(flags[missing]), collapse = ", "), call. = FALSE)
  }
  hex <- c(
    expected_analysis_plan_sha256 = 64L,
    expected_analysis_amendment_sha256 = 64L,
    expected_analysis_amendment2_sha256 = 64L,
    expected_analysis_amendment3_sha256 = 64L,
    expected_package_commit = 40L, expected_cache_commit = 40L,
    expected_cache_sha256 = 64L, expected_cache_manifest_sha256 = 64L,
    expected_cache_digest_sha256 = 64L,
    expected_runtime_image_sha256 = 64L,
    expected_runtime_attestation_manifest_sha256 = 64L,
    expected_container_launcher_sha256 = 64L,
    expected_engine_manifest_sha256 = 64L,
    expected_environment_manifest_sha256 = 64L,
    expected_provenance_script_sha256 = 64L,
    expected_provenance_manifest_sha256 = 64L,
    expected_provenance_union_inventory_sha256 = 64L,
    expected_analysis_code_id = 16L
  )
  for (name in names(hex)) {
    if (!grepl(sprintf("^[0-9a-f]{%d}$", hex[[name]]), out[[name]])) {
      stop(flags[[name]], " must be an exact lowercase hexadecimal pin of ",
           hex[[name]], " characters.", call. = FALSE)
    }
  }
  selected <- strsplit(out$methods, ",", fixed = TRUE)[[1L]]
  if (!length(selected) || any(!nzchar(selected)) || anyDuplicated(selected) ||
      length(setdiff(selected, .ranktree_methods()))) {
    stop("--methods must be a unique comma-separated subset of the seven ",
         "registered sensitivity methods.", call. = FALSE)
  }
  out$methods <- .ranktree_methods()[.ranktree_methods() %in% selected]
  out
}

.ranktree_sha <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

.ranktree_assert_file_pin <- function(path, expected, label) {
  if (!file.exists(path)) stop(label, " does not exist: ", path, call. = FALSE)
  observed <- .ranktree_sha(path)
  if (!identical(observed, expected)) {
    stop(label, " does not match its expected SHA-256 pin.", call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

.ranktree_validate_analysis_plan <- function(path, paper_root,
                                             expected_sha256) {
  paper_root <- normalizePath(paper_root, mustWork = TRUE)
  canonical <- normalizePath(file.path(
    paper_root, "plan", "ranktree_external_competitor_sensitivity.md"
  ), mustWork = TRUE)
  path <- .ranktree_assert_file_pin(path, expected_sha256,
                                    "Approved rank-tree analysis plan")
  if (!identical(path, canonical)) {
    stop("Rank-tree analysis plan must be the exact canonical path under the ",
         "declared paper root.", call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  approved <- grepl(
    "^Status:[[:space:]]+(?:`|\\*\\*)?APPROVED\\b", lines, perl = TRUE
  )
  if (sum(approved) != 1L) {
    stop("Rank-tree analysis plan must contain exactly one explicit ",
         "'Status: APPROVED' marker.", call. = FALSE)
  }
  path
}

.ranktree_validate_analysis_amendment <- function(path, paper_root,
                                                  expected_sha256,
                                                  amendment_number = 1L) {
  amendment_number <- as.integer(amendment_number)
  if (length(amendment_number) != 1L || is.na(amendment_number) ||
      !amendment_number %in% c(1L, 2L, 3L)) {
    stop("Rank-tree amendment number must be 1, 2, or 3.", call. = FALSE)
  }
  paper_root <- normalizePath(paper_root, mustWork = TRUE)
  canonical <- normalizePath(file.path(
    paper_root, "plan", sprintf(
      "ranktree_external_competitor_amendment%d.md",
      amendment_number
    )
  ), mustWork = TRUE)
  path <- .ranktree_assert_file_pin(
    path, expected_sha256, "Approved rank-tree analysis amendment"
  )
  if (!identical(path, canonical)) {
    stop("Rank-tree analysis amendment must be the exact canonical path under ",
         "the declared paper root.", call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  approved <- grepl(
    "^Status:[[:space:]]+(?:`|\\*\\*)?APPROVED\\b", lines, perl = TRUE
  )
  if (sum(approved) != 1L) {
    stop("Rank-tree analysis amendment must contain exactly one explicit ",
         "'Status: APPROVED' marker.", call. = FALSE)
  }
  path
}

.ranktree_git_provenance <- function(package_root) {
  package_root <- normalizePath(package_root, mustWork = TRUE)
  git <- Sys.which("git")
  if (!nzchar(git)) stop("A git executable is required.", call. = FALSE)
  commit <- system2(git, c("-C", package_root, "rev-parse", "HEAD"),
                    stdout = TRUE, stderr = TRUE)
  status <- system2(git, c("-C", package_root, "status", "--porcelain=v1",
                           "--untracked-files=all"),
                    stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(commit, "status")) || length(commit) != 1L ||
      !grepl("^[0-9a-f]{40}$", commit) || !is.null(attr(status, "status"))) {
    stop("Could not resolve exact OmicSelector git provenance.", call. = FALSE)
  }
  list(root = package_root, commit = unname(commit[[1L]]),
       dirty = length(status) > 0L)
}

.ranktree_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
  }
  normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
}

.ranktree_validate_engine_manifest <- function(path, expected_sha256) {
  path <- .ranktree_assert_file_pin(path, expected_sha256, "Engine manifest")
  manifest <- data.table::fread(path)
  required <- c("role", "path", "sha256")
  roles <- c("benchmark_engine", "paired_auc_engine", "fold_engine",
             "public_design")
  if (!all(required %in% names(manifest)) || nrow(manifest) != length(roles) ||
      anyNA(manifest[, ..required]) || anyDuplicated(manifest$role) ||
      !setequal(manifest$role, roles) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    stop("Engine manifest must enumerate exactly the four frozen benchmark ",
         "source roles with complete SHA-256 pins.", call. = FALSE)
  }
  manifest[, path := vapply(seq_len(.N), function(i) {
    .ranktree_assert_file_pin(path[[i]], sha256[[i]],
                              paste0("Benchmark source ", role[[i]]))
  }, character(1L))]
  manifest <- manifest[match(roles, role)]
  manifest
}

.ranktree_installed_db <- function() {
  raw <- utils::installed.packages()
  packages <- unique(as.character(raw[, "Package"]))
  lib_paths <- vapply(unique(raw[, "LibPath"]), normalizePath, character(1L),
                      mustWork = TRUE)
  names(lib_paths) <- unique(raw[, "LibPath"])
  selected <- vapply(packages, function(package) {
    # Restrict resolution to the installed libraries. pkgload::load_all() makes
    # find.package() without lib.loc prefer the clean source checkout, which is
    # intentionally not an installed.packages() record and is pinned separately
    # by Git provenance.
    path <- find.package(package, lib.loc = .libPaths(), quiet = TRUE)
    if (!nzchar(path)) return(NA_integer_)
    library <- normalizePath(dirname(path), mustWork = TRUE)
    matches <- which(raw[, "Package"] == package &
                       unname(lib_paths[raw[, "LibPath"]]) == library)
    if (length(matches) == 1L) matches else NA_integer_
  }, integer(1L))
  if (anyNA(selected)) {
    stop("Could not resolve one active installed-library record per package.",
         call. = FALSE)
  }
  out <- raw[selected, , drop = FALSE]
  rownames(out) <- out[, "Package"]
  out
}

.ranktree_runtime_closure <- function(package = "ranktreeEnsemble") {
  installed <- .ranktree_installed_db()
  if (!package %in% rownames(installed)) {
    stop(package, " is not installed.", call. = FALSE)
  }
  dependencies <- tools::package_dependencies(
    package, db = installed, which = c("Depends", "Imports", "LinkingTo"),
    recursive = TRUE
  )[[package]]
  dependencies <- setdiff(unique(c(package, dependencies)), "R")
  missing <- setdiff(dependencies, rownames(installed))
  if (length(missing)) {
    stop("The installed ranktreeEnsemble dependency closure is incomplete: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  data.table::data.table(
    package = sort(dependencies),
    version = as.character(installed[sort(dependencies), "Version"])
  )
}

.ranktree_validate_environment_manifest <- function(
    path, expected_sha256, closure = .ranktree_runtime_closure()) {
  path <- .ranktree_assert_file_pin(path, expected_sha256,
                                    "Dependency environment manifest")
  manifest <- data.table::fread(path)
  required <- c("package", "version", "source_path", "source_sha256")
  if (!all(required %in% names(manifest)) || anyNA(manifest[, ..required]) ||
      anyDuplicated(manifest$package) ||
      any(!nzchar(manifest$package) | !nzchar(manifest$version) |
            !nzchar(manifest$source_path)) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$source_sha256))) {
    stop("Dependency environment manifest has an invalid schema or pin.",
         call. = FALSE)
  }
  closure <- data.table::as.data.table(closure)
  if (!all(c("package", "version") %in% names(closure)) ||
      anyDuplicated(closure$package) ||
      !setequal(manifest$package, closure$package)) {
    stop("Dependency environment manifest does not enumerate the exact ",
         "recursive ranktreeEnsemble closure.", call. = FALSE)
  }
  observed <- closure$version[match(manifest$package, closure$package)]
  if (any(manifest$version != observed)) {
    stop("Installed dependency versions differ from the frozen environment ",
         "manifest.", call. = FALSE)
  }
  if (!identical(manifest[package == "ranktreeEnsemble", version], "0.24")) {
    stop("The external sensitivity requires ranktreeEnsemble 0.24.", call. = FALSE)
  }
  manifest[, source_path := vapply(
    seq_len(.N), function(i) .ranktree_assert_file_pin(
      source_path[[i]], source_sha256[[i]],
      paste0("Dependency source archive for ", package[[i]])), character(1L)
  )]
  data.table::setorder(manifest, package)
  manifest
}

.ranktree_attestation_closure <- function(
    roots = c("data.table", "digest", "metafor", "pkgload",
              "ranktreeEnsemble")) {
  installed <- .ranktree_installed_db()
  missing <- setdiff(roots, rownames(installed))
  if (length(missing)) {
    stop("Runtime attestation is missing required package(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  dependencies <- unique(unlist(tools::package_dependencies(
    roots, db = installed, which = c("Depends", "Imports", "LinkingTo"),
    recursive = TRUE
  ), use.names = FALSE))
  packages <- sort(setdiff(unique(c(roots, dependencies)), "R"))
  data.table::data.table(
    package = packages,
    version = as.character(installed[packages, "Version"]),
    root_path = vapply(packages, function(package) {
      normalizePath(file.path(installed[package, "LibPath"], package),
                    mustWork = TRUE)
    }, character(1L))
  )
}

.ranktree_attestation_files <- function(root, scope, owner, version) {
  paths <- list.files(root, recursive = TRUE, all.files = TRUE,
                      full.names = TRUE, include.dirs = FALSE,
                      no.. = TRUE)
  info <- file.info(paths)
  paths <- paths[!is.na(info$isdir) & !info$isdir]
  info <- file.info(paths)
  if (!length(paths) || anyNA(info$size)) {
    stop("Runtime attestation found an empty or unreadable root: ", root,
         call. = FALSE)
  }
  root_prefix <- paste0(normalizePath(root, mustWork = TRUE), .Platform$file.sep)
  normalized <- normalizePath(paths, mustWork = TRUE)
  relative <- substring(normalized, nchar(root_prefix) + 1L)
  if (any(!nzchar(relative)) || any(startsWith(relative, "../"))) {
    stop("Runtime attestation could not derive safe relative paths.",
         call. = FALSE)
  }
  data.table::data.table(
    scope = scope, owner = owner, version = version,
    root_path = normalizePath(root, mustWork = TRUE),
    relative_path = relative, bytes = as.numeric(info$size),
    sha256 = vapply(normalized, .ranktree_sha, character(1L))
  )
}

.ranktree_capture_runtime_attestation <- function() {
  image_sha <- Sys.getenv("OMICSELECTOR_RUNTIME_IMAGE_SHA256", unset = "")
  launcher_sha <- Sys.getenv("OMICSELECTOR_CONTAINER_LAUNCHER_SHA256",
                             unset = "")
  container_marker <- Sys.getenv("SINGULARITY_CONTAINER", unset = "")
  container_runtime <- "singularity"
  if (!nzchar(container_marker)) {
    container_marker <- Sys.getenv("APPTAINER_CONTAINER", unset = "")
    container_runtime <- "apptainer"
  }
  if (!grepl("^[0-9a-f]{64}$", image_sha) ||
      !grepl("^[0-9a-f]{64}$", launcher_sha) ||
      !nzchar(container_marker)) {
    stop("Runtime attestation must run inside the package launcher and an ",
         "active Singularity/Apptainer container.", call. = FALSE)
  }
  runtime_roots <- unique(normalizePath(
    c(R.home("bin"), R.home("etc"), R.home("lib"), R.home("modules")),
    mustWork = TRUE
  ))
  runtime <- data.table::rbindlist(lapply(runtime_roots, function(root) {
    .ranktree_attestation_files(
      root, "r_runtime", basename(root), R.version.string
    )
  }))
  closure <- .ranktree_attestation_closure()
  libraries <- data.table::rbindlist(lapply(seq_len(nrow(closure)), function(i) {
    .ranktree_attestation_files(
      closure$root_path[[i]], "r_library", closure$package[[i]],
      closure$version[[i]]
    )
  }))
  out <- data.table::rbindlist(list(runtime, libraries), use.names = TRUE)
  out[, `:=`(
    r_version = R.version.string,
    r_home = normalizePath(R.home(), mustWork = TRUE),
    container_runtime = container_runtime,
    container_marker = container_marker,
    runtime_image_sha256 = image_sha,
    container_launcher_sha256 = launcher_sha
  )]
  data.table::setorder(out, scope, owner, relative_path)
  if (anyDuplicated(out[, .(scope, owner, relative_path)]) ||
      any(!grepl("^[0-9a-f]{64}$", out$sha256))) {
    stop("Runtime attestation contains duplicate or invalid file records.",
         call. = FALSE)
  }
  out
}

.ranktree_write_runtime_attestation <- function(path) {
  parent <- normalizePath(dirname(path), mustWork = TRUE)
  target <- file.path(parent, basename(path))
  if (file.exists(target)) {
    stop("Runtime attestation target already exists: ", target, call. = FALSE)
  }
  stage <- tempfile(".runtime-attestation-", tmpdir = parent)
  on.exit(if (file.exists(stage)) unlink(stage), add = TRUE)
  attestation <- .ranktree_capture_runtime_attestation()
  .public_fwrite(attestation, stage,
                     sep = "\t", quote = FALSE)
  .ranktree_validate_runtime_attestation(
    stage, .ranktree_sha(stage), unique(attestation$r_version),
    unique(attestation$runtime_image_sha256),
    unique(attestation$container_launcher_sha256), current = attestation
  )
  if (!file.rename(stage, target)) {
    stop("Could not atomically write runtime attestation.", call. = FALSE)
  }
  cat("Created runtime attestation:", target, "\n")
  invisible(target)
}

.ranktree_validate_runtime_attestation <- function(
    path, expected_sha256, expected_r_version, expected_image_sha256,
    expected_launcher_sha256,
    current = .ranktree_capture_runtime_attestation()) {
  path <- .ranktree_assert_file_pin(path, expected_sha256,
                                    "Runtime closure attestation")
  manifest <- data.table::fread(path, colClasses = "character")
  required <- c(
    "scope", "owner", "version", "root_path", "relative_path", "bytes",
    "sha256", "r_version", "r_home", "container_runtime",
    "container_marker", "runtime_image_sha256", "container_launcher_sha256"
  )
  if (!identical(names(manifest), required) || !nrow(manifest) ||
      anyNA(manifest) || any(!nzchar(as.matrix(manifest))) ||
      anyDuplicated(manifest[, .(scope, owner, relative_path)]) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    stop("Runtime closure attestation has an invalid fail-closed schema.",
         call. = FALSE)
  }
  scalar <- function(column, expected, label) {
    observed <- unique(manifest[[column]])
    if (length(observed) != 1L || !identical(observed, expected)) {
      stop(label, " differs from the exact runtime attestation.",
           call. = FALSE)
    }
  }
  scalar("r_version", expected_r_version, "R version")
  scalar("runtime_image_sha256", expected_image_sha256,
         "Runtime image SHA-256")
  scalar("container_launcher_sha256", expected_launcher_sha256,
         "Container launcher SHA-256")
  current <- data.table::as.data.table(current)
  current[, bytes := as.character(bytes)]
  manifest[, bytes := as.character(bytes)]
  if (!identical(current, manifest)) {
    stop("In-process R runtime/library bytes differ from the exact frozen ",
         "attestation.", call. = FALSE)
  }
  list(path = path, sha256 = expected_sha256, manifest = manifest)
}

.ranktree_run_bound_provenance <- function(
    cohorts, paper_root, provenance_script, provenance_manifest,
    union_inventory, expected_script_sha256, expected_manifest_sha256,
    expected_union_sha256, log_path, preflight_fn) {
  paper_root <- normalizePath(paper_root, mustWork = TRUE)
  canonical <- c(
    script = file.path(paper_root, "code", "check_provenance.R"),
    manifest = file.path(paper_root, "data_bucket",
                         "_PROVENANCE_MANIFEST.md"),
    union = file.path(paper_root, "knowledge", "union_inventory.tsv")
  )
  supplied <- c(script = provenance_script, manifest = provenance_manifest,
                union = union_inventory)
  supplied <- vapply(supplied, normalizePath, character(1L), mustWork = TRUE)
  canonical <- vapply(canonical, normalizePath, character(1L), mustWork = TRUE)
  if (!identical(supplied, canonical)) {
    stop("Provenance inputs must be the exact canonical paths under the ",
         "declared paper root.", call. = FALSE)
  }
  expected <- c(expected_script_sha256, expected_manifest_sha256,
                expected_union_sha256)
  observed_before <- vapply(supplied, .ranktree_sha, character(1L))
  if (!identical(unname(observed_before), expected)) {
    stop("A canonical provenance input differs from its exact pin.",
         call. = FALSE)
  }
  inventory <- data.table::fread(supplied[["union"]], colClasses = "character")
  if (!"accession" %in% names(inventory) || !nrow(inventory) ||
      anyNA(inventory$accession) || any(!nzchar(inventory$accession)) ||
      anyDuplicated(toupper(inventory$accession))) {
    stop("Canonical provenance union inventory has an invalid accession key.",
         call. = FALSE)
  }
  source_accessions <- unique(toupper(unlist(lapply(cohorts, function(cohort) {
    value <- cohort$source_accessions
    if (is.null(value) || !length(value)) value <- cohort$preflight_accessions
    if (is.null(value) || !length(value)) value <- cohort$accession
    as.character(value)
  }), use.names = FALSE)))
  if (!length(source_accessions) || anyNA(source_accessions) ||
      any(!nzchar(source_accessions)) ||
      length(setdiff(source_accessions, toupper(inventory$accession)))) {
    stop("Every benchmark source accession must resolve in the exact union ",
         "inventory.", call. = FALSE)
  }
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(paper_root)
  result <- preflight_fn(cohorts, paper_root, log_path)
  observed_after <- vapply(supplied, .ranktree_sha, character(1L))
  if (!identical(observed_after, observed_before)) {
    stop("A canonical provenance input mutated during preflight.",
         call. = FALSE)
  }
  result
}

.ranktree_structural_reason <- function() {
  paste0("n_valid_folds=0 < 3; per-fold: ", paste(rep(
    "method/baseline constant or all-NA on test fold", 3L), collapse = "; "))
}

.ranktree_abort_engine <- function(stage, method, cohort, seed, fold,
                                   parent = NULL, detail = NULL) {
  message <- paste0(
    "Strict structural-proof engine failure [", stage, "] for ", method,
    " / ", cohort, " / seed ", seed, " / fold ", fold,
    if (!is.null(detail)) paste0(": ", detail) else if (!is.null(parent)) {
      paste0(": ", conditionMessage(parent))
    } else ""
  )
  condition <- structure(
    list(message = message, call = NULL, stage = stage, method = method,
         cohort = cohort, seed = as.integer(seed), fold = as.integer(fold),
         parent = parent),
    class = c("ranktree_engine_error", "error", "condition")
  )
  stop(condition)
}

.ranktree_strict_engine_call <- function(fn, args, stage, method, cohort, seed,
                                         fold) {
  tryCatch(
    do.call(fn, args),
    error = function(error) {
      .ranktree_abort_engine(stage, method, cohort, seed, fold, parent = error)
    }
  )
}

.ranktree_strict_score_vector <- function(value, n, stage, method, cohort,
                                          seed, fold) {
  if (!is.numeric(value) || length(value) != n || anyNA(value) ||
      any(!is.finite(value))) {
    .ranktree_abort_engine(
      stage, method, cohort, seed, fold,
      detail = "score must be a complete finite numeric held-out vector"
    )
  }
  as.numeric(value)
}

.ranktree_prove_structural_exception <- function(
    method, cohort_id, cohort, roster, seed, outer_k,
    grouped_fold_fn = get(".grouped_folds", envir = .GlobalEnv),
    stratified_fold_fn = get(".stratified_folds", envir = .GlobalEnv),
    baseline_fn = get(".srb_rclr_panel_score", envir = .GlobalEnv),
    fit_dispatch = get(".srb_call_fit", envir = .GlobalEnv),
    score_dispatch = OmicSelector::singlesample_score_call,
    direction_fn = get(".srb_training_direction", envir = .GlobalEnv)) {
  if (!identical(cohort_id, "GSE83977") || !identical(as.integer(outer_k), 3L) ||
      !as.integer(seed) %in% .ranktree_seeds()) {
    stop("Structural proof is restricted to the registered GSE83977/3-fold ",
         "exception.", call. = FALSE)
  }
  X <- cohort$expr_per_sample
  y <- as.integer(cohort$y_bin)
  group_id <- cohort$group_id
  if (!is.matrix(X) || nrow(X) != length(y) ||
      any(!y %in% c(0L, 1L)) || length(unique(y)) != 2L) {
    stop("Structural-proof cohort inputs are invalid.", call. = FALSE)
  }
  has_groups <- !is.null(group_id) && length(group_id) == nrow(X)
  folds <- if (has_groups) {
    grouped_fold_fn(y, group_id, k = outer_k, seed = seed)
  } else {
    stratified_fold_fn(y, k = outer_k, seed = seed)
  }
  if (!is.list(folds) || length(folds) != outer_k) {
    stop("Structural proof did not reproduce the exact outer-fold count.",
         call. = FALSE)
  }
  row <- roster[roster$method_id == method, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Structural proof did not resolve one registered method.", call. = FALSE)
  }
  proof <- vector("list", outer_k)
  for (fold in seq_len(outer_k)) {
    test <- as.integer(folds[[fold]])
    train <- setdiff(seq_len(nrow(X)), test)
    if (!length(test) || !length(train) || anyDuplicated(test) ||
        any(test < 1L | test > nrow(X)) ||
        (has_groups && length(intersect(unique(group_id[test]),
                                        unique(group_id[train]))))) {
      stop("Structural proof found an invalid or group-leaking outer fold.",
           call. = FALSE)
    }
    y_train <- y[train]
    y_test <- y[test]
    if (length(unique(y_train)) != 2L || length(unique(y_test)) != 2L ||
        sum(y_test == 0L) < cohort$min_per_fold ||
        sum(y_test == 1L) < cohort$min_per_fold) {
      stop("Structural exception is not attributable solely to constant ",
           "score support.", call. = FALSE)
    }
    X_train <- X[train, , drop = FALSE]
    X_test <- X[test, , drop = FALSE]
    baseline <- .ranktree_strict_engine_call(
      baseline_fn, list(X_train, y_train, X_test, panel_size = 20L),
      "baseline", method, cohort_id, seed, fold
    )
    baseline <- .ranktree_strict_score_vector(
      baseline, length(test), "baseline", method, cohort_id, seed, fold
    )
    if (identical(method, "ws-rclr-trimmed")) {
      method_score <- baseline
    } else {
      fit_name <- as.character(row$fit_fn[[1L]])
      if (is.na(fit_name) || !nzchar(fit_name)) {
        .ranktree_abort_engine(
          "fit", method, cohort_id, seed, fold,
          detail = "registered fit function is absent"
        )
      }
      model <- .ranktree_strict_engine_call(
        fit_dispatch, list(fit_name, X_train, y_train, NULL, list()),
        "fit", method, cohort_id, seed, fold
      )
      training_score <- .ranktree_strict_engine_call(
        score_dispatch,
        list(method, model, X_train, NULL, roster = roster),
        "training_prediction", method, cohort_id, seed, fold
      )
      training_score <- .ranktree_strict_score_vector(
        training_score, length(train), "training_prediction", method,
        cohort_id, seed, fold
      )
      direction <- .ranktree_strict_engine_call(
        direction_fn, list(training_score, y_train), "direction", method,
        cohort_id, seed, fold
      )
      if (!is.numeric(direction) || length(direction) != 1L ||
          !is.finite(direction) || !direction %in% c(-1, 1)) {
        .ranktree_abort_engine(
          "direction", method, cohort_id, seed, fold,
          detail = "training-only direction is unavailable"
        )
      }
      method_score <- .ranktree_strict_engine_call(
        score_dispatch, list(method, model, X_test, NULL, roster = roster),
        "heldout_prediction", method, cohort_id, seed, fold
      )
      method_score <- direction * .ranktree_strict_score_vector(
        method_score, length(test), "heldout_prediction", method, cohort_id,
        seed, fold
      )
    }
    baseline_constant <- data.table::uniqueN(baseline) < 2L
    method_constant <- data.table::uniqueN(method_score) < 2L
    if (!baseline_constant && !method_constant) {
      stop("Strict replay did not prove constant held-out support in every ",
           "outer fold.", call. = FALSE)
    }
    proof[[fold]] <- data.table::data.table(
      fold = as.integer(fold), n_train = length(train), n_test = length(test),
      n_positive_test = sum(y_test == 1L),
      n_negative_test = sum(y_test == 0L),
      baseline_constant = baseline_constant, method_constant = method_constant,
      test_index_sha256 = digest::digest(test, algo = "sha256", serialize = TRUE),
      baseline_score_sha256 = digest::digest(
        baseline, algo = "sha256", serialize = TRUE),
      method_score_sha256 = digest::digest(
        method_score, algo = "sha256", serialize = TRUE)
    )
  }
  proof <- data.table::rbindlist(proof)
  list(
    proven = TRUE,
    basis = paste0(
      "complete finite replay: method or rCLR baseline constant on every ",
      "exact group-safe held-out fold"
    ),
    proof = proof,
    proof_sha256 = digest::digest(proof, algo = "sha256", serialize = TRUE)
  )
}

.ranktree_validate_cell_result <- function(cell) {
  required <- c("method", "cohort", "seed", "eligible", "ineligible_reason",
                "n_eff", "n_valid_folds", "outer_k",
                "n_group_split_violations", "structural_exception_proven",
                "structural_exception_proof_sha256")
  if (!is.data.frame(cell) || nrow(cell) != 1L ||
      !all(required %in% names(cell)) || !is.logical(cell$eligible) ||
      is.na(cell$eligible)) {
    stop("Benchmark cell result has an invalid schema.", call. = FALSE)
  }
  if (isTRUE(cell$eligible)) {
    if (is.na(cell$n_eff) || cell$n_eff < 1L ||
        is.na(cell$n_valid_folds) || is.na(cell$outer_k) ||
        cell$n_valid_folds != cell$outer_k ||
        is.na(cell$n_group_split_violations) ||
        cell$n_group_split_violations != 0L ||
        (!is.na(cell$ineligible_reason) && nzchar(cell$ineligible_reason)) ||
        !identical(cell$structural_exception_proven, FALSE) ||
        !is.na(cell$structural_exception_proof_sha256)) {
      stop("Eligible benchmark cell violates fold/group accounting.", call. = FALSE)
    }
    return(invisible(TRUE))
  }
  accepted <- identical(as.character(cell$cohort), "GSE83977") &&
    as.integer(cell$seed) %in% .ranktree_seeds() &&
    identical(as.integer(cell$n_valid_folds), 0L) &&
    identical(as.integer(cell$outer_k), 3L) &&
    identical(as.integer(cell$n_group_split_violations), 0L) &&
    identical(as.character(cell$ineligible_reason), .ranktree_structural_reason()) &&
    identical(cell$structural_exception_proven, TRUE) &&
    is.character(cell$structural_exception_proof_sha256) &&
    length(cell$structural_exception_proof_sha256) == 1L &&
    grepl("^[0-9a-f]{64}$", cell$structural_exception_proof_sha256)
  if (!accepted) {
    stop("Unreviewed ineligibility, fit/score failure, or dependency/runtime ",
         "failure: ", as.character(cell$ineligible_reason), call. = FALSE)
  }
  invisible(TRUE)
}

.ranktree_cell_row <- function(result, seed, package_pins) {
  `%rt||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
  data.table::data.table(
    method = result$method_id, cohort = result$cohort, seed = as.integer(seed),
    eligible = isTRUE(result$eligible),
    ineligible_reason = result$ineligible_reason %rt||% NA_character_,
    auc_method = result$auc_method %rt||% NA_real_,
    auc_baseline = result$auc_baseline %rt||% NA_real_,
    lift = result$est %rt||% NA_real_, lift_se = result$se %rt||% NA_real_,
    n_eff = result$n_eff %rt||% NA_integer_,
    n_valid_folds = result$n_valid_folds %rt||% NA_integer_,
    min_pos = result$min_pos %rt||% NA_integer_,
    min_neg = result$min_neg %rt||% NA_integer_,
    outer_k = result$outer_k %rt||% NA_integer_,
    grouping_policy = result$grouping_policy %rt||% NA_character_,
    n_group_split_violations =
      result$n_group_split_violations %rt||% NA_integer_,
    structural_exception_proven = FALSE,
    structural_exception_basis = NA_character_,
    structural_exception_proof_sha256 = NA_character_,
    package_version = package_pins$version,
    package_commit = package_pins$commit,
    cache_package_commit = package_pins$cache_commit,
    cache_sha256 = package_pins$cache_sha256,
    engine_manifest_sha256 = package_pins$engine_manifest_sha256,
    environment_manifest_sha256 = package_pins$environment_manifest_sha256,
    runtime_image_sha256 = package_pins$runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      package_pins$runtime_attestation_manifest_sha256,
    container_launcher_sha256 = package_pins$container_launcher_sha256,
    r_version = package_pins$r_version,
    analysis_plan_sha256 = package_pins$analysis_plan_sha256,
    analysis_amendment_sha256 = package_pins$analysis_amendment_sha256,
    analysis_amendment2_sha256 = package_pins$analysis_amendment2_sha256,
    analysis_amendment3_sha256 = package_pins$analysis_amendment3_sha256,
    provenance_manifest_sha256 = package_pins$provenance_manifest_sha256,
    provenance_union_inventory_sha256 =
      package_pins$provenance_union_inventory_sha256,
    analysis_code_id = package_pins$analysis_code_id
  )
}

.ranktree_add_split_ids <- function(predictions) {
  predictions <- data.table::as.data.table(predictions)
  if (!nrow(predictions)) return(predictions[, split_id := character()])
  required <- c("method", "cohort", "seed", "fold", "sample_id", "group_id",
                "y", "score_m", "n_train", "n_test")
  if (!all(required %in% names(predictions)) ||
      anyNA(predictions[, ..required]) || any(!is.finite(predictions$score_m)) ||
      any(!predictions$y %in% c(0L, 1L)) ||
      anyDuplicated(predictions[, .(method, cohort, seed, sample_id)])) {
    stop("Held-out predictions are incomplete, invalid, or duplicated.",
         call. = FALSE)
  }
  group_fold <- predictions[, .(n_fold = data.table::uniqueN(fold),
                                n_y = data.table::uniqueN(y)),
                            by = .(method, cohort, seed, group_id)]
  if (any(group_fold$n_fold != 1L) || any(group_fold$n_y != 1L)) {
    stop("A provenance group crosses a held-out fold or outcome label.",
         call. = FALSE)
  }
  split_map <- predictions[, {
    reference <- .SD[, .(fold = as.integer(fold), sample_id = as.character(sample_id),
                         group_id = as.character(group_id), y = as.integer(y),
                         n_train = as.integer(n_train), n_test = as.integer(n_test))]
    data.table::setorder(reference, fold, sample_id, group_id, y)
    .(split_id = digest::digest(reference, algo = "sha256", serialize = TRUE))
  }, by = .(method, cohort, seed)]
  out <- merge(predictions, split_map,
               by = c("method", "cohort", "seed"), all.x = TRUE, sort = FALSE)
  cross_method <- unique(out[, .(method, cohort, seed, split_id)])[
    , .(n_split = data.table::uniqueN(split_id)), by = .(cohort, seed)]
  if (any(cross_method$n_split != 1L)) {
    stop("Selected methods do not share byte-identical held-out splits.",
         call. = FALSE)
  }
  out
}

.ranktree_units <- function(cohorts) {
  data.table::rbindlist(lapply(names(cohorts), function(id) {
    cohort <- cohorts[[id]]
    data.table::data.table(
      unit_id = id, disease = as.character(cohort$disease),
      modality = as.character(cohort$modality),
      biospecimen = as.character(cohort$biospecimen),
      provenance_block = as.character(cohort$provenance_block),
      source_accessions = paste(as.character(cohort$source_accessions),
                                collapse = ";"),
      outer_k = as.integer(cohort$outer_k),
      primary_unit = !grepl("whole[ -]?blood", cohort$biospecimen,
                            ignore.case = TRUE),
      n_deposited_profiles = nrow(cohort$expr_per_sample),
      n_unique_groups = data.table::uniqueN(as.character(cohort$group_id))
    )
  }), fill = TRUE)
}

.ranktree_validate_runtime_bundle <- function(cells, predictions, units,
                                               methods,
                                               runtime_image_sha256 = NULL,
                                               runtime_attestation_sha256 = NULL,
                                               container_launcher_sha256 = NULL,
                                               r_version = NULL) {
  expected <- data.table::CJ(method = methods, cohort = units$unit_id,
                             seed = .ranktree_seeds(), unique = TRUE)
  if (nrow(units) != 34L || sum(units$primary_unit) != 33L ||
      anyDuplicated(units$unit_id) ||
      nrow(cells) != nrow(expected) ||
      anyDuplicated(cells[, .(method, cohort, seed)]) ||
      !setequal(do.call(paste, cells[, .(method, cohort, seed)]),
                do.call(paste, expected))) {
    stop("Runtime bundle is not the exact method-by-34-unit-by-five-seed grid.",
         call. = FALSE)
  }
  invisible(lapply(seq_len(nrow(cells)), function(i) {
    .ranktree_validate_cell_result(cells[i])
  }))
  if (any(!cells$eligible)) {
    stop("Amendment 2 requires every rank-tree sensitivity cell to satisfy ",
         "the complete-finite held-out support contract.", call. = FALSE)
  }
  if (!is.null(runtime_image_sha256) &&
      (!"runtime_image_sha256" %in% names(cells) ||
       !"runtime_image_sha256" %in% names(predictions) ||
       anyNA(cells$runtime_image_sha256) ||
       anyNA(predictions$runtime_image_sha256) ||
       any(cells$runtime_image_sha256 != runtime_image_sha256) ||
       any(predictions$runtime_image_sha256 != runtime_image_sha256))) {
    stop("Runtime image pin is absent or inconsistent in cells/predictions.",
         call. = FALSE)
  }
  runtime_columns <- c(
    runtime_attestation_manifest_sha256 = runtime_attestation_sha256,
    container_launcher_sha256 = container_launcher_sha256,
    r_version = r_version
  )
  for (column in names(runtime_columns)) {
    expected_value <- runtime_columns[[column]]
    if (!is.null(expected_value) &&
        (!column %in% names(cells) || !column %in% names(predictions) ||
         anyNA(cells[[column]]) || anyNA(predictions[[column]]) ||
         any(cells[[column]] != expected_value) ||
         any(predictions[[column]] != expected_value))) {
      stop("Runtime attestation field is absent or inconsistent: ", column,
           call. = FALSE)
    }
  }
  eligible <- cells[eligible == TRUE, .(method, cohort, seed, n_eff, outer_k)]
  pred_cells <- unique(predictions[, .(method, cohort, seed)])
  if (!setequal(do.call(paste, eligible[, .(method, cohort, seed)]),
                do.call(paste, pred_cells))) {
    stop("Eligible cells and held-out prediction cells disagree.", call. = FALSE)
  }
  accounting <- merge(
    eligible, predictions[, .(n_prediction = .N,
                              n_fold = data.table::uniqueN(fold)),
                          by = .(method, cohort, seed)],
    by = c("method", "cohort", "seed"), all.x = TRUE
  )
  if (anyNA(accounting$n_prediction) ||
      any(accounting$n_prediction != accounting$n_eff) ||
      any(accounting$n_fold != accounting$outer_k)) {
    stop("Eligible prediction cardinality/fold accounting is incomplete.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.ranktree_copy <- function(from, to, label) {
  if (!file.copy(from, to, overwrite = FALSE, copy.mode = TRUE)) {
    stop("Could not copy ", label, ".", call. = FALSE)
  }
  invisible(to)
}

.ranktree_checkpoint_contract <- function(output_dir, methods, units, pins,
                                          runner_sha, r_version) {
  data.table::as.data.table(data.frame(
    key = c(
      "schema_version", "output_dir", "methods", "seeds", "units",
      "package_version", "package_commit", "cache_package_commit",
      "cache_sha256", "cache_manifest_sha256", "cache_digest_sha256",
      "runner_script_sha256", "engine_manifest_sha256",
      "environment_manifest_sha256", "runtime_image_sha256",
      "runtime_attestation_manifest_sha256",
      "container_launcher_sha256",
      "analysis_plan_sha256", "analysis_amendment_sha256",
      "analysis_amendment2_sha256",
      "analysis_amendment3_sha256",
      "provenance_script_sha256",
      "provenance_manifest_sha256", "provenance_union_inventory_sha256",
      "analysis_code_id", "r_version"
    ),
    value = c(
      "omicselector-paper3-ranktree-checkpoint-v2",
      normalizePath(output_dir, mustWork = FALSE),
      paste(methods, collapse = ","), paste(.ranktree_seeds(), collapse = ","),
      paste(units$unit_id, collapse = ","), pins$version, pins$commit,
      pins$cache_commit, pins$cache_sha256, pins$cache_manifest_sha256,
      pins$cache_digest_sha256, runner_sha, pins$engine_manifest_sha256,
      pins$environment_manifest_sha256, pins$runtime_image_sha256,
      pins$runtime_attestation_manifest_sha256,
      pins$container_launcher_sha256,
      pins$analysis_plan_sha256,
      pins$analysis_amendment_sha256,
      pins$analysis_amendment2_sha256,
      pins$analysis_amendment3_sha256,
      pins$provenance_script_sha256,
      pins$provenance_manifest_sha256,
      pins$provenance_union_inventory_sha256,
      pins$analysis_code_id, r_version
    ), stringsAsFactors = FALSE
  ))
}

.ranktree_validate_checkpoint_contract <- function(checkpoint_dir, expected) {
  if (!dir.exists(checkpoint_dir) || nzchar(Sys.readlink(checkpoint_dir))) {
    stop("Persistent checkpoint path is absent, not a directory, or a symlink.",
         call. = FALSE)
  }
  path <- file.path(checkpoint_dir, "checkpoint_contract.tsv")
  if (!file.exists(path)) {
    stop("Persistent checkpoint directory has an incomplete contract.",
         call. = FALSE)
  }
  observed <- tryCatch(data.table::fread(
    path, colClasses = c(key = "character", value = "character")),
    error = identity)
  if (inherits(observed, "error") ||
      !identical(names(observed), c("key", "value")) ||
      nrow(observed) != nrow(expected) || anyNA(observed) ||
      anyDuplicated(observed$key) || !identical(observed, expected)) {
    stop("Persistent checkpoint contract is mismatched or incomplete.",
         call. = FALSE)
  }
  invisible(path)
}

.ranktree_open_checkpoints <- function(output_dir, expected) {
  checkpoint_dir <- paste0(normalizePath(output_dir, mustWork = FALSE),
                           ".checkpoints")
  if (file.exists(checkpoint_dir)) {
    .ranktree_validate_checkpoint_contract(checkpoint_dir, expected)
    return(checkpoint_dir)
  }
  parent <- normalizePath(dirname(checkpoint_dir), mustWork = TRUE)
  stage <- tempfile(".ranktree-checkpoint-init-", tmpdir = parent)
  if (!dir.create(stage, recursive = FALSE)) {
    stop("Could not create checkpoint-contract staging directory.",
         call. = FALSE)
  }
  promoted <- FALSE
  on.exit(if (!promoted && dir.exists(stage)) unlink(stage, recursive = TRUE),
          add = TRUE)
  .public_fwrite(expected, file.path(stage, "checkpoint_contract.tsv"),
                     sep = "\t", quote = FALSE)
  if (!file.rename(stage, checkpoint_dir)) {
    stop("Could not atomically create the persistent checkpoint directory.",
         call. = FALSE)
  }
  promoted <- TRUE
  .ranktree_validate_checkpoint_contract(checkpoint_dir, expected)
  checkpoint_dir
}

.ranktree_bundle_artifacts <- function() {
  c("analysis_plan.md", "analysis_amendment.md", "analysis_amendment2.md",
    "analysis_amendment3.md",
    "engine_manifest.tsv", "environment_manifest.tsv",
    "methods.tsv",
    "provenance_preflight.log", "provenance_resolution.tsv", "report.md",
    "runtime_attestation.tsv",
    "units.tsv", "within_bundle.rds", "within_cells.tsv",
    "within_predictions.tsv.gz")
}

.ranktree_validate_completed_output <- function(output_dir) {
  if (!dir.exists(output_dir) || nzchar(Sys.readlink(output_dir))) {
    stop("Final output is absent, not a directory, or a symlink.",
         call. = FALSE)
  }
  manifest_path <- file.path(output_dir, "output_manifest.tsv")
  if (!file.exists(manifest_path)) {
    stop("Final output lacks its completion manifest.", call. = FALSE)
  }
  manifest <- data.table::fread(manifest_path)
  required <- .ranktree_bundle_artifacts()
  if (!all(c("file", "bytes", "sha256") %in% names(manifest)) ||
      anyNA(manifest[, .(file, bytes, sha256)]) || anyDuplicated(manifest$file) ||
      !setequal(manifest$file, required) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    stop("Final output completion manifest is malformed or incomplete.",
         call. = FALSE)
  }
  paths <- file.path(output_dir, manifest$file)
  if (any(!file.exists(paths)) || any(file.info(paths)$size != manifest$bytes) ||
      !identical(unname(vapply(paths, .ranktree_sha, character(1L))),
                 as.character(manifest$sha256))) {
    stop("Final output bytes differ from the completion manifest.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.ranktree_remove_completed_checkpoints <- function(output_dir, expected) {
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  checkpoint_dir <- paste0(output_dir, ".checkpoints")
  .ranktree_validate_completed_output(output_dir)
  .ranktree_validate_checkpoint_contract(checkpoint_dir, expected)
  status <- unlink(checkpoint_dir, recursive = TRUE, force = FALSE)
  if (!identical(status, 0L) || dir.exists(checkpoint_dir)) {
    stop("Final output is complete but its checkpoint directory could not be ",
         "removed.", call. = FALSE)
  }
  invisible(TRUE)
}

.ranktree_strict_exception_cell <- function(method, cohort_id, cohort, roster,
                                            out_dir, seed, outer_k, ...) {
  if (!identical(cohort_id, "GSE83977")) {
    stop("Strict exception cell is restricted to GSE83977.", call. = FALSE)
  }
  checkpoint <- file.path(
    out_dir,
    paste0("strict__", cohort_id, "__", gsub("[^A-Za-z0-9_.-]", "_", method),
           ".rds")
  )
  if (file.exists(checkpoint)) {
    cached <- .public_readRDS(checkpoint)
    proof <- cached$structural_exception_proof
    if (!is.list(cached) || !identical(cached$method_id, method) ||
        !identical(cached$cohort, cohort_id) ||
        !identical(as.integer(cached$seed), as.integer(seed)) ||
        !identical(as.integer(cached$outer_k), as.integer(outer_k)) ||
        !is.data.frame(proof) || nrow(proof) != outer_k ||
        !identical(
          digest::digest(proof, algo = "sha256", serialize = TRUE),
          cached$structural_exception_proof_sha256
        )) {
      stop("Strict structural-exception checkpoint is invalid.", call. = FALSE)
    }
    return(cached)
  }
  replay <- .ranktree_prove_structural_exception(
    method, cohort_id, cohort, roster, seed, outer_k, ...
  )
  result <- list(
    method_id = method, cohort = cohort_id, seed = as.integer(seed),
    eligible = FALSE, ineligible_reason = .ranktree_structural_reason(),
    n_eff = NA_integer_, n_valid_folds = 0L,
    min_pos = NA_integer_, min_neg = NA_integer_, outer_k = as.integer(outer_k),
    grouping_policy = if (!is.null(cohort$group_id)) {
      "groupKFold"
    } else "ungrouped_stratified",
    n_group_split_violations = 0L, fold_scores = NULL,
    structural_exception_proven = replay$proven,
    structural_exception_basis = replay$basis,
    structural_exception_proof = replay$proof,
    structural_exception_proof_sha256 = replay$proof_sha256
  )
  stage <- tempfile(".strict-cell-", tmpdir = out_dir, fileext = ".rds")
  on.exit(if (file.exists(stage)) unlink(stage), add = TRUE)
  saveRDS(result, stage, version = 3L)
  if (!file.rename(stage, checkpoint)) {
    stop("Could not atomically persist strict structural-exception cell.",
         call. = FALSE)
  }
  result
}

.ranktree_call_benchmark_cell <- function(
    benchmark_fn, method, cohort_id, cohort, roster, out_dir, seed, outer_k,
    strict_exception_fn = .ranktree_strict_exception_cell) {
  # Amendment 2 applies one uniform complete-finite score contract. Constant
  # held-out scores are valid midrank-AUC evidence, including for GSE83977; the
  # reference strict-exception helper is retained only for audit replay.
  benchmark_fn(
    method, cohort_id, cohort, roster, out_dir = out_dir, seed = seed,
    outer_k = outer_k, min_valid_folds = outer_k,
    min_pos_per_fold = cohort$min_per_fold,
    min_neg_per_fold = cohort$min_per_fold,
    baseline_mode = "rclr_panel", panel_size = 20L,
    overwrite = FALSE, rerun_present = TRUE, orient = "fixed", hp = list()
  )
}

.ranktree_runner_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- .ranktree_runner_args(args)
  code_start <- .ranktree_git_provenance(opt$package_root)
  if (code_start$dirty) {
    stop("Rank-tree sensitivity execution requires a clean OmicSelector ",
         "worktree.", call. = FALSE)
  }
  if (!identical(code_start$commit, opt$expected_package_commit)) {
    stop("OmicSelector checkout does not match the expected package commit.",
         call. = FALSE)
  }
  description <- read.dcf(file.path(code_start$root, "DESCRIPTION"))
  package_version <- unname(description[, "Version"])
  if (!identical(package_version, opt$expected_package_version)) {
    stop("OmicSelector DESCRIPTION version does not match its expected pin.",
         call. = FALSE)
  }
  script_path <- .ranktree_script_path()
  expected_script <- normalizePath(file.path(
    code_start$root, "inst", "scripts", basename(script_path)), mustWork = TRUE)
  if (!identical(script_path, expected_script)) {
    stop("Runner must execute from the pinned OmicSelector checkout.",
         call. = FALSE)
  }
  runner_sha <- .ranktree_sha(script_path)
  launcher_path <- .ranktree_assert_file_pin(
    file.path(code_start$root, "inst", "scripts",
              "paper3_run_ranktree_sensitivity_container.sh"),
    opt$expected_container_launcher_sha256, "Package container launcher"
  )

  engine_manifest <- .ranktree_validate_engine_manifest(
    opt$engine_manifest, opt$expected_engine_manifest_sha256)
  environment_manifest <- .ranktree_validate_environment_manifest(
    opt$environment_manifest, opt$expected_environment_manifest_sha256)
  runtime_image <- .ranktree_assert_file_pin(
    opt$runtime_image, opt$expected_runtime_image_sha256,
    "Singularity runtime image")
  if (!identical(R.version.string, opt$expected_r_version)) {
    stop("In-process R version differs from its exact pre-execution pin.",
         call. = FALSE)
  }
  runtime_attestation <- .ranktree_validate_runtime_attestation(
    opt$runtime_attestation_manifest,
    opt$expected_runtime_attestation_manifest_sha256,
    opt$expected_r_version, opt$expected_runtime_image_sha256,
    opt$expected_container_launcher_sha256
  )

  paper_root <- normalizePath(opt$paper_root, mustWork = TRUE)
  analysis_plan <- .ranktree_validate_analysis_plan(
    opt$analysis_plan, paper_root, opt$expected_analysis_plan_sha256
  )
  analysis_amendment <- .ranktree_validate_analysis_amendment(
    opt$analysis_amendment, paper_root,
    opt$expected_analysis_amendment_sha256, amendment_number = 1L
  )
  analysis_amendment2 <- .ranktree_validate_analysis_amendment(
    opt$analysis_amendment2, paper_root,
    opt$expected_analysis_amendment2_sha256, amendment_number = 2L
  )
  analysis_amendment3 <- .ranktree_validate_analysis_amendment(
    opt$analysis_amendment3, paper_root,
    opt$expected_analysis_amendment3_sha256, amendment_number = 3L
  )
  cache <- .ranktree_assert_file_pin(opt$cache, opt$expected_cache_sha256,
                                    "Cohort cache")
  cache_manifest_path <- file.path(dirname(cache), "cohort_cache_manifest.tsv")
  cache_digest_path <- file.path(dirname(cache), "cohort_cache_digest.tsv")
  .ranktree_assert_file_pin(cache_manifest_path,
                            opt$expected_cache_manifest_sha256,
                            "Cohort cache manifest")
  .ranktree_assert_file_pin(cache_digest_path,
                            opt$expected_cache_digest_sha256,
                            "Cohort cache digest")
  provenance_script <- .ranktree_assert_file_pin(
    file.path(paper_root, "code", "check_provenance.R"),
    opt$expected_provenance_script_sha256, "Provenance preflight script")
  provenance_manifest <- .ranktree_assert_file_pin(
    opt$provenance_manifest,
    opt$expected_provenance_manifest_sha256, "Canonical provenance manifest")
  provenance_union_inventory <- .ranktree_assert_file_pin(
    opt$provenance_union_inventory,
    opt$expected_provenance_union_inventory_sha256,
    "Canonical provenance union inventory")

  analysis_code_id <- substr(digest::digest(paste(
    runner_sha, opt$expected_package_commit, opt$expected_cache_sha256,
    opt$expected_analysis_plan_sha256,
    opt$expected_analysis_amendment_sha256,
    opt$expected_analysis_amendment2_sha256,
    opt$expected_analysis_amendment3_sha256,
    opt$expected_engine_manifest_sha256,
    opt$expected_environment_manifest_sha256,
    opt$expected_runtime_image_sha256,
    opt$expected_runtime_attestation_manifest_sha256,
    opt$expected_container_launcher_sha256, opt$expected_r_version,
    opt$expected_provenance_script_sha256,
    opt$expected_provenance_manifest_sha256,
    opt$expected_provenance_union_inventory_sha256, sep = "|"),
    algo = "sha256"),
    1L, 16L)
  if (!identical(analysis_code_id, opt$expected_analysis_code_id)) {
    stop("Computed analysis code id does not match its pre-execution pin.",
         call. = FALSE)
  }

  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The exact sensitivity runner requires pkgload.", call. = FALSE)
  }
  suppressMessages(suppressWarnings(
    pkgload::load_all(code_start$root, quiet = TRUE, export_all = FALSE,
                      helpers = FALSE)
  ))
  if (!identical(as.character(utils::packageVersion("OmicSelector")),
                 opt$expected_package_version)) {
    stop("Loaded OmicSelector version differs from the pinned checkout.",
         call. = FALSE)
  }

  old_engine_root <- Sys.getenv("OMICSELECTOR_PAPER_ROOT", unset = NA_character_)
  on.exit({
    if (is.na(old_engine_root)) Sys.unsetenv("OMICSELECTOR_PAPER_ROOT")
    else Sys.setenv(OMICSELECTOR_PAPER_ROOT = old_engine_root)
  }, add = TRUE)
  engine_root <- dirname(dirname(dirname(
    engine_manifest[role == "benchmark_engine", path]
  )))
  Sys.setenv(OMICSELECTOR_PAPER_ROOT = engine_root)
  for (role_value in c("paired_auc_engine", "fold_engine", "benchmark_engine",
                       "public_design")) {
    source_path <- engine_manifest[role == role_value, path]
    if (length(source_path) != 1L) {
      stop("Frozen engine role did not resolve uniquely: ", role_value,
           call. = FALSE)
    }
    sys.source(source_path, envir = .GlobalEnv)
  }
  required_engine_fns <- c("srb_verify_cache_bundle",
                           "srb_run_provenance_preflight",
                           "srb_benchmark_within_cell", "srb_expected_outer_k",
                           "srb_sanitize_character_fields")
  if (!all(vapply(required_engine_fns, exists, logical(1L), mode = "function",
                  envir = .GlobalEnv, inherits = FALSE))) {
    stop("Frozen benchmark sources did not expose their required API.",
         call. = FALSE)
  }
  srb_verify_cache_bundle(cache, cache_manifest_path, cache_digest_path)
  cache_manifest <- data.table::fread(cache_manifest_path)
  if (!"package_commit" %in% names(cache_manifest) ||
      data.table::uniqueN(cache_manifest$package_commit) != 1L ||
      !identical(unique(cache_manifest$package_commit),
                 opt$expected_cache_commit)) {
    stop("Cohort cache manifest does not match the expected cache-building ",
         "commit.", call. = FALSE)
  }
  cohorts <- .public_readRDS(cache)
  if (!is.list(cohorts) || length(cohorts) != 34L ||
      anyNA(names(cohorts)) || any(!nzchar(names(cohorts))) ||
      anyDuplicated(names(cohorts))) {
    stop("Cohort cache is not the exact 34-unit public benchmark corpus.",
         call. = FALSE)
  }
  units <- .ranktree_units(cohorts)
  expected_k <- as.integer(srb_expected_outer_k(units$unit_id))
  if (any(units$outer_k != expected_k) || sum(units$primary_unit) != 33L) {
    stop("Cohort cache violates frozen outer-fold or primary-unit contracts.",
         call. = FALSE)
  }

  base_roster <- OmicSelector::singlesample_method_roster()
  anchor <- base_roster[match(.ranktree_methods()[1:5], base_roster$method_id), ]
  external <- OmicSelector::singlesample_external_competitor_roster()
  roster <- data.table::rbindlist(list(anchor, external), fill = TRUE)
  if (nrow(roster) != 7L || anyNA(roster$method_id) ||
      !identical(roster$method_id, .ranktree_methods()) ||
      any(roster$estimand != "within")) {
    stop("Installed package does not expose the frozen seven-method ",
         "sensitivity registry.", call. = FALSE)
  }

  output_parent <- normalizePath(dirname(opt$output_dir), mustWork = TRUE)
  output_dir <- file.path(output_parent, basename(opt$output_dir))
  if (file.exists(output_dir)) stop("Output bundle already exists: ", output_dir)
  stage <- tempfile(".ranktree-run-", tmpdir = output_parent)
  dir.create(stage)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
  provenance_log <- file.path(stage, "provenance_preflight.log")
  provenance <- .ranktree_run_bound_provenance(
    cohorts, paper_root, provenance_script, provenance_manifest,
    provenance_union_inventory, opt$expected_provenance_script_sha256,
    opt$expected_provenance_manifest_sha256,
    opt$expected_provenance_union_inventory_sha256, provenance_log,
    srb_run_provenance_preflight
  )
  if (!file.exists(provenance$resolution_path)) {
    stop("Provenance preflight did not emit its resolution table.",
         call. = FALSE)
  }

  pins <- list(
    version = package_version, commit = code_start$commit,
    cache_commit = opt$expected_cache_commit,
    cache_sha256 = opt$expected_cache_sha256,
    cache_manifest_sha256 = opt$expected_cache_manifest_sha256,
    cache_digest_sha256 = opt$expected_cache_digest_sha256,
    engine_manifest_sha256 = opt$expected_engine_manifest_sha256,
    environment_manifest_sha256 = opt$expected_environment_manifest_sha256,
    runtime_image_sha256 = opt$expected_runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      opt$expected_runtime_attestation_manifest_sha256,
    container_launcher_sha256 = opt$expected_container_launcher_sha256,
    r_version = opt$expected_r_version,
    analysis_plan_sha256 = opt$expected_analysis_plan_sha256,
    analysis_amendment_sha256 =
      opt$expected_analysis_amendment_sha256,
    analysis_amendment2_sha256 =
      opt$expected_analysis_amendment2_sha256,
    analysis_amendment3_sha256 =
      opt$expected_analysis_amendment3_sha256,
    provenance_script_sha256 = opt$expected_provenance_script_sha256,
    provenance_manifest_sha256 = opt$expected_provenance_manifest_sha256,
    provenance_union_inventory_sha256 =
      opt$expected_provenance_union_inventory_sha256,
    analysis_code_id = analysis_code_id
  )
  checkpoint_contract <- .ranktree_checkpoint_contract(
    output_dir, opt$methods, units, pins, runner_sha, opt$expected_r_version)
  checkpoint <- .ranktree_open_checkpoints(output_dir, checkpoint_contract)
  cell_rows <- prediction_rows <- list()
  index <- 0L
  for (method in opt$methods) {
    method_roster <- as.data.frame(roster)
    for (seed in .ranktree_seeds()) {
      seed_dir <- file.path(checkpoint, sprintf("seed_%03d", seed))
      dir.create(seed_dir, recursive = TRUE, showWarnings = FALSE)
      for (cohort_id in names(cohorts)) {
        message(sprintf("[%d] %s / %s", seed, method, cohort_id))
        outer_k <- as.integer(srb_expected_outer_k(cohort_id))
        result <- .ranktree_call_benchmark_cell(
          srb_benchmark_within_cell, method, cohort_id,
          cohorts[[cohort_id]], method_roster, seed_dir, seed, outer_k
        )
        result$outer_k <- outer_k
        cell <- .ranktree_cell_row(result, seed, pins)
        if (!isTRUE(result$eligible)) {
          stop(
            "Rank-tree sensitivity cell failed the Amendment 2 complete-finite ",
            "support contract: ", method, "/", cohort_id, "/seed ", seed,
            ": ", result$ineligible_reason, call. = FALSE
          )
        }
        .ranktree_validate_cell_result(cell)
        index <- index + 1L
        cell_rows[[index]] <- cell
        if (isTRUE(result$eligible)) {
          if (is.null(result$fold_scores)) {
            stop("Eligible cell lacks held-out fold scores.", call. = FALSE)
          }
          prediction <- data.table::as.data.table(result$fold_scores)
          prediction[, `:=`(
            method = method, cohort = cohort_id, seed = as.integer(seed),
            package_version = package_version,
            package_commit = code_start$commit,
            cache_package_commit = opt$expected_cache_commit,
            runtime_image_sha256 = opt$expected_runtime_image_sha256,
            runtime_attestation_manifest_sha256 =
              opt$expected_runtime_attestation_manifest_sha256,
            container_launcher_sha256 =
              opt$expected_container_launcher_sha256,
            r_version = opt$expected_r_version,
            analysis_plan_sha256 = opt$expected_analysis_plan_sha256,
            analysis_amendment_sha256 =
              opt$expected_analysis_amendment_sha256,
            analysis_amendment2_sha256 =
              opt$expected_analysis_amendment2_sha256,
            analysis_amendment3_sha256 =
              opt$expected_analysis_amendment3_sha256,
            provenance_manifest_sha256 =
              opt$expected_provenance_manifest_sha256,
            provenance_union_inventory_sha256 =
              opt$expected_provenance_union_inventory_sha256,
            analysis_code_id = analysis_code_id
          )]
          prediction_rows[[length(prediction_rows) + 1L]] <- prediction
        }
      }
    }
  }
  cells <- data.table::rbindlist(cell_rows, fill = TRUE)
  predictions <- data.table::rbindlist(prediction_rows, fill = TRUE)
  cells <- srb_sanitize_character_fields(cells)
  predictions <- srb_sanitize_character_fields(predictions)
  predictions <- .ranktree_add_split_ids(predictions)
  data.table::setorder(cells, method, cohort, seed)
  data.table::setorder(predictions, method, cohort, seed, fold, sample_id)
  .ranktree_validate_runtime_bundle(
    cells, predictions, units, opt$methods,
    runtime_image_sha256 = opt$expected_runtime_image_sha256,
    runtime_attestation_sha256 =
      opt$expected_runtime_attestation_manifest_sha256,
    container_launcher_sha256 = opt$expected_container_launcher_sha256,
    r_version = opt$expected_r_version)

  bundle <- list(
    cells = cells, predictions = predictions, units = units,
    methods = opt$methods, method_registry = roster[method_id %in% opt$methods],
    seeds = .ranktree_seeds(), cohorts = units$unit_id,
    package_version = package_version, package_commit = code_start$commit,
    cache_package_commit = opt$expected_cache_commit,
    cache_sha256 = opt$expected_cache_sha256,
    cache_manifest_sha256 = opt$expected_cache_manifest_sha256,
    cache_digest_sha256 = opt$expected_cache_digest_sha256,
    runner_script_sha256 = runner_sha,
    engine_manifest_sha256 = opt$expected_engine_manifest_sha256,
    environment_manifest_sha256 = opt$expected_environment_manifest_sha256,
    runtime_image_path = runtime_image,
    runtime_image_sha256 = opt$expected_runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      opt$expected_runtime_attestation_manifest_sha256,
    container_launcher_path = launcher_path,
    container_launcher_sha256 = opt$expected_container_launcher_sha256,
    analysis_plan_path = analysis_plan,
    analysis_plan_sha256 = opt$expected_analysis_plan_sha256,
    analysis_amendment_path = analysis_amendment,
    analysis_amendment_sha256 = opt$expected_analysis_amendment_sha256,
    analysis_amendment2_path = analysis_amendment2,
    analysis_amendment2_sha256 = opt$expected_analysis_amendment2_sha256,
    analysis_amendment3_path = analysis_amendment3,
    analysis_amendment3_sha256 = opt$expected_analysis_amendment3_sha256,
    provenance_script_sha256 = opt$expected_provenance_script_sha256,
    provenance_manifest_sha256 = opt$expected_provenance_manifest_sha256,
    provenance_union_inventory_path = provenance_union_inventory,
    provenance_union_inventory_sha256 =
      opt$expected_provenance_union_inventory_sha256,
    analysis_code_id = analysis_code_id,
    r_version = opt$expected_r_version,
    ranktree_dependency_closure = environment_manifest[, .(package, version,
                                                           source_sha256)],
    provenance_exit_status = provenance$raw_status,
    provenance_overlap_pairs = provenance$n_overlap_pairs
  )
  saveRDS(bundle, file.path(stage, "within_bundle.rds"), version = 3L)
  .public_fwrite(cells, file.path(stage, "within_cells.tsv"), sep = "\t",
                     quote = FALSE, na = "")
  .public_fwrite(predictions,
                     file.path(stage, "within_predictions.tsv.gz"), sep = "\t",
                     quote = FALSE, na = "", compress = "gzip")
  .public_fwrite(units, file.path(stage, "units.tsv"), sep = "\t",
                     quote = FALSE, na = "")
  .public_fwrite(roster[method_id %in% opt$methods],
                     file.path(stage, "methods.tsv"), sep = "\t",
                     quote = FALSE, na = "")
  .ranktree_copy(opt$engine_manifest,
                 file.path(stage, "engine_manifest.tsv"), "engine manifest")
  .ranktree_copy(opt$environment_manifest,
                 file.path(stage, "environment_manifest.tsv"),
                 "environment manifest")
  .ranktree_copy(runtime_attestation$path,
                 file.path(stage, "runtime_attestation.tsv"),
                 "runtime closure attestation")
  .ranktree_copy(analysis_plan, file.path(stage, "analysis_plan.md"),
                 "approved analysis plan")
  .ranktree_copy(analysis_amendment, file.path(stage, "analysis_amendment.md"),
                 "approved analysis amendment 1")
  .ranktree_copy(analysis_amendment2,
                 file.path(stage, "analysis_amendment2.md"),
                 "approved analysis amendment 2")
  .ranktree_copy(analysis_amendment3,
                 file.path(stage, "analysis_amendment3.md"),
                 "approved analysis amendment 3")
  resolution_target <- file.path(stage, "provenance_resolution.tsv")
  if (!identical(normalizePath(provenance$resolution_path, mustWork = TRUE),
                 normalizePath(resolution_target, mustWork = FALSE))) {
    .ranktree_copy(provenance$resolution_path, resolution_target,
                   "provenance resolution")
    unlink(provenance$resolution_path)
  }
  writeLines(c(
    "# External rank-tree sensitivity method bundle", "",
    paste0("Methods: ", paste(opt$methods, collapse = ", ")),
    paste0("Analysis units: ", nrow(units), " (",
           sum(units$primary_unit), " primary)"),
    paste0("Deposited profiles represented: ", sum(units$n_deposited_profiles)),
    paste0("Deduplicated biological/provenance groups: ",
           sum(units$n_unique_groups)),
    paste0("Singularity runtime image SHA-256: ",
           opt$expected_runtime_image_sha256),
    paste0("Exact in-process R runtime/library attestation SHA-256: ",
           opt$expected_runtime_attestation_manifest_sha256),
    paste0("Package container launcher SHA-256: ",
           opt$expected_container_launcher_sha256),
    paste0("Canonical provenance union inventory SHA-256: ",
           opt$expected_provenance_union_inventory_sha256),
    paste0("Approved analysis plan SHA-256: ",
           opt$expected_analysis_plan_sha256),
    paste0("Approved analysis amendment 1 SHA-256: ",
           opt$expected_analysis_amendment_sha256),
    paste0("Approved analysis amendment 2 SHA-256: ",
           opt$expected_analysis_amendment2_sha256),
    paste0("Approved analysis amendment 3 SHA-256: ",
           opt$expected_analysis_amendment3_sha256),
    "Rank forest contract: unweighted preliminary screen; inverse-frequency ",
    "pair-forest case weights; numeric y=1 (`case`) is the positive target.",
    "Performance is computed as group-collapsed held-out AUC: one mean score ",
    "per provenance group within each exact test fold.",
    "", "This is an unreviewed execution bundle. It is not reportable until ",
    "assembly, independent results review, plausibility assessment and promotion."
  ), file.path(stage, "report.md"))

  # Recheck all mutable external inputs and package state before promotion.
  .ranktree_validate_runtime_attestation(
    runtime_attestation$path,
    opt$expected_runtime_attestation_manifest_sha256,
    opt$expected_r_version, opt$expected_runtime_image_sha256,
    opt$expected_container_launcher_sha256
  )
  end_pins <- c(
    .ranktree_sha(cache), .ranktree_sha(cache_manifest_path),
    .ranktree_sha(cache_digest_path), .ranktree_sha(opt$engine_manifest),
    .ranktree_sha(analysis_plan),
    .ranktree_sha(analysis_amendment),
    .ranktree_sha(analysis_amendment2),
    .ranktree_sha(analysis_amendment3),
    .ranktree_sha(opt$environment_manifest), .ranktree_sha(runtime_image),
    .ranktree_sha(runtime_attestation$path), .ranktree_sha(launcher_path),
    .ranktree_sha(provenance_script),
    .ranktree_sha(provenance_manifest),
    .ranktree_sha(provenance_union_inventory)
  )
  expected_end <- c(
    opt$expected_cache_sha256, opt$expected_cache_manifest_sha256,
    opt$expected_cache_digest_sha256, opt$expected_engine_manifest_sha256,
    opt$expected_analysis_plan_sha256,
    opt$expected_analysis_amendment_sha256,
    opt$expected_analysis_amendment2_sha256,
    opt$expected_analysis_amendment3_sha256,
    opt$expected_environment_manifest_sha256,
    opt$expected_runtime_image_sha256,
    opt$expected_runtime_attestation_manifest_sha256,
    opt$expected_container_launcher_sha256,
    opt$expected_provenance_script_sha256,
    opt$expected_provenance_manifest_sha256,
    opt$expected_provenance_union_inventory_sha256
  )
  code_end <- .ranktree_git_provenance(opt$package_root)
  if (!identical(end_pins, expected_end) || !identical(code_end, code_start) ||
      !identical(.ranktree_sha(script_path), runner_sha)) {
    stop("A pinned input or OmicSelector code provenance changed during ",
         "execution.", call. = FALSE)
  }
  artifacts <- list.files(stage, full.names = TRUE)
  manifest <- data.table::data.table(
    file = basename(artifacts), bytes = file.info(artifacts)$size,
    sha256 = vapply(artifacts, .ranktree_sha, character(1L)),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_version = package_version, package_commit = code_start$commit,
    runner_script_sha256 = runner_sha, analysis_code_id = analysis_code_id,
    engine_manifest_sha256 = opt$expected_engine_manifest_sha256,
    environment_manifest_sha256 = opt$expected_environment_manifest_sha256,
    analysis_plan_sha256 = opt$expected_analysis_plan_sha256,
    analysis_amendment_sha256 = opt$expected_analysis_amendment_sha256,
    analysis_amendment2_sha256 = opt$expected_analysis_amendment2_sha256,
    analysis_amendment3_sha256 = opt$expected_analysis_amendment3_sha256,
    runtime_image_sha256 = opt$expected_runtime_image_sha256,
    runtime_attestation_manifest_sha256 =
      opt$expected_runtime_attestation_manifest_sha256,
    container_launcher_sha256 = opt$expected_container_launcher_sha256,
    r_version = opt$expected_r_version,
    provenance_manifest_sha256 = opt$expected_provenance_manifest_sha256,
    provenance_union_inventory_sha256 =
      opt$expected_provenance_union_inventory_sha256
  )
  .public_fwrite(manifest, file.path(stage, "output_manifest.tsv"),
                     sep = "\t", quote = FALSE)
  if (!file.rename(stage, output_dir)) {
    stop("Could not atomically promote rank-tree method bundle.", call. = FALSE)
  }
  .ranktree_remove_completed_checkpoints(output_dir, checkpoint_contract)
  cat("Created external-competitor method bundle:", output_dir, "\n")
  invisible(output_dir)
}

if (sys.nframe() == 0L) {
  trailing <- commandArgs(trailingOnly = TRUE)
  capture_index <- match("--write-runtime-attestation", trailing)
  if (!is.na(capture_index)) {
    if (capture_index == length(trailing) ||
        !nzchar(trailing[[capture_index + 1L]])) {
      stop("--write-runtime-attestation requires a target path.",
           call. = FALSE)
    }
    .ranktree_write_runtime_attestation(trailing[[capture_index + 1L]])
  } else {
    .ranktree_runner_main(trailing)
  }
}
