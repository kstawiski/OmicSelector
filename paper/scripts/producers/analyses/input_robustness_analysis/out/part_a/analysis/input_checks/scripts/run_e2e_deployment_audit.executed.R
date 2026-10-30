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
# End-to-end single-specimen deployment audit (strengthening).
# Retrospectively specified extension; exploratory. Frozen parameters live in
# frozen_parameters.tsv (sha256 checked below) and were written before any
# held-out specimen was preprocessed or scored.
#
# For each cohort, split, and A/K route, one fitted object is frozen. Held-out
# specimens then pass through preprocessing + score_specimen() along several
# frozen training-state preprocessing). Every path is compared with the
# benchmark batch path (B_ref). Required-input removal scenarios test whether
# the scorer rejects explicitly or returns a silent score.
#
# Run from this directory: TMPDIR="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript --vanilla run_e2e_deployment_audit.R
# Writes only inside this directory.

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
.libPaths(c("./runtime/reference-home/R/omicselector", "./runtime/reference-home/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  source('./analyses/input_robustness_analysis/out/part_a/environment.R')
})
data.table::setDTthreads(1L)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
ARGS <- commandArgs(trailingOnly = TRUE)
SMOKE <- "--smoke" %in% ARGS
RESUME <- "--resume" %in% ARGS
ROOT <- "."
AUDIT_DIR <- file.path(ROOT, "analysis/opportunity/strengthening/e2e_deployment_audit")
OUT <- './analyses/input_robustness_analysis/out/part_a/analysis/input_checks'
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PKG <- "./software/OmicSelector-267"
FROZEN_SHA <- "f8fba6b1833a1d8016a6dcaa08b2fddc56934303089f243254c51548c2e60371"
t_start <- Sys.time()

fp <- file.path(AUDIT_DIR, "frozen_parameters.tsv")
if (!identical(digest(fp, file = TRUE, algo = "sha256"), FROZEN_SHA)) {
  stop("frozen_parameters.tsv hash mismatch; tolerances must not change after the freeze")
}
frozen <- fread(fp, sep = "\t", quote = "")
getp <- function(k) frozen$value[frozen$parameter == k]
TOL <- as.numeric(getp("tol_score_abs"))
N_SPLITS <- as.integer(getp("splits"))
N_SHUF <- as.integer(getp("shuffles_per_split"))
ROUTES <- strsplit(getp("routes"), ";", fixed = TRUE)[[1]]
NONDEP_A <- strsplit(getp("nondeployable_A_routes_probe"), ";", fixed = TRUE)[[1]]
BASE_SEED <- 20260930L
if (SMOKE) { N_SPLITS <- 1L; N_SHUF <- 1L }

log_con <- file(file.path(OUT, "run_log.txt"), open = "at")
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), " ", ...)
  cat(msg, "\n"); writeLines(msg, log_con); flush(log_con)
}
say("Sequential cohort analysis; pkgload source load; OMP_NUM_THREADS=", Sys.getenv("OMP_NUM_THREADS"))
say("frozen_parameters sha256 verified: ", FROZEN_SHA)

pkg_desc <- packageDescription("OmicSelector")
say("OmicSelector ", pkg_desc$Version, " RemoteSha=", pkg_desc$RemoteSha %||% "NA",
    " lib=", find.package("OmicSelector"))
src_head <- trimws(system2("git", c("-C", PKG, "rev-parse", "HEAD"), stdout = TRUE))
say("OmicSelector source HEAD: ", src_head)

# ---------------------------------------------------------------------------
# Benchmark helpers (read-only sources)
# ---------------------------------------------------------------------------
setwd(ROOT)
source(file.path(ROOT, "code/methods/mirna_alias_resolver.R"))
# Cache only the immutable feature-name lookup, never expression values or scores.
resolve_mirna_aliases <- local({
  original <- resolve_mirna_aliases
  memo <- new.env(parent = emptyenv())
  function(...) {
    args <- list(...)
    key <- digest(args, algo = "sha256")
    if (exists(key, envir = memo, inherits = FALSE)) return(get(key, envir = memo))
    value <- do.call(original, args)
    stopifnot(identical(value, do.call(original, args)))
    assign(key, value, envir = memo)
    value
  }
})
say("Feature-name lookup memoized; each new key checked against repeated original lookup; no expression or score caching")

source(file.path(ROOT, "code/methods/outcome_dictionary.R"))
inv <- fread(file.path(ROOT, "knowledge/union_inventory.tsv"), sep = "\t", quote = "", fill = TRUE)
dict <- load_outcome_dictionary_v0_4()
source('./analyses/input_robustness_analysis/out/part_a/scripts/cohort_loader.R')
builder <- file.path(ROOT, "code/analyses/build_singlesample_public_cache.R")
say("cache builder sha256: ", digest(builder, file = TRUE, algo = "sha256"))
cache_path <- file.path(ROOT, "analysis/singlesample_public_benchmark/cohort_cache.rds")
say("cohort cache sha256: ", digest(cache_path, file = TRUE, algo = "sha256"))
cache <- .public_readRDS(cache_path)

# Verbatim copy of canonical_abundance() from the cache builder (lines 66-85).
canonical_abundance <- function(X, scale_state = c("abundance", "log"),
                                log_base = 2, ct = FALSE,
                                nondetect_ct = 40) {
  scale_state <- match.arg(scale_state)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  colnames(X) <- .normalize_feature_ids(colnames(X))
  X <- resolve_and_collapse_features(X, scale_state = scale_state)
  X <- normalize_to_canonical_mirna_namespace(X, scale_state = scale_state)
  X <- X[, colSums(is.finite(X)) > 0L, drop = FALSE]
  if (ct) {
    X[!is.finite(X)] <- nondetect_ct
    X <- log_base ^ (40 - X)
  } else if (scale_state == "log") {
    X <- log_base ^ X
  }
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0
  X
}

# Verbatim copy of feature_by_sample() from the cache builder.
feature_by_sample <- function(d, feature_col = 1L, sample_cols = NULL) {
  if (is.null(sample_cols)) sample_cols <- setdiff(seq_len(ncol(d)), feature_col)
  feature <- as.character(d[[feature_col]])
  X <- t(as.matrix(d[, ..sample_cols]))
  storage.mode(X) <- "double"
  rownames(X) <- names(d)[sample_cols]
  colnames(X) <- feature
  X
}

as_double <- function(X) { X <- as.matrix(X); suppressWarnings(storage.mode(X) <- "double"); X }

# Per-specimen canonicalisation shared by benchmark and frozen paths.
canon_rows <- function(X, scale_state) {
  X <- as_double(X)
  colnames(X) <- .normalize_feature_ids(colnames(X))
  X <- resolve_and_collapse_features(X, scale_state = scale_state)
  X <- normalize_to_canonical_mirna_namespace(X, scale_state = scale_state)
  X
}
restrict_universe <- function(X, universe) {
  miss <- setdiff(universe, colnames(X))
  if (length(miss)) {
    stop("frozen preprocessing: ", length(miss),
         " frozen universe feature(s) absent from the specimen record")
  }
  X[, universe, drop = FALSE]
}
ct_to_abundance <- function(X) {
  X[!is.finite(X)] <- 40
  X <- 2 ^ (40 - X)
  X[!is.finite(X)] <- 0
  X[X < 0] <- 0
  X
}

# ---------------------------------------------------------------------------
# Cohort definitions: raw start, benchmark preprocessing, frozen preprocessing
# ---------------------------------------------------------------------------
make_ct40 <- function(id, raw, y, start, note) {
  list(
    id = id, kind = "ct40", raw = raw, y = y, start = start, note = note,
    bench_pre = function(R) canonical_abundance(R, "log", log_base = 2, ct = TRUE),
    state_from = function(R) NULL,
    frozen_pre = function(R, universe, state) {
      ct_to_abundance(restrict_universe(canon_rows(R, "log"), universe))
    },
    companion = "ct"
  )
}

read_gse64591 <- function() {
  expr_path <- file.path(ROOT, "data_bucket/GSE64591/raw/suppl/GSE64591_non-normalized.txt.gz")
  meta_path <- file.path(ROOT, "data_bucket/GSE64591/metadata.tsv")
  d <- fread(expr_path, skip = 5L, check.names = FALSE)
  feature <- d[[1L]]
  X <- t(as.matrix(d[, -1L]))
  colnames(X) <- feature
  rownames(X) <- colnames(d)[-1L]
  meta <- fread(meta_path)
  idx <- match(rownames(X), meta$title)
  if (anyNA(idx)) stop("GSE64591: title-to-GSM join failed")
  y <- as.integer(meta$sample_type[idx] == "cancer")
  rownames(X) <- meta$gsm_id[idx]
  make_ct40("GSE64591", as_double(X), setNames(y, rownames(X)), "raw_ct_non_normalized",
            "TaqMan TLDA raw Ct; non-detect -> Ct 40 per specimen")
}

read_gse41922 <- function() {
  p <- file.path(ROOT, "data_bucket/GSE41922/GSE41922_non-normalized.txt.gz")
  d <- fread(p, skip = 5L, check.names = FALSE)
  X <- feature_by_sample(d)
  y <- as.integer(startsWith(rownames(X), "s"))
  make_ct40("GSE41922", X, setNames(y, rownames(X)), "raw_ct_non_normalized",
            "Exiqon raw Ct; non-detect -> Ct 40 per specimen")
}

read_gse65071 <- function() {
  p <- file.path(ROOT, "data_bucket/GSE65071/GSE65071_human_raw_data.txt.gz")
  soft <- file.path(ROOT, "data_bucket/GSE65071/GSE65071_family.soft.gz")
  d <- fread(p, check.names = FALSE)
  lines <- readLines(gzfile(soft), warn = FALSE)
  b <- grep("^!platform_table_begin", lines)
  e <- grep("^!platform_table_end", lines)
  platform <- fread(text = paste(lines[(b + 1L):(e - 1L)], collapse = "\n"),
                    sep = "\t", quote = "", fill = TRUE)
  idx <- match(d[[1L]], platform$REF_ID)
  keep <- !is.na(idx) & nzchar(platform$miRNA_ID[idx])
  d <- d[keep]
  X <- feature_by_sample(d)
  colnames(X) <- platform$miRNA_ID[idx[keep]]
  y <- as.integer(startsWith(rownames(X), "Disease"))
  make_ct40("GSE65071", X, setNames(y, rownames(X)), "raw_ct",
            "Exiqon raw Ct; non-detect -> Ct 40 per specimen")
}

read_generic <- function(acc, modality, start, note) {
  ref <- .load_cohort(acc, modality)
  if (is.null(ref)) stop(acc, ": generic loader failed")
  inv_row <- inv[accession == acc][1L]
  path <- ref$expression_path
  raw <- .read_raw_expression_path(path, acc)
  scale_attr <- tolower(trimws(as.character(inv_row$scale %||% "")))
  scale_state <- if (scale_attr %in% c("count", "rpm", "abundance")) "abundance" else "log"
  keep_attrs <- function(X, like) {
    for (nm in c("scale", "log_base", "transform_direction", "expression_path")) {
      attr(X, nm) <- attr(like, nm)
    }
    X
  }
  template <- .postprocess_expression(raw[1:2, , drop = FALSE], inv_row, path)
  list(
    id = acc, kind = "median_ct", raw = as_double(raw), y = setNames(ref$y, rownames(ref$X)),
    start = start, note = note, expression_path = path,
    bench_pre = function(R) .to_abundance(.postprocess_expression(R, inv_row, path)),
    state_from = function(R) {
      C <- canon_rows(R, scale_state)
      med <- vapply(colnames(C), function(j) {
        m <- suppressWarnings(stats::median(C[, j], na.rm = TRUE))
        if (is.finite(m)) m else NA_real_
      }, numeric(1L))
      med
    },
    frozen_pre = function(R, universe, state) {
      X <- restrict_universe(canon_rows(R, scale_state), universe)
      for (j in universe) {
        bad <- !is.finite(X[, j])
        if (any(bad)) {
          v <- state[[j]]
          X[bad, j] <- if (is.null(v) || !is.finite(v)) 40 else v
        }
      }
      .to_abundance(keep_attrs(X, template))
    },
    companion = "ct"
  )
}

read_gse70080_well <- function(y70080) {
  base <- file.path(ROOT, "data_bucket/GSE70080/raw/suppl")
  rd <- function(f, card) {
    L <- readLines(gzfile(file.path(base, f)), warn = FALSE)
    h <- grep("^ID_REF\t", L)[1L]
    d <- fread(text = paste(L[h:length(L)], collapse = "\n"), header = TRUE,
               na.strings = c("null", ""))
    d <- d[, grepl("^ID_REF$|^SAMPLE", names(d)), with = FALSE]
    m <- t(as.matrix(d[, setdiff(names(d), "ID_REF"), with = FALSE]))
    storage.mode(m) <- "double"
    colnames(m) <- paste0(card, "_", d$ID_REF)
    m
  }
  R <- cbind(rd("GSE70080_cardA_non-normalized.txt.gz", "A"),
             rd("GSE70080_cardB_non-normalized.txt.gz", "B"))
  sm <- readLines(gzfile(file.path(ROOT, "data_bucket/GSE70080/raw/matrix/GSE70080_series_matrix.txt.gz")),
                  warn = FALSE)
  gsm <- gsub('"', "", strsplit(sm[startsWith(sm, "!Sample_geo_accession")], "\t")[[1]][-1])
  tit <- gsub('"', "", strsplit(sm[startsWith(sm, "!Sample_title")], "\t")[[1]][-1])
  # Raw 'SAMPLE n' labels: 1-16 and 18-71 align in order with the 70 non-calibrator
  # GEO samples (verified value-by-value between the deposited normalized file and
  # the series matrix); 'SAMPLE 17' is absent from the normalized file and is
  # treated as the presumed calibrator, not analysed.
  lab <- paste("SAMPLE", c(1:16, 18:71))
  noncal <- gsm[tit != "kalibrator"]
  if (length(noncal) != 70L || !all(lab %in% rownames(R))) stop("GSE70080 well map failed")
  R <- R[c(lab, "SAMPLE 17"), , drop = FALSE]
  rownames(R) <- c(noncal, "PRESUMED_CALIBRATOR_SAMPLE_17")
  y <- y70080[intersect(names(y70080), rownames(R))]
  list(
    id = "GSE70080_RAWCT_WELL", kind = "ct40_well", raw = R, y = y,
    start = "raw_ct_non_normalized_well_keyed",
    note = "Raw card A/B Ct keyed by card well (no well-to-assay map deposited); benchmark ct40 convention restarted from raw Ct",
    bench_pre = function(X) {
      X <- X[, colSums(is.finite(X)) > 0L, drop = FALSE]
      ct_to_abundance(X)
    },
    state_from = function(R) NULL,
    frozen_pre = function(R, universe, state) ct_to_abundance(restrict_universe(R, universe)),
    companion = "ct"
  )
}

read_cache_abundance <- function(acc, note) {
  e <- cache[[acc]]
  X <- e$expr_per_sample
  list(
    id = acc, kind = "cache_abundance", raw = X, y = setNames(e$y_bin, rownames(X)),
    start = "benchmark_cache_after_submitter_normalization", note = note,
    bench_pre = function(R) R,
    state_from = function(R) NULL,
    frozen_pre = function(R, universe, state) restrict_universe(R, universe),
    companion = "abundance"
  )
}

say("reading cohorts")
cohorts <- list(
  GSE64591 = read_gse64591(),
  GSE41922 = read_gse41922(),
  GSE65071 = read_gse65071(),
  GSE279209 = read_generic("GSE279209", "OpenArray qPCR", "raw_crt_non_normalized",
                           "OpenArray raw Crt; benchmark imputes non-detects with cohort medians"),
  GSE70080 = read_generic("GSE70080", "TLDA qPCR", "submitter_normalized_ddct_series_matrix",
                          "Series-matrix values: ExpressionSuite global-normalized delta-delta Ct; benchmark imputes nulls with cohort medians"),
  GSE188232 = read_cache_abundance("GSE188232",
                                   "Submitter-normalized log2 values converted to abundance in the benchmark cache; raw counts not deposited")
)
cohorts$GSE70080_RAWCT_WELL <- read_gse70080_well(cohorts$GSE70080$y)

# ---------------------------------------------------------------------------
# Reproduce the benchmark cache (batch path) before any held-out scoring
# ---------------------------------------------------------------------------
repro <- list()
for (cid in names(cohorts)) {
  co <- cohorts[[cid]]
  B <- co$bench_pre(co$raw)
  ids <- names(co$y)
  B <- B[ids, , drop = FALSE]
  co$B_all <- B
  co$universe <- colnames(B)
  if (cid %in% names(cache)) {
    ce <- cache[[cid]]
    same_rows <- identical(sort(rownames(ce$expr_per_sample)), sort(ids))
    same_cols <- identical(sort(colnames(ce$expr_per_sample)), sort(colnames(B)))
    maxd <- if (same_rows && same_cols) {
      max(abs(B[rownames(ce$expr_per_sample), colnames(ce$expr_per_sample)] - ce$expr_per_sample) /
            pmax(1, abs(ce$expr_per_sample)))
    } else NA_real_
    same_y <- identical(as.integer(co$y[rownames(ce$expr_per_sample)]), as.integer(ce$y_bin))
    repro[[cid]] <- data.table(cohort = cid, rows_match = same_rows, cols_match = same_cols,
                               labels_match = same_y, max_rel_diff_vs_cache = maxd)
    if (!isTRUE(same_rows && same_cols && same_y && maxd == 0)) {
      stop(cid, ": benchmark batch path does not reproduce the cohort cache exactly")
    }
  } else {
    repro[[cid]] <- data.table(cohort = cid, rows_match = NA, cols_match = NA,
                               labels_match = NA, max_rel_diff_vs_cache = NA_real_)
  }
  cohorts[[cid]] <- co
  say(cid, ": deposited rows in preprocessing=", nrow(co$raw), " analysed=", length(ids),
      " universe=", ncol(B), " cases/controls=", sum(co$y == 1L), "/", sum(co$y == 0L))
}
repro <- rbindlist(repro)
say("cache reproduction:\n", paste(capture.output(print(repro)), collapse = "\n"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
score_try <- function(dep, X) {
  warn <- character(0)
  res <- withCallingHandlers(
    tryCatch(list(s = score_specimen(dep, X), err = NA_character_),
             error = function(e) list(s = rep(NA_real_, NROW(X)), err = conditionMessage(e))),
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  res$warn <- if (length(warn)) warn[[1L]] else NA_character_
  res
}

youden_cutoff <- function(s, y) {
  u <- sort(unique(s))
  cand <- c(min(u) - 1, if (length(u) > 1L) (u[-1L] + u[-length(u)]) / 2, max(u) + 1)
  J <- vapply(cand, function(t) mean(s[y == 1L] >= t) + mean(s[y == 0L] < t) - 1, numeric(1L))
  best <- which(J >= max(J) - 1e-12)
  cand[best[ceiling(length(best) / 2)]]
}

required_features <- function(dep) {
  m <- dep$model
  switch(dep$method_id,
    "ws-rclr-trimmed" = m$panel,
    "ws-balance-ilr" = intersect(m$sbp_canonical_features, m$feature_universe),
    "reo-ktsp" = unique(c(m$ktsp$pairs$feature_a, m$ktsp$pairs$feature_b)),
    "reo-pairratio" = unique(c(m$selected_pairs$feature_a, m$selected_pairs$feature_b)),
    character(0))
}

companions_for <- function(co, target_raw, c4_raw) {
  if (co$companion == "ct") {
    C1 <- target_raw; C1[] <- NA_real_
    C2 <- target_raw - 20
    C3 <- target_raw + 10
  } else {
    C1 <- target_raw; C1[] <- 0
    C2 <- target_raw * 2^20
    C3 <- target_raw * 2^-10
  }
  out <- rbind(C1, C2, C3, c4_raw)
  rownames(out) <- c("C1_all_nondetect", "C2_high", "C3_low", "C4_real_most_missing")
  out
}

missing_count <- function(co, R) {
  if (co$companion == "ct") rowSums(!is.finite(R)) else rowSums(!is.finite(R) | R <= 0)
}

stratified_split <- function(y, seed) {
  set.seed(seed)
  test <- unlist(lapply(split(names(y), y), function(ids) {
    sample(ids, round(0.3 * length(ids)))
  }), use.names = FALSE)
  list(train = setdiff(names(y), test), test = sort(test))
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
# Cohorts run sequentially in this R process. Every random draw retains its
# frozen split/shuffle seed. Progress is written and flushed immediately.
run_cohort <- function(cid) {
  rec <- list()        # specimen-level path records
  rj <- list()         # rejection-scenario records
  fits <- list()       # fit-level records
  naud <- list()       # preprocessing dependence numbers
  msgs <- character(0)
  say <- function(...) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), " ", ...)
    cat(msg, "\n"); writeLines(msg, log_con); flush(log_con)
    msgs <<- c(msgs, msg)
  }
  # Single-row benchmark preprocessing depends only on the specimen: cache it.
  S1_cache <- list()
  co <- cohorts[[cid]]
  y <- co$y
  routes <- ROUTES
  if (co$kind == "ct40_well") routes <- setdiff(routes, "ws-balance-ilr")
  S1_cache[[cid]] <- list()
  state_all <- co$state_from(co$raw)                     # replayed benchmark state (all deposited rows)

  for (r in seq_len(N_SPLITS)) {
    seed <- BASE_SEED + r
    sp <- stratified_split(y, seed)
    Xtr <- co$B_all[sp$train, , drop = FALSE]
    ytr <- as.integer(y[sp$train])
    test <- sp$test
    if (SMOKE) test <- test[seq_len(min(4L, length(test)))]
    state_tr <- co$state_from(co$raw[sp$train, , drop = FALSE])

    # preprocessing-dependence audit numbers for this split
    tr_raw_can <- if (co$kind %in% c("ct40", "median_ct")) {
      canon_rows(co$raw[sp$train, , drop = FALSE], "log")
    } else co$raw[sp$train, , drop = FALSE]
    tr_det <- if (co$companion == "ct") colSums(is.finite(tr_raw_can)) else colSums(tr_raw_can > 0)
    only_other <- sum(co$universe %in% names(tr_det)[tr_det == 0] | !(co$universe %in% names(tr_det)))
    te_can <- if (co$kind %in% c("ct40", "median_ct")) canon_rows(co$raw[test, , drop = FALSE], "log") else co$raw[test, , drop = FALSE]
    te_u <- te_can[, intersect(co$universe, colnames(te_can)), drop = FALSE]
    te_nd <- if (co$companion == "ct") !is.finite(te_u) else (te_u <= 0)
    naud[[paste(cid, r)]] <- data.table(
      cohort = cid, split = r, n_train = length(sp$train), n_test = length(test),
      n_deposited_rows_in_batch_preprocessing = nrow(co$raw),
      n_rows_outside_analysis_set = nrow(co$raw) - length(y),
      universe_size = length(co$universe),
      universe_features_undetected_in_all_training_rows = only_other,
      heldout_cells_nondetect_in_universe = sum(te_nd),
      heldout_rows_with_any_nondetect = sum(rowSums(te_nd) > 0),
      nondetect_fill = switch(co$kind, ct40 = "fixed Ct 40 (per specimen)",
                              ct40_well = "fixed Ct 40 (per specimen)",
                              median_ct = "cohort median over all deposited rows (cross-specimen)",
                              cache_abundance = "none (submitter-normalized start)")
    )

    # companions and C4 choice (specified per split)
    mc <- missing_count(co, co$raw[test, , drop = FALSE])
    c4_order <- test[order(-mc, test)]

    # B_shuffled preprocessing (route-independent)
    shuf <- lapply(seq_len(N_SHUF), function(s) {
      set.seed(seed * 10L + s)
      perm <- sample.int(nrow(co$raw))
      Bp <- co$bench_pre(co$raw[perm, , drop = FALSE])
      ord <- rownames(co$raw)[perm]
      ord <- ord[ord %in% test]
      Bp[ord, , drop = FALSE]
    })
    B_tb <- co$bench_pre(co$raw[test, , drop = FALSE])

    # single-row benchmark preprocessing (cached by specimen)
    for (id in test) {
      if (is.null(S1_cache[[cid]][[id]])) {
        S1_cache[[cid]][[id]] <- tryCatch(co$bench_pre(co$raw[id, , drop = FALSE]),
                                          error = function(e) e)
      }
    }
    # frozen training-state and replayed-state single-row preprocessing
    S2_pre <- lapply(test, function(id) tryCatch(co$frozen_pre(co$raw[id, , drop = FALSE], co$universe, state_tr), error = function(e) e))
    S3_pre <- lapply(test, function(id) tryCatch(co$frozen_pre(co$raw[id, , drop = FALSE], co$universe, state_all), error = function(e) e))
    names(S2_pre) <- names(S3_pre) <- test
    # companion batches
    Xb_pre <- list(); Xf_pre <- list()
    for (id in test) {
      c4 <- setdiff(c4_order, id)[1L]
      batch <- rbind(co$raw[id, , drop = FALSE], companions_for(co, co$raw[id, , drop = FALSE], co$raw[c4, , drop = FALSE]))
      rownames(batch)[1L] <- id
      Xb_pre[[id]] <- tryCatch(co$bench_pre(batch), error = function(e) e)
      Xf_pre[[id]] <- tryCatch(co$frozen_pre(batch, co$universe, state_tr), error = function(e) e)
    }

    for (m in routes) {
      t0 <- Sys.time()
      dep <- tryCatch(suppressWarnings(deploy_singlesample(Xtr, ytr, method = m, verify = TRUE)),
                      error = function(e) e)
      if (inherits(dep, "error")) {
        fits[[paste(cid, r, m)]] <- data.table(cohort = cid, split = r, route = m, fit_ok = FALSE,
                                               fit_error = conditionMessage(dep))
        say(cid, " split ", r, " ", m, ": FIT FAILED: ", conditionMessage(dep))
        next
      }
      dig0 <- digest(serialize(dep$model, NULL), algo = "sha256")
      s_tr <- score_specimen(dep, Xtr)
      cut <- youden_cutoff(s_tr, ytr)
      gate <- is_singlesample_deployable(dep, X_probe = co$B_all[test, , drop = FALSE])
      req <- required_features(dep)

      ref <- score_try(dep, co$B_all[test, , drop = FALSE])
      if (!is.na(ref$err)) stop(cid, " ", m, ": reference batch scoring failed: ", ref$err)
      s_ref <- setNames(ref$s, test)

      add <- function(path, comparator, ids, s, err, warn, s_cmp = s_ref[ids], extra = NA_character_) {
        rec[[length(rec) + 1L]] <<- data.table(
          cohort = cid, split = r, route = m, path = path, comparator = comparator,
          specimen = ids, y = as.integer(y[ids]), s_cmp = s_cmp, s = s, err = err, warn = warn,
          cutoff = cut, detail = extra)
      }

      # B_shuffled
      for (s in seq_len(N_SHUF)) {
        Xs <- shuf[[s]]
        o <- score_try(dep, Xs)
        add("B_shuffled", "B_ref", rownames(Xs), o$s, o$err, o$warn, extra = paste0("shuffle_", s))
      }
      # B_testbatch
      o <- score_try(dep, B_tb[test, , drop = FALSE])
      add("B_testbatch", "B_ref", test, o$s, o$err, o$warn)
      # matched_baseline / companion_dependence / reassay_design single-row paths
      for (path in c("S1_single_benchmark_preproc", "S2_single_frozen_training_state",
                     "S3_single_frozen_replayed_state")) {
        pre <- switch(path,
                      S1_single_benchmark_preproc = S1_cache[[cid]][test],
                      S2_single_frozen_training_state = S2_pre,
                      S3_single_frozen_replayed_state = S3_pre)
        for (id in test) {
          p <- pre[[id]]
          if (inherits(p, "error")) {
            add(path, "B_ref", id, NA_real_, paste0("preprocessing: ", conditionMessage(p)), NA_character_)
            next
          }
          miss_req <- length(setdiff(req, colnames(p)))
          o <- score_try(dep, p)
          add(path, "B_ref", id, o$s, o$err, o$warn,
              extra = paste0("features_in=", ncol(p), ";required_missing=", miss_req))
        }
      }
      # companion paths (score the whole companion batch; read the target row)
      for (id in test) {
        for (path in c("X_companions_benchmark_preproc", "X_companions_frozen_preproc")) {
          p <- if (path == "X_companions_benchmark_preproc") Xb_pre[[id]] else Xf_pre[[id]]
          if (inherits(p, "error")) {
            add(path, "B_ref", id, NA_real_, paste0("preprocessing: ", conditionMessage(p)), NA_character_)
            next
          }
          o <- score_try(dep, p)
          st <- if (is.na(o$err)) o$s[match(id, rownames(p))] else NA_real_
          add(path, "B_ref", id, st, o$err, o$warn)
          if (path == "X_companions_frozen_preproc") {
            add("X_companions_frozen_preproc", "S2_single_frozen_training_state", id, st, o$err, o$warn,
                s_cmp = NA_real_)
          }
        }
      }

      # ---------------- rejection scenarios (baseline = B_ref row) ----------------
      addrj <- function(scn, id, x_or_err, base_s) {
        if (inherits(x_or_err, "error")) {
          o <- list(s = NA_real_, err = paste0("input construction: ", conditionMessage(x_or_err)), warn = NA_character_)
        } else o <- score_try(dep, x_or_err)
        rj[[length(rj) + 1L]] <<- data.table(
          cohort = cid, split = r, route = m, scenario = scn, specimen = id,
          n_required = length(req), s_base = base_s, s = o$s[1L], err = o$err, warn = o$warn, cutoff = cut)
      }
      nonreq <- setdiff(co$universe, req)
      drop_half <- nonreq[seq(1L, length(nonreq), by = 2L)]
      for (id in test) {
        x0 <- co$B_all[id, , drop = FALSE]
        b0 <- s_ref[[id]]
        if (length(req)) {
          addrj("RJ1_all_required_absent", id, x0[, setdiff(colnames(x0), req), drop = FALSE], b0)
          addrj("RJ2_one_required_absent", id, x0[, setdiff(colnames(x0), req[1L]), drop = FALSE], b0)
          x3 <- x0; x3[, req[1L]] <- NA_real_
          addrj("RJ3_one_required_nonfinite", id, x3, b0)
        }
        if (length(drop_half)) {
          addrj("RJ4_half_nonrequired_absent", id, x0[, setdiff(colnames(x0), drop_half), drop = FALSE], b0)
        }
        ext <- matrix(stats::median(x0) * seq_len(10L), nrow = 1L,
                      dimnames = list(id, paste0("UNREGISTERED_FEATURE_", seq_len(10L))))
        addrj("RJ5_extra_unregistered_features", id, cbind(x0, ext), b0)
        blank <- co$raw[id, , drop = FALSE]
        blank[] <- if (co$companion == "ct") NA_real_ else 0
        addrj("RJ6_all_nondetect_specimen", id,
              tryCatch(co$frozen_pre(blank, co$universe, state_tr), error = function(e) e), b0)
        addrj("RJ7_unnamed_vector", id, as.numeric(x0), b0)
      }

      dig1 <- digest(serialize(dep$model, NULL), algo = "sha256")
      fits[[paste(cid, r, m)]] <- data.table(
        cohort = cid, split = r, route = m, fit_ok = TRUE, fit_error = NA_character_,
        n_train = nrow(Xtr), n_test = length(test), deploy_verified = isTRUE(dep$verified),
        verify_n = dep$verify_n, pkg_gate_heldout = isTRUE(gate), model_state_unchanged = identical(dig0, dig1),
        n_required = length(req), cutoff = cut,
        fit_status = if (!is.null(dep$model$fit_status)) dep$model$fit_status else NA_character_,
        seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
      say(cid, " split ", r, " ", m, ": gate=", isTRUE(gate), " required=", length(req),
          " cut=", signif(cut, 6), " t=", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s")
    }

    # non-deployable A routes: explicit rejection at deploy time?
    for (m in NONDEP_A) {
      o <- tryCatch({ suppressWarnings(deploy_singlesample(Xtr, ytr, method = m, verify = TRUE)); NA_character_ },
                    error = function(e) conditionMessage(e))
      rj[[length(rj) + 1L]] <- data.table(
        cohort = cid, split = r, route = m, scenario = "RJ9_nondeployable_A_route_deploy", specimen = NA_character_,
        n_required = NA_integer_, s_base = NA_real_, s = if (is.na(o)) 0 else NA_real_, err = o, warn = NA_character_,
        cutoff = NA_real_)
    }
  }
  list(rec = rbindlist(rec, fill = TRUE), rj = rbindlist(rj, fill = TRUE),
       fits = rbindlist(fits, fill = TRUE), naud = rbindlist(naud, fill = TRUE), msgs = msgs)
}

say("running ", length(cohorts), " cohorts sequentially in one R process")
res <- setNames(vector("list", length(cohorts)), names(cohorts))
for (cid in names(cohorts)) {
  checkpoint <- file.path(OUT, switch(as.character(cid),'single_benchmark_preproc'='checkpoint_single_benchmark_preproc.rds','single_frozen_training_state'='checkpoint_single_frozen_training_state.rds','reassay_design_single_frozen_replayed_state'='checkpoint_reassay_design_single_frozen_replayed_state.rds',paste0("checkpoint_", cid, ".rds")))
  if (RESUME && file.exists(checkpoint)) {
    done <- .public_readRDS(checkpoint)
    expected_routes <- if (cohorts[[cid]]$kind == "ct40_well") setdiff(ROUTES, "ws-balance-ilr") else ROUTES
    stopifnot(identical(unique(done$naud$cohort), cid),
              identical(sort(unique(done$fits$split)), seq_len(N_SPLITS)),
              setequal(done$fits$route, expected_routes),
              nrow(done$fits) == N_SPLITS * length(expected_routes))
    res[[cid]] <- done
    say(cid, ": RESUMED completed checkpoint sha256=", digest(checkpoint, file = TRUE, algo = "sha256"))
    preview <- done$rec[path == "S1_single_benchmark_preproc", .(
      n_evaluations = .N, n_rejected = sum(!is.na(err)),
      max_abs_diff = if (any(is.finite(s))) max(abs(s - s_cmp), na.rm = TRUE) else NA_real_,
      n_discordant = sum((s >= cutoff) != (s_cmp >= cutoff), na.rm = TRUE)), by = route]
    say(cid, ": completed S1 checkpoint summary:\n", paste(capture.output(print(preview)), collapse = "\n"))
    next
  }
  say(cid, ": START")
  res[[cid]] <- tryCatch(run_cohort(cid), error = function(e) e)
  if (inherits(res[[cid]], "error") || is.null(res[[cid]]$rec)) {
    msg <- if (inherits(res[[cid]], "error")) conditionMessage(res[[cid]]) else "no result"
    say(cid, ": HOLD: ", msg)
    stop(cid, ": sequential cohort failed: ", msg)
  }
  saveRDS(res[[cid]], file.path(OUT, switch(as.character(cid),'single_benchmark_preproc'='checkpoint_single_benchmark_preproc.rds','single_frozen_training_state'='checkpoint_single_frozen_training_state.rds','reassay_design_single_frozen_replayed_state'='checkpoint_reassay_design_single_frozen_replayed_state.rds',paste0("checkpoint_", cid, ".rds"))))
  writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  say(cid, ": COMPLETE; checkpoint written")
}
rec <- rbindlist(lapply(res, `[[`, "rec"), fill = TRUE)
rj <- rbindlist(lapply(res, `[[`, "rj"), fill = TRUE)
fits <- rbindlist(lapply(res, `[[`, "fits"), fill = TRUE)
naud <- rbindlist(lapply(res, `[[`, "naud"), fill = TRUE)

# X_companions_frozen vs companion_dependence: fill comparator scores
s2 <- rec[path == "S2_single_frozen_training_state", .(cohort, split, route, specimen, s2 = s)]
rec[s2, on = .(cohort, split, route, specimen), s2_join := i.s2]
rec[comparator == "S2_single_frozen_training_state", s_cmp := s2_join]
rec[, s2_join := NULL]

# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------
rec[, rejected := !is.na(err)]
rec[, abs_diff := abs(s - s_cmp)]
rec[, discordant := !rejected & is.finite(s_cmp) & ((s >= cutoff) != (s_cmp >= cutoff))]
rec[, sem_detail := ifelse(grepl("required_missing=[1-9]", detail), "required_missing", NA_character_)]

path_meta <- data.table(
  path = c("B_shuffled", "B_testbatch", "S1_single_benchmark_preproc", "S2_single_frozen_training_state",
           "S3_single_frozen_replayed_state", "X_companions_benchmark_preproc", "X_companions_frozen_preproc"),
  path_description = c(
    'Full deposited cohort rows permuted, benchmark preprocessing, held-out rows scored together in permuted order.',
    'Benchmark preprocessing applied to the held-out rows only (a new batch), scored together.',
    'Each held-out specimen alone through the benchmark preprocessing code, scored alone.',
    'Each held-out specimen alone through per-specimen preprocessing with the frozen feature universe and training-only fill constants, scored alone.',
    'Control. As S2 but fill constants replayed from the full deposited cohort (the benchmark state).',
    'Held-out specimen plus four extreme companions through the benchmark preprocessing, whole batch scored.',
    'Held-out specimen plus four extreme companions through the frozen per-specimen preprocessing, whole batch scored.'),
  uses_other_specimens_in_preprocessing = c("yes_same_as_benchmark", "yes_heldout_batch", "no", "training_constants_only",
                                            "yes_replayed_cohort_constants", "yes_companion_batch", "training_constants_only"))

summ <- rec[, .(
  n_evaluations = .N,
  n_scored = sum(!rejected),
  n_rejected = sum(rejected),
  rejection_rate = mean(rejected),
  n_scored_with_required_missing = sum(!rejected & sem_detail %in% "required_missing"),
  max_abs_diff = if (any(!rejected & is.finite(abs_diff))) max(abs_diff[!rejected], na.rm = TRUE) else NA_real_,
  median_abs_diff = if (any(!rejected & is.finite(abs_diff))) stats::median(abs_diff[!rejected], na.rm = TRUE) else NA_real_,
  n_diff_gt_tol = sum(!rejected & abs_diff > TOL, na.rm = TRUE),
  n_discordant = sum(discordant, na.rm = TRUE),
  discordance_rate_scored = if (sum(!rejected)) sum(discordant, na.rm = TRUE) / sum(!rejected) else NA_real_,
  n_warned = sum(!is.na(warn)),
  example_rejection = if (any(rejected)) substr(names(sort(table(err[rejected]), decreasing = TRUE))[1L], 1L, 160L) else NA_character_
), by = .(cohort, route, path, comparator)]
summ[, verdict := fifelse(n_rejected > 0L, "REJECTED_ANY",
                  fifelse(n_discordant > 0L, "NONEQUIVALENT_DECISION",
                  fifelse(n_diff_gt_tol > 0L, "NONEQUIVALENT_SCORE_ONLY", "EQUIVALENT")))]
# a path can both reject and change decisions; report the decision flag separately
summ[, decision_change_among_scored := n_discordant > 0L]
fit_s <- fits[fit_ok == TRUE, .(n_fits = .N, deploy_verified_n = sum(deploy_verified),
                                pkg_gate_heldout_pass_n = sum(pkg_gate_heldout),
                                model_state_unchanged_n = sum(model_state_unchanged),
                                median_n_test = as.numeric(stats::median(n_test))), by = .(cohort, route)]
summ <- merge(summ, fit_s, by = c("cohort", "route"), all.x = TRUE)
summ <- merge(summ, path_meta, by = "path", all.x = TRUE)
cmeta <- rbindlist(lapply(cohorts, function(co) data.table(cohort = co$id, start_point = co$start,
                                                         preprocessing_kind = co$kind)))
summ <- merge(summ, cmeta, by = "cohort", all.x = TRUE)
summ[, tol_score_abs := TOL]
path_order <- c(path_meta$path)
setorder(summ, cohort, route)
summ[, path := factor(path, levels = path_order)]
setorder(summ, cohort, route, path, comparator)
summ[, path := as.character(path)]
setcolorder(summ, c("cohort", "start_point", "preprocessing_kind", "route", "path", "comparator",
                    "uses_other_specimens_in_preprocessing", "n_fits", "median_n_test", "n_evaluations",
                    "n_scored", "n_rejected", "rejection_rate", "n_scored_with_required_missing",
                    "max_abs_diff", "median_abs_diff", "n_diff_gt_tol", "tol_score_abs", "n_discordant",
                    "discordance_rate_scored", "decision_change_among_scored", "verdict", "n_warned",
                    "deploy_verified_n", "pkg_gate_heldout_pass_n", "model_state_unchanged_n",
                    "example_rejection", "path_description"))
.public_fwrite(summ, file.path(OUT, "e2e_equivalence.tsv"), sep = "\t", quote = FALSE, na = "NA")

# Rejection audit
rj[, rejected := !is.na(err)]
rj[, abs_diff := abs(s - s_base)]
rj[, flip := !rejected & is.finite(s_base) & is.finite(s) & ((s >= cutoff) != (s_base >= cutoff))]
scn_meta <- data.table(
  scenario = c("RJ1_all_required_absent", "RJ2_one_required_absent", "RJ3_one_required_nonfinite",
               "RJ4_half_nonrequired_absent", "RJ5_extra_unregistered_features", "RJ6_all_nondetect_specimen",
               "RJ7_unnamed_vector", "RJ8_single_row_benchmark_preproc_required_missing",
               "RJ9_nondeployable_A_route_deploy"),
  required_inputs_absent = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, NA),
  expected = c("explicit_error", "explicit_error", "explicit_error", "score_unchanged_or_error",
               "score_unchanged_or_error", "explicit_error", "explicit_error", "explicit_error",
               "explicit_error_at_deploy"),
  description = c(
    "Every frozen required feature (panel, pair, or SBP feature) removed from the specimen",
    "The first frozen required feature removed",
    "The first frozen required feature set to NA",
    "Every second non-required universe feature removed",
    "Ten features outside the training universe appended",
    "A specimen with no detected target (all Ct undetermined; all-zero for abundance start) through frozen preprocessing",
    "Feature values supplied without names",
    "Single-specimen benchmark preprocessing dropped at least one required feature (from path S1)",
    "deploy_singlesample() called for a registered A-family route that is not single-sample deployable"))
# RJ8 from matched_baseline records
s1 <- rec[path == "S1_single_benchmark_preproc" & grepl("required_missing=[1-9]", detail)]
if (nrow(s1)) {
  rj <- rbind(rj, s1[, .(cohort, split, route, scenario = "RJ8_single_row_benchmark_preproc_required_missing",
                         specimen, n_required = NA_integer_, s_base = s_cmp, s, err, warn, cutoff,
                         rejected, abs_diff = abs(s - s_cmp),
                         flip = !rejected & is.finite(s) & ((s >= cutoff) != (s_cmp >= cutoff)))], fill = TRUE)
}
agg_rj <- function(d) d[, .(
  n_attempts = .N,
  n_explicit_rejection = sum(rejected),
  n_silent_score = sum(!rejected),
  explicit_rejection_rate = mean(rejected),
  n_silent_changed = sum(!rejected & is.finite(abs_diff) & abs_diff > TOL),
  max_abs_diff_silent = if (any(!rejected & is.finite(abs_diff))) max(abs_diff[!rejected], na.rm = TRUE) else NA_real_,
  n_silent_decision_flip = sum(flip, na.rm = TRUE),
  example_message = if (any(rejected)) substr(names(sort(table(err[rejected]), decreasing = TRUE))[1L], 1L, 160L) else NA_character_
), by = .(cohort, route, scenario)]
rja <- rbind(agg_rj(rj), agg_rj(copy(rj)[, cohort := "ALL"]))
rja <- merge(rja, scn_meta, by = "scenario", all.x = TRUE)
rja[, verdict := fifelse(expected %in% c("explicit_error", "explicit_error_at_deploy"),
                         fifelse(explicit_rejection_rate >= 1, "MEETS_EXPECTATION", "SILENT_DEGRADATION"),
                         fifelse(n_silent_changed == 0L, "MEETS_EXPECTATION", "SILENT_SCORE_SHIFT"))]
rja[, tol_score_abs := TOL]
rja[, cohort_order := cohort == "ALL"]
setorder(rja, scenario, -cohort_order, route, cohort)
rja[, cohort_order := NULL]
setcolorder(rja, c("scenario", "required_inputs_absent", "expected", "route", "cohort", "n_attempts",
                   "n_explicit_rejection", "n_silent_score", "explicit_rejection_rate", "n_silent_changed",
                   "max_abs_diff_silent", "n_silent_decision_flip", "tol_score_abs", "verdict",
                   "example_message", "description"))
.public_fwrite(rja, file.path(OUT, "rejection_audit.tsv"), sep = "\t", quote = FALSE, na = "NA")

# Supporting evidence (specimen level) and audit numbers
.public_fwrite(rec, file.path(OUT, "specimen_level_paths.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
.public_fwrite(fits, file.path(OUT, "fit_level_records.tsv"), sep = "\t", quote = FALSE, na = "NA")
.public_fwrite(naud, file.path(OUT, "preprocessing_dependence_by_split.tsv"), sep = "\t", quote = FALSE, na = "NA")
.public_fwrite(repro, file.path(OUT, "cache_reproduction.tsv"), sep = "\t", quote = FALSE, na = "NA")

say("rows: e2e_equivalence=", nrow(summ), " rejection_audit=", nrow(rja), " specimen records=", nrow(rec))
say("elapsed minutes: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
say("done")
close(log_con)
