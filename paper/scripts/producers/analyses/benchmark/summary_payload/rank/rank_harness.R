
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
# =============================================================================
# singlesample_roster_benchmark.R  —  benchmark / W1
#
# Per-cohort grouped nested-CV benchmark harness for the frozen 74-method
# single-sample roster. For each cohort x method it
#   (1) builds grouped outer folds ONCE per cohort (so EVERY method AND the
#       ws-rclr-trimmed signed-panel reference are scored on byte-identical
#       folds), reusing the fold construction +
#       group-split assertion + eligibility gates from
#       code/methods/matched_null_benchmark.R;
#   (2) per outer fold fits the method on the TRAINING fold only
#       (model <- <fit_fn>(X_train, y_train, meta_train, hp = <frozen analysis_3_m
#       default = list()>)) and scores the held-out TEST fold via
#       OmicSelector::singlesample_score_call(); collects out-of-fold scores for
#       the method AND for the rCLR baseline on the same folds;
#   (3) preserves paired per-fold AUC differences for corrected repeated-CV
#       uncertainty and cross-cohort synthesis in the downstream summarizer;
#   (4) checkpoints analysis/roster_benchmark/<cohort>__<method>.rds (analysis-safe);
#   (5) records ineligible cohort x method cells as NA lift + reason (never a
#       silent drop, never a crash);
#   (6) for transfer methods (estimand == transfer) fits on the OTHER cohorts
#       (provenance-block-aware split; check_provenance run first) and scores the
#       held-out cohort -> table_S6x_roster_transfer.tsv.
#
# This is the W1 of plan/benchmark_execution_plan.md (triple-consensus companion_dependence). It is a
# analysis to be §7-verified + gated before any production run. It does NOT modify
# the committed pipeline and does NOT launch the full benchmark.
#
# --- DESIGN DECISIONS FLAGGED FOR THE GATE (see analysis/benchmark/W1_SEVEN_evidence.md
#     §"Flagged design decisions") -------------------------------------------
#   decision_rules. The locked production comparator is the canonical OmicSelector
#       trimmed-rCLR signed-panel route. On every outer-training fold it ranks
#       features by training-only absolute standardized mean difference, freezes
#       the top-20 panel and case-minus-control signs, and applies the frozen
#       one-score-per-specimen scorer to the identical held-out rows. The reference
#       full-footprint diagnostic remains switchable but is not the production
#       reference.
#   benchmark. Present-method score adapters. OmicSelector::.onLoad does NOT register
#       the reference ("present") score adapters, so singlesample_score_call() errors
#       for every present-status method (the baseline and all present transfer
#       methods). The harness registers the adapters it needs at runtime via
#       singlesample_register_score_adapter() (.register_present_adapters()).
#       Flagged as a package gap to be specified in the package, not the harness.
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("singlesample_roster_benchmark requires pROC for AUC.")
})

# Reuse the gated fold/group machinery + the DeLong paired-diff SE. These two
# sources are READ-ONLY here: the harness reuses .grouped_folds /
# .stratified_folds and delong_lift_se, it does not redefine them.
.srb_source_once <- function(path) {
  if (!file.exists(path)) stop("singlesample_roster_benchmark: missing dependency ", path)
  sys.source(path, envir = topenv())
}
local({
  here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  root_candidates <- c(
    Sys.getenv("OMICSELECTOR_PAPER_ROOT", unset = ""),
    ".",
    getwd())
  root <- root_candidates[nzchar(root_candidates) &
                            file.exists(file.path(root_candidates,
                                                  "code/methods/paired_auc_diff_se.R"))][1L]
  if (is.na(root) || !length(root)) root <- "."
  assign(".SRB_ROOT", root, envir = topenv())
})

if (!exists("delong_lift_se", mode = "function"))
  .srb_source_once(file.path(.SRB_ROOT, "code/methods/paired_auc_diff_se.R"))
if (!exists(".grouped_folds", mode = "function") ||
    !exists(".stratified_folds", mode = "function"))
  .srb_source_once(file.path(.SRB_ROOT, "code/methods/matched_null_benchmark.R"))

`%srb||%` <- function(a, b) if (is.null(a) || length(a) == 0L || (length(a) == 1L && is.na(a))) b else a

srb_verify_cache_bundle <- function(cache_path,
                                    manifest_path = file.path(dirname(cache_path),
                                                              "cohort_cache_manifest.tsv"),
                                    digest_path = file.path(dirname(cache_path),
                                                            "cohort_cache_digest.tsv")) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Cache verification requires the digest package")
  }
  if (!file.exists(cache_path) || !file.exists(manifest_path) || !file.exists(digest_path)) {
    stop("Cache verification requires cache, source manifest, and cache digest files")
  }
  cache_record <- utils::read.delim(digest_path, stringsAsFactors = FALSE,
                                    check.names = FALSE)
  if (nrow(cache_record) != 1L || !all(c("path", "sha256") %in% names(cache_record))) {
    stop("Cache digest record is malformed")
  }
  observed_cache_sha <- digest::digest(cache_path, file = TRUE, algo = "sha256")
  if (!identical(observed_cache_sha, cache_record$sha256[[1L]])) {
    stop("Cohort cache bytes do not match cohort_cache_digest.tsv; rebuild the cache")
  }
  manifest <- utils::read.delim(manifest_path, stringsAsFactors = FALSE,
                                check.names = FALSE)
  if (!all(c("path", "sha256") %in% names(manifest)) || !nrow(manifest) ||
      anyDuplicated(manifest[, c("accession", "path")])) {
    stop("Cohort source manifest is malformed")
  }
  missing <- !file.exists(manifest$path)
  if (any(missing)) {
    stop("Cohort source file is missing: ", paste(manifest$path[missing], collapse = ", "))
  }
  observed_source_sha <- vapply(manifest$path, digest::digest, character(1L),
                                file = TRUE, algo = "sha256")
  mismatch <- observed_source_sha != manifest$sha256
  if (any(mismatch)) {
    stop("Cohort source bytes changed after cache construction: ",
         paste(manifest$path[mismatch], collapse = ", "))
  }
  invisible(TRUE)
}

srb_run_provenance_preflight <- function(cohorts, root, log_path) {
  if (!is.list(cohorts) || !length(cohorts)) stop("Provenance preflight requires cohorts")
  root <- normalizePath(root, mustWork = TRUE)
  check_script <- file.path(root, "code/check_provenance.R")
  if (!file.exists(check_script)) stop("Missing provenance preflight script: ", check_script)
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)

  source_map <- do.call(rbind, lapply(names(cohorts), function(unit) {
    z <- cohorts[[unit]]
    accessions <- unique(as.character(z$source_accessions %srb||%
                                       z$preflight_accessions %srb||% z$accession))
    if (!length(accessions) || anyNA(accessions) || any(!nzchar(accessions))) {
      stop("Missing source accession for analysis unit ", unit)
    }
    data.frame(accession = toupper(accessions), analysis_unit = unit,
               provenance_block = as.character(z$provenance_block %srb||% NA_character_),
               stringsAsFactors = FALSE)
  }))
  if (anyDuplicated(source_map$accession)) {
    dup <- unique(source_map$accession[duplicated(source_map$accession) |
                                         duplicated(source_map$accession, fromLast = TRUE)])
    stop("Source accession assigned to multiple analysis units: ", paste(dup, collapse = ", "))
  }

  status <- system2("/usr/bin/Rscript", c(check_script, sort(source_map$accession)),
                    stdout = log_path, stderr = log_path,
                    env = c("LD_LIBRARY_PATH="))
  status <- as.integer(status)
  if (!(status %in% c(0L, 2L))) {
    stop("Provenance preflight found an unknown accession or execution error (status ",
         status, "); see ", log_path)
  }

  lines <- readLines(log_path, warn = FALSE)
  pair_lines <- grep("^   [A-Z0-9_-]+ <-> [A-Z0-9_-]+[[:space:]]*$", lines,
                     value = TRUE)
  overlap_pairs <- if (length(pair_lines)) {
    matched <- regexec("^   ([A-Z0-9_-]+) <-> ([A-Z0-9_-]+)[[:space:]]*$",
                       pair_lines)
    values <- regmatches(pair_lines, matched)
    unique(data.frame(
      accession_a = vapply(values, `[[`, character(1L), 2L),
      accession_b = vapply(values, `[[`, character(1L), 3L),
      stringsAsFactors = FALSE
    ))
  } else {
    data.frame(accession_a = character(), accession_b = character(),
               stringsAsFactors = FALSE)
  }
  if (status == 2L && !nrow(overlap_pairs)) {
    stop("Provenance preflight reported overlap but no overlap pairs could be parsed; see ",
         log_path)
  }
  hit_line <- grep("^ Overlap hits:[[:space:]]+[0-9]+[[:space:]]*$", lines, value = TRUE)
  reported_hits <- if (length(hit_line) == 1L) {
    as.integer(sub("^ Overlap hits:[[:space:]]+([0-9]+)[[:space:]]*$", "\\1", hit_line))
  } else NA_integer_
  if (is.na(reported_hits) || reported_hits != nrow(overlap_pairs)) {
    stop("Parsed provenance-pair count does not match the preflight summary; see ", log_path)
  }

  if (nrow(overlap_pairs)) {
    unit_for <- stats::setNames(source_map$analysis_unit, source_map$accession)
    overlap_pairs$unit_a <- unname(unit_for[overlap_pairs$accession_a])
    overlap_pairs$unit_b <- unname(unit_for[overlap_pairs$accession_b])
    if (anyNA(overlap_pairs$unit_a) || anyNA(overlap_pairs$unit_b)) {
      stop("Parsed provenance overlap includes an accession outside the cache source map")
    }
    cross_unit <- overlap_pairs$unit_a != overlap_pairs$unit_b
    if (any(cross_unit)) {
      bad <- overlap_pairs[cross_unit, , drop = FALSE]
      stop("Known specimen/cohort reuse crosses analysis units: ",
           paste(paste0(bad$accession_a, "(", bad$unit_a, ") <-> ",
                        bad$accession_b, "(", bad$unit_b, ")"), collapse = "; "))
    }
  }

  composite_units <- unique(source_map$analysis_unit[duplicated(source_map$analysis_unit) |
                                                       duplicated(source_map$analysis_unit,
                                                                  fromLast = TRUE)])
  for (unit in composite_units) {
    z <- cohorts[[unit]]
    groups <- as.character(z$group_id)
    n_profiles <- nrow(z$expr_per_sample)
    if (!length(groups) || length(groups) != n_profiles || anyNA(groups) ||
        any(!nzchar(groups)) || length(unique(groups)) > n_profiles) {
      stop("Grouped provenance unit is missing a valid profile-to-group map: ", unit)
    }
    if (grepl("^audited_independent", as.character(z$provenance_block))) {
      stop("Multi-accession analysis unit cannot be labelled audited_independent: ", unit)
    }
  }

  resolution_path <- sub("\\.log$", "_resolution.tsv", log_path)
  if (identical(resolution_path, log_path)) resolution_path <- switch(as.character(log_path),'fixed'='specified.resolution.tsv',paste0('',log_path,'.resolution.tsv'))
  resolved <- if (nrow(overlap_pairs)) overlap_pairs else data.frame(
    accession_a = character(), accession_b = character(),
    unit_a = character(), unit_b = character(), stringsAsFactors = FALSE)
  .public_write_table(resolved, resolution_path, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")
  invisible(list(raw_status = status, n_overlap_pairs = nrow(overlap_pairs),
                 n_source_accessions = nrow(source_map),
                 n_analysis_units = length(cohorts),
                 resolution_path = resolution_path))
}

# --- the rCLR baseline scorer (decision_rules) -------------------------------------------
# Per-fold, single-sample-exact trimmed-rCLR summary over the test rows. Uses the
# package transform OmicSelector::ws_rclr_trimmed (matrix in -> matrix out),
# summed across features. NA columns from the transform contribute 0 (na.rm).
.srb_rclr_baseline_score <- function(X_test) {
  if (!requireNamespace("OmicSelector", quietly = TRUE))
    stop("OmicSelector must be loaded (devtools::load_all) for the rCLR baseline.")
  Z <- tryCatch(OmicSelector::ws_rclr_trimmed(X_test), error = function(e) NULL)
  if (is.null(Z)) return(rep(NA_real_, nrow(X_test)))
  if (!is.matrix(Z)) Z <- matrix(Z, nrow = nrow(X_test))
  s <- rowSums(Z, na.rm = TRUE)
  s[!is.finite(s)] <- NA_real_
  s
}

# Canonical training-frozen rCLR signed-panel reference. This calls the public
# OmicSelector fit/score route directly so the benchmark and shipped software
# cannot drift into two implementations of the comparator.
.srb_rclr_panel_score <- function(X_train, y_train, X_test, panel_size = 20L) {
  if (!requireNamespace("OmicSelector", quietly = TRUE))
    stop("OmicSelector must be loaded for the panel baseline.")
  model <- OmicSelector::fit_ws_rclr_panel(
    X_train, y_train, hp = list(panel_size = as.integer(panel_size)))
  OmicSelector::score_ws_rclr_panel(model, X_test)
}

.srb_training_direction <- function(score, y) {
  ok <- is.finite(score) & !is.na(y)
  if (sum(ok) < 4L || length(unique(y[ok])) < 2L) return(NA_real_)
  d <- stats::median(score[ok & y == 1L], na.rm = TRUE) -
    stats::median(score[ok & y == 0L], na.rm = TRUE)
  if (!is.finite(d)) NA_real_ else if (d < 0) -1 else 1
}

.srb_signed_panel <- function(Z_train, y_train, Z_test, panel_size = 20L) {
  if (!is.matrix(Z_train)) Z_train <- as.matrix(Z_train)
  if (!is.matrix(Z_test)) Z_test <- as.matrix(Z_test)
  common <- intersect(colnames(Z_train), colnames(Z_test))
  if (length(common) < 2L) stop("fewer than two common transformed features")
  effect <- vapply(common, function(j) {
    x <- Z_train[, j]; ok <- is.finite(x) & !is.na(y_train)
    den <- stats::sd(x[ok])
    if (sum(ok) < 5L || !is.finite(den) || den <= 0) return(0)
    abs(mean(x[ok & y_train == 1L]) - mean(x[ok & y_train == 0L])) / den
  }, numeric(1L))
  panel <- common[order(effect, decreasing = TRUE)[seq_len(min(panel_size, length(common)))]]
  weights <- vapply(panel, function(j) {
    d <- stats::median(Z_train[y_train == 1L, j], na.rm = TRUE) -
      stats::median(Z_train[y_train == 0L, j], na.rm = TRUE)
    if (is.finite(d) && d < 0) -1 else 1
  }, numeric(1L))
  list(train = as.numeric(Z_train[, panel, drop = FALSE] %*% weights),
       test = as.numeric(Z_test[, panel, drop = FALSE] %*% weights),
       panel = panel, weights = weights)
}

.srb_present_within_fold <- function(method_id, X_train, y_train, X_test,
                                     panel_size = 20L) {
  ns <- asNamespace("OmicSelector")
  call_ns <- function(fn, ...) get(fn, envir = ns, inherits = FALSE)(...)
  if (identical(method_id, "ws-dominance")) {
    tr <- call_ns("ws_dominance_score", X_train)
    te <- call_ns("ws_dominance_score", X_test)
    direction <- .srb_training_direction(tr, y_train)
    if (!is.finite(direction)) stop("dominance score direction unavailable")
    return(list(train = direction * tr, test = direction * te,
                panel = character(0), weights = direction))
  }
  if (identical(method_id, "ws-alr-pivot")) {
    Z_train <- call_ns("ws_alr_pivot", X_train, allow_global_fallback = TRUE)
    Z_test <- call_ns("ws_alr_pivot", X_test, allow_global_fallback = TRUE)
  } else if (identical(method_id, "ws-mad-logratio")) {
    Z_train <- call_ns("ws_mad_logratio", X_train)
    Z_test <- call_ns("ws_mad_logratio", X_test)
  } else if (identical(method_id, "frozen-ruv")) {
    fit <- call_ns("fit_frozen_ruv", X_train)
    Z_train <- call_ns("apply_frozen_ruv", X_train, fit)
    Z_test <- call_ns("apply_frozen_ruv", X_test, fit)
  } else if (identical(method_id, "robust-pca-resid")) {
    fit <- call_ns("fit_robust_pca_residual", X_train)
    Z_train <- call_ns("apply_robust_pca_residual", fit, X_train)
    Z_test <- call_ns("apply_robust_pca_residual", fit, X_test)
  } else if (identical(method_id, "hemolysis-resid")) {
    h_train <- call_ns("hemolysis_index_blondal", log(X_train + 0.5))
    h_test <- call_ns("hemolysis_index_blondal", log(X_test + 0.5))
    if (is.null(h_train) || is.null(h_test)) stop("required hemolysis anchors absent")
    fit <- call_ns("fit_hemolysis_rr", X_train, h_train)
    Z_train <- call_ns("apply_hemolysis_rr", X_train, fit, h_train)
    Z_test <- call_ns("apply_hemolysis_rr", X_test, fit, h_test)
  } else {
    stop("no present-method benchmark adapter for ", method_id)
  }
  .srb_signed_panel(Z_train, y_train, Z_test, panel_size = panel_size)
}

# --- present-method score adapters (benchmark) --------------------------------------
# Register canonical (model, X, meta) shims for the reference present-status scorers
# singlesample_score_call and need no adapter. Idempotent.
.srb_present_adapters_registered <- new.env(parent = emptyenv())
.register_present_adapters <- function(roster) {
  if (!requireNamespace("OmicSelector", quietly = TRUE)) return(invisible(FALSE))
  ns <- asNamespace("OmicSelector")
  present <- roster[roster$pkg_status == "present", , drop = FALSE]
  for (i in seq_len(nrow(present))) {
    mid <- present$method_id[i]; sfn <- present$score_fn[i]
    if (is.na(sfn) || !nzchar(sfn)) next
    if (exists(mid, envir = .srb_present_adapters_registered, inherits = FALSE)) next
    if (!exists(sfn, envir = ns, mode = "function", inherits = FALSE)) next
    fn <- get(sfn, envir = ns, mode = "function", inherits = FALSE)
    # Name-aware shim. Present scorers do NOT share one positional convention:
    #   - canonical (model, X, meta):          score_ecod_copod, score_sinkhorn_single
    #   - fit-first  (fit, x_test, ...):         apply_compositional_mahalanobis,
    #                                            apply_isolation_forest_logratio
    #   - DATA-FIRST (x, fit, ...) [reversed!]:  os_conformal_anomaly
    # A purely positional shim mis-binds the data-first scorer (passes the feature
    # matrix as `fit`) and coerces the spurious 3rd arg (return_pvalue/alpha) from a
    # data.frame -> "list cannot be coerced to double". So bind the fit-object and the
    # (e.g. return_pvalue/alpha/fpr_target). Pass meta only if the scorer names it.
    # (keeps any present DISCRIMINATOR adapter that relied on positional dispatch
    # working; those are skipped in W1 under option B but still registered here).
    # (companion_dependence-gate HOLD fix: lets the 3 present negative controls score correctly.)
    shim <- local({
      f <- fn
      fmls <- names(formals(f))
      fit_arg  <- intersect(c("fit", "model"), fmls)
      data_arg <- intersect(c("x_test", "X", "x", "newx", "x_new", "X_new", "X_test"), fmls)
      function(model, X, meta) {
        out <- if (length(fit_arg) && length(data_arg)) {
          args <- list(); args[[fit_arg[1]]] <- model; args[[data_arg[1]]] <- X
          if ("meta" %in% fmls) args[["meta"]] <- meta
          do.call(f, args)
        } else {
          o <- tryCatch(f(model, X, meta), error = function(e) NULL)
          if (is.null(o)) f(model, X) else o
        }
        # Present one-class scorers return a LIST carrying the per-sample atypicality
        # vector plus auxiliary fields (mahalanobis: distance/distance_sq/pvalue_chisq;
        # conformal: score/p_value/is_anomaly; isoforest: score/pvalue_emp/flagged).
        # Extract the canonical per-sample score (length == nrow(X)); fall back to the
        # first numeric field of the right length. A bare-numeric return passes through.
        if (is.list(out) && !is.data.frame(out)) {
          pick <- NULL
          for (nm in c("score", "distance", "distance_sq", "anomaly_score", "nonconformity")) {
            v <- out[[nm]]
            if (!is.null(v) && is.numeric(v) && length(v) == nrow(X)) { pick <- v; break }
          }
          if (is.null(pick))
            for (v in out) if (is.numeric(v) && length(v) == nrow(X)) { pick <- v; break }
          out <- pick
        }
        as.numeric(out)
      }
    })
    OmicSelector::singlesample_register_score_adapter(mid, shim)
    assign(mid, TRUE, envir = .srb_present_adapters_registered)
  }
  invisible(TRUE)
}

# --- per cohort x method evaluation (within estimand) ------------------------
# Returns a per-fold long data.frame of paired scores (method + baseline on the
# SAME test rows/folds) plus a fold audit, or an ineligible record.
.srb_validate_heldout_score_support <- function(score_m, score_b, y,
                                                expected_n = length(y)) {
  complete_numeric <- is.numeric(score_m) && is.numeric(score_b) &&
    length(score_m) == expected_n && length(score_b) == expected_n
  if (!complete_numeric) {
    return(list(
      valid = FALSE,
      reason = "method/baseline score is not a complete numeric held-out vector"
    ))
  }
  if (expected_n < 4L || length(y) != expected_n || anyNA(y) ||
      length(unique(y)) < 2L) {
    return(list(
      valid = FALSE,
      reason = "held-out support has fewer than four paired rows or one outcome class"
    ))
  }
  if (any(!is.finite(score_m)) || any(!is.finite(score_b))) {
    return(list(
      valid = FALSE,
      reason = "method/baseline score is not a complete finite held-out vector"
    ))
  }
  # A complete constant score is valid no-discrimination evidence. Midrank AUC
  # assigns every positive-negative tie 0.5, and paired DeLong placements remain
  # defined (possibly with exact zero variance).
  list(valid = TRUE, reason = NA_character_)
}

.srb_eval_within_cell <- function(method_id, fit_fn, X, y, group_id, roster,
                                  outer_k = 5L, seed = 42L,
                                  min_valid_folds = outer_k,
                                  min_pos_per_fold = 5L, min_neg_per_fold = 5L,
                                  baseline_mode = "rclr_full",
                                  baseline_method_id = "ws-rclr-trimmed",
                                  panel_size = 20L,
                                  assay_panel_size = Inf,
                                  hp = list()) {
  n <- nrow(X)
  has_groups <- !is.null(group_id) && length(group_id) == n
  grouping_policy <- if (has_groups) "groupKFold" else "ungrouped_stratified"

  folds <- if (has_groups) .grouped_folds(y, group_id, k = outer_k, seed = seed)
           else .stratified_folds(y, k = outer_k, seed = seed)

  # Hard group-split assertion (reused contract from matched_null_benchmark.R).
  n_group_split_violations <- 0L
  if (has_groups) {
    for (fi in seq_along(folds)) {
      train_idx <- setdiff(seq_len(n), folds[[fi]])
      viol <- length(intersect(unique(group_id[folds[[fi]]]),
                                unique(group_id[train_idx])))
      n_group_split_violations <- n_group_split_violations + viol
      if (viol > 0L)
        stop(sprintf("singlesample_roster_benchmark: %d group(s) split across train/test in fold %d (%s).",
                     viol, fi, method_id))
    }
  }

  is_baseline_cell <- identical(method_id, baseline_method_id)
  rows <- list()
  fold_audit <- vector("list", outer_k)
  valid_fold_mask <- logical(outer_k)

  for (fi in seq_len(outer_k)) {
    test_idx  <- folds[[fi]]
    train_idx <- setdiff(seq_len(n), test_idx)
    y_test <- y[test_idx]; y_train <- y[train_idx]
    n_pos_test <- sum(y_test == 1L); n_neg_test <- sum(y_test == 0L)
    fa <- list(fold = fi, n_train = length(train_idx), n_test = length(test_idx),
               n_pos_test = n_pos_test, n_neg_test = n_neg_test)
    if (n_pos_test < min_pos_per_fold || n_neg_test < min_neg_per_fold) {
      fa$skip_reason <- sprintf("test n_pos=%d<%d or n_neg=%d<%d",
                                n_pos_test, min_pos_per_fold, n_neg_test, min_neg_per_fold)
      fold_audit[[fi]] <- fa; next
    }
    if (length(unique(y_train)) < 2L) {
      fa$skip_reason <- "y_train one class"; fold_audit[[fi]] <- fa; next
    }
    X_train <- X[train_idx, , drop = FALSE]; X_test <- X[test_idx, , drop = FALSE]
    zz<-m18_fold_prepare(X_train,X_test);X_train<-zz$train;X_test<-zz$test
    if (is.finite(assay_panel_size) && ncol(X_train) > assay_panel_size) {
      prevalence <- colMeans(is.finite(X_train) & X_train > 0)
      dispersion <- apply(log1p(X_train), 2L, stats::mad, na.rm = TRUE)
      dispersion[!is.finite(dispersion)] <- -Inf
      ranked <- order(prevalence, dispersion, colnames(X_train), decreasing = TRUE)
      assay_features <- colnames(X_train)[head(ranked, as.integer(assay_panel_size))]
      X_train <- X_train[, assay_features, drop = FALSE]
      X_test <- X_test[, assay_features, drop = FALSE]
    } else {
      assay_features <- colnames(X_train)
    }
    fa$n_assay_features <- length(assay_features)
    fa$assay_features <- assay_features
    meta_train <- NULL; meta_test <- NULL

    # Baseline score on this fold's TEST rows (train-free for rclr_full).
    s_base <- tryCatch({
      if (identical(baseline_mode, "rclr_panel"))
        .srb_rclr_panel_score(X_train, y_train, X_test, panel_size = panel_size)
      else
        .srb_rclr_baseline_score(X_test)
    }, error = function(e) rep(NA_real_, length(test_idx)))

    # Method score. For the baseline cell, method == baseline (lift 0 by
    # construction). Otherwise fit on train fold ONLY, score test fold.
    if (is_baseline_cell) {
      s_meth <- s_base
    } else if (method_id %in% c("ws-alr-pivot", "ws-mad-logratio", "ws-dominance",
                                "frozen-ruv", "robust-pca-resid", "hemolysis-resid")) {
      present_err <- NULL
      present <- tryCatch(.srb_present_within_fold(method_id, X_train, y_train, X_test,
                                                   panel_size = panel_size),
                          error = function(e) {
                            present_err <<- conditionMessage(e)
                            NULL
                          })
      if (is.null(present)) {
        fa$skip_reason <- paste0("present adapter failed: ", present_err %srb||% "unknown error")
        fold_audit[[fi]] <- fa; next
      }
      s_meth <- present$test
    } else {
      model <- tryCatch(
        if (is.na(fit_fn) || !nzchar(fit_fn)) NULL
        else .srb_call_fit(fit_fn, X_train, y_train, meta_train, hp),
        error = function(e) structure(NA, srb_err = conditionMessage(e)))
      fit_err <- attr(model, "srb_err")
      if (is.null(model) || !is.null(fit_err)) {
        fa$skip_reason <- paste0("fit failed: ", fit_err %srb||% "no fit_fn")
        fold_audit[[fi]] <- fa; next
      }
      s_train <- tryCatch(
        OmicSelector::singlesample_score_call(method_id, model, X_train, meta_train, roster = roster),
        error = function(e) rep(NA_real_, nrow(X_train)))
      direction <- .srb_training_direction(s_train, y_train)
      if (!is.finite(direction)) {
        fa$skip_reason <- "training-only score direction unavailable"
        fold_audit[[fi]] <- fa; next
      }
      s_meth <- direction * tryCatch(
        OmicSelector::singlesample_score_call(method_id, model, X_test, meta_test, roster = roster),
        error = function(e) { attr(fa, "score_err") <- conditionMessage(e); rep(NA_real_, length(test_idx)) })
    }

    support <- .srb_validate_heldout_score_support(
      s_meth, s_base, y_test, expected_n = length(test_idx)
    )
    if (!isTRUE(support$valid)) {
      fa$skip_reason <- support$reason
      fold_audit[[fi]] <- fa; next
    }
    ok <- rep(TRUE, length(test_idx))
    valid_fold_mask[fi] <- TRUE
    fold_audit[[fi]] <- fa
    rows[[length(rows) + 1L]] <- data.frame(
      fold = fi,
      n_train = length(train_idx),
      n_test = length(test_idx),
      sample_idx = test_idx,
      sample_id = rownames(X)[test_idx] %srb||% as.character(test_idx),
      group_id = if (has_groups) as.character(group_id[test_idx]) else rownames(X)[test_idx] %srb||% as.character(test_idx),
      y = as.integer(y_test),
      score_m = as.numeric(s_meth),
      score_b = as.numeric(s_base),
      stringsAsFactors = FALSE)
  }

  n_valid_folds <- sum(valid_fold_mask)
  pos_v <- vapply(which(valid_fold_mask), function(fi) fold_audit[[fi]]$n_pos_test %srb||% NA_integer_, integer(1L))
  neg_v <- vapply(which(valid_fold_mask), function(fi) fold_audit[[fi]]$n_neg_test %srb||% NA_integer_, integer(1L))
  skip_summary <- paste(vapply(fold_audit, function(fa)
    (fa$skip_reason %srb||% "ok"), character(1L)), collapse = "; ")

  if (n_valid_folds < min_valid_folds) {
    return(list(
      eligible = FALSE,
      ineligible_reason = sprintf("n_valid_folds=%d < %d; per-fold: %s",
                                  n_valid_folds, min_valid_folds, skip_summary),
      fold_scores = NULL,
      n_valid_folds = n_valid_folds,
      grouping_policy = grouping_policy,
      n_group_split_violations = n_group_split_violations,
      min_pos_per_fold = if (length(pos_v)) min(pos_v) else NA_integer_,
      min_neg_per_fold = if (length(neg_v)) min(neg_v) else NA_integer_,
      fold_audit = fold_audit))
  }
  list(
    eligible = TRUE,
    ineligible_reason = NA_character_,
    fold_scores = do.call(rbind, rows),
    n_valid_folds = n_valid_folds,
    grouping_policy = grouping_policy,
    n_group_split_violations = n_group_split_violations,
    min_pos_per_fold = if (length(pos_v)) min(pos_v) else NA_integer_,
    min_neg_per_fold = if (length(neg_v)) min(neg_v) else NA_integer_,
    fold_audit = fold_audit)
}

# Resolve + call a fit_fn from the OmicSelector namespace. Two fit conventions
# coexist in the roster and are dispatched by signature:
#     or `...`): (X_train, y_train, meta_train, hp = hp).
#   - reference one-class fitters (the 3 present negative controls
#     fit_compositional_mahalanobis / fit_conformal_anomaly /
#     fit_isolation_forest_logratio): take ONLY the feature matrix as their first
#     positional arg (named x_train / X_ref) with their own defaults for the rest, and
#     ignore the label by construction -- an UNSUPERVISED atypicality scorer cannot use
#     y, which is exactly why it is a negative control. Calling them with the canonical
#     4-arg signature errors ("unused argument (hp = hp)").
# Dispatch to the canonical signature whenever the fitter accepts the benchmark args
# (hp / ... / y_train+meta_train); otherwise call it train-only as a one-class fitter.
# (companion_dependence-gate HOLD fix: lets the present NCs actually run instead of "fit failed".)
.srb_call_fit <- function(fit_fn, X_train, y_train, meta_train, hp) {
  ns <- asNamespace("OmicSelector")
  if (!exists(fit_fn, envir = ns, mode = "function", inherits = FALSE))
    stop(sprintf("fit_fn '%s' not found in OmicSelector namespace.", fit_fn))
  f <- get(fit_fn, envir = ns, mode = "function", inherits = FALSE)
  fmls <- names(formals(f))
  if (("hp" %in% fmls) || ("..." %in% fmls) ||
      all(c("y_train", "meta_train") %in% fmls)) {
    f(X_train, y_train, meta_train, hp = hp)
  } else {
    f(X_train)
  }
}

#' Benchmark ONE within-estimand cohort x method cell (analysis-safe).
#'
#' @return invisibly the per-cell result list (also written to checkpoint).
srb_benchmark_within_cell <- function(method_id, cohort_id, cohort, roster,
                                      out_dir, seed = 42L, outer_k = 5L,
                                      min_valid_folds = 4L,
                                      min_pos_per_fold = 5L, min_neg_per_fold = 5L,
                                      baseline_mode = "rclr_full",
                                      panel_size = 20L,
                                      assay_panel_size = Inf,
                                      hp = list(), overwrite = FALSE,
                                      rerun_present = FALSE,
                                      orient = "fixed") {
  ckpt <- file.path(out_dir, sprintf("%s__%s.rds", .srb_safe(cohort_id), .srb_safe(method_id)))
  if (!overwrite && file.exists(ckpt)) {
    cached <- .srb_read_ckpt(ckpt)
    if (!is.null(cached)) return(invisible(cached))   # unreadable -> fall through, recompute
  }

  row <- roster[roster$method_id == method_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("method_id not unique in roster: ", method_id)
  fit_fn <- row$fit_fn[1]

  # (method_expansion_prespecification.md:81/103/479 -- the existing methods' lifts
  # are "already computed" in reassay_design/method_cohort_performance-sensitivity_effects; their verdicts are recomputed under the
  # present-status method is SKIPPED here with an explicit reason -- never the
  # misleading "fit failed: no fit_fn" drop, and never a fabricated lift=0. W6
  # sources its lift from the existing table. (benchmark_cells-gate HOLD fix; see
  # analysis/benchmark/W1_design_fork_present_methods.md option B.)
  # NB: the skip is keyed on role=="discriminator", NOT merely pkg_status=="present".
  # Two present non-discriminator roles MUST run in W1 (they have NO existing reassay_design/method_cohort_performance-sensitivity_effects
  # lift row to reuse in W6):
  #   (i)  the baseline ws-rclr-trimmed (role=="baseline") -- the rCLR comparator /
  #        self-lift-0 anchor; and
  #   (ii) the 3 present negative controls nc-mahalanobis / nc-conformal-anom /
  #        nc-isoforest (role=="negative_control", pkg_status=="present") -- the prespec
  #        (method_expansion_prespecification.md:229-236,401-404) mandates ALL 5 NCs be
  #        "run and reported in full" as the expected-null sanity tier; only the 20
  #        present DISCRIMINATORS carry published S3b/method_cohort_performance-sensitivity_effects lifts to reuse.
  # Keying on pkg_status alone would skip the 3 present NCs to a non-existent table join
  # -- the same silent-drop failure class as the benchmark_cells HOLD, relocated to the NC tier
  # (companion_dependence-gate HOLD fix; verified: S3b holds the 6 present discriminators + baseline but
  # no mahalanobis/conformal/isoforest row).
  if (!isTRUE(rerun_present) && identical(row$pkg_status[1], "present") &&
      identical(row$role[1], "discriminator")) {
    res <- list(method_id = method_id, cohort = cohort_id, estimand = "within",
                est = NA_real_, se = NA_real_, df = NA_real_, n_eff = NA_integer_,
                k_valid_folds = 0L, n_valid_folds = 0L,
                min_pos = NA_integer_, min_neg = NA_integer_,
                grouping_policy = NA_character_, n_group_split_violations = NA_integer_,
                eligible = FALSE, skipped_present = TRUE,
                ineligible_reason = "present-status: lift from existing S3/S6-S9 table; not re-run in W1 (W6 join)",
                baseline_mode = baseline_mode, seed = seed)
    .srb_save_ckpt(res, ckpt); return(invisible(res))
  }

  X <- cohort$expr_per_sample
  y <- as.integer(cohort$y_bin)
  group_id <- cohort$group_id

  ev <- tryCatch(
    .srb_eval_within_cell(method_id, fit_fn, X, y, group_id, roster,
                          outer_k = outer_k, seed = seed,
                          min_valid_folds = min_valid_folds,
                          min_pos_per_fold = min_pos_per_fold,
                          min_neg_per_fold = min_neg_per_fold,
                          baseline_mode = baseline_mode,
                          panel_size = panel_size,
                          assay_panel_size = assay_panel_size,
                          hp = hp),
    error = function(e) list(eligible = FALSE,
                             ineligible_reason = paste0("cell error: ", conditionMessage(e)),
                             fold_scores = NULL, n_valid_folds = 0L,
                             grouping_policy = NA_character_,
                             n_group_split_violations = NA_integer_,
                             min_pos_per_fold = NA_integer_, min_neg_per_fold = NA_integer_))

  if (!isTRUE(ev$eligible)) {
    res <- list(method_id = method_id, cohort = cohort_id, estimand = "within",
                est = NA_real_, se = NA_real_, df = NA_real_, n_eff = NA_integer_,
                k_valid_folds = ev$n_valid_folds %srb||% 0L,
                n_valid_folds = ev$n_valid_folds %srb||% 0L,
                min_pos = ev$min_pos_per_fold, min_neg = ev$min_neg_per_fold,
                grouping_policy = ev$grouping_policy,
                n_group_split_violations = ev$n_group_split_violations,
                eligible = FALSE, ineligible_reason = ev$ineligible_reason,
                baseline_mode = baseline_mode, panel_size = panel_size,
                assay_panel_size = assay_panel_size, seed = seed)
    .srb_save_ckpt(res, ckpt); return(invisible(res))
  }

  fold_dt <- ev$fold_scores
  ls <- delong_lift_se(fold_dt, orient = orient)
  fold_auc <- lapply(split(fold_dt, fold_dt$fold), function(d) {
    c(method = .delong_placements(d$y, d$score_m, orient = orient)$auc,
      baseline = .delong_placements(d$y, d$score_b, orient = orient)$auc)
  })
  fold_auc <- do.call(rbind, fold_auc)
  n_eff <- nrow(fold_dt)
  res <- list(
    method_id = method_id, cohort = cohort_id, estimand = "within",
    est = ls$lift, se = ls$se, df = as.double(max(0L, ls$k - 1L)), n_eff = n_eff,
    auc_method = mean(fold_auc[, "method"], na.rm = TRUE),
    auc_baseline = mean(fold_auc[, "baseline"], na.rm = TRUE),
    k_valid_folds = ls$k, n_valid_folds = ev$n_valid_folds,
    min_pos = ev$min_pos_per_fold, min_neg = ev$min_neg_per_fold,
    grouping_policy = ev$grouping_policy,
    n_group_split_violations = ev$n_group_split_violations,
    eligible = TRUE, ineligible_reason = NA_character_,
    baseline_mode = baseline_mode, panel_size = panel_size,
    assay_panel_size = assay_panel_size, seed = seed,
    fold_scores = fold_dt)
  .srb_save_ckpt(res, ckpt)
  invisible(res)
}

# --- transfer estimand: fit on OTHER provenance blocks, score held-out cohort -
#' Benchmark ONE transfer cohort x method cell (held-out = `cohort_id`).
#'
#' Provenance-block-aware: training cohorts are every cohort whose
#' provenance_block differs from the held-out cohort's block (NA block ==
#' audited_independent singleton, kept). check_provenance() must have been run by
#' the driver before assembling the training accessions.
srb_benchmark_transfer_cell <- function(method_id, held_out_id, cohorts, roster,
                                        out_dir, seed = 42L, overwrite = FALSE,
                                        baseline_mode = "rclr_full",
                                        hp = list(), orient = "fixed") {
  ckpt <- file.path(out_dir, sprintf("transfer__%s__%s.rds",
                                     .srb_safe(held_out_id), .srb_safe(method_id)))
  if (!overwrite && file.exists(ckpt)) {
    cached <- .srb_read_ckpt(ckpt)
    if (!is.null(cached)) return(invisible(cached))   # unreadable -> fall through, recompute
  }

  row <- roster[roster$method_id == method_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("method_id not unique: ", method_id)
  fit_fn <- row$fit_fn[1]

  # from the existing method_cohort_performance-sensitivity_effects/S47/S48/S65 tables (W6 join), verdicts recomputed under
  # d1-d4. Skip-with-reason here -- NEVER the prior fabricated lift=0 (baseline
  # substitution). Keyed on role=="discriminator" (symmetry with the within cell): all
  # 5 negative controls are within-estimand, so no present NC reaches the transfer cell
  # today, but keying on role keeps a future present-NC-transfer from being mis-skipped
  # to a non-existent table join. (benchmark_cells-gate HOLD fix + companion_dependence NC fix;
  # analysis/benchmark/W1_design_fork_present_methods.md option B.)
  if (identical(row$pkg_status[1], "present") && identical(row$role[1], "discriminator")) {
    res <- .srb_transfer_ineligible(method_id, held_out_id,
            "present-status: lift from existing S6-S9 table; not re-run in W1 (W6 join)",
            baseline_mode, seed)
    res$skipped_present <- TRUE
    .srb_save_ckpt(res, ckpt); return(invisible(res))
  }

  ho <- cohorts[[held_out_id]]
  ho_block <- ho$provenance_block
  train_ids <- setdiff(names(cohorts), held_out_id)
  # Drop any cohort sharing the held-out cohort's provenance block (specimen
  # leakage guard). NA block = independent singleton -> never matches a real block.
  if (!is.null(ho_block) && !is.na(ho_block)) {
    train_ids <- train_ids[vapply(train_ids, function(id) {
      b <- cohorts[[id]]$provenance_block
      is.null(b) || is.na(b) || !identical(b, ho_block)
    }, logical(1L))]
  }
  if (length(train_ids) < 1L) {
    res <- .srb_transfer_ineligible(method_id, held_out_id,
            sprintf("no provenance-clean training cohort for block '%s'", ho_block %srb||% "NA"),
            baseline_mode, seed)
    .srb_save_ckpt(res, ckpt); return(invisible(res))
  }

  # Build a pooled candidate-feature matrix across clean training cohorts.
  # Canonicalisation is transfer-local: never mutate the shared cohort cache.
  if (!exists(".normalize_feature_ids", mode = "function") ||
      !exists("normalize_to_canonical_mirna_namespace", mode = "function")) {
    .srb_source_once(file.path(.SRB_ROOT, "code/methods/cohort_loader.R"))
  }
  .srb_transfer_canonicalize <- function(X, label) {
    if (!is.matrix(X)) X <- as.matrix(X)
    if (is.null(colnames(X)))
      stop("transfer canonicalization: ", label, " has no feature names")
    X <- X[, , drop = FALSE]
    colnames(X) <- .normalize_feature_ids(colnames(X))
    X <- normalize_to_canonical_mirna_namespace(X, scale_state = "abundance")
    if (any(!is.finite(X)))
      stop("transfer canonicalization: ", label, " has non-finite abundance values")
    if (any(X < 0, na.rm = TRUE))
      stop("transfer canonicalization: ", label, " has negative abundance values")
    X
  }
  .srb_zero_fill_candidates <- function(X, features, label) {
    out <- matrix(0, nrow = nrow(X), ncol = length(features),
                  dimnames = list(rownames(X), features))
    present <- intersect(features, colnames(X))
    if (length(present) > 0L)
      out[, present] <- X[, present, drop = FALSE]
    storage.mode(out) <- "double"
    if (any(!is.finite(out)) || any(out < 0, na.rm = TRUE))
      stop("transfer zero-fill: ", label, " produced non-finite or negative values")
    out
  }

  Xho <- .srb_transfer_canonicalize(ho$expr_per_sample, held_out_id)
  Xtr <- lapply(train_ids, function(id)
    .srb_transfer_canonicalize(cohorts[[id]]$expr_per_sample, id))
  names(Xtr) <- train_ids
  feature_counts <- table(unlist(lapply(Xtr, colnames), use.names = FALSE))
  candidates <- sort(names(feature_counts)[feature_counts >= 3L])
  n_held_candidate_features <- length(intersect(candidates, colnames(Xho)))
  fill_policy <- "zero_fill_absent_candidate_features_full_candidate_universe"
  if (length(candidates) < 8L) {
    res <- .srb_transfer_ineligible(method_id, held_out_id,
            sprintf("only %d candidate features across training cohorts (min 8; fill_policy=%s)",
                    length(candidates), fill_policy),
            baseline_mode, seed)
    res$n_train_cohorts <- length(train_ids); res$train_cohorts <- train_ids
    res$n_candidate_features <- length(candidates)
    res$n_held_candidate_features <- n_held_candidate_features
    res$candidate_fill_policy <- fill_policy
    .srb_save_ckpt(res, ckpt); return(invisible(res))
  }
  if (n_held_candidate_features < 8L) {
    res <- .srb_transfer_ineligible(method_id, held_out_id,
            sprintf("only %d held-out candidate features present (min 8; n_candidate_features=%d; fill_policy=%s)",
                    n_held_candidate_features, length(candidates), fill_policy),
            baseline_mode, seed)
    res$n_train_cohorts <- length(train_ids); res$train_cohorts <- train_ids
    res$n_candidate_features <- length(candidates)
    res$n_held_candidate_features <- n_held_candidate_features
    res$candidate_fill_policy <- fill_policy
    .srb_save_ckpt(res, ckpt); return(invisible(res))
  }
  X_train <- do.call(rbind, lapply(train_ids, function(id)
    .srb_zero_fill_candidates(Xtr[[id]], candidates, id)))
  y_train <- unlist(lapply(train_ids, function(id) as.integer(cohorts[[id]]$y_bin)), use.names = FALSE)
  meta_train <- data.frame(
    accession = rep(train_ids, vapply(train_ids, function(id) nrow(Xtr[[id]]), integer(1L))),
    stringsAsFactors = FALSE)
  X_test <- .srb_zero_fill_candidates(Xho, candidates, held_out_id)
  y_test <- as.integer(ho$y_bin)
  meta_test <- data.frame(accession = rep(held_out_id, nrow(X_test)), stringsAsFactors = FALSE)

  s_base <- tryCatch(.srb_rclr_baseline_score(X_test), error = function(e) rep(NA_real_, length(y_test)))
  s_meth <- tryCatch({
    # A NA / empty fit_fn here would be a stateless transfer scorer. We do NOT
    # substitute s_base (that fabricates lift=0 -- the benchmark_cells transfer HOLD). Instead
    # NA-out with a reason so the downstream eligibility guard records an honest
    # and the present NCs all carry a fit_fn), so this branch is unreachable today --
    # it exists to prevent silent re-introduction of the fabricated-0 bug.
    if (is.na(fit_fn) || !nzchar(fit_fn))
      structure(rep(NA_real_, length(y_test)),
                srb_err = "no fit_fn (transfer): stateless transfer scorer cannot be benchmarked train-free; NA-out (never fabricate baseline lift=0)")
    else {
      hp_fit <- hp
      if (is.list(hp_fit) &&
          method_id %in% c("cvae", "dann", "dg-fishr", "dg-ibirm",
                           "dro-group", "dro-vrex", "icp",
                           "moe-gated", "sel-stablemate") &&
          is.null(hp_fit$seed)) {
        hp_fit$seed <- as.integer(seed)
      }
      set.seed(as.integer(seed))
      model <- .srb_call_fit(fit_fn, X_train, y_train, meta_train, hp_fit)
      OmicSelector::singlesample_score_call(method_id, model, X_test, meta_test, roster = roster)
    }
  }, error = function(e) structure(rep(NA_real_, length(y_test)), srb_err = conditionMessage(e)))

  err <- attr(s_meth, "srb_err")
  ok <- is.finite(s_meth) & is.finite(s_base) & !is.na(y_test)
  if (sum(ok) < 4L || length(unique(s_meth[ok])) < 2L ||
      length(unique(s_base[ok])) < 2L || length(unique(y_test[ok])) < 2L) {
    res <- .srb_transfer_ineligible(method_id, held_out_id,
            paste0("held-out scores constant/all-NA",
                   if (!is.null(err)) paste0(" (", err, ")") else ""),
            baseline_mode, seed)
    res$n_train_cohorts <- length(train_ids); res$train_cohorts <- train_ids
    res$n_candidate_features <- length(candidates)
    res$n_held_candidate_features <- n_held_candidate_features
    res$candidate_fill_policy <- fill_policy
    .srb_save_ckpt(res, ckpt); return(invisible(res))
  }

  # Single held-out cohort: lift = oriented AUC_method - AUC_baseline on the
  # held-out samples; lift_se = sqrt(DeLong paired-diff variance) on that cohort.
  fold_dt <- data.frame(fold = 1L, y = as.integer(y_test[ok]),
                        score_m = as.numeric(s_meth[ok]), score_b = as.numeric(s_base[ok]))
  vd <- .delong_var_diff_fold(fold_dt$y, fold_dt$score_m, fold_dt$score_b, orient = orient)
  res <- list(
    method_id = method_id, cohort = held_out_id, estimand = "transfer",
    est = vd$auc_m - vd$auc_b, se = sqrt(max(0, vd$var)), df = NA_real_,
    n_eff = sum(ok), auc_method = vd$auc_m, auc_baseline = vd$auc_b,
    n_train_cohorts = length(train_ids), train_cohorts = train_ids,
    n_candidate_features = length(candidates),
    n_held_candidate_features = n_held_candidate_features,
    candidate_fill_policy = fill_policy,
    held_out_auc_gt_0_95 = isTRUE(vd$auc_m > 0.95 || vd$auc_b > 0.95),
    held_out_block = ho_block, eligible = TRUE, ineligible_reason = NA_character_,
    baseline_mode = baseline_mode, seed = seed)
  .srb_save_ckpt(res, ckpt)
  invisible(res)
}

.srb_transfer_ineligible <- function(method_id, held_out_id, reason, baseline_mode, seed) {
  list(method_id = method_id, cohort = held_out_id, estimand = "transfer",
       est = NA_real_, se = NA_real_, df = NA_real_, n_eff = NA_integer_,
       auc_method = NA_real_, auc_baseline = NA_real_,
       n_train_cohorts = NA_integer_, train_cohorts = character(0),
       held_out_block = NA_character_, eligible = FALSE,
       ineligible_reason = reason, baseline_mode = baseline_mode, seed = seed)
}

.srb_safe <- function(s) gsub("[^A-Za-z0-9._-]", "-", as.character(s))

# --- analysis-safe checkpoint I/O (atomic write + fail-open read) ---------------
# Atomic checkpoint write: serialise to a unique temp file in the authorized panel_summary/tmp directory,
# so a SIGKILL (e.g. an OOM kill) DURING the write leaves EITHER the prior file or no
# file -- never a truncated/corrupt .rds that a later readRDS would choke on. The temp
# name carries the PID, so concurrent disjoint-cell tranches never collide and a crash
# leaves only an orphan "<ckpt>.tmp-<pid>" (ignored by the "*.rds" analysis scan).
.srb_save_ckpt <- function(res, ckpt) {
  tmp <- tempfile(pattern = "rank_checkpoint_", tmpdir = Sys.getenv("TMPDIR"))
  saveRDS(res, tmp)
  if (!file.rename(tmp, ckpt)) {
    # rename failed (rare: cross-device tmp, permissions). Fall back to a direct
    # write so the cell is not lost, then drop the orphan temp.
    saveRDS(res, ckpt)
    unlink(tmp)
  }
  invisible(ckpt)
}

# analysis-guard read that FAILS OPEN to recompute (never swallows a corrupt cell).
# A checkpoint that exists but is unreadable (e.g. a pre-fix truncated file from an
# OOM-killed run) returns NULL, so the caller recomputes and atomically rewrites it --
# instead of the readRDS error propagating to the driver's per-cell tryCatch and being
# recorded as a SILENT eligible=FALSE "cell error" row in the assembled tables. No
# real per-cell result is ever NULL (always a list), so NULL unambiguously == unreadable.
.srb_read_ckpt <- function(ckpt) {
  tryCatch(.public_readRDS(ckpt), error = function(e) NULL)
}

# --- table assembly ----------------------------------------------------------
# reassay_design within-cohort column contract (matches supplement/tables/
# table_reassay_design_per_cell_results.tsv key fields the W6 ALL_COHORT_ROWS bind needs:
# Accession, Scoring method, est(lift), se, df, n_eff + audit fields). We emit a
# normalised contract with the canonical W6 keys plus audit columns.
SRB_S3_COLUMNS <- c("method", "cohort", "estimand", "est", "se", "df", "n_eff",
                    "n_valid_folds", "min_pos", "min_neg", "grouping_policy",
                    "eligible", "ineligible_reason", "baseline_mode")

srb_assemble_within_table <- function(results) {
  rows <- lapply(results, function(r) data.frame(
    method = r$method_id, cohort = r$cohort, estimand = "within",
    est = r$est %srb||% NA_real_, se = r$se %srb||% NA_real_,
    df = r$df %srb||% NA_real_, n_eff = r$n_eff %srb||% NA_integer_,
    n_valid_folds = r$n_valid_folds %srb||% NA_integer_,
    min_pos = r$min_pos %srb||% NA_integer_, min_neg = r$min_neg %srb||% NA_integer_,
    grouping_policy = r$grouping_policy %srb||% NA_character_,
    eligible = isTRUE(r$eligible),
    ineligible_reason = r$ineligible_reason %srb||% NA_character_,
    baseline_mode = r$baseline_mode %srb||% NA_character_,
    stringsAsFactors = FALSE))
  out <- do.call(rbind, rows)
  out[, SRB_S3_COLUMNS, drop = FALSE]
}

SRB_S6X_COLUMNS <- c("method", "cohort", "estimand", "est", "se", "df", "n_eff",
                     "auc_method", "auc_baseline", "n_train_cohorts",
                     "held_out_block", "eligible", "ineligible_reason", "baseline_mode")

srb_assemble_transfer_table <- function(results) {
  rows <- lapply(results, function(r) data.frame(
    method = r$method_id, cohort = r$cohort, estimand = "transfer",
    est = r$est %srb||% NA_real_, se = r$se %srb||% NA_real_,
    df = r$df %srb||% NA_real_, n_eff = r$n_eff %srb||% NA_integer_,
    auc_method = r$auc_method %srb||% NA_real_, auc_baseline = r$auc_baseline %srb||% NA_real_,
    n_train_cohorts = r$n_train_cohorts %srb||% NA_integer_,
    held_out_block = r$held_out_block %srb||% NA_character_,
    eligible = isTRUE(r$eligible),
    ineligible_reason = r$ineligible_reason %srb||% NA_character_,
    baseline_mode = r$baseline_mode %srb||% NA_character_,
    stringsAsFactors = FALSE))
  out <- do.call(rbind, rows)
  out[, SRB_S6X_COLUMNS, drop = FALSE]
}

srb_sanitize_character_fields <- function(df) {
  out <- if (inherits(df, "data.table")) data.table::copy(df) else df
  # Sanitise character columns: free-text ineligible_reason values can carry embedded
  # newlines/tabs/CRs (e.g. multi-line fit-error messages, "; "-joined per-fold skip
  # summaries). With quote=FALSE those would split a logical row across multiple physical
  # lines and corrupt the field count, making the TSV unparseable by read_tsv/awk. Collapse
  # any run of \r\n\t into a single space so every row is exactly one physical line.
  char_cols <- vapply(out, is.character, logical(1L))
  for (j in which(char_cols))
    out[[j]] <- gsub("[\r\n\t]+", " ", out[[j]], perl = TRUE)
  out
}

srb_write_tsv <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  df <- srb_sanitize_character_fields(df)
  .public_write_table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}
