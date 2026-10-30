# Shared disease-matched transfer design for the public single-sample benchmark.
# Keeping this definition in one reviewed source prevents the runner and shard
# reconciler from silently accepting different direction sets.

srb_transfer_disease_sets <- function() {
  list(
    colorectal_cancer = c("GSE188627", "GSE198834"),
    nsclc = c("GSE188232", "GSE204951", "GSE279209", "GSE64591", "GSE70080"),
    gastric_cancer = c("GSE130654", "GSE227778"),
    hepatocellular_carcinoma = c("GSE50013", "GSE83977"),
    pancreatic_ductal_adenocarcinoma = c(
      "GSE259327_DIAGNOSTIC", "GSE304572", "GSE85589"
    ),
    breast_cancer = c("GSE210329", "GSE270497", "GSE320031", "GSE41922")
  )
}
# Frozen outer-fold contract for the 34 public analysis units. These eight
# units require K=3 to retain both outcome classes and the prespecified minimum
# per-class test count; every other included unit uses K=5.
srb_expected_outer_k <- function(cohort_id) {
  k3 <- c("GSE117063", "GSE130654", "GSE200489", "GSE204951",
          "GSE210546", "GSE227778", "GSE266859", "GSE83977")
  cohort_id <- as.character(cohort_id)
  stats::setNames(ifelse(cohort_id %in% k3, 3L, 5L), cohort_id)
}
