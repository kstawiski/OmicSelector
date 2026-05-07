#' @title Paper 3 specimen-overlap provenance pre-flight gate
#'
#' @description
#' Audits a candidate accession list against a manifest of known specimen-
#' overlap pairs and returns one of three exit conditions:
#' \describe{
#'   \item{\code{NO_OVERLAP}}{No accession in the input list overlaps with any
#'     other accession in the input list according to the manifest, and every
#'     supplied accession appears (as either side of a pair) somewhere in the
#'     manifest.}
#'   \item{\code{KNOWN_OVERLAP}}{At least one pair of input accessions is
#'     listed in the manifest. The cross-cohort analysis must apply
#'     block-aware multiplicity correction (\code{\link{paper3_bh_fdr_correct_blocked}}).}
#'   \item{\code{UNKNOWN_ACCESSION}}{At least one input accession is not in
#'     the manifest. Audit required before combining with anything else.
#'     Takes precedence over \code{KNOWN_OVERLAP}.}
#' }
#'
#' Intended as the first stage in the Paper-3 pipeline (Methods section
#' "Provenance auditing and specimen-overlap handling").
#'
#' @param accessions Character vector of GEO/ArrayExpress accessions
#'   (e.g., \code{c("GSE211692", "GSE106817")}).
#' @param manifest_path Path to a TSV file with columns
#'   \code{acc1, acc2, block_id} (and optionally \code{fraction},
#'   \code{evidence_source}). If \code{NULL}, uses the bundled manifest at
#'   \code{system.file("extdata", "provenance_manifest.tsv", package = "OmicSelector")}.
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{status}}{One of \code{"NO_OVERLAP"}, \code{"KNOWN_OVERLAP"},
#'       \code{"UNKNOWN_ACCESSION"}.}
#'     \item{\code{overlapping_pairs}}{\code{data.frame} subset of the manifest
#'       with rows whose \code{acc1} and \code{acc2} are both in the input
#'       accession list.}
#'     \item{\code{unknown_accessions}}{Character vector of input accessions
#'       not present in the manifest.}
#'     \item{\code{manifest_path}}{Resolved path of the manifest used.}
#'   }
#'
#' @examples
#' res <- os_provenance_preflight(c("GSE211692", "GSE106817"))
#' res$status
#'
#' @export
os_provenance_preflight <- function(accessions, manifest_path = NULL) {
  if (is.null(manifest_path)) {
    manifest_path <- system.file("extdata", "provenance_manifest.tsv",
                                  package = "OmicSelector")
  }
  if (!nzchar(manifest_path) || !file.exists(manifest_path)) {
    stop("Provenance manifest not found: ",
         if (nzchar(manifest_path)) manifest_path else "(empty path)")
  }
  manifest <- utils::read.delim(manifest_path, stringsAsFactors = FALSE)
  required_cols <- c("acc1", "acc2", "block_id")
  missing_cols <- setdiff(required_cols, names(manifest))
  if (length(missing_cols) > 0L) {
    stop("Manifest missing columns: ", paste(missing_cols, collapse = ", "))
  }

  accessions <- as.character(accessions)
  known_accs <- unique(c(manifest$acc1, manifest$acc2))
  unknown <- setdiff(accessions, known_accs)

  # Self-pair rows (acc1 == acc2) register an accession as "audited
  # independent" without implying overlap. Treat them as registration
  # records, not overlap evidence.
  is_self_pair <- manifest$acc1 == manifest$acc2
  in_input <- manifest$acc1 %in% accessions & manifest$acc2 %in% accessions
  overlapping_pairs <- manifest[in_input & !is_self_pair, , drop = FALSE]
  rownames(overlapping_pairs) <- NULL

  status <- if (length(unknown) > 0L) {
    "UNKNOWN_ACCESSION"
  } else if (nrow(overlapping_pairs) > 0L) {
    "KNOWN_OVERLAP"
  } else {
    "NO_OVERLAP"
  }

  list(
    status = status,
    overlapping_pairs = overlapping_pairs,
    unknown_accessions = unknown,
    manifest_path = manifest_path
  )
}
