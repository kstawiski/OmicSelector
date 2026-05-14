#' @title MIMAT accession \eqn{\leftrightarrow} hsa-miR name lookup
#'
#' @description
#' Functions to translate between MIMAT accession identifiers (e.g.
#' \code{MIMAT0000062}) and canonical \emph{hsa-miR-*} mature miRNA names
#' (e.g. \code{hsa-let-7a-5p}) using the bundled miRBase 22.1 lookup table.
#'
#' The lookup table is shipped at
#' \code{inst/extdata/mimat_hsa_lookup_v22_1.tsv} and covers all 2 656 human
#' mature miRNAs in miRBase release 22.1.  Each MIMAT ID maps to exactly one
#' primary name (1-to-1; no ambiguity in the source data).
#'
#' @name mimat_hsa_lookup
NULL

# Package-level lazy singleton so the TSV is read only once per session.
.mimat_env <- new.env(parent = emptyenv())

.get_lookup <- function() {
  if (!exists("tbl", envir = .mimat_env)) {
    path <- system.file(
      "extdata", "mimat_hsa_lookup_v22_1.tsv",
      package = "OmicSelector",
      mustWork = TRUE
    )
    tbl <- utils::read.delim(path, comment.char = "#", stringsAsFactors = FALSE)
    assign("tbl", tbl, envir = .mimat_env)
  }
  get("tbl", envir = .mimat_env)
}

#' Resolve MIMAT accession IDs to canonical hsa-miR-* names
#'
#' @param ids Character vector of MIMAT accession identifiers, e.g.
#'   \code{c("MIMAT0000062", "MIMAT0004481")}.
#' @return Character vector (same length as \code{ids}) of canonical
#'   hsa-miR-* names.  Entries not found in the lookup table are returned as
#'   \code{NA_character_}.
#' @details
#' Matching is case-insensitive on input; the returned value is always in
#' canonical case (MIMAT IDs uppercase, e.g. \code{MIMAT0000062}; hsa-miR
#' names lowercase, e.g. \code{hsa-let-7a-5p}). Unknown inputs return
#' \code{NA_character_}.
#' @examples
#' resolve_mimat_to_hsa(c("MIMAT0000062", "MIMAT9999999"))
#' # [1] "hsa-let-7a-5p" NA
#' @export
resolve_mimat_to_hsa <- function(ids) {
  stopifnot(is.character(ids))
  tbl <- .get_lookup()
  idx <- match(toupper(ids), toupper(tbl$mimat_id))
  tbl$hsa_mir_name[idx]
}

#' Resolve canonical hsa-miR-* names to MIMAT accession IDs
#'
#' @param names Character vector of canonical hsa-miR-* mature miRNA names,
#'   e.g. \code{c("hsa-let-7a-5p", "hsa-miR-21-5p")}.
#' @return Character vector (same length as \code{names}) of MIMAT accession
#'   IDs.  Entries not found in the lookup table are returned as
#'   \code{NA_character_}.
#' @details
#' Matching is case-insensitive on input; the returned value is always in
#' canonical case (MIMAT IDs uppercase, e.g. \code{MIMAT0000062}; hsa-miR
#' names lowercase, e.g. \code{hsa-let-7a-5p}). Unknown inputs return
#' \code{NA_character_}.
#' @examples
#' resolve_hsa_to_mimat(c("hsa-let-7a-5p", "hsa-not-real"))
#' # [1] "MIMAT0000062" NA
#' @export
resolve_hsa_to_mimat <- function(names) {
  stopifnot(is.character(names))
  tbl <- .get_lookup()
  idx <- match(tolower(names), tolower(tbl$hsa_mir_name))
  tbl$mimat_id[idx]
}
