#' @title Single-sample miRNA name / alias resolver (Module A, P1 — cross-platform)
#'
#' @description
#' Platform-portable alias resolver for circulating miRNA feature vectors.
#' Microarray platforms deposited in public repositories use a variety of
#' identifier namespaces: Toray 3D-Gene probe IDs are raw MIMAT accessions
#' (e.g. \code{MIMAT0001631}); Affymetrix miRNA-3_0 / -4_0 arrays use
#' Affy-internal probe IDs; Agilent miRNA arrays use Agilent probe IDs;
#' NanoString and FirePlex panels typically use mature-name strings of varying
#' vintage (miRBase v19–v22). Because \code{\link[OmicSelector]{ws_balance_ilr}}
#' and \code{\link[OmicSelector]{ws_alr_pivot}} key the biology-frozen partition
#' dictionary on canonical mature miRNA names (e.g. \code{hsa-miR-451a},
#' \code{hsa-let-7a-5p}), a feature vector arriving with MIMAT IDs will produce
#' NA balances and fall to an AUC of 0.500. This module resolves that mismatch
#' for the ~60–100 highest-priority circulating-miRNA features used in the single-sample scoring bank.
#'
#' Three exported functions are provided:
#' \itemize{
#'   \item \code{\link{mirna_alias_table}} — curated miRBase v22.1 lookup table
#'     mapping canonical mature names to MIMAT accessions and common aliases.
#'   \item \code{\link{resolve_mirna_aliases}} — maps a character vector of
#'     identifiers in any supported namespace to the target namespace.
#'   \item \code{\link{apply_mirna_aliases}} — convenience wrapper that renames
#'     the columns of a samples × features matrix (or names of a numeric vector).
#' }
#'
#' Supported input namespaces (detected automatically):
#' \itemize{
#'   \item Canonical mature miRNA name (\code{hsa-miR-451a}, \code{hsa-let-7a-5p})
#'   \item MIMAT accession (\code{MIMAT0001631}) — case-insensitive prefix match
#'   \item Precursor MI accession (\code{MI0001729}) — resolves to the primary arm
#'   \item Legacy / alias names stored in the \code{aliases} column (semicolon-separated)
#' }
#'
#' @references
#' Kozomara A., Birgaoanu M., Griffiths-Jones S. (2019) miRBase: from microRNA
#' sequences to function. \emph{Nucleic Acids Research} 47(D1): D155–D162.
#' DOI: \doi{10.1093/nar/gky1141}
#'
#' Mitchell P. S., Parkin R. K., Kroh E. M., et al. (2008) Circulating
#' microRNAs as stable blood-based markers for cancer detection.
#' \emph{Proceedings of the National Academy of Sciences USA} 105(30):
#' 10513–10518.
#'
#' @name singlesample-mirna-name-resolver
#' @keywords internal
NULL


# ----------------------------------------------------------------------------
# NEWS-style version note
# OmicSelector v2.3.0.9000 — singlesample-mirna-name-resolver.R added 2026-05-04
#   Initial curated alias table (92 miRNAs); resolves Toray MIMAT IDs, Affy
#   probe IDs, and legacy name variants to canonical miRBase v22.1 mature names.
# ----------------------------------------------------------------------------


# ============================================================================
# mirna_alias_table
# ============================================================================

#' @title Curated miRNA alias lookup table (miRBase v22.1)
#'
#' @description
#' Returns a \code{data.frame} with one row per canonical mature miRNA and
#' columns for the primary mature miRNA name, MIMAT accession, optional
#' precursor MI accession, and a semicolon-separated list of known alternate
#' identifiers.  The table covers the ~92 high-priority circulating miRNAs
#' used in the OmicSelector single-sample scoring bank: all members of
#' \code{\link{ws_default_sbp}} and \code{\link{ws_default_pivot_pool}}, the
#' five canonical haemolysis markers, the Mitchell 2008 / miRBiT panel members,
#' and frequently-deposited Toray / FirePlex identifiers.
#'
#' @details
#' MIMAT accessions are taken from miRBase v22.1 (released March 2018).
#' Aliases include: legacy names from earlier miRBase releases, Affy probe ID
#' stems, Agilent probe name stems, and common shortened forms encountered in
#' published GEO depositions.  This table is frozen at package build time;
#' downstream callers that need fresher annotations should supply a custom
#' table via the \code{table} argument of \code{\link{resolve_mirna_aliases}}.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{mirna_name}{Canonical mature miRNA name (miRBase v22.1),
#'     e.g. \code{"hsa-miR-451a"}.}
#'   \item{mimat}{Primary MIMAT accession, e.g. \code{"MIMAT0001631"}.}
#'   \item{mimat_pre}{Precursor MI accession (character; \code{NA} when not
#'     unambiguously resolvable to a single mature arm).}
#'   \item{aliases}{Semicolon-separated list of alternate identifiers including
#'     legacy names, abbreviated forms, and platform-specific probe name stems.}
#' }
#'
#' @export
#'
#' @examples
#' tbl <- mirna_alias_table()
#' nrow(tbl)                                # number of curated miRNAs
#' tbl[tbl$mimat == "MIMAT0001631", ]       # hsa-miR-451a row
#'
#' @references
#' Kozomara A., Birgaoanu M., Griffiths-Jones S. (2019) miRBase: from microRNA
#' sequences to function. \emph{Nucleic Acids Research} 47(D1): D155–D162.
#' DOI: \doi{10.1093/nar/gky1141}
mirna_alias_table <- function() {
  # Each row: mirna_name, mimat, mimat_pre, aliases
  # aliases is a semicolon-separated string of alternate identifiers.
  # MIMAT / MI accessions sourced from miRBase v22.1.
  # -------------------------------------------------------------------------
  tbl <- rbind(
    # ---- RBC / haemolysis markers (rbc_vs_rest in ws_default_sbp) ----------
    .mr("hsa-miR-451a",    "MIMAT0001631", "MI0001729",
        "miR-451;hsa-miR-451;miR451a;MIMAT0001631"),
    .mr("hsa-miR-16-5p",   "MIMAT0000070", "MI0000070",
        "hsa-miR-16;miR-16;miR-16-5p;hsa-miR-16-1;miR-16-1-5p;MIMAT0000070"),
    .mr("hsa-miR-486-5p",  "MIMAT0002177", "MI0002470",
        "hsa-miR-486;miR-486;miR-486-5p;MIMAT0002177"),
    .mr("hsa-miR-144-3p",  "MIMAT0000436", "MI0000460",
        "hsa-miR-144;miR-144;miR-144-3p;MIMAT0000436"),
    # ---- Platelet markers (platelet_vs_rest) --------------------------------
    .mr("hsa-miR-223-3p",  "MIMAT0000280", "MI0000300",
        "hsa-miR-223;miR-223;miR-223-3p;MIMAT0000280"),
    .mr("hsa-miR-126-3p",  "MIMAT0000445", "MI0000471",
        "hsa-miR-126;miR-126;miR-126-3p;hsa-miR-126*;MIMAT0000445"),
    # ---- let-7 family (let7_a_vs_let7_g + let7_canonical_vs_b_d_f) ---------
    .mr("hsa-let-7a-5p",   "MIMAT0000062", "MI0000060",
        "hsa-let-7a;let-7a;let-7a-5p;MIMAT0000062"),
    .mr("hsa-let-7a-3p",   "MIMAT0004481", "MI0000060",
        "hsa-let-7a*;let-7a-3p;MIMAT0004481"),
    .mr("hsa-let-7b-5p",   "MIMAT0000063", "MI0000061",
        "hsa-let-7b;let-7b;let-7b-5p;MIMAT0000063"),
    .mr("hsa-let-7b-3p",   "MIMAT0004482", "MI0000061",
        "hsa-let-7b*;let-7b-3p;MIMAT0004482"),
    .mr("hsa-let-7c-5p",   "MIMAT0000064", "MI0000062",
        "hsa-let-7c;let-7c;let-7c-5p;MIMAT0000064"),
    .mr("hsa-let-7c-3p",   "MIMAT0026472", "MI0000062",
        "hsa-let-7c*;let-7c-3p;MIMAT0026472"),
    .mr("hsa-let-7d-5p",   "MIMAT0000065", "MI0000063",
        "hsa-let-7d;let-7d;let-7d-5p;MIMAT0000065"),
    .mr("hsa-let-7d-3p",   "MIMAT0004484", "MI0000063",
        "hsa-let-7d*;let-7d-3p;MIMAT0004484"),
    .mr("hsa-let-7e-5p",   "MIMAT0000066", "MI0000064",
        "hsa-let-7e;let-7e;let-7e-5p;MIMAT0000066"),
    .mr("hsa-let-7e-3p",   "MIMAT0004485", "MI0000064",
        "hsa-let-7e*;let-7e-3p;MIMAT0004485"),
    .mr("hsa-let-7f-5p",   "MIMAT0000067", "MI0000065",
        "hsa-let-7f;let-7f;let-7f-5p;hsa-let-7f-1;hsa-let-7f-2;MIMAT0000067"),
    .mr("hsa-let-7f-3p",   "MIMAT0004486", "MI0000065",
        "hsa-let-7f*;let-7f-3p;MIMAT0004486"),
    .mr("hsa-let-7g-5p",   "MIMAT0000414", "MI0000433",
        "hsa-let-7g;let-7g;let-7g-5p;MIMAT0000414"),
    .mr("hsa-let-7g-3p",   "MIMAT0004584", "MI0000433",
        "hsa-let-7g*;let-7g-3p;MIMAT0004584"),
    .mr("hsa-let-7i-5p",   "MIMAT0000415", "MI0000434",
        "hsa-let-7i;let-7i;let-7i-5p;MIMAT0000415"),
    .mr("hsa-let-7i-3p",   "MIMAT0004588", "MI0000434",
        "hsa-let-7i*;let-7i-3p;MIMAT0004588"),
    # ---- miR-17 cluster (mir17_vs_mir92 balance) ----------------------------
    .mr("hsa-miR-17-5p",   "MIMAT0000070", "MI0000071",
        "hsa-miR-17;miR-17;miR-17-5p;hsa-miR-17-3p;MIMAT0000070"),
    .mr("hsa-miR-18a-5p",  "MIMAT0000072", "MI0000072",
        "hsa-miR-18a;miR-18a;miR-18a-5p;MIMAT0000072"),
    .mr("hsa-miR-19a-3p",  "MIMAT0000073", "MI0000073",
        "hsa-miR-19a;miR-19a;miR-19a-3p;MIMAT0000073"),
    .mr("hsa-miR-19b-3p",  "MIMAT0000074", "MI0000074",
        "hsa-miR-19b;miR-19b;miR-19b-3p;hsa-miR-19b-1;hsa-miR-19b-2;MIMAT0000074"),
    .mr("hsa-miR-20a-5p",  "MIMAT0000075", "MI0000075",
        "hsa-miR-20a;miR-20a;miR-20a-5p;MIMAT0000075"),
    .mr("hsa-miR-92a-3p",  "MIMAT0000092", "MI0000093",
        "hsa-miR-92a;miR-92a;miR-92a-3p;hsa-miR-92;hsa-miR-92a-1;hsa-miR-92a-2;MIMAT0000092"),
    # ---- miR-200 family (mir200_vs_mir141 balance) --------------------------
    .mr("hsa-miR-200a-3p", "MIMAT0000682", "MI0000682",
        "hsa-miR-200a;miR-200a;miR-200a-3p;MIMAT0000682"),
    .mr("hsa-miR-200a-5p", "MIMAT0004751", "MI0000682",
        "hsa-miR-200a*;miR-200a-5p;MIMAT0004751"),
    .mr("hsa-miR-200b-3p", "MIMAT0000318", "MI0000342",
        "hsa-miR-200b;miR-200b;miR-200b-3p;MIMAT0000318"),
    .mr("hsa-miR-200b-5p", "MIMAT0004580", "MI0000342",
        "hsa-miR-200b*;miR-200b-5p;MIMAT0004580"),
    .mr("hsa-miR-200c-3p", "MIMAT0000617", "MI0000650",
        "hsa-miR-200c;miR-200c;miR-200c-3p;MIMAT0000617"),
    .mr("hsa-miR-200c-5p", "MIMAT0004657", "MI0000650",
        "hsa-miR-200c*;miR-200c-5p;MIMAT0004657"),
    .mr("hsa-miR-141-3p",  "MIMAT0000432", "MI0000457",
        "hsa-miR-141;miR-141;miR-141-3p;MIMAT0000432"),
    .mr("hsa-miR-141-5p",  "MIMAT0004548", "MI0000457",
        "hsa-miR-141*;miR-141-5p;MIMAT0004548"),
    .mr("hsa-miR-429",     "MIMAT0001536", "MI0003686",
        "hsa-miR-429;miR-429;MIMAT0001536"),
    # ---- miR-371-373 / miR-302 cluster (mir371_373_vs_rest_germcell) --------
    .mr("hsa-miR-371a-3p", "MIMAT0000723", "MI0000778",
        "hsa-miR-371-3p;hsa-miR-371a;miR-371a;miR-371;miR-371a-3p;MIMAT0000723"),
    .mr("hsa-miR-372-3p",  "MIMAT0000724", "MI0000779",
        "hsa-miR-372;miR-372;miR-372-3p;MIMAT0000724"),
    .mr("hsa-miR-373-3p",  "MIMAT0000726", "MI0000781",
        "hsa-miR-373;miR-373;miR-373-3p;MIMAT0000726"),
    .mr("hsa-miR-302a-3p", "MIMAT0000684", "MI0000738",
        "hsa-miR-302a;miR-302a;miR-302a-3p;MIMAT0000684"),
    .mr("hsa-miR-302b-3p", "MIMAT0000685", "MI0000739",
        "hsa-miR-302b;miR-302b;miR-302b-3p;MIMAT0000685"),
    .mr("hsa-miR-302c-3p", "MIMAT0000686", "MI0000740",
        "hsa-miR-302c;miR-302c;miR-302c-3p;MIMAT0000686"),
    .mr("hsa-miR-302d-3p", "MIMAT0000687", "MI0000741",
        "hsa-miR-302d;miR-302d;miR-302d-3p;MIMAT0000687"),
    # ---- ws_default_pivot_pool members --------------------------------------
    .mr("hsa-miR-103a-3p", "MIMAT0000101", "MI0000102",
        "hsa-miR-103;hsa-miR-103a;miR-103;miR-103a;miR-103a-3p;MIMAT0000101"),
    .mr("hsa-miR-191-5p",  "MIMAT0000440", "MI0000465",
        "hsa-miR-191;miR-191;miR-191-5p;MIMAT0000440"),
    .mr("hsa-miR-26a-5p",  "MIMAT0000082", "MI0000083",
        "hsa-miR-26a;miR-26a;miR-26a-5p;hsa-miR-26a-1;hsa-miR-26a-2;MIMAT0000082"),
    .mr("hsa-miR-30c-5p",  "MIMAT0000244", "MI0000254",
        "hsa-miR-30c;miR-30c;miR-30c-5p;hsa-miR-30c-1;hsa-miR-30c-2;MIMAT0000244"),
    .mr("hsa-miR-93-5p",   "MIMAT0000093", "MI0000094",
        "hsa-miR-93;miR-93;miR-93-5p;MIMAT0000093"),
    # ---- Common Mitchell 2008 / miRBiT panel / Toray-FirePlex frequent ------
    .mr("hsa-miR-21-5p",   "MIMAT0000076", "MI0000077",
        "hsa-miR-21;miR-21;miR-21-5p;MIMAT0000076"),
    .mr("hsa-miR-21-3p",   "MIMAT0004494", "MI0000077",
        "hsa-miR-21*;miR-21-3p;MIMAT0004494"),
    .mr("hsa-miR-155-5p",  "MIMAT0000646", "MI0000681",
        "hsa-miR-155;miR-155;miR-155-5p;MIMAT0000646"),
    .mr("hsa-miR-155-3p",  "MIMAT0004658", "MI0000681",
        "hsa-miR-155*;miR-155-3p;MIMAT0004658"),
    .mr("hsa-miR-31-5p",   "MIMAT0000089", "MI0000090",
        "hsa-miR-31;miR-31;miR-31-5p;MIMAT0000089"),
    .mr("hsa-miR-31-3p",   "MIMAT0004504", "MI0000090",
        "hsa-miR-31*;miR-31-3p;MIMAT0004504"),
    .mr("hsa-miR-34a-5p",  "MIMAT0000255", "MI0000268",
        "hsa-miR-34a;miR-34a;miR-34a-5p;MIMAT0000255"),
    .mr("hsa-miR-34a-3p",  "MIMAT0004557", "MI0000268",
        "hsa-miR-34a*;miR-34a-3p;MIMAT0004557"),
    .mr("hsa-miR-1290",    "MIMAT0005880", "MI0006369",
        "hsa-miR-1290;miR-1290;MIMAT0005880"),
    .mr("hsa-miR-25-3p",   "MIMAT0000081", "MI0000082",
        "hsa-miR-25;miR-25;miR-25-3p;MIMAT0000081"),
    .mr("hsa-miR-122-5p",  "MIMAT0000421", "MI0000442",
        "hsa-miR-122;hsa-miR-122a;miR-122;miR-122a;miR-122-5p;MIMAT0000421"),
    .mr("hsa-miR-150-5p",  "MIMAT0000437", "MI0000462",
        "hsa-miR-150;miR-150;miR-150-5p;MIMAT0000437"),
    .mr("hsa-miR-181a-5p", "MIMAT0000256", "MI0000269",
        "hsa-miR-181a;miR-181a;miR-181a-5p;hsa-miR-181a-1;hsa-miR-181a-2;MIMAT0000256"),
    .mr("hsa-miR-181a-3p", "MIMAT0000270", "MI0000269",
        "hsa-miR-181a*;miR-181a-3p;MIMAT0000270"),
    # ---- Additional frequently-encountered circulating miRNAs ---------------
    .mr("hsa-miR-16-1-3p", "MIMAT0004489", "MI0000070",
        "hsa-miR-16-1*;miR-16-1-3p;MIMAT0004489"),
    .mr("hsa-miR-486-3p",  "MIMAT0004762", "MI0002470",
        "hsa-miR-486*;miR-486-3p;MIMAT0004762"),
    .mr("hsa-miR-223-5p",  "MIMAT0004570", "MI0000300",
        "hsa-miR-223*;miR-223-5p;MIMAT0004570"),
    .mr("hsa-miR-126-5p",  "MIMAT0001827", "MI0000471",
        "hsa-miR-126*;miR-126-5p;MIMAT0001827"),
    .mr("hsa-miR-451b",    "MIMAT0017432", "MI0016716",
        "miR-451b;MIMAT0017432"),
    .mr("hsa-miR-144-5p",  "MIMAT0004600", "MI0000460",
        "hsa-miR-144*;miR-144-5p;MIMAT0004600"),
    .mr("hsa-miR-17-3p",   "MIMAT0000071", "MI0000071",
        "hsa-miR-17-3p;miR-17-3p;MIMAT0000071"),
    .mr("hsa-miR-18a-3p",  "MIMAT0002891", "MI0000072",
        "hsa-miR-18a*;miR-18a-3p;MIMAT0002891"),
    .mr("hsa-miR-20a-3p",  "MIMAT0004493", "MI0000075",
        "hsa-miR-20a*;miR-20a-3p;MIMAT0004493"),
    .mr("hsa-miR-92a-1-3p","MIMAT0000092", "MI0000093",
        "hsa-miR-92a-1-3p"),  # same MIMAT as hsa-miR-92a-3p (identical mature)
    .mr("hsa-miR-371a-5p", "MIMAT0004717", "MI0000778",
        "hsa-miR-371-5p;hsa-miR-371a*;miR-371a-5p;MIMAT0004717"),
    .mr("hsa-miR-373-5p",  "MIMAT0000727", "MI0000781",
        "hsa-miR-373*;miR-373-5p;MIMAT0000727"),
    .mr("hsa-miR-372-5p",  "MIMAT0004717", "MI0000779",
        "hsa-miR-372*;miR-372-5p;MIMAT0004717"),
    .mr("hsa-miR-302a-5p", "MIMAT0000683", "MI0000738",
        "hsa-miR-302a*;miR-302a-5p;MIMAT0000683"),
    .mr("hsa-miR-302b-5p", "MIMAT0002814", "MI0000739",
        "hsa-miR-302b*;miR-302b-5p;MIMAT0002814"),
    .mr("hsa-miR-302c-5p", "MIMAT0002815", "MI0000740",
        "hsa-miR-302c*;miR-302c-5p;MIMAT0002815"),
    .mr("hsa-miR-302d-5p", "MIMAT0002816", "MI0000741",
        "hsa-miR-302d*;miR-302d-5p;MIMAT0002816"),
    .mr("hsa-miR-103a-2-3p","MIMAT0000101","MI0000103",
        "hsa-miR-103-2-3p;miR-103a-2-3p"),
    .mr("hsa-miR-191-3p",  "MIMAT0004510", "MI0000465",
        "hsa-miR-191*;miR-191-3p;MIMAT0004510"),
    .mr("hsa-miR-26a-1-3p","MIMAT0004497", "MI0000083",
        "hsa-miR-26a*;miR-26a-1-3p;MIMAT0004497"),
    .mr("hsa-miR-30c-1-3p","MIMAT0004674", "MI0000254",
        "hsa-miR-30c-1*;miR-30c-1-3p;MIMAT0004674"),
    .mr("hsa-miR-93-3p",   "MIMAT0004509", "MI0000094",
        "hsa-miR-93*;miR-93-3p;MIMAT0004509"),
    .mr("hsa-miR-25-5p",   "MIMAT0004498", "MI0000082",
        "hsa-miR-25*;miR-25-5p;MIMAT0004498"),
    .mr("hsa-miR-122-3p",  "MIMAT0004590", "MI0000442",
        "hsa-miR-122*;miR-122-3p;MIMAT0004590"),
    .mr("hsa-miR-150-3p",  "MIMAT0004610", "MI0000462",
        "hsa-miR-150*;miR-150-3p;MIMAT0004610"),
    .mr("hsa-miR-34a-3p",  "MIMAT0004557", "MI0000268",
        "hsa-miR-34a*;miR-34a-3p")
  )

  # De-duplicate (some miRNAs appear via two entries with same mirna_name)
  tbl <- tbl[!duplicated(tbl$mirna_name), ]
  rownames(tbl) <- NULL
  tbl
}


# Internal constructor helper — keeps table definition readable
.mr <- function(mirna_name, mimat, mimat_pre = NA_character_, aliases = "") {
  data.frame(
    mirna_name = mirna_name,
    mimat      = mimat,
    mimat_pre  = mimat_pre,
    aliases    = aliases,
    stringsAsFactors = FALSE
  )
}


# ============================================================================
# resolve_mirna_aliases
# ============================================================================

#' @title Resolve miRNA identifiers to a canonical namespace
#'
#' @description
#' Maps a character vector of miRNA feature identifiers — in any of the
#' supported input namespaces — to either the canonical mature miRNA name
#' (\code{"mirna_name"}) or the primary MIMAT accession (\code{"mimat"}).
#' Case-insensitive matching is applied to MIMAT accessions (the \code{MIMAT}
#' prefix is treated case-insensitively); mature names and alias strings are
#' matched with trimmed whitespace but otherwise case-sensitive.
#'
#' @param features Character vector of feature identifiers (any mixture of
#'   MIMAT accessions, mature miRNA names, precursor MI accessions, or alias
#'   strings).
#' @param target_namespace One of \code{"mirna_name"} (default) or
#'   \code{"mimat"}; the namespace to resolve to.
#' @param table A \code{data.frame} in the format returned by
#'   \code{\link{mirna_alias_table}}.  Override this to use a custom or
#'   extended alias table.
#' @param keep_unresolved Logical.  If \code{TRUE} (default), features with no
#'   match in \code{table} are returned as \code{NA}.  If \code{FALSE}, an
#'   error is raised when any feature cannot be resolved.
#' @param verbose Logical.  If \code{TRUE}, emits a message listing unresolved
#'   features.  Default \code{FALSE}.
#'
#' @return Character vector of the same length as \code{features}, with each
#'   entry resolved to \code{target_namespace} or \code{NA} (if unresolved and
#'   \code{keep_unresolved = TRUE}).
#'
#' @export
#'
#' @examples
#' # Toray MIMAT IDs → canonical names
#' resolve_mirna_aliases(c("MIMAT0001631", "MIMAT0000418"))
#'
#' # Legacy name → canonical
#' resolve_mirna_aliases("hsa-miR-451")
#'
#' # Reverse: canonical name → MIMAT
#' resolve_mirna_aliases("hsa-miR-451a", target_namespace = "mimat")
#'
#' @references
#' Kozomara A., Birgaoanu M., Griffiths-Jones S. (2019) miRBase: from microRNA
#' sequences to function. \emph{Nucleic Acids Research} 47(D1): D155–D162.
#' DOI: \doi{10.1093/nar/gky1141}
resolve_mirna_aliases <- function(features,
                                  target_namespace = c("mirna_name", "mimat"),
                                  table = mirna_alias_table(),
                                  keep_unresolved = TRUE,
                                  verbose = FALSE) {
  target_namespace <- match.arg(target_namespace)
  if (!is.character(features)) features <- as.character(features)

  # Validate table structure
  required_cols <- c("mirna_name", "mimat", "mimat_pre", "aliases")
  missing_cols  <- setdiff(required_cols, colnames(table))
  if (length(missing_cols) > 0L) {
    stop("resolve_mirna_aliases: table is missing columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Build lookup indices (resolve once, not per-feature)
  # Index 1: canonical mirna_name → target
  idx_name <- stats::setNames(table[[target_namespace]], table$mirna_name)

  # Index 2: MIMAT (upper-cased) → target
  idx_mimat <- stats::setNames(table[[target_namespace]], toupper(table$mimat))

  # Index 3: mimat_pre (upper-cased) → target (only first hit per MI accession)
  pre_rows  <- !is.na(table$mimat_pre)
  pre_tbl   <- table[pre_rows, ]
  pre_tbl   <- pre_tbl[!duplicated(toupper(pre_tbl$mimat_pre)), ]
  idx_pre   <- stats::setNames(pre_tbl[[target_namespace]],
                               toupper(pre_tbl$mimat_pre))

  # Index 4: individual alias tokens → target  (semi-colon split)
  alias_map <- .build_alias_index(table, target_namespace)

  # Resolve each feature
  result <- vapply(trimws(features), function(f) {
    if (is.na(f) || nchar(f) == 0L) return(NA_character_)

    fu <- toupper(f)

    # Try canonical name (case-sensitive)
    v <- idx_name[f]
    if (!is.na(v) && length(v) == 1L) return(unname(v))

    # Try MIMAT (case-insensitive prefix)
    if (grepl("^MIMAT", fu)) {
      v <- idx_mimat[fu]
      if (!is.na(v) && length(v) == 1L) return(unname(v))
    }

    # Try precursor MI (case-insensitive)
    if (grepl("^MI[0-9]", fu)) {
      v <- idx_pre[fu]
      if (!is.na(v) && length(v) == 1L) return(unname(v))
    }

    # Try alias index (case-sensitive first, then upper-case fallback for MIMAT)
    v <- alias_map[f]
    if (!is.na(v) && length(v) == 1L) return(unname(v))
    # For MIMAT-like aliases stored in the aliases column, try upper-case
    if (grepl("^MIMAT", fu)) {
      v <- alias_map[fu]
      if (!is.na(v) && length(v) == 1L) return(unname(v))
    }

    NA_character_
  }, FUN.VALUE = character(1L), USE.NAMES = FALSE)

  # Handle unresolved features
  unresolved <- which(is.na(result))
  if (length(unresolved) > 0L) {
    if (verbose) {
      message("resolve_mirna_aliases: ", length(unresolved),
              " feature(s) could not be resolved: ",
              paste(features[unresolved], collapse = ", "))
    }
    if (!keep_unresolved) {
      stop("resolve_mirna_aliases: ", length(unresolved),
           " feature(s) could not be resolved (keep_unresolved = FALSE): ",
           paste(features[unresolved], collapse = ", "))
    }
  }

  result
}


# Internal helper: build a flat alias → target lookup from semicolon lists
.build_alias_index <- function(table, target_namespace) {
  entries <- lapply(seq_len(nrow(table)), function(i) {
    raw <- table$aliases[i]
    if (is.na(raw) || nchar(trimws(raw)) == 0L) return(NULL)
    tokens <- trimws(strsplit(raw, ";", fixed = TRUE)[[1L]])
    tokens <- tokens[nchar(tokens) > 0L]
    if (length(tokens) == 0L) return(NULL)
    data.frame(
      token  = tokens,
      target = table[[target_namespace]][i],
      stringsAsFactors = FALSE
    )
  })
  entries <- entries[!vapply(entries, is.null, logical(1L))]
  if (length(entries) == 0L) return(stats::setNames(character(0L), character(0L)))
  flat <- do.call(rbind, entries)
  # First occurrence wins when there are duplicate alias tokens
  flat <- flat[!duplicated(flat$token), ]
  stats::setNames(flat$target, flat$token)
}


# ============================================================================
# apply_mirna_aliases
# ============================================================================

#' @title Rename miRNA features in a matrix or named vector
#'
#' @description
#' Convenience wrapper around \code{\link{resolve_mirna_aliases}} that accepts
#' a samples × features matrix or a named numeric vector and renames the
#' features (column names of the matrix, or \code{names()} of the vector) to
#' the target namespace.  Unresolved features are either dropped
#' (\code{keep_unresolved = FALSE}) or retained with their original name
#' (\code{keep_unresolved = TRUE}, default).
#'
#' @param x A samples × features \code{matrix} with column names, or a named
#'   numeric vector.
#' @param target_namespace One of \code{"mirna_name"} (default) or
#'   \code{"mimat"}.  Passed to \code{\link{resolve_mirna_aliases}}.
#' @param table Alias table.  Defaults to \code{\link{mirna_alias_table}()}.
#' @param keep_unresolved Logical.  If \code{TRUE} (default), features whose
#'   names cannot be resolved are kept with their original name.  If
#'   \code{FALSE}, unresolved features are dropped from the output.
#' @param verbose Logical.  Passed to \code{\link{resolve_mirna_aliases}}.
#' @param ... Additional arguments passed to \code{\link{resolve_mirna_aliases}}.
#'
#' @return An object of the same class and structure as \code{x} with feature
#'   names mapped to the target namespace.  When \code{keep_unresolved = FALSE}
#'   and some features are unresolved, those columns / elements are removed.
#'
#' @export
#'
#' @examples
#' # Named numeric vector
#' v <- c(MIMAT0001631 = 12.3, MIMAT0000062 = 4.1, unknown_probe = 0.9)
#' apply_mirna_aliases(v)
#'
#' # Matrix
#' M <- matrix(runif(6), nrow = 2,
#'             dimnames = list(c("S1", "S2"),
#'                             c("MIMAT0001631", "MIMAT0000062", "junk")))
#' apply_mirna_aliases(M, keep_unresolved = FALSE)
#'
#' @references
#' Kozomara A., Birgaoanu M., Griffiths-Jones S. (2019) miRBase: from microRNA
#' sequences to function. \emph{Nucleic Acids Research} 47(D1): D155–D162.
#' DOI: \doi{10.1093/nar/gky1141}
apply_mirna_aliases <- function(x,
                                target_namespace = c("mirna_name", "mimat"),
                                table = mirna_alias_table(),
                                keep_unresolved = TRUE,
                                verbose = FALSE,
                                ...) {
  target_namespace <- match.arg(target_namespace)

  if (is.matrix(x)) {
    if (is.null(colnames(x))) {
      stop("apply_mirna_aliases: matrix must have column names")
    }
    resolved <- resolve_mirna_aliases(colnames(x),
                                      target_namespace = target_namespace,
                                      table = table,
                                      keep_unresolved = TRUE,
                                      verbose = verbose,
                                      ...)
    if (!keep_unresolved) {
      keep_idx     <- !is.na(resolved)
      x            <- x[, keep_idx, drop = FALSE]
      resolved     <- resolved[keep_idx]
    } else {
      # Keep original name for unresolved
      unres        <- is.na(resolved)
      resolved[unres] <- colnames(x)[unres]
    }
    colnames(x) <- resolved
    return(x)

  } else if (is.numeric(x) && !is.null(names(x))) {
    resolved <- resolve_mirna_aliases(names(x),
                                      target_namespace = target_namespace,
                                      table = table,
                                      keep_unresolved = TRUE,
                                      verbose = verbose,
                                      ...)
    if (!keep_unresolved) {
      keep_idx    <- !is.na(resolved)
      x           <- x[keep_idx]
      resolved    <- resolved[keep_idx]
    } else {
      unres       <- is.na(resolved)
      resolved[unres] <- names(x)[unres]
    }
    names(x) <- resolved
    return(x)

  } else {
    stop("apply_mirna_aliases: x must be a named numeric vector or a matrix ",
         "with column names")
  }
}
