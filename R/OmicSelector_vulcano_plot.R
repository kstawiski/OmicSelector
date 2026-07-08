#' OmicSelector_vulcano_plot
#'
#' Draw a volcano plot of selected features from a differential expression table.
#'
#' @param selected_miRNAs Character vector of selected feature names.
#' @param DE Differential expression data.frame containing columns for log2FC and p-values.
#' @param only_label Optional character vector: label only these features (must be in DE$miR).
#' @param take_adjusted_p Use adjusted p-values (BH) if TRUE.
#' @param log2fc_col Column name for log2 fold-change (default: "log2FC").
#' @param p_col Column name for raw p-values (default: "p-value").
#' @param p_adj_col Column name for adjusted p-values (default: "p-value BH").
#' @param title Optional plot title.
#'
#' @return A ggplot object.
#'
#' @export
OmicSelector_vulcano_plot <- function(selected_miRNAs,
                                      DE,
                                      only_label = NULL,
                                      take_adjusted_p = FALSE,
                                      log2fc_col = "log2FC",
                                      p_col = "p-value",
                                      p_adj_col = "p-value BH",
                                      title = NULL) {
  if (missing(DE) || is.null(DE)) {
    stop("DE must be provided (data.frame with log2FC and p-value columns).", call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for OmicSelector_vulcano_plot().", call. = FALSE)
  }

  if (!log2fc_col %in% names(DE)) {
    stop(sprintf("Column '%s' not found in DE.", log2fc_col), call. = FALSE)
  }

  pcol <- if (isTRUE(take_adjusted_p)) p_adj_col else p_col
  if (!pcol %in% names(DE)) {
    stop(sprintf("Column '%s' not found in DE.", pcol), call. = FALSE)
  }

  if (!"miR" %in% names(DE)) {
    stop("DE must contain a 'miR' column with feature names.", call. = FALSE)
  }

  idx <- match(selected_miRNAs, DE$miR)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) {
    stop("No selected features found in DE$miR.", call. = FALSE)
  }

  temp <- DE[idx, , drop = FALSE]
  log2fc <- temp[[log2fc_col]]
  pvals <- temp[[pcol]]

  label_vec <- gsub("\\.", "-", temp$miR)
  if (!is.null(only_label)) {
    label_vec[!temp$miR %in% only_label] <- NA_character_
  }

  df <- data.frame(
    log2fc = log2fc,
    pval = pvals,
    label = label_vec,
    stringsAsFactors = FALSE
  )

  df$pval <- pmax(df$pval, .Machine$double.xmin)
  df$neglog10p <- -log10(df$pval)

  thres <- -log10(0.05)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2fc, y = neglog10p, label = label)) +
    ggplot2::geom_point(color = "black") +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::labs(
      x = "Log2(FC)",
      y = if (take_adjusted_p) "-Log10(adjusted p-value)" else "-Log10(p-value)",
      title = title
    ) +
    ggplot2::geom_hline(yintercept = thres, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "blue", size = 0.5)

  if (any(!is.na(df$label))) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        arrow = ggplot2::arrow(length = grid::unit(0.01, "npc"),
                               type = "closed",
                               ends = "last"),
        force = 50
      )
    } else {
      p <- p + ggplot2::geom_text(check_overlap = TRUE, size = 3)
    }
  }

  p
}
