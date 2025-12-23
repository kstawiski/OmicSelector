#' OmicSelector_profileplot
#'
#' Profile plot for comparing score patterns across methods or samples.
#'
#' @param form Numeric matrix/data.frame of scores (rows = methods, cols = subscores).
#' @param Method.id Optional vector of method identifiers (length = nrow(form)).
#' @param standardize Logical; if TRUE, z-score each column.
#' @param interval Number of intervals on y-axis for the pattern view.
#' @param by.pattern Logical; if TRUE, draw pattern view with ggplot2.
#' @param original.names Logical; if TRUE, use column names as labels.
#'
#' @return A ggplot object (by.pattern = TRUE) or NULL (base plot).
#'
#' @export
OmicSelector_profileplot <- function(form,
                                     Method.id,
                                     standardize = TRUE,
                                     interval = 10,
                                     by.pattern = TRUE,
                                     original.names = TRUE) {
  form <- as.matrix(form)
  if (!is.numeric(form)) {
    stop("'form' must be numeric.", call. = FALSE)
  }

  if (standardize) {
    form <- scale(form)
  }

  n <- ncol(form)
  k <- nrow(form)
  labels <- if (isTRUE(original.names) && !is.null(colnames(form))) {
    colnames(form)
  } else {
    paste0("v", seq_len(n))
  }

  id <- if (missing(Method.id)) {
    seq_len(k)
  } else {
    if (length(Method.id) != k) {
      stop("Method.id must have length equal to nrow(form).", call. = FALSE)
    }
    Method.id
  }

  if (isTRUE(by.pattern)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for by.pattern plots.", call. = FALSE)
    }

    level <- rowMeans(form, na.rm = TRUE)
    pattern <- form - level

    subscore <- as.vector(t(form))
    pattern_vec <- as.vector(t(pattern))

    df <- data.frame(
      id = factor(rep(id, each = n)),
      subscore.id = rep(labels, times = k),
      subscore = subscore,
      pattern = pattern_vec,
      level = rep(level, each = n),
      stringsAsFactors = FALSE
    )

    max_s <- max(df$subscore, na.rm = TRUE)
    min_s <- min(df$subscore, na.rm = TRUE)
    int <- if (interval <= 0) 1 else (max_s - min_s) / interval

    p <- ggplot2::ggplot(df,
                         ggplot2::aes(x = subscore.id, y = subscore,
                                      group = id, color = id)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2, fill = "white") +
      ggplot2::scale_colour_hue(name = "Method", l = 30) +
      ggplot2::theme_classic() +
      ggplot2::scale_x_discrete(name = " ", labels = labels) +
      ggplot2::scale_y_continuous(
        name = "Scores",
        limits = c(min_s, max_s),
        breaks = seq(min_s, max_s, int)
      )

    return(p)
  }

  colours <- grDevices::hcl.colors(n, "Set2")
  mymin <- min(form, na.rm = TRUE)
  mymax <- max(form, na.rm = TRUE)

  for (i in seq_len(n)) {
    Scoresi <- form[, i]
    namei <- labels[i]
    colouri <- colours[i]
    if (i == 1) {
      graphics::plot(Scoresi, col = colouri, type = "l",
                     ylim = c(mymin, mymax), ylab = "Score", xlab = "Index")
    } else {
      graphics::points(Scoresi, col = colouri, type = "l")
    }
    lastxval <- length(Scoresi)
    lastyval <- Scoresi[length(Scoresi)]
    graphics::text(lastxval, lastyval, namei, col = "black", cex = 0.7)
  }

  invisible(NULL)
}
