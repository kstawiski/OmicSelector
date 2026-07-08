#' OmicSelector_correlation_plot
#'
#' Draw a correlation scatterplot with optional y = x reference and summary stats.
#'
#' @param var1 First numeric vector.
#' @param var2 Second numeric vector.
#' @param labvar1 Label for x-axis.
#' @param labvar2 Label for y-axis.
#' @param title Optional plot title.
#' @param yx Logical, draw y = x reference line (useful for calibration plots).
#' @param metoda Correlation method: "pearson" or "spearman".
#' @param gdzie_legenda Legend position passed to `legend()`.
#'
#' @export
OmicSelector_correlation_plot <- function(var1,
                                          var2,
                                          labvar1,
                                          labvar2,
                                          title = NULL,
                                          yx = TRUE,
                                          metoda = "pearson",
                                          gdzie_legenda = "topleft") {
  metoda <- match.arg(metoda, c("pearson", "spearman"))

  x <- as.numeric(var1)
  y <- as.numeric(var2)
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 2) {
    stop("Not enough complete observations for correlation plot.", call. = FALSE)
  }

  graphics::plot(x, y,
                 pch = 19,
                 cex = 0.6,
                 xlab = labvar1,
                 ylab = labvar2,
                 main = title)

  if (isTRUE(yx)) {
    graphics::abline(0, 1, col = "gray")
  }

  fit <- stats::lm(y ~ x)
  graphics::abline(fit, col = "black", lty = 2)

  cor_val <- stats::cor(x, y, use = "complete.obs", method = metoda)
  cor_test <- stats::cor.test(x, y, method = metoda)

  if (metoda == "pearson") {
    label <- if (cor_test$p.value < 1e-4) {
      sprintf("r = %.2f, p < 0.0001\nadj. R2 = %.4f",
              cor_val,
              summary(fit)$adj.r.squared)
    } else {
      sprintf("r = %.2f, p = %.4f\nadj. R2 = %.4f",
              cor_val,
              cor_test$p.value,
              summary(fit)$adj.r.squared)
    }
  } else {
    label <- if (cor_test$p.value < 1e-4) {
      sprintf("rho = %.2f, p < 0.0001", cor_val)
    } else {
      sprintf("rho = %.2f, p = %.4f", cor_val, cor_test$p.value)
    }
  }

  graphics::legend(gdzie_legenda, bty = "n", legend = label)

  invisible(list(correlation = cor_val, p_value = cor_test$p.value, model = fit))
}
