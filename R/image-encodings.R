#' @title Image Encoding Methods for Biomarker Panels
#'
#' @description
#' Converts tabular biomarker expression data into 2D image representations
#' suitable for CNN classification. Four encoding methods are supported:
#'
#' \enumerate{
#'   \item \code{encode_simple_grid}: Fixed row-major layout in rows×cols grid
#'   \item \code{encode_corr_grid}: Features ordered by hierarchical correlation clustering
#'   \item \code{encode_deepinsight}: Feature positions learned via t-SNE on training data
#'     (Sharma et al., 2019, Scientific Reports)
#'   \item \code{encode_ratio_image}: Pairwise log-ratio matrix (see ratio-image-cnn.R)
#' }
#'
#' These complement the within-sample normalization approaches in
#' \code{within-sample.R} and can be combined with \code{train_ratio_cnn}
#' and related CNN training utilities.
#'
#' @references
#' Sharma A et al. (2019). DeepInsight: A methodology to transform a
#' non-image data to an image for convolution neural network architecture.
#' Scientific Reports 9:11399.
#'
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Quinn TP, Erb I, Richardson MF, Crowley TM. (2020).
#' Understanding sequencing data as compositions: an outlook and review.
#' \emph{Bioinformatics}, 36(16), 4424-4432.
#'
#' @name image-encodings
NULL


#' @title Simple Row-Major Grid Encoding
#'
#' @description
#' Places feature values in a rectangular grid in row-major order. Simple
#' reshape with zero-padding if (rows × cols) > p.
#'
#' @param x Numeric vector of feature values (length p)
#' @param rows Number of rows in the grid (default: ceiling(sqrt(p)))
#' @param cols Number of columns (default: ceiling(p / rows))
#'
#' @return A rows × cols numeric matrix
#'
#' @examples
#' encode_simple_grid(1:6, rows = 2, cols = 3)
#'
#' @references
#' Sharma A et al. (2019). DeepInsight: A methodology to transform a
#' non-image data to an image for convolution neural network architecture.
#' \emph{Scientific Reports}, 9, 11399.
#'
#' @export
encode_simple_grid <- function(x, rows = NULL, cols = NULL) {
  p <- length(x)
  if (is.null(rows)) rows <- ceiling(sqrt(p))
  if (is.null(cols)) cols <- ceiling(p / rows)

  img <- matrix(0, nrow = rows, ncol = cols)
  for (i in seq_len(p)) {
    if (i > rows * cols) break
    r <- ((i - 1) %/% cols) + 1
    c <- ((i - 1) %% cols) + 1
    img[r, c] <- x[i]
  }
  img
}


#' @title Correlation-Ordered Grid Encoding
#'
#' @description
#' Reorders features by hierarchical clustering on feature-feature correlation
#' distance, then fills a grid in row-major order. Features that co-vary end up
#' as neighbors in the image, helping the CNN's local receptive fields exploit
#' correlation structure.
#'
#' Must be FIT on training data only to avoid leakage.
#'
#' @param X_train Numeric matrix (n_train × p features) for fitting the order
#' @param method Linkage method for \code{hclust} (default: "average")
#'
#' @return Integer vector of length p, the feature order to use for encoding
#'
#' @examples
#' X_train <- matrix(rnorm(30), nrow = 6, ncol = 5)
#' fit_corr_order(X_train)
#'
#' @references
#' Sharma A et al. (2019). DeepInsight: A methodology to transform a
#' non-image data to an image for convolution neural network architecture.
#' \emph{Scientific Reports}, 9, 11399.
#'
#' @export
fit_corr_order <- function(X_train, method = "average") {
  corr <- stats::cor(X_train, use = "pairwise.complete.obs")
  corr[is.na(corr)] <- 0
  d <- stats::as.dist(1 - abs(corr))
  tree <- stats::hclust(d, method = method)
  tree$order
}


#' @title Encode with Correlation-Ordered Grid
#'
#' @param x Numeric feature vector (length p)
#' @param order Integer vector from \code{fit_corr_order}
#' @param rows Number of rows in the output grid. Defaults to
#'   \code{ceiling(sqrt(length(x)))}.
#' @param cols Number of columns in the output grid. Defaults to
#'   \code{ceiling(length(x) / rows)}.
#'
#' @return A rows × cols numeric matrix
#'
#' @examples
#' encode_corr_grid(1:6, order = c(6, 5, 4, 3, 2, 1), rows = 2, cols = 3)
#'
#' @references
#' Sharma A et al. (2019). DeepInsight: A methodology to transform a
#' non-image data to an image for convolution neural network architecture.
#' \emph{Scientific Reports}, 9, 11399.
#'
#' @export
encode_corr_grid <- function(x, order, rows = NULL, cols = NULL) {
  x_ordered <- x[order]
  encode_simple_grid(x_ordered, rows = rows, cols = cols)
}


#' @title Fit DeepInsight Feature Layout
#'
#' @description
#' Uses t-SNE on feature profiles (each feature viewed across samples) to
#' project features to 2D positions. Discretizes to a grid. Resolves grid
#' collisions by nearest-free-cell search.
#'
#' Must be FIT on training data only to avoid leakage.
#'
#' Requires \code{Rtsne} package.
#'
#' @param X_train Numeric matrix (n_train × p features)
#' @param grid_size Integer, size of the square output grid (default: 14)
#' @param perplexity t-SNE perplexity; if NULL, auto-chosen as min(5, p-1)
#' @param seed Random seed for reproducibility
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{positions}: p × 2 integer matrix of (row, col) per feature
#'     \item \code{grid_size}: the grid dimension
#'   }
#'
#' @references Sharma A et al. (2019). Scientific Reports 9:11399.
#'
#' @examples
#' \dontrun{
#' X_train <- matrix(rnorm(70), nrow = 10, ncol = 7)
#' layout <- fit_deepinsight(X_train, grid_size = 7, seed = 1)
#' }
#'
#' @export
fit_deepinsight <- function(X_train, grid_size = 14L, perplexity = NULL, seed = 42L) {
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("Package 'Rtsne' is required for DeepInsight. Install with: install.packages('Rtsne')")
  }

  p <- ncol(X_train)
  if (is.null(perplexity)) {
    perplexity <- max(1, min(5, floor((p - 1L) / 3L)))
  }

  # Each feature's profile across training samples = one "point" in t-SNE
  feat_profiles <- t(X_train)  # p features × n_train samples

  set.seed(seed)
  tsne_fit <- Rtsne::Rtsne(feat_profiles, dims = 2, perplexity = perplexity,
                             pca = TRUE, check_duplicates = FALSE)
  coords <- tsne_fit$Y

  # Normalize to [0, grid_size - 1] and discretize
  coords_n <- sweep(coords, 2, apply(coords, 2, min), "-")
  coords_n <- sweep(coords_n, 2, pmax(apply(coords_n, 2, max), 1e-10), "/")
  coords_g <- pmin(floor(coords_n * grid_size), grid_size - 1L)
  mode(coords_g) <- "integer"

  # Resolve collisions
  occupied <- matrix(FALSE, nrow = grid_size, ncol = grid_size)
  positions <- matrix(0L, nrow = p, ncol = 2)
  for (i in seq_len(p)) {
    r <- coords_g[i, 1] + 1L
    c <- coords_g[i, 2] + 1L
    if (occupied[r, c]) {
      # Find nearest empty cell
      done <- FALSE
      for (d in 1:grid_size) {
        for (dr in -d:d) {
          for (dc in -d:d) {
            nr <- r + dr; nc <- c + dc
            if (nr >= 1 && nr <= grid_size && nc >= 1 && nc <= grid_size &&
                !occupied[nr, nc]) {
              positions[i, ] <- c(nr, nc)
              occupied[nr, nc] <- TRUE
              done <- TRUE; break
            }
          }
          if (done) break
        }
        if (done) break
      }
    } else {
      positions[i, ] <- c(r, c)
      occupied[r, c] <- TRUE
    }
  }

  list(positions = positions, grid_size = grid_size)
}


#' @title Encode with DeepInsight Layout
#'
#' @param x Numeric feature vector (length p)
#' @param layout Output of \code{fit_deepinsight}
#'
#' @return A grid_size × grid_size numeric matrix
#'
#' @examples
#' layout <- list(positions = matrix(c(1, 1, 1, 2, 2, 1), ncol = 2, byrow = TRUE), grid_size = 2L)
#' encode_deepinsight(1:3, layout)
#'
#' @references
#' Sharma A et al. (2019). DeepInsight: A methodology to transform a
#' non-image data to an image for convolution neural network architecture.
#' \emph{Scientific Reports}, 9, 11399.
#'
#' @export
encode_deepinsight <- function(x, layout) {
  positions <- layout$positions
  g <- layout$grid_size
  img <- matrix(0, nrow = g, ncol = g)
  for (i in seq_along(x)) {
    r <- positions[i, 1]; c <- positions[i, 2]
    img[r, c] <- x[i]
  }
  img
}


#' @title Batch-Encode an Expression Matrix with Any Encoding
#'
#' @description
#' Convenience wrapper that applies an encoding method to every row of an
#' expression matrix, returning a 3D array suitable as CNN input.
#'
#' @param X Numeric matrix (n × p)
#' @param method One of "simple", "corr", "deepinsight", "ratio"
#' @param fit_params For "corr" or "deepinsight": the fitted order/layout.
#'   For "simple" or "ratio": ignored.
#' @param ... Additional arguments passed to the encoding function
#'
#' @return A 3D array (n × height × width)
#'
#' @examples
#' X <- matrix(rnorm(28), nrow = 4, ncol = 7)
#' encode_batch(X, method = "simple", rows = 1, cols = 7)
#'
#' @references
#' Sharma A et al. (2019). DeepInsight: A methodology to transform a
#' non-image data to an image for convolution neural network architecture.
#' \emph{Scientific Reports}, 9, 11399.
#'
#' Aitchison J. (1986). \emph{The Statistical Analysis of Compositional Data}.
#' Chapman and Hall.
#'
#' Quinn TP, Erb I, Richardson MF, Crowley TM. (2020).
#' Understanding sequencing data as compositions: an outlook and review.
#' \emph{Bioinformatics}, 36(16), 4424-4432.
#'
#' @export
encode_batch <- function(X, method = c("simple", "corr", "deepinsight", "ratio"),
                         fit_params = NULL, ...) {
  method <- match.arg(method)
  n <- nrow(X)
  args <- list(...)

  # Determine output dimensions from first sample
  first_img <- switch(method,
    simple = encode_simple_grid(X[1, ], rows = args$rows, cols = args$cols),
    corr = encode_corr_grid(X[1, ], order = fit_params, rows = args$rows, cols = args$cols),
    deepinsight = encode_deepinsight(X[1, ], layout = fit_params),
    ratio = encode_ratio_image(X[1, ])
  )
  h <- nrow(first_img); w <- ncol(first_img)
  imgs <- array(0, dim = c(n, h, w))
  imgs[1, , ] <- first_img

  for (i in 2:n) {
    imgs[i, , ] <- switch(method,
      simple = encode_simple_grid(X[i, ], rows = args$rows, cols = args$cols),
      corr = encode_corr_grid(X[i, ], order = fit_params, rows = args$rows, cols = args$cols),
      deepinsight = encode_deepinsight(X[i, ], layout = fit_params),
      ratio = encode_ratio_image(X[i, ])
    )
  }
  imgs
}
