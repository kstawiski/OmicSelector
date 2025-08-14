#' Create Diverging Color Palette for Heatmaps
#'
#' @description
#' Modern helper function to create diverging color palettes centered on a specific value.
#' Useful for heatmap visualization where you want to emphasize deviations from a baseline.
#'
#' @param data Numeric vector or matrix of values
#' @param centered_on Numeric value to center the color scale on (default: 0)
#' @param n_colors Integer number of colors to generate (default: 100)
#' @param colors Character vector of colors for the palette (default: blue-white-red)
#' @param symmetric Logical indicating whether to make the scale symmetric around center (default: TRUE)
#'
#' @return List containing:
#'   \item{breaks}{Numeric vector of break points for the color scale}
#'   \item{colors}{Character vector of colors corresponding to breaks}
#'   \item{cuts}{Cut intervals for mapping data to colors}
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' data <- rnorm(100)
#' palette <- create_diverging_palette(data)
#' 
#' # Custom center and colors
#' palette <- create_diverging_palette(data, centered_on = 1, 
#'                                   colors = c("green", "white", "purple"))
#' }
#'
#' @export
create_diverging_palette <- function(data, 
                                   centered_on = 0, 
                                   n_colors = 100,
                                   colors = c("blue", "white", "red"),
                                   symmetric = TRUE) {
  
  # Input validation
  if (!is.numeric(data)) {
    stop("'data' must be numeric", call. = FALSE)
  }
  
  if (length(colors) < 2) {
    stop("'colors' must contain at least 2 colors", call. = FALSE)
  }
  
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Package 'grDevices' is required", call. = FALSE)
  }
  
  # Remove missing values for range calculation
  data_clean <- data[!is.na(data)]
  
  if (length(data_clean) == 0) {
    stop("No non-missing data values found", call. = FALSE)
  }
  
  # Calculate data range
  data_min <- min(data_clean)
  data_max <- max(data_clean)
  
  # Ensure center is within data range
  if (centered_on < data_min || centered_on > data_max) {
    warning("centered_on is outside data range; may produce unexpected results", 
            call. = FALSE)
  }
  
  # Create color palette
  color_palette <- grDevices::colorRampPalette(colors)(n_colors)
  
  if (symmetric) {
    # Create symmetric breaks around center
    max_deviation <- max(abs(data_max - centered_on), abs(data_min - centered_on))
    breaks <- seq(centered_on - max_deviation, 
                  centered_on + max_deviation, 
                  length.out = n_colors + 1)
  } else {
    # Create asymmetric breaks
    n_half <- floor(n_colors / 2)
    breaks_low <- seq(data_min, centered_on, length.out = n_half + 1)
    breaks_high <- seq(centered_on, data_max, length.out = n_half + 1)[-1]
    breaks <- c(breaks_low, breaks_high)
  }
  
  # Create cut intervals
  cuts <- cut(data_clean, breaks = breaks, include.lowest = TRUE)
  
  return(list(
    breaks = breaks,
    colors = color_palette,
    cuts = cuts,
    data_range = c(data_min, data_max),
    center = centered_on
  ))
}

#' @rdname create_diverging_palette
#' @export
OmicSelector_diverge_color <- function(data, centeredOn = 0) {
  # Backward compatibility wrapper
  .Deprecated("create_diverging_palette", 
              msg = "OmicSelector_diverge_color is deprecated. Use create_diverging_palette instead.")
  
  if (!requireNamespace("classInt", quietly = TRUE)) {
    stop("Package 'classInt' is required for legacy function", call. = FALSE)
  }
  
  # Legacy implementation for exact compatibility
  n_half <- 50
  data_min <- min(data, na.rm = TRUE)
  data_max <- max(data, na.rm = TRUE)
  
  # Create legacy color palette
  pal <- grDevices::colorRampPalette(c("blue", "white", "red"))(n = 11)
  rc1 <- grDevices::colorRampPalette(colors = c(pal[1], pal[2]), space = "Lab")(10)
  
  for (i in 2:10) {
    tmp <- grDevices::colorRampPalette(colors = c(pal[i], pal[i + 1]), space = "Lab")(10)
    rc1 <- c(rc1, tmp)
  }
  
  # Create breaks
  rb1 <- seq(data_min, centeredOn, length.out = n_half + 1)
  rb2 <- seq(centeredOn, data_max, length.out = n_half + 1)[-1]
  rampbreaks <- c(rb1, rb2)
  
  cuts <- classInt::classIntervals(data, style = "fixed", fixedBreaks = rampbreaks)
  
  return(list(cuts, rc1))
}
