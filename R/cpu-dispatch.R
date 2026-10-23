# Set before any OmicSelector backend imports torch or NumPy. Explicit caller
# settings, including empty strings, take precedence. AVX2 names are x86-only.
.ss_cpu_has_avx512 <- function() {
  if (!file.exists("/proc/cpuinfo")) return(FALSE)
  info <- tryCatch(readLines("/proc/cpuinfo", n = 100L, warn = FALSE),
                   error = function(e) character())
  any(grepl("\\bavx512f\\b", info, ignore.case = TRUE))
}

.ss_configure_cpu_dispatch <- function(avx512 = .ss_cpu_has_avx512()) {
  if (!tolower(Sys.info()[["machine"]]) %in% c("x86_64", "amd64")) {
    return(invisible(FALSE))
  }
  defaults <- c(ATEN_CPU_CAPABILITY = "avx2", MKL_ENABLE_INSTRUCTIONS = "AVX2",
    ONEDNN_MAX_CPU_ISA = "AVX2", DNNL_MAX_CPU_ISA = "AVX2")
  # Native AVX2 OpenBLAS kernels must remain unchanged on AVX2 hosts. On
  # Linux AVX-512 hosts, Zen selects the validated AVX2 implementation.
  if (isTRUE(avx512)) defaults <- c(defaults, OPENBLAS_CORETYPE = "Zen")
  missing <- is.na(Sys.getenv(names(defaults), unset = NA_character_))
  if (any(missing)) do.call(Sys.setenv, as.list(defaults[missing]))
  invisible(TRUE)
}

.ss_cpu_dispatch_provenance <- function(capability = NULL) {
  if (is.null(capability)) {
    torch <- reticulate::import("torch", convert = FALSE)
    capability <- reticulate::py_to_r(torch$backends$cpu$get_cpu_capability())
  }
  list(capability = capability, settings = Sys.getenv(c("ATEN_CPU_CAPABILITY",
    "MKL_ENABLE_INSTRUCTIONS", "ONEDNN_MAX_CPU_ISA", "DNNL_MAX_CPU_ISA",
    "OPENBLAS_CORETYPE"), unset = NA_character_))
}

.ss_check_cpu_dispatch <- function(model, active = NULL,
    strict = getOption("OmicSelector.cpu_dispatch_strict", FALSE)) {
  .ss_check_strict(strict, "foundation CPU dispatch")
  fitted <- model$provenance$cpu_dispatch
  # Legacy models have no recorded capability. Strict mode requires evidence.
  if (is.null(fitted) && !strict) return(invisible(TRUE))
  if (is.null(active)) active <- .ss_cpu_dispatch_provenance()
  if (identical(fitted, active)) return(invisible(TRUE))
  describe <- function(x) {
    if (is.null(x)) return("unrecorded")
    paste(c(paste0("capability=", x$capability),
      paste0(names(x$settings), "=", x$settings)), collapse = ", ")
  }
  msg <- paste0("Foundation-model CPU dispatch differs from fitted provenance. ",
    "Fitted: ", describe(fitted), "; active: ", describe(active), ". ",
    "Scores may differ. Restart R with matching dispatch settings before ",
    "Python imports torch; changing environment variables cannot reset an ",
    "already initialized backend.")
  if (strict) stop(msg, call. = FALSE)
  warning(msg, call. = FALSE)
  invisible(FALSE)
}
