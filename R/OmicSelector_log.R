#' OmicSelector_log
#'
#' Logging to sink.
#'
#' @param message Message to log.
#' @param logfile Path to logfile.
#'
#'
#' @export
OmicSelector_log = function(message_to_log = "Working...", logfile = "temp.log") {
  hostname = as.character(R.utils::GString$getBuiltinHostname())
  timestamp = as.character(Sys.time())
  pid = as.character(R.utils::getBuiltinPid.GString())
  sink(logfile, append=TRUE)
  to_write = paste0("\n[",timestamp," | pid:", pid,"] ", message_to_log)
cat(to_write)
sink()
  message(to_write)
}
