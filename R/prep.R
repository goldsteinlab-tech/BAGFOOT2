
ds239487_ptm <- proc.time();  # initial time
#' Prep
#'
#' @param output Message to print.
#' @param ... something
#'
#' @return something
#' @export
#'
dlog <- function(output, ...) {
  etime = proc.time()-ds239487_ptm;
  cat(sprintf("[%g sec] %s", etime[["elapsed"]], output, ...));
  ds239487_ptm <- proc.time();
}
