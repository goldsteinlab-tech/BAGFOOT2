#' Deprecated aliases kept from the "tiger" naming
#'
#' The package was renamed from \code{tiger} to \code{bagfoot2}. These wrappers
#' keep older scripts working; they emit a deprecation warning and forward all
#' arguments to the new function. They will be removed in a future release.
#'
#' @param ... Passed on to the corresponding \code{*_bagfoot2()} function.
#' @name tiger-deprecated
NULL

#' @rdname tiger-deprecated
#' @export
gen_bagplot_tiger <- function(...) {
  .Deprecated("gen_bagplot_bagfoot2")
  gen_bagplot_bagfoot2(...)
}

#' @rdname tiger-deprecated
#' @export
gen_evt_barplot_tiger <- function(...) {
  .Deprecated("gen_evt_barplot_bagfoot2")
  gen_evt_barplot_bagfoot2(...)
}

#' @rdname tiger-deprecated
#' @export
download_tiger_files <- function(...) {
  .Deprecated("download_bagfoot2_files")
  download_bagfoot2_files(...)
}
