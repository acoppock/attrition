#' Print bounds
#'
#' Prints the named vector of results on its own. `print.default()` would
#' otherwise append the class vector to every result, which is noise.
#'
#' @param x An object of class `"attrition_bounds"`, produced by
#'   [estimator_ev()], [estimator_ds()], or [estimator_ds_sens()].
#' @param ... Passed to `print.default()`.
#'
#' @return `x`, invisibly.
#' @export
print.attrition_bounds <- function(x, ...) {
  print(setNames(as.numeric(x), names(x)), ...)
  invisible(x)
}

#' Print trimming bounds
#'
#' @param x An object of class `"attrition_trim"`, produced by [estimator_trim()].
#' @param ... Passed to `print.default()`.
#'
#' @return `x`, invisibly.
#' @export
print.attrition_trim <- function(x, ...) {
  print(setNames(as.numeric(x), names(x)), ...)
  invisible(x)
}
