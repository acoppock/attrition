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

#' Print a sensitivity analysis
#'
#' Reports delta*, the smallest value of the sensitivity parameter at which the
#' confidence interval reaches zero, and names the components of the object.
#' Printing the list itself would draw the plot, which is not what a glance at
#' the result should do.
#'
#' @param x An object of class `"attrition_sensitivity"`, produced by
#'   [sensitivity_ds()].
#' @param ... Unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#' @export
print.attrition_sensitivity <- function(x, ...) {
  level <- format_pct(1 - (attr(x, "alpha") %||% 0.05), digits = 0)
  cat("Sensitivity analysis on ", attr(x, "outcome") %||% "the outcome", "\n", sep = "")
  if (is.na(x$delta_star)) {
    writeLines(strwrap(paste0(
      "delta*: none. The ", level, " confidence interval contains zero at every ",
      "delta, including delta = 0, where ignorability is assumed for all follow-up ",
      "nonrespondents."), width = 76, exdent = 2))
  } else {
    writeLines(strwrap(paste0(
      "delta* = ", format(round(x$delta_star, 3), nsmall = 3), ": the ", level,
      " confidence interval first includes zero when ignorability is dropped for ",
      format_pct(x$delta_star), " of the follow-up nonrespondents."),
      width = 76, exdent = 2))
  }
  cat("Components: sensitivity_plot, sims_df (", nrow(x$sims_df), " values of delta), delta_star\n", sep = "")
  invisible(x)
}
