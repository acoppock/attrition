#' @importFrom generics tidy
#' @export
generics::tidy

#' Tidy an attrition bounds object
#'
#' Returns a three-row tibble. The `bounds` row is the estimator's own named
#' vector transcribed, so `estimate_lower`, `estimate_upper`,
#' `std.error_lower`, `std.error_upper`, `conf.low` and `conf.high` are the
#' same names in both places. `estimate` and `std.error` are `NA` on that row,
#' because bounds do not yield a single point estimate; the `lower_bound` and
#' `upper_bound` rows carry one endpoint each in broom's long form, which is
#' what DeclareDesign selects on with `term`.
#'
#' @param x An object of class `"attrition_bounds"` (produced by
#'   [estimator_ev()] or [estimator_ds()]).
#' @param ... Unused; included for S3 compatibility.
#'
#' @return A [tibble::tibble()] with columns `term`, `estimate`, `std.error`,
#'   `conf.low`, `conf.high`, `estimate_lower`, `estimate_upper`,
#'   `std.error_lower`, `std.error_upper`, `outcome`.
#' @export
tidy.attrition_bounds <- function(x, ...) {
  tibble::tibble(
    term            = c("bounds",  "lower_bound",                  "upper_bound"),
    estimate        = c(NA_real_,  unname(x["estimate_lower"]),    unname(x["estimate_upper"])),
    std.error       = c(NA_real_,  unname(x["std.error_lower"]),   unname(x["std.error_upper"])),
    conf.low        = c(unname(x["conf.low"]),         NA_real_,   NA_real_),
    conf.high       = c(unname(x["conf.high"]),        NA_real_,   NA_real_),
    estimate_lower  = c(unname(x["estimate_lower"]),   NA_real_,   NA_real_),
    estimate_upper  = c(unname(x["estimate_upper"]),   NA_real_,   NA_real_),
    std.error_lower = c(unname(x["std.error_lower"]),  NA_real_,   NA_real_),
    std.error_upper = c(unname(x["std.error_upper"]),  NA_real_,   NA_real_),
    outcome         = attr(x, "outcome") %||% NA_character_
  )
}

#' Tidy a trimming bounds object
#'
#' Returns a three-row tibble with the same columns as [tidy.attrition_bounds()].
#' Standard errors and the joint Imbens-Manski confidence interval are filled in
#' when [estimator_trim()] was called with `se = "analytic"` or `se = "bootstrap"`,
#' and are `NA` when it was called with `se = "none"` or when monotonicity failed.
#'
#' @param x An object of class `"attrition_trim"` (produced by
#'   [estimator_trim()]).
#' @param ... Unused; included for S3 compatibility.
#'
#' @return A [tibble::tibble()] with columns `term`, `estimate`, `std.error`,
#'   `conf.low`, `conf.high`, `estimate_lower`, `estimate_upper`,
#'   `std.error_lower`, `std.error_upper`, `outcome`.
#' @export
tidy.attrition_trim <- function(x, ...) {
  # Elements are absent, not NA, when se = "none" was never wired in
  element <- function(nm) if (nm %in% names(x)) unname(x[nm]) else NA_real_
  tibble::tibble(
    term            = c("bounds",  "lower_bound",             "upper_bound"),
    estimate        = c(NA_real_,  element("estimate_lower"), element("estimate_upper")),
    std.error       = c(NA_real_,  element("std.error_lower"), element("std.error_upper")),
    conf.low        = c(element("conf.low"),         NA_real_, NA_real_),
    conf.high       = c(element("conf.high"),        NA_real_, NA_real_),
    estimate_lower  = c(element("estimate_lower"),   NA_real_, NA_real_),
    estimate_upper  = c(element("estimate_upper"),   NA_real_, NA_real_),
    std.error_lower = c(element("std.error_lower"),  NA_real_, NA_real_),
    std.error_upper = c(element("std.error_upper"),  NA_real_, NA_real_),
    outcome         = attr(x, "outcome") %||% NA_character_
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x


#' Tidy a sensitivity analysis
#'
#' Returns the bounds and joint Imbens-Manski interval at every value of the
#' sensitivity parameter, one row per `delta`, under the same names the
#' estimators return.
#'
#' @param x An object of class `"attrition_sensitivity"` (produced by
#'   [sensitivity_ds()]).
#' @param ... Unused; included for S3 compatibility.
#'
#' @return A [tibble::tibble()] with columns `delta`, `estimate_lower`,
#'   `estimate_upper`, `std.error_lower`, `std.error_upper`, `conf.low`,
#'   `conf.high`, `outcome`.
#' @export
tidy.attrition_sensitivity <- function(x, ...) {
  tibble::tibble(
    delta = x$sims_df$delta,
    estimate_lower = x$sims_df$estimate_lower,
    estimate_upper = x$sims_df$estimate_upper,
    std.error_lower = x$sims_df$std.error_lower,
    std.error_upper = x$sims_df$std.error_upper,
    conf.low = x$sims_df$conf.low,
    conf.high = x$sims_df$conf.high,
    outcome = attr(x, "outcome") %||% NA_character_
  )
}
