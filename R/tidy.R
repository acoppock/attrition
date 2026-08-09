#' @importFrom generics tidy
#' @export
generics::tidy

#' Tidy an attrition bounds object
#'
#' Returns a one-row tibble representing the identified interval. Designed for
#' use with DeclareDesign: `estimate` and `std.error` are `NA` because bounds
#' do not yield a single point estimate, while `conf.low`/`conf.high` hold the
#' Imbens-Manski joint confidence interval and `estimate.low`/`estimate.high`
#' hold the bound point estimates.
#'
#' @param x An object of class `"attrition_bounds"` (produced by
#'   [estimator_ev()] or [estimator_ds()]).
#' @param ... Unused; included for S3 compatibility.
#'
#' @return A [tibble::tibble()] with columns `term`, `estimate`, `std.error`,
#'   `conf.low`, `conf.high`, `estimate.low`, `estimate.high`.
#' @export
tidy.attrition_bounds <- function(x, ...) {
  tibble::tibble(
    term          = c("bounds",               "lower_bound",            "upper_bound"),
    estimate      = c(NA_real_,               unname(x["low_est"]),     unname(x["upp_est"])),
    std.error     = c(NA_real_,               sqrt(unname(x["low_var"])), sqrt(unname(x["upp_var"]))),
    conf.low      = c(unname(x["ci_lower"]),  NA_real_,                 NA_real_),
    conf.high     = c(unname(x["ci_upper"]),  NA_real_,                 NA_real_),
    estimate.low  = c(unname(x["low_est"]),   NA_real_,                 NA_real_),
    estimate.high = c(unname(x["upp_est"]),   NA_real_,                 NA_real_),
    outcome       = attr(x, "outcome") %||% NA_character_,
    nobs          = attr(x, "nobs") %||% NA_integer_
  )
}

#' Tidy a trimming bounds object
#'
#' Returns a three-row tibble matching the structure of [tidy.attrition_bounds()].
#' Standard errors and the joint Imbens-Manski confidence interval are filled in
#' when [estimator_trim()] was called with `se = "analytic"` or `se = "bootstrap"`,
#' and are `NA` when it was called with `se = "none"` or when monotonicity failed.
#'
#' @param x An object of class `"attrition_trim"` (produced by
#'   [estimator_trim()]).
#' @param ... Unused; included for S3 compatibility.
#'
#' @return A [tibble::tibble()] with columns `term`, `estimate`, `std.error`,
#'   `conf.low`, `conf.high`, `estimate.low`, `estimate.high`.
#' @export
tidy.attrition_trim <- function(x, ...) {
  # Elements are absent, not NA, when se = "none" was never wired in
  element <- function(nm) if (nm %in% names(x)) unname(x[nm]) else NA_real_
  tibble::tibble(
    term          = c("bounds",                 "lower_bound",            "upper_bound"),
    estimate      = c(NA_real_,                 element("lower_bound"),   element("upper_bound")),
    std.error     = c(NA_real_,                 element("lower_se"),      element("upper_se")),
    conf.low      = c(element("ci_lower"),      NA_real_,                 NA_real_),
    conf.high     = c(element("ci_upper"),      NA_real_,                 NA_real_),
    estimate.low  = c(element("lower_bound"),   NA_real_,                 NA_real_),
    estimate.high = c(element("upper_bound"),   NA_real_,                 NA_real_),
    outcome       = attr(x, "outcome") %||% NA_character_,
    nobs          = attr(x, "nobs") %||% NA_integer_
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Number of observations
#'
#' The number of units the estimator was given, counting those whose outcome is
#' missing. Every variance in the package is an asymptotic result indexed to the
#' full randomized sample rather than to the respondents: the worst-case bounds
#' divide by the number assigned to each arm, and Lee (2009) Proposition 3 is a
#' root-n result in which the response and trimming rates appear as proportions
#' of that same n. Counts of respondents, of retained observations after
#' trimming, and of the groups behind them are returned as named elements of the
#' estimator output.
#'
#' @param object An object of class `"attrition_bounds"` or `"attrition_trim"`.
#' @param ... Unused; included for S3 compatibility.
#'
#' @return An integer.
#' @export
nobs.attrition_bounds <- function(object, ...) attr(object, "nobs")

#' @rdname nobs.attrition_bounds
#' @export
nobs.attrition_trim <- function(object, ...) attr(object, "nobs")
