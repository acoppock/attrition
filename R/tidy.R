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
    estimate.high = c(unname(x["upp_est"]),   NA_real_,                 NA_real_)
  )
}

#' Tidy a trimming bounds object
#'
#' Returns a three-row tibble matching the structure of [tidy.attrition_bounds()].
#' Standard errors and confidence intervals are not available for trimming
#' bounds, so those columns are `NA` throughout.
#'
#' @param x An object of class `"attrition_trim"` (produced by
#'   [estimator_trim()]).
#' @param ... Unused; included for S3 compatibility.
#'
#' @return A [tibble::tibble()] with columns `term`, `estimate`, `std.error`,
#'   `conf.low`, `conf.high`, `estimate.low`, `estimate.high`.
#' @export
tidy.attrition_trim <- function(x, ...) {
  tibble::tibble(
    term          = c("bounds",                      "lower_bound",              "upper_bound"),
    estimate      = c(NA_real_,                      unname(x["lower_bound"]),   unname(x["upper_bound"])),
    std.error     = c(NA_real_,                      NA_real_,                   NA_real_),
    conf.low      = c(NA_real_,                      NA_real_,                   NA_real_),
    conf.high     = c(NA_real_,                      NA_real_,                   NA_real_),
    estimate.low  = c(unname(x["lower_bound"]),      NA_real_,                   NA_real_),
    estimate.high = c(unname(x["upper_bound"]),      NA_real_,                   NA_real_)
  )
}
