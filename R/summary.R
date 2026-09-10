#' Summarize bounds
#'
#' Reports the identification region, its standard errors and the joint
#' Imbens-Manski interval with a line naming the estimand and the assumptions
#' that produced them. Where [print()] shows the returned vector, `summary()`
#' says what the numbers are of.
#'
#' @param object An object of class `"attrition_bounds"`, produced by
#'   [estimator_ev()], [estimator_ds()], or [estimator_ds_sens()].
#' @param ... Unused; included for S3 compatibility.
#'
#' @return The [tidy()] data frame, invisibly. The report is printed.
#' @export
#'
#' @examples
#' summary(estimator_ev(Y = Y_polarization_w2, Z = Z, R = R1,
#'                      minY = 0, maxY = 6, data = levendusky_replication))
summary.attrition_bounds <- function(object, ...) {
  estimator <- switch(class(object)[1],
                      attrition_ev = "Extreme value (Manski) bounds",
                      attrition_ds = "Extreme value bounds with double sampling",
                      attrition_ds_sens = "Extreme value bounds with double sampling and sensitivity",
                      "Bounds")
  delta <- attr(object, "delta")

  assumptions <- c(
    "the outcome lies within the stated minimum and maximum",
    if (!identical(class(object)[1], "attrition_ev"))
      "the follow-up sample was drawn at random from the nonrespondents",
    if (!is.null(delta))
      paste0("ignorability for ", format_pct(1 - delta),
             " of the follow-up nonrespondents (delta = ", delta, ")")
  )

  cat(estimator, " on ", attr(object, "outcome") %||% "the outcome", "\n", sep = "")
  cat("Estimand: the average treatment effect among all subjects\n\n")
  report_bounds(object)
  cat("\n")
  writeLines(strwrap(paste0("Assuming ", paste(assumptions, collapse = "; "), "."),
                     width = 76, exdent = 2))
  if (isTRUE(attr(object, "strata"))) {
    cat("Poststratified, which sharpens the estimate without changing the estimand.\n")
  }
  invisible(tidy(object))
}

#' Summarize trimming bounds
#'
#' Reports the identification region, its standard errors and the joint
#' Imbens-Manski interval, together with the design and the selection assumption
#' that produced them. The two are separate choices in [estimator_trim()] and
#' both change what the numbers mean, so both are named here.
#'
#' @param object An object of class `"attrition_trim"`, produced by
#'   [estimator_trim()].
#' @param ... Unused; included for S3 compatibility.
#'
#' @return The [tidy()] data frame, invisibly. The report is printed.
#' @export
#'
#' @examples
#' summary(estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
#'                        monotonicity = "treatment_decreases_response",
#'                        data = levendusky_replication))
summary.attrition_trim <- function(object, ...) {
  monotonicity <- attr(object, "monotonicity") %||% "treatment_increases_response"
  single_stage <- attr(object, "single_stage") %||% TRUE

  assumption <- switch(
    monotonicity,
    treatment_increases_response =
      "treatment never lowered the chance of responding, so the control respondents are the always-reporters and the treatment group is trimmed",
    treatment_decreases_response =
      "treatment never raised the chance of responding, so the treated respondents are the always-reporters and the control group is trimmed",
    none =
      "random assignment alone, so both groups are trimmed by the largest share of each that could fail to be an always-reporter")

  cat("Trimming bounds on ", attr(object, "outcome") %||% "the outcome", "\n", sep = "")
  cat("Estimand: the average treatment effect among always-reporters\n\n")

  if (is.na(object["estimate_lower"])) {
    cat("  No bounds: monotonicity is violated in the direction assumed.\n\n")
  } else {
    report_bounds(object)
    cat("\n")
  }

  cat("Design:     ", if (single_stage) "single sample" else "double sampling", "\n", sep = "")
  writeLines(strwrap(paste0("Assumption: ", assumption), width = 76, exdent = 12))
  if ("Q" %in% names(object)) {
    trimmed_group <- if (monotonicity == "treatment_decreases_response") "control" else "treatment"
    cat("Trimmed:    ", format_pct(unname(object["Q"])), " of the ", trimmed_group,
        " group\n", sep = "")
  }
  if (all(c("trim0", "trim1") %in% names(object))) {
    cat("Trimmed:    ", format_pct(unname(object["trim1"])), " of the treatment group and ",
        format_pct(unname(object["trim0"])), " of the control group\n", sep = "")
  }
  se_method <- attr(object, "se_method") %||% "none"
  cat("Std errors: ",
      switch(se_method,
             analytic = "analytic (Lee 2009, Proposition 3)",
             bootstrap = "bootstrap, resampled within treatment arm",
             none = "not computed"),
      "\n", sep = "")
  invisible(tidy(object))
}

# The three rows every summary shares, aligned and rounded for reading rather
# than for further computation.
report_bounds <- function(object) {
  level <- format_pct(1 - (attr(object, "alpha") %||% 0.05), digits = 0)
  fmt <- function(a, b) paste0("[", format(round(a, 3), nsmall = 3), ", ",
                               format(round(b, 3), nsmall = 3), "]")
  line <- function(label, value) cat("  ", formatC(label, width = -22), value, "\n", sep = "")

  line("Identification region",
       fmt(unname(object["estimate_lower"]), unname(object["estimate_upper"])))
  if (!is.na(object["std.error_lower"])) {
    line("Standard errors",
         paste0(format(round(unname(object["std.error_lower"]), 3), nsmall = 3), ", ",
                format(round(unname(object["std.error_upper"]), 3), nsmall = 3)))
  }
  if (!is.na(object["conf.low"])) {
    line(paste0(level, " Imbens-Manski CI"),
         fmt(unname(object["conf.low"]), unname(object["conf.high"])))
  }
}

format_pct <- function(x, digits = 1) {
  paste0(format(round(100 * x, digits), nsmall = digits, trim = TRUE), "%")
}
