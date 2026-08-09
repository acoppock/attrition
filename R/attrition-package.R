#' attrition: bounds for experiments with missing outcomes
#'
#' When subjects go missing from an experiment and the reason they went missing
#' is related to what their outcome would have been, no amount of covariate
#' adjustment will fix the problem. The package takes the other route: rather
#' than assume the missingness away, it reports the range of average treatment
#' effects consistent with the data, and it offers a research design that makes
#' that range small enough to be useful.
#'
#' @section The design:
#' Worst-case bounds fill in every missing outcome with the smallest and largest
#' values the outcome could take. The resulting interval is honest, and it is
#' usually far too wide to settle anything. Double sampling narrows it. After the
#' first round of data collection, draw a random sample of the nonrespondents and
#' pursue them harder: pay more, call again, send an interviewer. Because those
#' subjects are a random sample of the nonrespondents, their recovered outcomes
#' stand in for all of them, and only the residual group who refuse twice needs
#' worst-case treatment. In the application shipped with the package, chasing 100
#' of 536 nonrespondents cut the width of the 95 percent confidence interval from
#' 3.50 to 1.23.
#'
#' @section The estimators:
#' \describe{
#'   \item{\code{\link{estimator_ev}}}{Worst-case (Manski) bounds from a single
#'     round of data collection.}
#'   \item{\code{\link{estimator_ds}}}{Double-sampling bounds, with analytic
#'     variances and Imbens-Manski confidence intervals. The estimator of
#'     Coppock, Gerber, Green, and Kern (2017).}
#'   \item{\code{\link{estimator_ds_sens}}}{Double-sampling bounds at a chosen
#'     value of delta, the fraction of follow-up nonrespondents for whom
#'     ignorability is allowed to fail.}
#'   \item{\code{\link{sensitivity_ds}}}{A search over delta for the point at
#'     which the confidence interval starts to include zero.}
#'   \item{\code{\link{estimator_trim}}}{Lee (2009) trimming bounds, which
#'     assume monotone selection instead of a bounded outcome.}
#' }
#'
#' \code{estimator_ev}, \code{estimator_ds}, and \code{estimator_ds_sens} accept
#' a \code{strata} argument for poststratification on a discrete covariate. The
#' identified set is the same either way; poststratification estimates it more
#' precisely, and by the law of total variance the asymptotic variance is no
#' larger. Every estimator has a \code{\link[=tidy.attrition_bounds]{tidy()}} method
#' and a formula interface for use with \pkg{DeclareDesign}.
#'
#' @section Where to start:
#' \code{vignette("attrition")} walks through the design and all five estimators
#' on the replication data in \code{\link{levendusky}}, reproducing the published
#' table as it goes.
#'
#' @references
#' Coppock, Alexander, Alan S. Gerber, Donald P. Green, and Holger L. Kern (2017).
#' Combining Double Sampling and Bounds to Address Nonignorable Missing Outcomes
#' in Randomized Experiments. \emph{Political Analysis} 25(2):188-206.
#' \doi{10.1017/pan.2016.6}
#'
#' Imbens, Guido W., and Charles F. Manski (2004). Confidence Intervals for
#' Partially Identified Parameters. \emph{Econometrica} 72(6):1845-1857.
#'
#' Lee, David S. (2009). Training, Wages, and Sample Selection: Estimating Sharp
#' Bounds on Treatment Effects. \emph{Review of Economic Studies} 76(3):1071-1102.
#'
#' @keywords internal
"_PACKAGE"
