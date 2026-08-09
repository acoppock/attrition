#' Extreme Value Bounds with Double Sampling
#'
#' @param Y The (unquoted) outcome variable, or a formula \code{outcome ~ treatment}
#'   for use with \code{declare_estimator(.method = estimator_ds)}. Must be numeric.
#' @param Z The (unquoted) assignment indicator variable. Must be numeric and take values 0 or 1.
#'   Ignored when \code{Y} is a formula.
#' @param R1 The initial sample response indicator: unquoted column name, or a quoted string column
#'   name when using the formula interface. Must be numeric and take values 0 or 1.
#' @param Attempt The follow-up attempt indicator: unquoted column name, or quoted string.
#'   Must be numeric and take values 0 or 1.
#' @param R2 The follow-up response indicator: unquoted column name, or quoted string.
#'   Must be numeric and take values 0 or 1.
#' @param minY The minimum possible value of the outcome (Y) variable.
#' @param maxY The maximum possible value of the outcome (Y) variable.
#' @param strata Stratification variable: unquoted column name or a quoted string column name.
#' @param alpha The desired significance level. 0.05 by default.
#' @param data A dataframe. Must be given by name: \code{data} is the last
#'   argument, so passing it positionally assigns it to another argument.
#'
#' @return A named numeric vector with elements \code{ci_lower} and \code{ci_upper},
#'   the joint Imbens-Manski confidence interval; \code{low_est} and \code{upp_est},
#'   the bound point estimates; and \code{low_var} and \code{upp_var}, their variances.
#'   Pass to \code{\link[=tidy.attrition_bounds]{tidy()}} for a data frame.
#' @export
#'
#' @examples
#' set.seed(343) # For reproducibility
#' N <- 1000
#'
#' # Potential Outcomes
#' Y_0 <- sample(1:5, N, replace=TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
#' Y_1 <- sample(1:5, N, replace=TRUE, prob = c(0.1, 0.1, 0.4, 0.3, 0.1))
#'
#' R1_0 <- rbinom(N, 1, prob = 0.7)
#' R1_1 <- rbinom(N, 1, prob = 0.8)
#'
#' R2_0 <- rbinom(N, 1, prob = 0.9)
#' R2_1 <- rbinom(N, 1, prob = 0.95)
#'
#' # Covariate
#' strata <- as.numeric(Y_0 > 2)
#'
#' # Random Assignment
#' Z <- rbinom(N, 1, .5)
#'
#' # Reveal Initial Sample Outcomes
#' R1 <- Z*R1_1 + (1-Z)*R1_0 # Initial sample response
#' Y_star <- Z*Y_1 + (1-Z)*Y_0 # True outcomes
#' Y <- Y_star
#' Y[R1==0] <- NA # Mask outcome of non-responders
#'
#' # Conduct Double Sampling
#' Attempt <- rep(0, N)
#' Attempt[is.na(Y)] <- rbinom(sum(is.na(Y)), 1, .5)
#'
#' R2 <- rep(0, N)
#' R2[Attempt==1] <-  (Z*R2_1 + (1-Z)*R2_0)[Attempt==1]
#'
#' Y[R2==1 & Attempt==1] <- Y_star[R2==1 & Attempt==1]
#'
#' df <- data.frame(Y, Z, R1, Attempt, R2, strata)
#'
#' # Without post-stratification
#' estimator_ds(Y, Z, R1, Attempt, R2, minY=1, maxY=5, data=df)
#'
#' # With post-stratification
#' estimator_ds(Y, Z, R1, Attempt, R2, minY=1, maxY=5, strata=strata, data=df)
#'
estimator_ds <- function(Y, Z, R1, Attempt, R2, minY, maxY, strata = NULL, alpha = 0.05, data){
  if (missing(data)) require_data("estimator_ds")
  # Formula interface: estimator_ds(outcome ~ treatment, R1 = "R1", Attempt = "Attempt", R2 = "R2", data = ., ...)
  # R1/Attempt/R2 may be unquoted column names (NSE) or quoted strings.
  yz <- resolve_yz(substitute(Y), substitute(Z), data, parent.frame())
  Y  <- yz$Y
  Z  <- yz$Z
  R1      <- resolve_column(eval(substitute(R1),      data, parent.frame()), data)
  Attempt <- resolve_column(eval(substitute(Attempt), data, parent.frame()), data)
  R2      <- resolve_column(eval(substitute(R2),      data, parent.frame()), data)
  if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
  if(!all(Z %in% c(0,1))){stop("The treatment variable (Z) must be numeric and take values zero or one.")}
  if(!all(R1 %in% c(0,1))){stop("The initial sample response variable (R1) must be numeric and take values zero or one.")}
  if(!all(R2 %in% c(0,1))){stop("The follow-up sample response variable (R2) must be numeric and take values zero or one.")}
  if(!all(Attempt %in% c(0,1))){stop("The follow-up sample attempt variable (Attempt) must be numeric and take values zero or one.")}

  validate_support(Y, minY, maxY, alpha)

  strata <- resolve_column(eval(substitute(strata), data, parent.frame()), data)
  if(is.null(strata)) {
    n1_c <- sum(Z==0)
    n1_t <- sum(Z==1)

    p1_c <- sum(R1==1 & Z==0)/n1_c
    p1_t <- sum(R1==1 & Z==1)/n1_t

    y1m_c <- mean(Y[R1==1 & Z==0])
    y1m_t <- mean(Y[R1==1 & Z==1])

    n2_c <- sum(Attempt ==1 & Z ==0)
    n2_t <- sum(Attempt ==1 & Z ==1)

    p2_c <- sum(R2==1 & Z==0)/n2_c
    p2_t <- sum(R2==1 & Z==1)/n2_t

    y2m_nm_c <- mean(Y[R2==1 & Z==0])
    y2m_nm_t <- mean(Y[R2==1 & Z==1])

    s1_c <- sd(Y[R1==1 & Z==0])
    s1_t <- sd(Y[R1==1 & Z==1])
    s2_nm_c <- sd(Y[R2==1 & Z==0])
    s2_nm_t <- sd(Y[R2==1 & Z==1])

    cis_out <- ds_manski_cis_2s(n1_t=n1_t,n2_t=n2_t,
                                n1_c=n1_c,n2_c=n2_c,
                                p1_t=p1_t,p2_t=p2_t,
                                s1_t=s1_t,s1_c=s1_c,
                                s2_nm_t=s2_nm_t,
                                s2_nm_c=s2_nm_c,
                                y1m_t=y1m_t,y1m_c=y1m_c,
                                y2m_nm_t=y2m_nm_t,
                                y2m_nm_c=y2m_nm_c,
                                p1_c=p1_c,p2_c=p2_c,
                                minY=minY,maxY=maxY,alpha=alpha)
    return(structure(cis_out, class = c("attrition_ds", "attrition_bounds", "numeric"),
                   outcome = yz$outcome))
  }else{
    # With a stratification variable, estimate within each stratum by calling
    # this function recursively, then poststratify.

    if(sum(is.na(strata))!=0){stop("The stratification variable (strata) must not contain any missing values.")}

    unique_strata <- unique(strata)
    ds_df <- data.frame(Y, R1, Z, Attempt, R2, strata)

    strata_ests <- vapply(unique_strata, function(s)
      estimator_ds(Y = Y, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                   minY = minY, maxY = maxY, alpha = alpha,
                   data = subset(ds_df, strata == s)),
      numeric(6))
    proportions <- vapply(unique_strata, function(s) mean(strata == s), numeric(1))

    out <- pool_strata(strata_ests, proportions, alpha)
    return(structure(out, class = c("attrition_ds", "attrition_bounds", "numeric"),
                   outcome = yz$outcome))
  }
}


#' Extreme Value (Manski) Bounds
#'
#' @param Y The (unquoted) outcome variable, or a formula \code{outcome ~ treatment}
#'   for use with \code{declare_estimator(.method = estimator_ev)}. Must be numeric.
#' @param Z The (unquoted) assignment indicator variable. Must be numeric and take values 0 or 1.
#'   Ignored when \code{Y} is a formula.
#' @param R The response indicator variable: unquoted column name, or a quoted string column name
#'   when using the formula interface. Must be numeric and take values 0 or 1.
#' @param minY The minimum possible value of the outcome (Y) variable.
#' @param maxY The maximum possible value of the outcome (Y) variable.
#' @param strata Stratification variable: unquoted column name or a quoted string column name.
#' @param alpha The desired significance level. 0.05 by default.
#' @param data A dataframe. Must be given by name: \code{data} is the last
#'   argument, so passing it positionally assigns it to another argument.
#'
#' @return A named numeric vector with elements \code{ci_lower} and \code{ci_upper},
#'   the joint Imbens-Manski confidence interval; \code{low_est} and \code{upp_est},
#'   the bound point estimates; and \code{low_var} and \code{upp_var}, their variances.
#'   Pass to \code{\link[=tidy.attrition_bounds]{tidy()}} for a data frame.
#' @export
#'
#' @examples
#' set.seed(343)
#' N <- 1000
#' Y_0 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
#' Y_1 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.1, 0.4, 0.3, 0.1))
#' Z <- rbinom(N, 1, 0.5)
#' Y_star <- Z * Y_1 + (1 - Z) * Y_0
#'
#' # Treated units respond at a higher rate, so the missingness is nonignorable
#' R <- rbinom(N, 1, prob = 0.7 + 0.1 * Z)
#' Y <- Y_star
#' Y[R == 0] <- NA
#' df <- data.frame(Y, Z, R)
#'
#' estimator_ev(Y, Z, R, minY = 1, maxY = 5, data = df)
#'
#' # Equivalently, via the formula interface
#' estimator_ev(Y ~ Z, R = "R", minY = 1, maxY = 5, data = df)
estimator_ev <- function(Y, Z, R, minY, maxY, strata = NULL, alpha = 0.05, data){
  if (missing(data)) require_data("estimator_ev")
  # Formula interface: estimator_ev(outcome ~ treatment, R = "col_name", data = ., ...)
  yz <- resolve_yz(substitute(Y), substitute(Z), data, parent.frame())
  Y  <- yz$Y
  Z  <- yz$Z
  R  <- resolve_column(eval(substitute(R), data, parent.frame()), data)
  if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
  if(!all(Z %in% c(0,1))){stop("The treatment variable (Z) must be numeric and take values zero or one.")}
  if(!all(R %in% c(0,1))){stop("The response variable (R) must be numeric and take values zero or one.")}

  validate_support(Y, minY, maxY, alpha)

  strata <- resolve_column(eval(substitute(strata), data, parent.frame()), data)
  if(is.null(strata)) {
    n1_c <- sum(Z==0)
    n1_t <- sum(Z==1)

    p1_c <- sum(R==1 & Z==0)/n1_c
    p1_t <- sum(R==1 & Z==1)/n1_t

    y1m_c <- mean(Y[R==1 & Z==0])
    y1m_t <- mean(Y[R==1 & Z==1])

    s1_c <- sd(Y[R==1 & Z==0])
    s1_t <- sd(Y[R==1 & Z==1])

    cis_out <- manski_cis(n1_t = n1_t, n1_c = n1_c,
                          p1_t = p1_t, p1_c = p1_c,
                          y1m_t = y1m_t, y1m_c = y1m_c,
                          s1_t = s1_t, s1_c = s1_c,
                          minY = minY, maxY = maxY, alpha = alpha)

    return(structure(cis_out, class = c("attrition_ev", "attrition_bounds", "numeric"),
                   outcome = yz$outcome))
  }else{
    # With a stratification variable, estimate within each stratum by calling
    # this function recursively, then poststratify.

    if(sum(is.na(strata))!=0){stop("The stratification variable (strata) must not contain any missing values.")}

    unique_strata <- unique(strata)
    ds_df <- data.frame(Y, Z, R, strata)

    strata_ests <- vapply(unique_strata, function(s)
      estimator_ev(Y = Y, Z = Z, R = R,
                   minY = minY, maxY = maxY, alpha = alpha,
                   data = subset(ds_df, strata == s)),
      numeric(6))
    proportions <- vapply(unique_strata, function(s) mean(strata == s), numeric(1))

    out <- pool_strata(strata_ests, proportions, alpha)
    return(structure(out, class = c("attrition_ev", "attrition_bounds", "numeric"),
                   outcome = yz$outcome))
  }

}

#' Trimming Bounds
#'
#' @param Y The (unquoted) outcome variable, or a formula \code{outcome ~ treatment}
#'   for use with \code{declare_estimator(.method = estimator_trim)}. Must be numeric.
#' @param Z The (unquoted) assignment indicator variable. Must be numeric and take values 0 or 1.
#'   Ignored when \code{Y} is a formula.
#' @param R The single-stage response indicator: unquoted column name, or a quoted string column
#'   name when using the formula interface. Must be numeric and take values 0 or 1.
#'   Supply either \code{R} (single-stage) or \code{R1}/\code{Attempt}/\code{R2} (double-sampling).
#' @param R1 The initial sample response indicator. Unquoted or quoted string column name.
#'   Must be numeric and take values 0 or 1.
#' @param Attempt The follow-up attempt indicator. Unquoted or quoted string column name.
#'   Must be numeric and take values 0 or 1.
#' @param R2 The follow-up response indicator. Unquoted or quoted string column name.
#'   Must be numeric and take values 0 or 1.
#' @param strata Not supported; supplying any value raises an error.
#' @param alpha The desired significance level. 0.05 by default.
#' @param se How to obtain standard errors. \code{"analytic"} (the default) uses the
#'   closed-form asymptotic variance of Lee (2009), Proposition 3, and is available for
#'   the single-stage \code{R} path only. \code{"bootstrap"} resamples units within
#'   treatment arm and works for both paths. \code{"none"} returns bounds alone.
#' @param sims Number of bootstrap replicates when \code{se = "bootstrap"}. 1000 by default.
#' @param data A dataframe. Must be given by name: \code{data} is the last
#'   argument, so passing it positionally assigns it to another argument.
#'
#' @return A named numeric vector containing \code{lower_bound} and \code{upper_bound},
#'   the trimming bound estimates; \code{lower_se} and \code{upper_se}, their standard
#'   errors; \code{ci_lower} and \code{ci_upper}, the joint Imbens-Manski confidence
#'   interval; and the intermediate quantities used to build them. All elements are
#'   \code{NA} when monotonicity is violated. Pass to
#'   \code{\link[=tidy.attrition_trim]{tidy()}} for a data frame.
#'
#' @details
#' The analytic variance has four contributions: the variance of the retained
#' (trimmed) outcomes, the variance from estimating the trimming threshold, the
#' variance from estimating the trimming proportion, and the variance of the
#' control-group respondent mean. The third of these is often the largest, so
#' treating the trimming proportion as known would understate the uncertainty
#' substantially.
#'
#' Lee's derivation assumes the bounds are interior points, which fails when the
#' two response rates are equal: the trimming proportion is then zero, the bounds
#' collapse to a point, and the standard errors are not trustworthy. This case
#' warns.
#'
#' @references
#' Lee, David S. (2009). Training, Wages, and Sample Selection: Estimating Sharp
#' Bounds on Treatment Effects. \emph{Review of Economic Studies} 76(3):1071-1102.
#'
#' Tauchmann, Harald (2014). Lee (2009) Treatment-Effect Bounds for Nonrandom
#' Sample Selection. \emph{Stata Journal} 14(4):884-894.
#' @export
#'
#' @examples
#' set.seed(343)
#' N <- 1000
#' Y_0 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
#' Y_1 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.1, 0.4, 0.3, 0.1))
#' Z <- rbinom(N, 1, 0.5)
#' Y_star <- Z * Y_1 + (1 - Z) * Y_0
#'
#' # Treated units respond at a higher rate, so the missingness is nonignorable
#' R <- rbinom(N, 1, prob = 0.7 + 0.1 * Z)
#' Y <- Y_star
#' Y[R == 0] <- NA
#' df <- data.frame(Y, Z, R)
#'
#' # Single-stage: trimming bounds under monotonicity, with Lee (2009) standard errors
#' estimator_trim(Y, Z, R = R, data = df)
#'
#' # Bootstrap standard errors instead
#' estimator_trim(Y, Z, R = R, se = "bootstrap", sims = 200, data = df)
estimator_trim <-
  function(Y, Z, R = NULL, R1 = NULL, Attempt = NULL, R2 = NULL, strata = NULL,
           alpha = 0.05, se = c("analytic", "bootstrap", "none"), sims = 1000, data){
    if (missing(data)) require_data("estimator_trim")
    # Formula interface: estimator_trim(outcome ~ treatment, R = "col" | R1/Attempt/R2 = "col", data = .)
    yz <- resolve_yz(substitute(Y), substitute(Z), data, parent.frame())
    Y  <- yz$Y
    Z  <- yz$Z

    # Resolve R / R1 / Attempt / R2 — each accepts NSE or quoted string
    R           <- resolve_column(eval(substitute(R),       data, parent.frame()), data)
    R1_val      <- resolve_column(eval(substitute(R1),      data, parent.frame()), data)
    Attempt_val <- resolve_column(eval(substitute(Attempt), data, parent.frame()), data)
    R2_val      <- resolve_column(eval(substitute(R2),      data, parent.frame()), data)

    strata <- resolve_column(eval(substitute(strata), data, parent.frame()), data)
    if (!is.null(strata)) stop("Stratification is not yet supported for trimming bounds.")
    if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
    if(!all(Z %in% c(0,1))){stop("The treatment variable (Z) must be numeric and take values zero or one.")}
    se <- match.arg(se)
    if(!is.numeric(alpha) | length(alpha) != 1L){stop("The significance level (alpha) must be a single number.")}
    if(alpha <= 0 | alpha >= 1){stop("The significance level (alpha) must be strictly between zero and one.")}
    if(se == "bootstrap" && (!is.numeric(sims) | length(sims) != 1L | any(sims < 2))){
      stop("The number of bootstrap replicates (sims) must be a single number of at least two.")
    }

    na_trim <- structure(
      c(lower_bound = NA_real_, upper_bound = NA_real_, Q = NA_real_,
        lower_se = NA_real_, upper_se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_),
      class = c("attrition_trim", "numeric")
    )

    # One estimation routine per path, taking row indices, so the point estimate
    # and each bootstrap replicate go through identical code. The double-sampling
    # weights are recomputed inside, since they depend on the resampled counts.
    if (!is.null(R)) {
      if(!all(R %in% c(0,1))){stop("The response variable (R) must be numeric and take values zero or one.")}
      single_stage <- TRUE
      estimate <- function(idx) {
        trimming_bounds(Out = Y[idx], Treat = Z[idx], Fail = as.numeric(R[idx] == 0),
                        Weight = rep(1, length(idx)), monotonicity = TRUE)
      }
    } else {
      if (is.null(R1_val) || is.null(Attempt_val) || is.null(R2_val)) {
        stop("Supply either R (single-stage) or R1, Attempt, and R2 (double-sampling).")
      }
      R1      <- R1_val
      Attempt <- Attempt_val
      R2      <- R2_val
      if(!all(R1 %in% c(0,1))){stop("The initial sample response variable (R1) must be numeric and take values zero or one.")}
      if(!all(R2 %in% c(0,1))){stop("The follow-up sample response variable (R2) must be numeric and take values zero or one.")}
      if(!all(Attempt %in% c(0,1))){stop("The follow-up sample attempt variable (Attempt) must be numeric and take values zero or one.")}
      if (se == "analytic") {
        stop("Analytic standard errors follow Lee (2009) Proposition 3, which covers the ",
             "single-stage, unweighted, monotonicity case only. The double-sampling ",
             "estimator trims both groups and carries sampling weights, so use ",
             "se = \"bootstrap\" (or se = \"none\").")
      }
      single_stage <- FALSE
      estimate <- function(idx) {
        Yi <- Y[idx]; Zi <- Z[idx]; R1i <- R1[idx]; Ai <- Attempt[idx]; R2i <- R2[idx]
        Weight <- rep(NA, length(idx))
        Weight[R1i==1] <- 1
        Weight[Ai==1 & Zi==1] <- sum(Zi == 1 & R1i == 0)/sum(Zi == 1 & Ai == 1)
        Weight[Ai==1 & Zi==0] <- sum(Zi == 0 & R1i == 0)/sum(Zi == 0 & Ai == 1)
        Fail <- as.numeric(R1i == 0 & R2i == 0)
        Keep <- (R1i == 1 | Ai == 1)
        trimming_bounds(Out = Yi[Keep], Treat = Zi[Keep],
                        Fail = Fail[Keep], Weight = Weight[Keep], monotonicity = FALSE)
      }
    }

    all_rows <- seq_along(Y)
    out <- tryCatch(estimate(all_rows),
                    attrition_monotonicity_violation = \(e) NULL)
    if (is.null(out)) {
      return(structure(na_trim, outcome = yz$outcome))
    }

    variances <- c(lower_var = NA_real_, upper_var = NA_real_)
    if (se == "analytic") {
      if (unname(out["Q"]) <= 0) {
        warning("The trimming proportion is zero, so the bounds collapse to a point and sit ",
                "on the boundary of the parameter space. Lee (2009) Proposition 3 assumes an ",
                "interior point; the standard errors below are not reliable here.", call. = FALSE)
      }
      variances <- lee_variance(out, n_treat = sum(Z == 1), n_control = sum(Z == 0))
    } else if (se == "bootstrap") {
      boot <- bootstrap_trim_variance(
        function(idx) estimate(idx)[c("lower_bound", "upper_bound")], Z, sims)
      variances <- boot[c("lower_var", "upper_var")]
    }

    ci <- c(ci_lower = NA_real_, ci_upper = NA_real_)
    if (se != "none") {
      sig <- im_critical_value(unname(out["lower_bound"]), unname(out["upper_bound"]),
                               unname(variances["lower_var"]), unname(variances["upper_var"]), alpha)
      ci <- c(ci_lower = unname(out["lower_bound"]) - sig*unname(variances["lower_var"])^.5,
              ci_upper = unname(out["upper_bound"]) + sig*unname(variances["upper_var"])^.5)
    }

    # Drop the quantities that exist only to feed lee_variance
    out <- out[!names(out) %in% c("var_keep_U", "var_keep_L", "n_keep_U", "n_keep_L", "var_control")]
    out <- c(out,
             lower_se = unname(variances["lower_var"])^.5,
             upper_se = unname(variances["upper_var"])^.5,
             ci)
    return(structure(out, class = c("attrition_trim", "numeric"),
                     se_method = se, single_stage = single_stage,
                     outcome = yz$outcome))
  }


#' Extreme Value Bounds with Double Sampling with Sensitivity
#'
#' This function yields extreme value bounds under the assumption that the outcomes of 1-delta of the missing second-round units are ignorable, that is, that they are drawn from an unknown distribution with mean and variance equal to the observed second-round groups.
#'
#' @param Y The (unquoted) outcome variable, or a formula \code{outcome ~ treatment}
#'   for use with \code{declare_estimator(.method = estimator_ds_sens)}. Must be numeric.
#' @param Z The (unquoted) assignment indicator variable. Must be numeric and take values 0 or 1.
#'   Ignored when \code{Y} is a formula.
#' @param R1 The initial sample response indicator: unquoted column name, or a quoted string column
#'   name when using the formula interface. Must be numeric and take values 0 or 1.
#' @param Attempt The follow-up attempt indicator: unquoted column name, or quoted string.
#'   Must be numeric and take values 0 or 1.
#' @param R2 The follow-up response indicator: unquoted column name, or quoted string.
#'   Must be numeric and take values 0 or 1.
#' @param minY The minimum possible value of the outcome (Y) variable.
#' @param maxY The maximum possible value of the outcome (Y) variable.
#' @param strata Stratification variable: unquoted column name or a quoted string column name.
#' @param alpha The desired significance level. 0.05 by default.
#' @param data A dataframe
#' @param delta Sensitivity parameter in [0, 1]. At delta = 1 (default) worst-case bounds apply; at delta = 0 ignorability holds for all follow-up non-responders.
#'
#' @return A named numeric vector with elements \code{ci_lower} and \code{ci_upper},
#'   the joint Imbens-Manski confidence interval; \code{low_est} and \code{upp_est},
#'   the bound point estimates; and \code{low_var} and \code{upp_var}, their variances.
#'   Pass to \code{\link[=tidy.attrition_bounds]{tidy()}} for a data frame.
#' @export
#'
#' @examples
#' set.seed(343)
#' N <- 1000
#' Y_0 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
#' Y_1 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.1, 0.4, 0.3, 0.1))
#' Z <- rbinom(N, 1, 0.5)
#' Y_star <- Z * Y_1 + (1 - Z) * Y_0
#' R1 <- rbinom(N, 1, prob = 0.7 + 0.1 * Z)
#' Y <- Y_star
#' Y[R1 == 0] <- NA
#'
#' # Follow up intensively with a random half of the initial non-responders
#' Attempt <- rep(0, N)
#' Attempt[R1 == 0] <- rbinom(sum(R1 == 0), 1, 0.5)
#' R2 <- rep(0, N)
#' R2[Attempt == 1] <- rbinom(sum(Attempt == 1), 1, 0.9)
#' Y[Attempt == 1 & R2 == 1] <- Y_star[Attempt == 1 & R2 == 1]
#' df <- data.frame(Y, Z, R1, Attempt, R2)
#'
#' # delta = 1 reproduces the worst-case double-sampling bounds
#' estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 1, data = df)
#'
#' # delta = 0 assumes ignorability among follow-up non-responders
#' estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 0, data = df)
estimator_ds_sens <- function(Y, Z, R1, Attempt, R2, minY, maxY, delta, strata = NULL, alpha = 0.05, data){
  if (missing(data)) require_data("estimator_ds_sens")
  # Formula interface: estimator_ds_sens(outcome ~ treatment, R1 = "R1", Attempt = "Attempt", R2 = "R2", data = ., ...)
  yz <- resolve_yz(substitute(Y), substitute(Z), data, parent.frame())
  Y  <- yz$Y
  Z  <- yz$Z
  R1      <- resolve_column(eval(substitute(R1),      data, parent.frame()), data)
  Attempt <- resolve_column(eval(substitute(Attempt), data, parent.frame()), data)
  R2      <- resolve_column(eval(substitute(R2),      data, parent.frame()), data)
  if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
  if(!all(Z %in% c(0,1))){stop("The treatment variable (Z) must be numeric and take values zero or one.")}
  if(!all(R1 %in% c(0,1))){stop("The initial sample response variable (R1) must be numeric and take values zero or one.")}
  if(!all(R2 %in% c(0,1))){stop("The follow-up sample response variable (R2) must be numeric and take values zero or one.")}
  if(!all(Attempt %in% c(0,1))){stop("The follow-up sample attempt variable (Attempt) must be numeric and take values zero or one.")}

  validate_support(Y, minY, maxY, alpha)
  if(!is.numeric(delta) | length(delta) != 1L){stop("The sensitivity parameter (delta) must be a single number.")}
  if(delta < 0 | delta > 1){stop("The sensitivity parameter (delta) must be between zero and one.")}

  strata <- resolve_column(eval(substitute(strata), data, parent.frame()), data)
  if(is.null(strata)) {
    n1_c <- sum(Z==0)
    n1_t <- sum(Z==1)

    p1_c <- sum(R1==1 & Z==0)/n1_c
    p1_t <- sum(R1==1 & Z==1)/n1_t

    y1m_c <- mean(Y[R1==1 & Z==0])
    y1m_t <- mean(Y[R1==1 & Z==1])

    n2_c <- sum(Attempt ==1 & Z ==0)
    n2_t <- sum(Attempt ==1 & Z ==1)

    p2_c <- sum(R2==1 & Z==0)/n2_c
    p2_t <- sum(R2==1 & Z==1)/n2_t

    y2m_nm_c <- mean(Y[R2==1 & Z==0])
    y2m_nm_t <- mean(Y[R2==1 & Z==1])

    s1_c <- sd(Y[R1==1 & Z==0])
    s1_t <- sd(Y[R1==1 & Z==1])
    s2_nm_c <- sd(Y[R2==1 & Z==0])
    s2_nm_t <- sd(Y[R2==1 & Z==1])

    cis_out <- ds_manski_cis_2s_sens(n1_t=n1_t,n2_t=n2_t,
                                     n1_c=n1_c,n2_c=n2_c,
                                     p1_t=p1_t,p2_t=p2_t,
                                     s1_t=s1_t,s1_c=s1_c,
                                     s2_nm_t=s2_nm_t,
                                     s2_nm_c=s2_nm_c,
                                     y1m_t=y1m_t,y1m_c=y1m_c,
                                     y2m_nm_t=y2m_nm_t,
                                     y2m_nm_c=y2m_nm_c,
                                     p1_c=p1_c,p2_c=p2_c,
                                     minY=minY,maxY=maxY,alpha=alpha, delta = delta)
    return(structure(cis_out, class = c("attrition_ds_sens", "attrition_bounds", "numeric"),
                   outcome = yz$outcome))
  }else{
    # With a stratification variable, estimate within each stratum by calling
    # this function recursively, then poststratify.

    if(sum(is.na(strata))!=0){stop("The stratification variable (strata) must not contain any missing values.")}

    unique_strata <- unique(strata)
    ds_df <- data.frame(Y, R1, Z, Attempt, R2, strata)

    strata_ests <- vapply(unique_strata, function(s)
      estimator_ds_sens(Y = Y, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                        minY = minY, maxY = maxY, alpha = alpha, delta = delta,
                        data = subset(ds_df, strata == s)),
      numeric(6))
    proportions <- vapply(unique_strata, function(s) mean(strata == s), numeric(1))

    out <- pool_strata(strata_ests, proportions, alpha)
    return(structure(out, class = c("attrition_ds_sens", "attrition_bounds", "numeric"),
                   outcome = yz$outcome))
  }
}


#' Sensitivity Analysis
#'
#' This function performs a line search over values of delta, the sensitivity parameter, in order to find (if it exists) delta*, the value of delta where the confidence interval no longer includes zero.
#'
#' @param Y The (unquoted) outcome variable, or a formula \code{outcome ~ treatment}.
#'   Must be numeric.
#' @param Z The (unquoted) assignment indicator variable. Must be numeric and take values 0 or 1.
#'   Ignored when \code{Y} is a formula.
#' @param R1 The initial sample response indicator: unquoted column name, or a quoted string column
#'   name when using the formula interface. Must be numeric and take values 0 or 1.
#' @param Attempt The follow-up attempt indicator: unquoted column name, or quoted string.
#'   Must be numeric and take values 0 or 1.
#' @param R2 The follow-up response indicator: unquoted column name, or quoted string.
#'   Must be numeric and take values 0 or 1.
#' @param minY The minimum possible value of the outcome (Y) variable.
#' @param maxY The maximum possible value of the outcome (Y) variable.
#' @param strata Stratification variable: unquoted column name or a quoted string column name.
#' @param alpha The desired significance level. 0.05 by default.
#' @param data A dataframe. Must be given by name: \code{data} is the last
#'   argument, so passing it positionally assigns it to another argument.
#' @param sims Number of values of delta at which to evaluate the bounds. Defaults to 100.
#'
#' @return A list with three elements: \code{sensitivity_plot}, a ggplot object;
#'   \code{sims_df}, a data frame of bounds and confidence intervals at each delta;
#'   and \code{p_star}, a one-row data frame giving delta* when it exists and an
#'   explanatory character string when it does not.
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point geom_text
#'   geom_hline xlab ylab theme_bw theme element_blank
#' @importFrom grid unit
#' @importFrom purrr map
#' @importFrom stats complete.cases pnorm qnorm sd setNames uniroot var weighted.mean
#' @export
#'
#' @examples
#' set.seed(343)
#' N <- 1000
#' Y_0 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
#' Y_1 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.1, 0.4, 0.3, 0.1))
#' Z <- rbinom(N, 1, 0.5)
#' Y_star <- Z * Y_1 + (1 - Z) * Y_0
#' R1 <- rbinom(N, 1, prob = 0.7 + 0.1 * Z)
#' Y <- Y_star
#' Y[R1 == 0] <- NA
#'
#' # Follow up intensively with a random half of the initial non-responders
#' Attempt <- rep(0, N)
#' Attempt[R1 == 0] <- rbinom(sum(R1 == 0), 1, 0.5)
#' R2 <- rep(0, N)
#' R2[Attempt == 1] <- rbinom(sum(Attempt == 1), 1, 0.9)
#' Y[Attempt == 1 & R2 == 1] <- Y_star[Attempt == 1 & R2 == 1]
#' df <- data.frame(Y, Z, R1, Attempt, R2)
#'
#' sens <- sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5,
#'                        sims = 20, data = df)
#' sens$sensitivity_plot
#' sens$p_star
sensitivity_ds <- function(Y, Z, R1, Attempt, R2, minY, maxY, sims = 100, strata = NULL, alpha = 0.05, data){
  if (missing(data)) require_data("sensitivity_ds")
  # Formula interface: sensitivity_ds(outcome ~ treatment, R1 = "R1", Attempt = "Attempt", R2 = "R2", data = ., ...)
  yz <- resolve_yz(substitute(Y), substitute(Z), data, parent.frame())
  Y  <- yz$Y
  Z  <- yz$Z
  R1      <- resolve_column(eval(substitute(R1),      data, parent.frame()), data)
  Attempt <- resolve_column(eval(substitute(Attempt), data, parent.frame()), data)
  R2      <- resolve_column(eval(substitute(R2),      data, parent.frame()), data)
  if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
  if(!all(Z %in% c(0,1))){stop("The treatment variable (Z) must be numeric and take values zero or one.")}
  if(!all(R1 %in% c(0,1))){stop("The initial sample response variable (R1) must be numeric and take values zero or one.")}
  if(!all(R2 %in% c(0,1))){stop("The follow-up sample response variable (R2) must be numeric and take values zero or one.")}
  if(!all(Attempt %in% c(0,1))){stop("The follow-up sample attempt variable (Attempt) must be numeric and take values zero or one.")}
  validate_support(Y, minY, maxY, alpha)
  if(!is.numeric(sims) | length(sims) != 1L | any(sims < 2)){stop("The number of simulations (sims) must be a single number of at least two.")}

  strata <- resolve_column(eval(substitute(strata), data, parent.frame()), data)

  ps <- seq(0, 1, length.out = sims)

  if (is.null(strata)) {
    df <- data.frame(Y, Z, R1, R2, Attempt)
    sims_df <-
      map(ps, \(d) estimator_ds_sens(Y = Y, Z = Z, R1 = R1, Attempt = Attempt, alpha = alpha,
                                     R2 = R2, minY = minY, maxY = maxY, data = df, delta = d)) |>
      (\(lst) do.call(rbind, lst))() |>
      as.data.frame() |>
      dplyr::mutate(p = ps,
                    change_lower = find_sign_changes(ci_lower),
                    change_upper = find_sign_changes(ci_upper),
                    change_any = change_lower | change_upper)
  } else {
    df <- data.frame(Y, Z, R1, R2, Attempt, strata)
    sims_df <-
      map(ps, \(d) estimator_ds_sens(Y = Y, Z = Z, R1 = R1, Attempt = Attempt,
                                     strata = strata, alpha = alpha,
                                     R2 = R2, minY = minY, maxY = maxY, data = df, delta = d)) |>
      (\(lst) do.call(rbind, lst))() |>
      as.data.frame() |>
      dplyr::mutate(p = ps,
                    change_lower = find_sign_changes(ci_lower),
                    change_upper = find_sign_changes(ci_upper),
                    change_any = change_lower | change_upper)
  }


  points_df <-
    data.frame(p = c(0, 1, 1),
               value = c(with(sims_df, low_est[p==0]),
                         with(sims_df, low_est[p==1]),
                         with(sims_df, upp_est[p==1])),
               hjust = c(-.3, 1.1, 1.1),
               vjust = c(NA, 1, -1),
               label = c("Naive Estimate", "Worst Case Lower Bound", "Worst Case Upper Bound"))

  g <-
    ggplot(sims_df, aes(x = p)) +
    geom_line(aes(y = upp_est), alpha = 0.5) +
    geom_line(aes(y = low_est), alpha = 0.5) +
    geom_ribbon(aes(ymax = ci_upper, ymin = ci_lower), alpha = 0.2) +
    geom_point(data = points_df, aes(y = value)) +
    geom_text(data = points_df, aes(y = value, label = label, hjust = hjust, vjust = vjust)) +
    ylab(paste0("Identification Regions and ", round((1-alpha)*100), "% Confidence Intervals")) +
    xlab(expression(paste("Sensitivity Parameter ", delta, " (0 = Ignorability)"))) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.key.width = unit(3, "lines"),
          legend.title = element_blank())


  p_star_df <- "No value of the sensitivity parameter yields a statistically significant result."

  p_star <- with(sims_df, p[change_any])
  if(length(p_star) >= 1){
    p_star <- min(p_star)
    p_star_df <- data.frame(p = p_star,
                            value = 0,
                            label = paste0("delta^'*' == ", round(p_star, 2)),
                            hjust = ifelse(p_star > 0.5, 1.1, -1.1),
                            vjust = ifelse(with(sims_df, low_est[p==0]) > 0, 1.3, -1.3))

    g <- g + geom_point(data = p_star_df, aes(y = value)) +
      geom_text(data = p_star_df, aes(label = label, y = value,  vjust = vjust), parse = TRUE)
  }

  return(list(sensitivity_plot = g, sims_df = sims_df, p_star = p_star_df))

}
