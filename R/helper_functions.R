# `data` is this package's last argument rather than its second, so the call a
# user reflexively types, estimator_ds(y ~ z, df, ...), binds df to the second
# formal and R reports a missing `data` without saying why. The order is kept
# for backward compatibility; this replaces the unhelpful message.
require_data <- function(fn) {
  stop("`data` must be given by name, as ", fn, "(y ~ z, ..., data = your_data).\n",
       "  It is this function's last argument, so passing a data frame ",
       "positionally assigns it to another argument instead.", call. = FALSE)
}

# A column argument may be given four ways: as a vector of the column contents,
# as a single string naming a column, as a one-sided formula naming a column, or
# as a bare column name under normal evaluation. The formula case is what
# DeclareDesign produces: declare_estimator(..., R1 = R1) hands the method ~R1
# rather than the symbol, so a method that does not understand formulas can only
# be driven with strings there.
resolve_column <- function(x, data) {
  if (inherits(x, "formula")) {
    vars <- all.vars(x)
    if (length(vars) != 1L) {
      stop("A column given as a formula must name exactly one column, as in R1 = ~R1.")
    }
    if (!vars %in% names(data)) stop("Column '", vars, "' was not found in the data.")
    return(data[[vars]])
  }
  if (is.character(x) && length(x) == 1L) return(data[[x]])
  x
}

# Resolve the outcome and treatment arguments. Y is either an outcome vector or
# a Y ~ Z formula; when it is a formula, Z is taken from the formula and the Z
# argument is never evaluated, so callers may leave it missing.
resolve_yz <- function(Y_expr, Z_expr, data, env) {
  Y_val <- eval(Y_expr, data, env)
  if (inherits(Y_val, "formula")) {
    yz <- parse_yz_formula(Y_val, data)
    yz$outcome <- all.vars(Y_val)[1L]
    return(yz)
  }
  list(Y = Y_val, Z = eval(Z_expr, data, env),
       outcome = if (is.symbol(Y_expr)) as.character(Y_expr) else NA_character_)
}

# Parse a Y ~ Z formula, returning c(outcome_col, treatment_col).
# Errors if the formula has more than two variables.
parse_yz_formula <- function(f, data) {
  if (!inherits(f, "formula")) stop("'Y' must be a formula or an unquoted column name.")
  vars <- all.vars(f)
  if (length(vars) != 2L) {
    stop(
      "Formula must be outcome ~ treatment with exactly two variables. ",
      "Additional variables (R, R1, Attempt, R2) are passed as quoted ",
      "column name strings, e.g. estimator_ev(Y ~ Z, R = \"R\", ...)."
    )
  }
  list(Y = data[[vars[1L]]], Z = data[[vars[2L]]])
}

# The three double-sampling response arguments, resolved the same way in every
# estimator that takes them. The expressions are the caller's substitute()d
# arguments, so bare column names, strings, and one-sided formulas all work.
resolve_ds_columns <- function(R1_expr, Attempt_expr, R2_expr, data, env) {
  list(R1 = resolve_column(eval(R1_expr, data, env), data),
       Attempt = resolve_column(eval(Attempt_expr, data, env), data),
       R2 = resolve_column(eval(R2_expr, data, env), data))
}

# Flag the first element whose sign differs from the first element's, treating a
# zero as a change. Position 1 can never be flagged; the last position can.
find_sign_changes <- function(x){
  first_change <- Position(function(xi) sign(xi) != sign(x[1]), x)
  out <- rep(FALSE, length(x))
  if(!is.na(first_change)){out[first_change] <- TRUE}
  return(out)
}

# Input checks ----

validate_indicator <- function(x, what) {
  if (!all(x %in% c(0, 1))) {
    stop("The ", what, " must be numeric and take values zero or one.")
  }
}

# Shared argument checks for the bounding estimators. The assumed support of Y
# has to cover the observed outcomes, or the bounds it produces are not bounds.
validate_support <- function(Y, minY, maxY, alpha){
  if(!is.numeric(minY) | !is.numeric(maxY)){stop("The minimum and maximum possible values of Y (minY and maxY) must be numeric")}
  if(minY > maxY){stop("The minimum possible value of Y (minY) must not be greater than the maximum (maxY).")}
  if(any(Y < minY | Y > maxY, na.rm = TRUE)){stop("Some observed outcomes fall outside the assumed support of Y: widen minY and maxY.")}
  if(!is.numeric(alpha) | length(alpha) != 1L){stop("The significance level (alpha) must be a single number.")}
  if(alpha <= 0 | alpha >= 1){stop("The significance level (alpha) must be strictly between zero and one.")}
}

# The checks every double-sampling estimator runs on its resolved arguments.
validate_ds_inputs <- function(Y, Z, R1, Attempt, R2, minY, maxY, alpha) {
  if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
  validate_indicator(Z, "treatment variable (Z)")
  validate_indicator(R1, "initial sample response variable (R1)")
  validate_indicator(R2, "follow-up sample response variable (R2)")
  validate_indicator(Attempt, "follow-up sample attempt variable (Attempt)")
  validate_support(Y, minY, maxY, alpha)
}

# The bounds are built from arm-level means, and a mean over an empty cell is
# NaN rather than a bound. Each arm needs at least one unit of each kind.
require_in_each_arm <- function(flag, Z, what) {
  for (z in c(1, 0)) {
    if (!any(flag & Z == z)) {
      stop("The ", if (z == 1) "treatment" else "control", " group has no ", what,
           ", so the bounds are undefined. When strata are supplied, every stratum ",
           "needs at least one in each group.", call. = FALSE)
    }
  }
}

# Moments ----

# The arm-level quantities every double-sampling bound is built from. For each
# arm: its size, the initial response rate, the mean and standard deviation of
# the initial respondents' outcomes, the number of follow-up attempts, the
# follow-up response rate, and the mean and standard deviation among the
# follow-up respondents.
ds_moments <- function(Y, Z, R1, Attempt, R2) {
  require_in_each_arm(R1 == 1, Z, "initial-sample respondents (R1 == 1)")
  require_in_each_arm(Attempt == 1, Z, "follow-up attempts (Attempt == 1)")
  require_in_each_arm(R2 == 1, Z, "follow-up respondents (R2 == 1)")
  arm <- function(z) {
    in_arm <- Z == z
    initial <- in_arm & R1 == 1
    followed <- in_arm & R2 == 1
    list(n1 = sum(in_arm),
         p1 = mean(R1[in_arm]),
         y1m = mean(Y[initial]),
         s1 = sd(Y[initial]),
         n2 = sum(in_arm & Attempt == 1),
         p2 = sum(followed)/sum(in_arm & Attempt == 1),
         y2m = mean(Y[followed]),
         s2 = sd(Y[followed]))
  }
  list(t = arm(1), c = arm(0))
}

# Mean outcome in one arm when a share 1 - p of it is unobserved. A fraction
# delta of that unobserved share is filled with the extreme value (minY for a
# lower bound, maxY for an upper bound) and the remaining 1 - delta with the
# observed mean. delta = 1 is the worst case; delta = 0 is ignorability.
gen_mean <- function(y_m, p, delta, lower_bound = TRUE, minY, maxY){
  const <- if (lower_bound) minY else maxY
  p*y_m + (1 - p)*delta*const + (1 - p)*(1 - delta)*y_m
}

# Variance of the outcome under the same imputation: a mixture of the observed
# distribution, with weight p + (1 - p)(1 - delta), and a point mass at the
# extreme value.
gen_var <- function(y_m, y_s, p, delta, lower_bound = TRUE, minY, maxY) {
  const <- if (lower_bound) minY else maxY
  mixture_weight <- p + (1 - p)*(1 - delta)
  mixture_weight*y_s^2 + mixture_weight*(1 - mixture_weight)*(y_m - const)^2
}

construct_manski_bounds <-
  function(p1_t, y1m_t,
           p1_c, y1m_c,
           y2m_t_L, y2m_t_U,
           y2m_c_L, y2m_c_U){
    lower_bound <- ((p1_t*y1m_t + (1-p1_t)*y2m_t_L)) - (p1_c*y1m_c + (1-p1_c)*y2m_c_U)
    upper_bound <- ((p1_t*y1m_t + (1-p1_t)*y2m_t_U)) - (p1_c*y1m_c + (1-p1_c)*y2m_c_L)
    return(c(lower_bound, upper_bound))
  }

ds_var <- function(n1,n2,p1,p2,s1,s2,y1m,y2m) {
  return_value <-
    p1*n1/n1^2 * s1^2 +
    ((1-p1)*n1)^2/(n2*n1^2)*s2^2 +
    ((1-p1)*n1)*p1*n1/n1^3*(y2m-y1m)^2
  return(return_value)
}

# Poststratification ----

# The estimator is run inside each stratum and the results combined by stratum
# share. `fit` takes a subset of `df` and returns the six-element vector the
# unstratified estimator returns.
poststratify <- function(df, strata, alpha, fit) {
  if (anyNA(strata)) stop("The stratification variable (strata) must not contain any missing values.")
  unique_strata <- unique(strata)
  strata_ests <- vapply(unique_strata, function(s) fit(df[strata == s, , drop = FALSE]), numeric(6))
  proportions <- vapply(unique_strata, function(s) mean(strata == s), numeric(1))
  pool_strata(strata_ests, proportions, alpha)
}

# Poststratified bounds: stratum bounds combined by stratum share, stratum
# variances by squared share, with a joint Imbens-Manski interval computed on
# the pooled quantities. Stratum shares are treated as fixed, not estimated.
# strata_ests is a 6-row matrix, one column per stratum, as returned by the
# unstratified estimators. Those carry standard errors, so they are squared back
# to variances here: it is variances that combine by squared stratum share.
pool_strata <- function(strata_ests, proportions, alpha) {
  im_interval(lower_est = sum(strata_ests["estimate_lower", ] * proportions),
              upper_est = sum(strata_ests["estimate_upper", ] * proportions),
              lower_var = sum(strata_ests["std.error_lower", ]^2 * proportions^2),
              upper_var = sum(strata_ests["std.error_upper", ]^2 * proportions^2),
              alpha = alpha)
}

# Imbens-Manski interval ----

# The six quantities every estimator reports: the two bound estimates, their
# standard errors, and the joint Imbens-Manski interval around the region.
im_interval <- function(lower_est, upper_est, lower_var, upper_var, alpha) {
  sig <- im_critical_value(lower_est, upper_est, lower_var, upper_var, alpha)
  c(estimate_lower = lower_est,
    estimate_upper = upper_est,
    std.error_lower = lower_var^0.5,
    std.error_upper = upper_var^0.5,
    conf.low = lower_est - sig*lower_var^0.5,
    conf.high = upper_est + sig*upper_var^0.5)
}

# Coverage of the Imbens-Manski interval at critical value ca, in excess of the
# 1 - alpha target: zero at the solution, negative below it, positive above.
im_crit <- function(ca, upper_bound_est, lower_bound_est, upper_bound_var_est, lower_bound_var_est, alpha) {
  return_value <-
    pnorm(ca + (upper_bound_est-lower_bound_est)/sqrt(max(upper_bound_var_est,lower_bound_var_est)))-pnorm(-ca)-(1-alpha)
  return(return_value)
}

# The Imbens-Manski critical value solves im_crit(ca) == 0. Its root always lies
# between z_(1-alpha/2), the value when the bounds coincide, and z_(1-alpha), the
# limit as the bounds separate, so the bracket is derived from alpha rather than
# fixed. im_crit is increasing in ca, so uniroot solves it to machine precision.
im_critical_value <- function(lower_bound_est, upper_bound_est,
                              lower_bound_var_est, upper_bound_var_est, alpha) {
  z_lower <- qnorm(1 - alpha)
  z_upper <- qnorm(1 - alpha/2)
  excess <- function(ca) im_crit(ca, upper_bound_est, lower_bound_est,
                                 upper_bound_var_est, lower_bound_var_est, alpha)
  if(is.na(excess(z_lower))){return(NA_real_)}
  if(excess(z_lower) >= 0){return(z_lower)}
  if(excess(z_upper) <= 0){return(z_upper)}
  uniroot(excess, lower = z_lower, upper = z_upper, tol = .Machine$double.eps^0.5)$root
}
