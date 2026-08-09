
# A column argument is either a vector or a single string naming a column of
# data. Anything else (a vector of column contents) is passed through as is.
resolve_column <- function(x, data) {
  if (is.character(x) && length(x) == 1L) data[[x]] else x
}

# Resolve the outcome and treatment arguments. Y is either an outcome vector or
# a Y ~ Z formula; when it is a formula, Z is taken from the formula and the Z
# argument is never evaluated, so callers may leave it missing.
resolve_yz <- function(Y_expr, Z_expr, data, env) {
  Y_val <- eval(Y_expr, data, env)
  if (inherits(Y_val, "formula")) return(parse_yz_formula(Y_val, data))
  list(Y = Y_val, Z = eval(Z_expr, data, env))
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

# Flag the first element whose sign differs from the first element's, treating a
# zero as a change. Position 1 can never be flagged; the last position can.
find_sign_changes <- function(x){
  first_change <- Position(function(xi) sign(xi) != sign(x[1]), x)
  out <- rep(FALSE, length(x))
  if(!is.na(first_change)){out[first_change] <- TRUE}
  return(out)
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

gen_mean <- function(y_m,p,lower_bound=TRUE,minY,maxY){
  if (lower_bound == TRUE){
    return(p*y_m + (1-p)*minY)
  }else{
    return(p*y_m+(1-p)*maxY)
  }
}

gen_var <- function(y_m , y_s, p, lower_bound=TRUE,minY,maxY) {
  if (lower_bound==TRUE){
    const <- minY
  }else{
    const <- maxY
  }
  wm <- gen_mean(y_m,p,lower_bound,minY,maxY)
  # formula for combined var
  return(p*y_s^2 + p*(y_m-wm)^2 + (1-p)*(const-wm)^2)
}

gen_mean_sens <- function(y_m, p, delta, lower_bound = TRUE, minY, maxY){
  if (lower_bound==TRUE){
    const <- minY
  }else{
    const <- maxY
  }
  return(p*y_m + (1-p)*delta*const + (1-p)*(1-delta)*y_m)
}

gen_var_sens <- function(y_m, y_s, p, delta, lower_bound = TRUE, minY, maxY) {

  if(lower_bound == TRUE){
    const <- minY
  }else{
    const <- maxY
  }
  mixture_weight <- p + (1-p)*(1-delta)

  var_sens <-
    mixture_weight*y_s^2 +
    mixture_weight*(1-mixture_weight)*(y_m - const)^2

  return(var_sens)
}

# Poststratified bounds: stratum bounds combined by stratum share, stratum
# variances by squared share, with a joint Imbens-Manski interval computed on
# the pooled quantities. Stratum shares are treated as fixed, not estimated.
# strata_ests is a 6-row matrix, one column per stratum, as returned by the
# unstratified estimators.
pool_strata <- function(strata_ests, proportions, alpha) {
  lower_bound_est <- sum(strata_ests["low_est", ] * proportions)
  upper_bound_est <- sum(strata_ests["upp_est", ] * proportions)
  lower_bound_var_est <- sum(strata_ests["low_var", ] * proportions^2)
  upper_bound_var_est <- sum(strata_ests["upp_var", ] * proportions^2)

  sig <- im_critical_value(lower_bound_est, upper_bound_est,
                           lower_bound_var_est, upper_bound_var_est, alpha)

  return(c(ci_lower = lower_bound_est - sig*lower_bound_var_est^.5,
           ci_upper = upper_bound_est + sig*upper_bound_var_est^.5,
           low_est = lower_bound_est,
           upp_est = upper_bound_est,
           low_var = lower_bound_var_est,
           upp_var = upper_bound_var_est))
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

ds_var_2s <- function(treatment_vec,control_vec) {
  ts_var <-
    ds_var(treatment_vec[1],treatment_vec[2],
           treatment_vec[3],treatment_vec[4],
           treatment_vec[5],treatment_vec[6],
           treatment_vec[7],treatment_vec[8]) +
    ds_var(control_vec[1],control_vec[2],
           control_vec[3],control_vec[4],
           control_vec[5],control_vec[6],
           control_vec[7],control_vec[8])
  return(ts_var)
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
