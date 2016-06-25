
find_sign_changes <- function(x){
  first_pos <- Position(function(x) x > 0, x )
  first_neg <- Position(function(x) x < 0, x )
  position_zero <- Position(function(x) x == 0,x)
  all_changes <- c(first_pos, first_neg, position_zero)
  all_changes <- all_changes[!all_changes %in% c(1, length(x))]

  return(1:length(x) %in% all_changes)

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

im_crit <- function(ca, upper_bound_est, lower_bound_est, upper_bound_var_est, lower_bound_var_est, alpha) {
  return_value <-
    abs(pnorm(ca + (upper_bound_est-lower_bound_est)/sqrt(max(upper_bound_var_est,lower_bound_var_est)))-pnorm(-ca)-(1-alpha))
  return(return_value)
}
