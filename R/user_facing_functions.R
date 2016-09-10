#' Extreme Value Bounds with Double Sampling
#'
#' @param Y The (unquoted) outcome variable. Must be numeric.
#' @param Z The (unquoted) assignment indicator variable. Must be numeric and take values 0 or 1.
#' @param R1 The (unquoted) initial sample respose indicator variable. Must be numeric and take values 0 or 1.
#' @param Attempt The (unquoted) follow-up sample attempt indicator variable. Must be numeric and take values 0 or 1.
#' @param R2 The (unquoted) follow-up sample respose indicator variable. Must be numeric and take values 0 or 1.
#' @param minY The minimum possible value of the outcome (Y) variable.
#' @param maxY The maximum possible value of the outcome (Y) variable.
#' @param strata A single (unquoted) variable that indicates which strata units are in.
#' @param alpha The desired significance level. 0.05 by default.
#' @param data A dataframe
#'
#' @return A results matrix
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
  Y <-  eval(substitute(Y), data)
  if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
  Z <-  eval(substitute(Z), data)
  if(!all(Z %in% c(0,1))){stop("The treatment variable (Z) must be numeric and take values zero or one.")}
  R1 <-  eval(substitute(R1), data)
  if(!all(R1 %in% c(0,1))){stop("The initial sample response variable (R1) must be numeric and take values zero or one.")}
  R2 <-  eval(substitute(R2), data)
  if(!all(R2 %in% c(0,1))){stop("The follow-up sample response variable (R2) must be numeric and take values zero or one.")}
  Attempt <-  eval(substitute(Attempt), data)
  if(!all(Attempt %in% c(0,1))){stop("The follow-up sample attempt variable (Attempt) must be numeric and take values zero or one.")}

  if(!is.numeric(minY) | !is.numeric(maxY)){stop("The minimum and maximum possible values of Y (minY and maxY) must be numeric")}

  strata <-  eval(substitute(strata), data)
  if(is.null(strata)) {
    n1_c_s <- sum(R1==1 & Z==0)
    n1_t_s <- sum(R1==1 & Z==1)
    n1_c <- sum(Z==0) # clean up
    n1_t <- sum(Z==1)

    p1_c <- n1_c_s/n1_c
    p1_t <- n1_t_s/n1_t

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

    # c1a_t <- c1r_t <- c2a_t <- c2r_t <- -99 # irrelevant

    cis_out <- ds_manski_cis_2s(n1_t=n1_t,n2_t=n2_t,
                                n1_c=n1_c,n2_c=n2_c,
                                p1_t=p1_t,p2_t=p2_t,
                                s1_t=s1_t,s1_c=s1_c,
                                s2_nm_t=s2_nm_t,
                                s2_nm_c=s2_nm_c,
                                y1m_t=y1m_t,y1m_c=y1m_c,
                                y2m_nm_t=y2m_nm_t,
                                y2m_nm_c=y2m_nm_c,
                                c1a_t=c1a_t,c1r_t=c1r_t,c2a_t=c2a_t,c2r_t=c2r_t,
                                p1_c=p1_c,p2_c=p2_c,
                                minY=minY,maxY=maxY,alpha=alpha)
    return(cis_out)
  }else{
    # If there is a stratification variable, call estimator_ds recursively.

    if(sum(is.na(strata))!=0){stop("The stratification variable (strata) must not contain any missing values.")}

    unique_strata <- unique(strata)
    n_strata <- length(unique_strata)
    m1_l_vec <- m1_u_vec <- v1_l_vec <- v1_u_vec <- proportions <- rep(NA, n_strata)

    ds_df <- data.frame(Y, R1, Z, Attempt, R2, strata)

    for(i in 1:n_strata){
      ests <- estimator_ds(Y = Y, R1 = R1,Z = Z,
                           Attempt = Attempt, R2 = R2,
                           minY = minY, maxY = maxY, alpha = alpha,
                           data=subset(ds_df, strata==unique_strata[i]))
      m1_l_vec[i] <- ests[3]
      m1_u_vec[i] <- ests[4]
      v1_l_vec[i] <- ests[5]
      v1_u_vec[i] <- ests[6]
      proportions[i] <- mean(strata == unique_strata[i])
    }

    lower_bound_est <-  sum(m1_l_vec *proportions)
    upper_bound_est <-  sum(m1_u_vec *proportions)

    lower_bound_var_est <- sum(v1_l_vec *proportions^2)
    upper_bound_var_est <- sum(v1_u_vec *proportions^2)

    sig <- optim(1.60,im_crit,method="Brent",lower=1,upper=2,
                 lower_bound_est = lower_bound_est,
                 upper_bound_est = upper_bound_est,
                 lower_bound_var_est = lower_bound_var_est,
                 upper_bound_var_est = upper_bound_var_est,
                 alpha = alpha)$par

    return(c(ci_lower=lower_bound_est - sig*lower_bound_var_est^.5,
             ci_upper=upper_bound_est + sig*upper_bound_var_est^.5,
             low_est=lower_bound_est,
             upp_est=upper_bound_est,
             low_var=lower_bound_var_est,
             upp_var=upper_bound_var_est))
  }
}


#' Extreme Value (Manski) Bounds
#'
#' @param Y The (unquoted) outcome variable. Must be numeric.
#' @param Z The (unquoted) assignment indicator variable. Must be numeric and take values 0 or 1.
#' @param R The (unquoted) respose indicator variable. Must be numeric and take values 0 or 1.
#' @param minY The minimum possible value of the outcome (Y) variable.
#' @param maxY The maximum possible value of the outcome (Y) variable.
#' @param strata A single (unquoted) variable that indicates which strata units are in.
#' @param alpha The desired significance level. 0.05 by default.
#' @param data A dataframe
#'
#' @return A results matrix
#' @export
estimator_ev <- function(Y, Z, R, minY, maxY, strata = NULL, alpha = 0.05, data){
  Y <-  eval(substitute(Y), data)
  if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
  Z <-  eval(substitute(Z), data)
  if(!all(Z %in% c(0,1))){stop("The treatment variable (Z) must be numeric and take values zero or one.")}
  R <-  eval(substitute(R), data)
  if(!all(R %in% c(0,1))){stop("The response variable (R) must be numeric and take values zero or one.")}

  if(!is.numeric(minY) | !is.numeric(maxY)){stop("The minimum and maximum possible values of Y (minY and maxY) must be numeric")}

  strata <-  eval(substitute(strata), data)
  if(is.null(strata)) {
    n1_c_s <- sum(R==1 & Z==0)
    n1_t_s <- sum(R==1 & Z==1)
    n1_c <- sum(Z==0) # clean up
    n1_t <- sum(Z==1)

    p1_c <- n1_c_s/n1_c
    p1_t <- n1_t_s/n1_t

    y1m_c <- mean(Y[R==1 & Z==0])
    y1m_t <- mean(Y[R==1 & Z==1])

    s1_c <- sd(Y[R==1 & Z==0])
    s1_t <- sd(Y[R==1 & Z==1])


    cis_out <- manski_cis(n1_t = n1t, n1_c = n1_c,
                          n1_t_s = n1_t_s, n1_c_s = n1_c_c,
                          p1_t = p1_t, p1_c = p1_c,
                          y1m_t = y1m_t, y1m_c = y1m_c,
                          s1_t = s1_t, s1_c = s1_c,
                          minY = minY, maxY = maxY, alpha = alpha)

    return(cis_out)
  }else{
    # If there is a stratification variable, call estimator_ds recursively.

    if(sum(is.na(strata))!=0){stop("The stratification variable (strata) must not contain any missing values.")}

    unique_strata <- unique(strata)
    n_strata <- length(unique_strata)
    m1_l_vec <- m1_u_vec <- v1_l_vec <- v1_u_vec <- proportions <- rep(NA, n_strata)

    ds_df <- data.frame(Y, Z, R, strata)

    for(i in 1:n_strata){
      ests <- estimator_ev(Y = Y, Z = Z, R = R,
                           minY = minY, maxY = maxY, alpha = alpha,
                           data=subset(ds_df, strata==unique_strata[i]))
      m1_l_vec[i] <- ests[3]
      m1_u_vec[i] <- ests[4]
      v1_l_vec[i] <- ests[5]
      v1_u_vec[i] <- ests[6]
      proportions[i] <- mean(strata == unique_strata[i])
    }

    lower_bound_est <-  sum(m1_l_vec *proportions)
    upper_bound_est <-  sum(m1_u_vec *proportions)

    lower_bound_var_est <- sum(v1_l_vec *proportions^2)
    upper_bound_var_est <- sum(v1_u_vec *proportions^2)

    sig <- optim(1.60,im_crit,method="Brent",lower=1,upper=2,
                 lower_bound_est = lower_bound_est,
                 upper_bound_est = upper_bound_est,
                 lower_bound_var_est = lower_bound_var_est,
                 upper_bound_var_est = upper_bound_var_est,
                 alpha = alpha)$par

    return(c(ci_lower=lower_bound_est - sig*lower_bound_var_est^.5,
             ci_upper=upper_bound_est + sig*upper_bound_var_est^.5,
             low_est=lower_bound_est,
             upp_est=upper_bound_est,
             low_var=lower_bound_var_est,
             upp_var=upper_bound_var_est))
  }

}

#' Trimming Bounds
#'
#' @param Y The (unquoted) outcome variable. Must be numeric.
#' @param Z The (unquoted) assignment indicator variable. Must be numeric and take values 0 or 1.
#' @param R The (unquoted) respose indicator variable. Must be numeric and take values 0 or 1.
#' @param R1 The (unquoted) initial sample respose indicator variable. Must be numeric and take values 0 or 1.
#' @param Attempt The (unquoted) follow-up sample attempt indicator variable. Must be numeric and take values 0 or 1.
#' @param R2 The (unquoted) follow-up sample respose indicator variable. Must be numeric and take values 0 or 1.
#' @param strata A single (unquoted) variable that indicates which strata units are in. DOES NOT WORK YET
#' @param alpha The desired significance level. 0.05 by default.
#' @param data A dataframe
#'
#' @return A results matrix
#' @export
estimator_trim <-
  function(Y, Z, R = NULL, R1 = NULL, Attempt = NULL, R2 = NULL, strata = NULL, alpha = 0.05, data){
    Y <-  eval(substitute(Y), data)
    if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
    Z <-  eval(substitute(Z), data)
    if(!all(Z %in% c(0,1))){stop("The treatment variable (Z) must be numeric and take values zero or one.")}

    R <-  eval(substitute(R), data)
    if(!is.null(R)){
      if(!all(R %in% c(0,1))){stop("The response variable (R) must be numeric and take values zero or one.")}

      out <- trimming_bounds(Out = Y, Treat = Z, Fail = as.numeric(R == 0), Weight = rep(1, length(Y)), monotonicity = TRUE)
    }else{
      R1 <-  eval(substitute(R1), data)
      if(!all(R1 %in% c(0,1))){stop("The initial sample response variable (R1) must be numeric and take values zero or one.")}
      R2 <-  eval(substitute(R2), data)
      if(!all(R2 %in% c(0,1))){stop("The follow-up sample response variable (R2) must be numeric and take values zero or one.")}
      Attempt <-  eval(substitute(Attempt), data)
      if(!all(Attempt %in% c(0,1))){stop("The follow-up sample attempt variable (Attempt) must be numeric and take values zero or one.")}

      Weight <- rep(NA, length(Y))

      Weight[R1==1] <- 1
      Weight[Attempt==1 & Z ==1] <- sum(Z == 1 & R1 == 0)/sum(Z == 1 & Attempt == 1)
      Weight[Attempt==1 & Z ==0] <- sum(Z == 0 & R1 == 0)/sum(Z == 0 & Attempt == 1)

      Fail <- as.numeric(R1 == 0 & R2 == 0)

      Keep <- (R1 == 1 | Attempt == 1)

      out <- trimming_bounds(Out = Y[Keep], Treat = Z[Keep],
                             Fail = Fail[Keep], Weight = Weight[Keep], monotonicity = FALSE)
    }

    return(out)
  }


#' Extreme Value Bounds with Double Sampling with Sensitivity
#'
#' This function yields extreme value bounds under the assumption that the outcomes of 1-delta of the missing second-round units are ignorable, that is, that they are drawn from an unknown distribution with mean and variance equal to the observed second-round groups.
#'
#' @param Y The (unquoted) outcome variable. Must be numeric.
#' @param Z The (unquoted) assignment indicator variable. Must be numeric and take values 0 or 1.
#' @param R1 The (unquoted) initial sample respose indicator variable. Must be numeric and take values 0 or 1.
#' @param Attempt The (unquoted) follow-up sample attempt indicator variable. Must be numeric and take values 0 or 1.
#' @param R2 The (unquoted) follow-up sample respose indicator variable. Must be numeric and take values 0 or 1.
#' @param minY The minimum possible value of the outcome (Y) variable.
#' @param maxY The maximum possible value of the outcome (Y) variable.
#' @param strata A single (unquoted) variable that indicates which strata units are in.
#' @param alpha The desired significance level. 0.05 by default.
#' @param data A dataframe
#' @param sims Number of points at which to evaluate sensitivity test. Defaults to 100
#'
#' @return A list containing a ggplot object, a dataframe of simulated bounds and cis, and a value of pstar, if it exists.
#' @export
#'
estimator_ds_sens <- function(Y, Z, R1, Attempt, R2, minY, maxY, delta, strata = NULL, alpha = 0.05, data){
  Y <-  eval(substitute(Y), data)
  if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
  Z <-  eval(substitute(Z), data)
  if(!all(Z %in% c(0,1))){stop("The treatment variable (Z) must be numeric and take values zero or one.")}
  R1 <-  eval(substitute(R1), data)
  if(!all(R1 %in% c(0,1))){stop("The initial sample response variable (R1) must be numeric and take values zero or one.")}
  R2 <-  eval(substitute(R2), data)
  if(!all(R2 %in% c(0,1))){stop("The follow-up sample response variable (R2) must be numeric and take values zero or one.")}
  Attempt <-  eval(substitute(Attempt), data)
  if(!all(Attempt %in% c(0,1))){stop("The follow-up sample attempt variable (Attempt) must be numeric and take values zero or one.")}

  if(!is.numeric(minY) | !is.numeric(maxY)){stop("The minimum and maximum possible values of Y (minY and maxY) must be numeric")}

  strata <-  eval(substitute(strata), data)
  if(is.null(strata)) {
    n1_c_s <- sum(R1==1 & Z==0)
    n1_t_s <- sum(R1==1 & Z==1)
    n1_c <- sum(Z==0) # clean up
    n1_t <- sum(Z==1)

    p1_c <- n1_c_s/n1_c
    p1_t <- n1_t_s/n1_t

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

    # c1a_t <- c1r_t <- c2a_t <- c2r_t <- -99 # irrelevant

    cis_out <- ds_manski_cis_2s_sens(n1_t=n1_t,n2_t=n2_t,
                                     n1_c=n1_c,n2_c=n2_c,
                                     p1_t=p1_t,p2_t=p2_t,
                                     s1_t=s1_t,s1_c=s1_c,
                                     s2_nm_t=s2_nm_t,
                                     s2_nm_c=s2_nm_c,
                                     y1m_t=y1m_t,y1m_c=y1m_c,
                                     y2m_nm_t=y2m_nm_t,
                                     y2m_nm_c=y2m_nm_c,
                                     c1a_t=c1a_t,c1r_t=c1r_t,c2a_t=c2a_t,c2r_t=c2r_t,
                                     p1_c=p1_c,p2_c=p2_c,
                                     minY=minY,maxY=maxY,alpha=alpha, delta = delta)
    return(cis_out)
  }else{
    # If there is a stratification variable, call estimator_ds recursively.

    if(sum(is.na(strata))!=0){stop("The stratification variable (strata) must not contain any missing values.")}

    unique_strata <- unique(strata)
    n_strata <- length(unique_strata)
    m1_l_vec <- m1_u_vec <- v1_l_vec <- v1_u_vec <- proportions <- rep(NA, n_strata)

    ds_df <- data.frame(Y, R1, Z, Attempt, R2, strata)

    for(i in 1:n_strata){
      ests <- estimator_ds_sens(Y = Y, R1 = R1,Z = Z,
                                Attempt = Attempt, R2 = R2,
                                minY = minY, maxY = maxY, alpha = alpha,
                                data=subset(ds_df, strata==unique_strata[i]), delta = delta)
      m1_l_vec[i] <- ests[3]
      m1_u_vec[i] <- ests[4]
      v1_l_vec[i] <- ests[5]
      v1_u_vec[i] <- ests[6]
      proportions[i] <- mean(strata == unique_strata[i])
    }

    lower_bound_est <-  sum(m1_l_vec *proportions)
    upper_bound_est <-  sum(m1_u_vec *proportions)

    lower_bound_var_est <- sum(v1_l_vec *proportions^2)
    upper_bound_var_est <- sum(v1_u_vec *proportions^2)

    sig <- optim(1.60,im_crit,method="Brent",lower=1,upper=2,
                 lower_bound_est = lower_bound_est,
                 upper_bound_est = upper_bound_est,
                 lower_bound_var_est = lower_bound_var_est,
                 upper_bound_var_est = upper_bound_var_est,
                 alpha = alpha)$par

    return(c(ci_lower=lower_bound_est - sig*lower_bound_var_est^.5,
             ci_upper=upper_bound_est + sig*upper_bound_var_est^.5,
             low_est=lower_bound_est,
             upp_est=upper_bound_est,
             low_var=lower_bound_var_est,
             upp_var=upper_bound_var_est))
  }
}


#' Sensitivity Analysis
#'
#' This function performs a line search over values of delta, the sensitivity parameter, in order to find (if it exists) delta*, the value of delta where the confidence interval no longer includes zero.
#'
#' @param Y The (unquoted) outcome variable. Must be numeric.
#' @param Z The (unquoted) assignment indicator variable. Must be numeric and take values 0 or 1.
#' @param R1 The (unquoted) initial sample respose indicator variable. Must be numeric and take values 0 or 1.
#' @param Attempt The (unquoted) follow-up sample attempt indicator variable. Must be numeric and take values 0 or 1.
#' @param R2 The (unquoted) follow-up sample respose indicator variable. Must be numeric and take values 0 or 1.
#' @param minY The minimum possible value of the outcome (Y) variable.
#' @param maxY The maximum possible value of the outcome (Y) variable.
#' @param strata A single (unquoted) variable that indicates which strata units are in.
#' @param alpha The desired significance level. 0.05 by default.
#' @param data A dataframe
#' @param sims Number of points at which to evaluate sensitivity test. Defaults to 100
#'
#' @return A list containing a ggplot object, a dataframe of simulated bounds and cis, and a value of pstar, if it exists.
#' @export
#'
sensitivity_ds <- function(Y, Z, R1, Attempt, R2, minY, maxY, sims = 100, strata = NULL, alpha = 0.05, data){
  require(ggplot2)
  require(dplyr)
  require(purrr)
  require(reshape2)

  Y <-  eval(substitute(Y), data)
  if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
  Z <-  eval(substitute(Z), data)
  if(!all(Z %in% c(0,1))){stop("The treatment variable (Z) must be numeric and take values zero or one.")}
  R1 <-  eval(substitute(R1), data)
  if(!all(R1 %in% c(0,1))){stop("The initial sample response variable (R1) must be numeric and take values zero or one.")}
  R2 <-  eval(substitute(R2), data)
  if(!all(R2 %in% c(0,1))){stop("The follow-up sample response variable (R2) must be numeric and take values zero or one.")}
  Attempt <-  eval(substitute(Attempt), data)
  if(!all(Attempt %in% c(0,1))){stop("The follow-up sample attempt variable (Attempt) must be numeric and take values zero or one.")}

  strata <-  eval(substitute(strata), data)

  if(is.null(strata)){
  df <- data.frame(Y, Z, R1, R2, Attempt)
  ps <- seq(0, 1, length.out = sims)

  sims_df <-
    map(ps, ~estimator_ds_sens(Y = Y, Z = Z, R1 = R1, Attempt = Attempt,alpha = alpha,
                               R2 = R2, minY=minY, maxY=maxY, data=df, delta = .x)) %>%
    do.call(rbind, .) %>%
    data.frame() %>%
    mutate(p = ps,
           change_lower = find_sign_changes(ci_lower),
           change_upper = find_sign_changes(ci_upper),
           change_any = change_lower | change_upper)
  }else{
    df <- data.frame(Y, Z, R1, R2, Attempt, strata)
    ps <- seq(0, 1, length.out = sims)

    sims_df <-
      map(ps, ~estimator_ds_sens(Y = Y, Z = Z, R1 = R1, Attempt = Attempt, strata = strata,alpha = alpha,
                                 R2 = R2, minY=minY, maxY=maxY, data=df, delta = .x)) %>%
      do.call(rbind, .) %>%
      data.frame() %>%
      mutate(p = ps,
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
    ylab("Identification Regions and 95% Confidence Intervals") +
    xlab(expression(paste("Sensitivity Parameter ", delta, " (0 = Ignorability)"))) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.key.width = unit(3, "lines"),
          legend.title = element_blank())


  p_star_df <- "No value of the sensitivity parameter yields a statistically significant result."

  p_star <- with(sims_df, p[change_any])
  if(length(p_star) == 1){
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
