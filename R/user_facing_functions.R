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
  strata <-  eval(substitute(strata), data)
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

  est_l <-  sum(m1_l_vec *proportions)
  est_u <-  sum(m1_u_vec *proportions)

  var_l <- sum(v1_l_vec *proportions^2)
  var_u <- sum(v1_u_vec *proportions^2)

  im_crit <- function(ca) abs(pnorm(ca + (est_u-est_l)/max(var_u,var_l))-pnorm(-ca)-(1-alpha))

  sig <- optim(1.60,im_crit,method="Brent",lower=1,upper=2)$par

  return(c(ci_lower=est_l - sig*var_l^.5,
           ci_upper=est_u + sig*var_u^.5,
           low_est=est_l,
           upp_est=est_u,
           low_var=var_l,
           upp_var=var_u))
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
estimator_ev <- function(Y, Z, R,  minY, maxY, alpha, data){
  Y <-  eval(substitute(Y), data)
  if(!is.numeric(Y)){stop("The outcome variable (Y) must be numeric.")}
  Z <-  eval(substitute(Z), data)
  if(!all(Z %in% c(0,1))){stop("The treatment variable (Z) must be numeric and take values zero or one.")}
  R <-  eval(substitute(R), data)
  if(!all(R %in% c(0,1))){stop("Theresponse variable (R) must be numeric and take values zero or one.")}

  if(!is.numeric(minY) | !is.numeric(maxY)){stop("The minimum and maximum possible values of Y (minY and maxY) must be numeric")}

  DV_l <- DV_u <- Y
  DV_l[R==0] <- minY
  DV_u[R==0] <- maxY

  est_u <- mean(DV_u[Z==1]) - mean(DV_l[Z==0])
  est_l <- mean(DV_l[Z==1]) - mean(DV_u[Z==0])

  var_u <- var(DV_u[Z==1])/sum(Z==1) + var(DV_l[Z==0])/sum(Z==0)
  var_l <- var(DV_l[Z==1])/sum(Z==1) + var(DV_u[Z==0])/sum(Z==0)

  im_crit <- function(ca) abs(pnorm(ca + (est_u-est_l)/max(var_u,var_l))-pnorm(-ca)-(1-alpha))
  sig <- optim(1.60,im_crit,method="Brent",lower=1,upper=2)$par

  return(c(ci_lower=est_l - sig*var_l^.5,
           ci_upper=est_u + sig*var_u^.5,
           low_est=est_l,
           upp_est=est_u,
           low_var=var_l,
           upp_var=var_u))
}


