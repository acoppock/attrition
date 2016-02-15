# calculate sampling variance, p2 is assumed to be 1 always in this equation
#' @export
ds_var <- function(n1,n2,p1,p2,s1,s2,y1m,y2m) p1*n1/n1^2 * s1^2 + ((1-p1)*n1)^2/(n2*n1^2)*s2^2 + ((1-p1)*n1)*p1*n1/n1^3*(y2m-y1m)^2

# calculate manski bounds
#' @export
ds_manski <- function(p1,p2,y1m,y2m_nm,minY,maxY) {
  const1 <- p1*y1m + (1-p1)*p2*y2m_nm
  const2 <- (1-p1)*(1-p2)
  return(c(const1+const2*minY,const1+const2*maxY))
}

# calculate two-sample sampling variance
#' @export
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

# estimate manski bounds, two-sample:
#' @export
ds_manski_2s <- function(p1_t,p2_t,y1m_t,y2m_t,p1_c,p2_c,y1m_c,y2m_c,minY,maxY) {
  const1_t <- p1_t*y1m_t + (1-p1_t)*p2_t*y2m_t
  const2_t <- (1-p1_t)*(1-p2_t)
  const1_c <- p1_c*y1m_c + (1-p1_c)*p2_c*y2m_c
  const2_c <- (1-p1_c)*(1-p2_c)
  return(c(const1_t+const2_t*minY - (const1_c+const2_c*maxY),const1_t+const2_t*maxY - (const1_c+const2_c*minY)))
}

#' @export
ds_manski_cis_2s <- function(n1_t,n2_t,n1_c,n2_c,
                             p1_t,p2_t,s1_t,s2_nm_t,
                             y1m_t,y2m_nm_t,
                             c1a_t,c1r_t,c2a_t,c2r_t,
                             p1_c,p2_c,
                             s1_c,s2_nm_c,
                             y1m_c,y2m_nm_c,
                             minY,maxY,alpha) {

  # compute 2nd round manski mean
  gen_mean <- function(y2m_nm,p2,lower_bound=TRUE,minY,maxY){
    if (lower_bound == TRUE){
      return(p2*y2m_nm + (1-p2)*minY)
    }else{
      return(p2*y2m_nm+(1-p2)*maxY)
    }
  }
  # compute 2nd round manski sd given prop, which bound, and mean
  gen_var <- function(y2m_nm,s2_nm,p2,lower_bound=TRUE,minY,maxY) {
    if (lower_bound==TRUE){
      const <- minY
    }else{
      const <- maxY
    }
    wm <- gen_mean(y2m_nm,p2,lower_bound,minY,maxY)
    # formula for combined var
    return(p2*s2_nm^2 + p2*(y2m_nm-wm)^2 + (1-p2)*(const-wm)^2)

  }

  y2m_t_L <- gen_mean(y2m_nm_t,p2_t,lower_bound=TRUE,minY,maxY)
  y2m_t_U <- gen_mean(y2m_nm_t,p2_t,lower_bound=FALSE,minY,maxY)
  y2m_c_L <- gen_mean(y2m_nm_c,p2_c,lower_bound=TRUE,minY,maxY)
  y2m_c_U <- gen_mean(y2m_nm_c,p2_c,lower_bound=FALSE,minY,maxY)

  s2_t_L <- gen_var(y2m_nm_t,s2_nm_t,p2_t,lower_bound=TRUE,minY,maxY)^.5
  s2_t_U <- gen_var(y2m_nm_t,s2_nm_t,p2_t,lower_bound=FALSE,minY,maxY)^.5
  s2_c_L <- gen_var(y2m_nm_c,s2_nm_c,p2_c,lower_bound=TRUE,minY,maxY)^.5
  s2_c_U <- gen_var(y2m_nm_c,s2_nm_c,p2_c,lower_bound=FALSE,minY,maxY)^.5

  manski_bounds_est <- ds_manski_2s(p1_t,p2_t,y1m_t,y2m_nm_t,p1_c,p2_c,y1m_c,y2m_nm_c,minY,maxY)
  lower_bound_est <- ds_manski_2s(p1_t,p2_t,y1m_t,y2m_nm_t,p1_c,p2_c,y1m_c,y2m_nm_c,minY,maxY)[1]
  upper_bound_est <- ds_manski_2s(p1_t,p2_t,y1m_t,y2m_nm_t,p1_c,p2_c,y1m_c,y2m_nm_c,minY,maxY)[2]

  lower_bound_var_est <- ds_var_2s(c(n1_t,n2_t,p1_t,p2_t,s1_t,s2_t_L,y1m_t,y2m_t_L),
                                   c(n1_c,n2_c,p1_c,p2_c,s1_c,s2_c_U,y1m_c,y2m_c_U))
  upper_bound_var_est <- ds_var_2s(c(n1_t,n2_t,p1_t,p2_t,s1_t,s2_t_U,y1m_t,y2m_t_U),
                                   c(n1_c,n2_c,p1_c,p2_c,s1_c,s2_c_L,y1m_c,y2m_c_L))


  im_crit <- function(ca) abs(pnorm(ca + (upper_bound_est-lower_bound_est)/
                                      sqrt(max(upper_bound_var_est,lower_bound_var_est))
                                    )-pnorm(-ca)-(1-alpha))

  sig <- optim(1.60,im_crit,method="Brent",lower=1,upper=2)$par

  return(c(ci_lower=lower_bound_est - sig*lower_bound_var_est^.5,
           ci_upper=upper_bound_est + sig*upper_bound_var_est^.5,
           low_est=lower_bound_est,upp_est=upper_bound_est,
           low_var=lower_bound_var_est,upp_var=upper_bound_var_est, sig = sig))
}

#' find the n1,n2 (number of attempts to measure) that minimize width of Manski CIs
#' @export
optim_manski_2s <- function(p1_t,p2_t,s1_t,s2_nm_t,
                            y1m_t,y2m_nm_t,
                            c1a_t,c1r_t,c2a_t,c2r_t,
                            p1_c,p2_c,s1_c,s2_nm_c,
                            y1m_c,y2m_nm_c,
                            c1a_c,c1r_c,c2a_c,c2r_c,
                            budget,minY,maxY,
                            mu=1e-04,niter=100,alpha=0.05) {


  matout <- matrix(NA,nrow=niter,ncol=5)
  colnames(matout) <- c("n1_t","n2_t","n1_c","n2_c","ci_width")
  iter <- 1

  while(iter <= niter) {

    #(c1a_t+c1r_t*p1_t)*n1_t + (c2a_t+c2r_t*p2_t)*n2_t + (c1a_c+c1r_c*p1_c)*n1_c + (c2a_c+c2r_c*p2_c)*n2_c = budget

    starting_values <- .5*runif(4)*budget/c((c1a_t+c1r_t*p1_t),(c2a_t+c2r_t*p2_t*(1-p1_t)),(c1a_c+c1r_c*p1_c),(c2a_c+c2r_c*p2_c*(1-p1_c)))

    n1g_t = starting_values[1]
    n2g_t = starting_values[2]
    n1g_c = starting_values[3]
    n2g_c = starting_values[4]


    # constraints
    #- (c1a_t+c1r_t*p1_t)*n1_t - (c2a_t+c2r_t*p2_t)*n2_t - (c1a_c+c1r_c*p1_c)*n1_c - (c2a_c+c2r_c*p2_c)*n2_c + budget >= 0
    # n1_t > 2
    # n2_t > 2
    # n1_c > 2
    # n2_c > 2
    # (1-p1_t)*n1_t - n2_t > 2
    # (1-p1_c)*n1_c - n2_c > 2

    constrMat <- rbind(c(-(c1a_t+c1r_t*p1_t),-(c2a_t+c2r_t*p2_t),-(c1a_c+c1r_c*p1_c),-(c2a_c+c2r_c*p2_c)),
                       c(1,0,0,0),
                       c(0,1,0,0),
                       c(0,0,1,0),
                       c(0,0,0,1),
                       c((1-p1_t),-1,0,0),
                       c(0,0,(1-p1_c),-1)
    ) # are these working?
    constrVec <- c(-budget,2,2,2,2,2,2)

    # compute width of 95% CIs
    ds_passthrough_2s <- function(ns,p1_t,p2_t,s1_t,s2_nm_t,
                                  y1m_t,y2m_nm_t,
                                  c1a_t,c1r_t,c2a_t,c2r_t,
                                  p1_c,p2_c,s1_c,s2_nm_c,
                                  y1m_c,y2m_nm_c,
                                  minY,maxY,alpha){
      diff(ds_manski_cis_2s(ns[1],ns[2],ns[3],ns[4],p1_t,p2_t,s1_t,s2_nm_t,y1m_t,y2m_nm_t,c1a_t,c1r_t,c2a_t,c2r_t,p1_c,p2_c,s1_c,s2_nm_c,y1m_c,y2m_nm_c,minY,maxY,alpha))
    }

    optimD <- try(constrOptim(c(n1g_t,n2g_t,n1g_c,n2g_c),
                              ds_passthrough_2s,grad=NULL,ui=constrMat,ci=constrVec,
                              p1_t=p1_t,p2_t=p2_t,
                              s1_t=s1_t,s2_nm_t=s2_nm_t,
                              y1m_t=y1m_t,y2m_nm_t=y2m_nm_t,
                              c1a_t=c1a_t,c1r_t=c1r_t,
                              c2a_t=c2a_t,c2r_t=c2r_t,
                              p1_c=p1_c,p2_c=p2_c,
                              s1_c=s1_c,s2_nm_c=s2_nm_c,
                              y1m_c=y1m_c,y2m_nm_c=y2m_nm_c,
                              minY=minY,maxY=maxY,alpha=alpha,
                              outer.iterations=1000,outer.eps = .Machine$double.eps^.5,mu=mu),silent=TRUE)

    if(!is.character(optimD)) {
      cat(iter,"")
      matout[iter,] <- c(n1_t = optimD$par[1],n2_t = optimD$par[2],n1_c = optimD$par[3],n2_c = optimD$par[4],ci_width=optimD$value)
      iter <- iter + 1
    }

  }
  return(matout[order(matout[,5]),])
}


trimming_bounds <-
  function(Out, Treat, Fail, Weight, monotonicity = FALSE) {

    dataf <- data.frame(Out, Treat, Fail, Weight)
    datafsort <- dataf[order(Out),]

    OutS0 <- datafsort[datafsort$Fail==0 & datafsort$Treat==0,]
    OutS1 <- datafsort[datafsort$Fail==0 & datafsort$Treat==1,]

    OutS0$Weight <- OutS0$Weight/sum(OutS0$Weight)
    OutS1$Weight <- OutS1$Weight/sum(OutS1$Weight)

    OutS0.CDF <- OutS0$Weight
    OutS1.CDF <- OutS1$Weight

    for(i in 2:length(OutS0.CDF)) OutS0.CDF[i] <- OutS0.CDF[i] + OutS0.CDF[i-1]
    for(i in 2:length(OutS1.CDF)) OutS1.CDF[i] <- OutS1.CDF[i] + OutS1.CDF[i-1]

    f0 <- sum(Weight[Fail==1 & Treat==0])/sum(Weight[Treat==0])
    f1 <- sum(Weight[Fail==1 & Treat==1])/sum(Weight[Treat==1])

    if(monotonicity){
      Q <- ((1 - f1) - (1 - f0))/(1-f1)
      if(Q < 0){stop("Monotonicity appears to be violated: The control group is more likely to be missing than the treatment group.")}

      Out0_mono <- weighted.mean(OutS0$Out, OutS0$Weight)

      Out1U_mono <- weighted.mean(OutS1$Out[OutS1.CDF>Q], OutS1$Weight[OutS1.CDF>Q])
      Out1L_mono <- weighted.mean(OutS1$Out[OutS1.CDF<(1-Q)], OutS1$Weight[OutS1.CDF<(1-Q)])

      upper_bound <- Out1U_mono - Out0_mono
      lower_bound <- Out1L_mono - Out0_mono

      return(c(upper_bound = upper_bound, lower_bound = lower_bound,
               Out0_mono = Out0_mono, Out1L_mono=Out1L_mono, Out1U_mono = Out1U_mono,
               control_group_N = nrow(OutS0), treat_group_N = nrow(OutS1), Q = Q, f1 = f1, f0 = f0, pi_r_1 = 1 - f1, pi_r_0 = 1 - f0))

    }else{

      trim0 <- (f1)/(1-f0)
      trim1 <- (f0)/(1-f1)

      Out0U <- weighted.mean(OutS0$Out[OutS0.CDF>trim0], OutS0$Weight[OutS0.CDF>trim0])
      Out0L <- weighted.mean(OutS0$Out[OutS0.CDF<(1-trim0)], OutS0$Weight[OutS0.CDF<(1-trim0)])

      Out1U <- weighted.mean(OutS1$Out[OutS1.CDF>trim1], OutS1$Weight[OutS1.CDF>trim1])
      Out1L <- weighted.mean(OutS1$Out[OutS1.CDF<(1-trim1)], OutS1$Weight[OutS1.CDF<(1-trim1)])

      upper_bound = Out1U - Out0L
      lower_bound = Out1L - Out0U

      return(c(upper_bound = upper_bound, lower_bound = lower_bound, Out0L=Out0L, Out0U=Out0U, Out1L=Out1L, Out1U=Out1U,
               control_group_N = nrow(OutS0), treat_group_N = nrow(OutS1), trim0 = trim0, trim1 = trim1))

    }
  }



#' @export
manski_cis <- function(n1_t, n1_c,
                       n1_t_s, n1_c_s,
                       p1_t, p1_c,
                       y1m_t, y1m_c,
                       s1_t, s1_c,
                       minY,maxY,alpha){

  gen_mean <- function(y1m, p, lower_bound=TRUE, minY, maxY){
    if (lower_bound == TRUE){
      return(p*y1m + (1-p)*minY)
    }else{
      return(p*y1m+ (1-p)*maxY)
    }
  }

  gen_var <- function(y1m, s1, p, lower_bound = TRUE, minY, maxY) {
    if (lower_bound==TRUE){
      const <- minY
    }else{
      const <- maxY
    }
    wm <- gen_mean(y1m,p,lower_bound,minY,maxY)
    # formula for combined var
    return(p*s1^2 + p*(y1m-wm)^2 + (1-p)*(const-wm)^2)
  }

  est_u <-
    gen_mean(y1m = y1m_t, p = p1_t, lower_bound = FALSE, minY = minY, maxY = maxY) -
    gen_mean(y1m = y1m_c, p = p1_c, lower_bound = TRUE, minY = minY, maxY = maxY)
  est_l <-
    gen_mean(y1m = y1m_t, p = p1_t, lower_bound = TRUE, minY = minY, maxY = maxY) -
    gen_mean(y1m = y1m_c, p = p1_c, lower_bound = FALSE, minY = minY, maxY = maxY)

  var_u <-
    gen_var(y1m = y1m_t, s1 = s1_t, p = p1_t, lower_bound = FALSE, minY = minY, maxY = maxY) +
    gen_var(y1m = y1m_c, s1 = s1_c, p = p1_c, lower_bound = TRUE, minY = minY, maxY = maxY)
  var_l <-
    gen_var(y1m = y1m_t, s1 = s1_t, p = p1_t, lower_bound = TRUE, minY = minY, maxY = maxY) +
    gen_var(y1m = y1m_c, s1 = s1_c, p = p1_c, lower_bound = FALSE, minY = minY, maxY = maxY)


  im_crit <- function(ca) abs(pnorm(ca + (est_u-est_l)/max(var_u,var_l))-pnorm(-ca)-(1-alpha))
  sig <- optim(1.60,im_crit,method="Brent",lower=1,upper=2)$par

  return(c(ci_lower=est_l - sig*var_l^.5,
           ci_upper=est_u + sig*var_u^.5,
           low_est=est_l,
           upp_est=est_u,
           low_var=var_l,
           upp_var=var_u))
}

