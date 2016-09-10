# Estimators for attrition package

manski_cis <- function(n1_t, n1_c,
                       n1_t_s, n1_c_s,
                       p1_t, p1_c,
                       y1m_t, y1m_c,
                       s1_t, s1_c,
                       minY,maxY,alpha){

  upper_bound_est <-
    gen_mean(y1m_t, p1_t, lower_bound = FALSE, minY = minY, maxY = maxY) -
    gen_mean(y1m_c, p1_c, lower_bound = TRUE, minY = minY, maxY = maxY)
  lower_bound_est <-
    gen_mean(y1m_t, p1_t, lower_bound = TRUE, minY = minY, maxY = maxY) -
    gen_mean(y1m_c, p1_c, lower_bound = FALSE, minY = minY, maxY = maxY)

  upper_bound_var_est <-
    gen_var(y1m_t, s1_t, p1_t, lower_bound = FALSE, minY = minY, maxY = maxY)/n1_t +
    gen_var(y1m_c, s1_c, p1_c, lower_bound = TRUE, minY = minY, maxY = maxY)/n1_c
  lower_bound_var_est <-
    gen_var(y1m_t, s1_t, p1_t, lower_bound = TRUE, minY = minY, maxY = maxY)/n1_t +
    gen_var(y1m_c, s1_c, p1_c, lower_bound = FALSE, minY = minY, maxY = maxY)/n1_c


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

ds_manski_cis_2s <- function(n1_t,n2_t,n1_c,n2_c,
                             p1_t,p2_t,s1_t,s2_nm_t,
                             y1m_t,y2m_nm_t,
                             c1a_t,c1r_t,c2a_t,c2r_t,
                             p1_c,p2_c,
                             s1_c,s2_nm_c,
                             y1m_c,y2m_nm_c,
                             minY,maxY,alpha) {

  y2m_t_L <- gen_mean(y2m_nm_t,p2_t,lower_bound=TRUE,minY,maxY)
  y2m_t_U <- gen_mean(y2m_nm_t,p2_t,lower_bound=FALSE,minY,maxY)
  y2m_c_L <- gen_mean(y2m_nm_c,p2_c,lower_bound=TRUE,minY,maxY)
  y2m_c_U <- gen_mean(y2m_nm_c,p2_c,lower_bound=FALSE,minY,maxY)

  s2_t_L <- gen_var(y2m_nm_t,s2_nm_t,p2_t,lower_bound=TRUE,minY,maxY)^.5
  s2_t_U <- gen_var(y2m_nm_t,s2_nm_t,p2_t,lower_bound=FALSE,minY,maxY)^.5
  s2_c_L <- gen_var(y2m_nm_c,s2_nm_c,p2_c,lower_bound=TRUE,minY,maxY)^.5
  s2_c_U <- gen_var(y2m_nm_c,s2_nm_c,p2_c,lower_bound=FALSE,minY,maxY)^.5

  manski_bounds_est <- construct_manski_bounds(p1_t, y1m_t,
                                               p1_c, y1m_c,
                                               y2m_t_L, y2m_t_U,
                                               y2m_c_L, y2m_c_U)

  lower_bound_est <-   manski_bounds_est[1]
  upper_bound_est <-   manski_bounds_est[2]

  lower_bound_var_est <- ds_var_2s(c(n1_t,n2_t,p1_t,p2_t,s1_t,s2_t_L,y1m_t,y2m_t_L),
                                   c(n1_c,n2_c,p1_c,p2_c,s1_c,s2_c_U,y1m_c,y2m_c_U))
  upper_bound_var_est <- ds_var_2s(c(n1_t,n2_t,p1_t,p2_t,s1_t,s2_t_U,y1m_t,y2m_t_U),
                                   c(n1_c,n2_c,p1_c,p2_c,s1_c,s2_c_L,y1m_c,y2m_c_L))

  sig <- optim(1.60,im_crit,method="Brent",lower=1,upper=2,
               lower_bound_est = lower_bound_est,
               upper_bound_est = upper_bound_est,
               lower_bound_var_est = lower_bound_var_est,
               upper_bound_var_est = upper_bound_var_est,
               alpha = alpha)$par

  return(c(ci_lower=lower_bound_est - sig*lower_bound_var_est^.5,
           ci_upper=upper_bound_est + sig*upper_bound_var_est^.5,
           low_est=lower_bound_est,upp_est=upper_bound_est,
           low_var=lower_bound_var_est,upp_var=upper_bound_var_est))
}

ds_manski_cis_2s_sens <- function(n1_t,n2_t,n1_c,n2_c,
                                  p1_t,p2_t,s1_t,s2_nm_t,
                                  y1m_t,y2m_nm_t,
                                  c1a_t,c1r_t,c2a_t,c2r_t,
                                  p1_c,p2_c,
                                  s1_c,s2_nm_c,
                                  y1m_c,y2m_nm_c,
                                  minY,maxY,alpha,delta){

  y2m_t_L <- gen_mean_sens(y2m_nm_t, p2_t, delta, lower_bound = TRUE, minY, maxY)
  y2m_t_U <- gen_mean_sens(y2m_nm_t, p2_t, delta, lower_bound = FALSE, minY, maxY)
  y2m_c_L <- gen_mean_sens(y2m_nm_c, p2_c, delta, lower_bound = TRUE, minY, maxY)
  y2m_c_U <- gen_mean_sens(y2m_nm_c, p2_c, delta, lower_bound = FALSE, minY, maxY)

  s2_t_L <- gen_var_sens(y2m_nm_t, s2_nm_t, p2_t, delta, lower_bound = TRUE, minY, maxY)^.5
  s2_t_U <- gen_var_sens(y2m_nm_t, s2_nm_t, p2_t, delta, lower_bound = FALSE, minY, maxY)^.5
  s2_c_L <- gen_var_sens(y2m_nm_c, s2_nm_c, p2_c, delta, lower_bound = TRUE, minY, maxY)^.5
  s2_c_U <- gen_var_sens(y2m_nm_c, s2_nm_c, p2_c, delta, lower_bound = FALSE, minY, maxY)^.5

  manski_bounds_est <- construct_manski_bounds(p1_t, y1m_t,
                                               p1_c, y1m_c,
                                               y2m_t_L, y2m_t_U,
                                               y2m_c_L, y2m_c_U)

  lower_bound_est <-   manski_bounds_est[1]
  upper_bound_est <-   manski_bounds_est[2]

  lower_bound_var_est <- ds_var_2s(c(n1_t,n2_t,p1_t,p2_t,s1_t,s2_t_L,y1m_t,y2m_t_L),
                                   c(n1_c,n2_c,p1_c,p2_c,s1_c,s2_c_U,y1m_c,y2m_c_U))
  upper_bound_var_est <- ds_var_2s(c(n1_t,n2_t,p1_t,p2_t,s1_t,s2_t_U,y1m_t,y2m_t_U),
                                   c(n1_c,n2_c,p1_c,p2_c,s1_c,s2_c_L,y1m_c,y2m_c_L))


  sig <- optim(1.60,im_crit,method="Brent",lower=1,upper=2,
               lower_bound_est = lower_bound_est,
               upper_bound_est = upper_bound_est,
               lower_bound_var_est = lower_bound_var_est,
               upper_bound_var_est = upper_bound_var_est,
               alpha = alpha)$par

  return(c(ci_lower= lower_bound_est - sig*lower_bound_var_est^.5,
           ci_upper= upper_bound_est + sig*upper_bound_var_est^.5,
           low_est= lower_bound_est,upp_est=upper_bound_est,
           low_var= lower_bound_var_est,upp_var=upper_bound_var_est))
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













#' # calculate manski bounds
##' @export
# ds_manski <- function(p1,p2,y1m,y2m_nm,minY,maxY) {
#   const1 <- p1*y1m + (1-p1)*p2*y2m_nm
#   const2 <- (1-p1)*(1-p2)
#   return(c(const1+const2*minY,const1+const2*maxY))
# }


#' #' find the n1,n2 (number of attempts to measure) that minimize width of Manski CIs
#' #' @export
#' optim_manski_2s <- function(p1_t,p2_t,s1_t,s2_nm_t,
#'                             y1m_t,y2m_nm_t,
#'                             c1a_t,c1r_t,c2a_t,c2r_t,
#'                             p1_c,p2_c,s1_c,s2_nm_c,
#'                             y1m_c,y2m_nm_c,
#'                             c1a_c,c1r_c,c2a_c,c2r_c,
#'                             budget,minY,maxY,
#'                             mu=1e-04,niter=100,alpha=0.05) {
#'
#'
#'   matout <- matrix(NA,nrow=niter,ncol=5)
#'   colnames(matout) <- c("n1_t","n2_t","n1_c","n2_c","ci_width")
#'   iter <- 1
#'
#'   while(iter <= niter) {
#'
#'     #(c1a_t+c1r_t*p1_t)*n1_t + (c2a_t+c2r_t*p2_t)*n2_t + (c1a_c+c1r_c*p1_c)*n1_c + (c2a_c+c2r_c*p2_c)*n2_c = budget
#'
#'     starting_values <- .5*runif(4)*budget/c((c1a_t+c1r_t*p1_t),(c2a_t+c2r_t*p2_t*(1-p1_t)),(c1a_c+c1r_c*p1_c),(c2a_c+c2r_c*p2_c*(1-p1_c)))
#'
#'     n1g_t = starting_values[1]
#'     n2g_t = starting_values[2]
#'     n1g_c = starting_values[3]
#'     n2g_c = starting_values[4]
#'
#'
#'     # constraints
#'     #- (c1a_t+c1r_t*p1_t)*n1_t - (c2a_t+c2r_t*p2_t)*n2_t - (c1a_c+c1r_c*p1_c)*n1_c - (c2a_c+c2r_c*p2_c)*n2_c + budget >= 0
#'     # n1_t > 2
#'     # n2_t > 2
#'     # n1_c > 2
#'     # n2_c > 2
#'     # (1-p1_t)*n1_t - n2_t > 2
#'     # (1-p1_c)*n1_c - n2_c > 2
#'
#'     constrMat <- rbind(c(-(c1a_t+c1r_t*p1_t),-(c2a_t+c2r_t*p2_t),-(c1a_c+c1r_c*p1_c),-(c2a_c+c2r_c*p2_c)),
#'                        c(1,0,0,0),
#'                        c(0,1,0,0),
#'                        c(0,0,1,0),
#'                        c(0,0,0,1),
#'                        c((1-p1_t),-1,0,0),
#'                        c(0,0,(1-p1_c),-1)
#'     ) # are these working?
#'     constrVec <- c(-budget,2,2,2,2,2,2)
#'
#'     # compute width of 95% CIs
#'     ds_passthrough_2s <- function(ns,p1_t,p2_t,s1_t,s2_nm_t,
#'                                   y1m_t,y2m_nm_t,
#'                                   c1a_t,c1r_t,c2a_t,c2r_t,
#'                                   p1_c,p2_c,s1_c,s2_nm_c,
#'                                   y1m_c,y2m_nm_c,
#'                                   minY,maxY,alpha){
#'       diff(ds_manski_cis_2s(ns[1],ns[2],ns[3],ns[4],p1_t,p2_t,s1_t,s2_nm_t,y1m_t,y2m_nm_t,c1a_t,c1r_t,c2a_t,c2r_t,p1_c,p2_c,s1_c,s2_nm_c,y1m_c,y2m_nm_c,minY,maxY,alpha))
#'     }
#'
#'     optimD <- try(constrOptim(c(n1g_t,n2g_t,n1g_c,n2g_c),
#'                               ds_passthrough_2s,grad=NULL,ui=constrMat,ci=constrVec,
#'                               p1_t=p1_t,p2_t=p2_t,
#'                               s1_t=s1_t,s2_nm_t=s2_nm_t,
#'                               y1m_t=y1m_t,y2m_nm_t=y2m_nm_t,
#'                               c1a_t=c1a_t,c1r_t=c1r_t,
#'                               c2a_t=c2a_t,c2r_t=c2r_t,
#'                               p1_c=p1_c,p2_c=p2_c,
#'                               s1_c=s1_c,s2_nm_c=s2_nm_c,
#'                               y1m_c=y1m_c,y2m_nm_c=y2m_nm_c,
#'                               minY=minY,maxY=maxY,alpha=alpha,
#'                               outer.iterations=1000,outer.eps = .Machine$double.eps^.5,mu=mu),silent=TRUE)
#'
#'     if(!is.character(optimD)) {
#'       cat(iter,"")
#'       matout[iter,] <- c(n1_t = optimD$par[1],n2_t = optimD$par[2],n1_c = optimD$par[3],n2_c = optimD$par[4],ci_width=optimD$value)
#'       iter <- iter + 1
#'     }
#'
#'   }
#'   return(matout[order(matout[,5]),])
#' }

# imputation_estimator_ds <-
#   function(p1_t,p2_t,y1m_t,y2m_t,p1_c,p2_c,y1m_c,y2m_c, minY, maxY, p) {
#
#     y_U_tilde_t <- (1-p)*y2m_t + p*maxY
#     y_L_tilde_t <- (1-p)*y2m_t + p*minY
#     y_U_tilde_c <- (1-p)*y2m_c + p*maxY
#     y_L_tilde_c <- (1-p)*y2m_c + p*minY
#
#     const1_t <- p1_t*y1m_t + (1-p1_t)*p2_t*y2m_t
#     const2_t <- (1-p1_t)*(1-p2_t)
#     const1_c <- p1_c*y1m_c + (1-p1_c)*p2_c*y2m_c
#     const2_c <- (1-p1_c)*(1-p2_c)
#     return(c(const1_t+const2_t*y_L_tilde_t - (const1_c+const2_c*y_U_tilde_c),
#              const1_t+const2_t*y_U_tilde_t - (const1_c+const2_c*y_L_tilde_c)))
#   }
