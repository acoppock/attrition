# Estimators for attrition package

# Extreme value (Manski) bounds from a single round of data collection, with a
# joint Imbens-Manski interval. Every nonrespondent gets the extreme value, so
# delta is fixed at one.
manski_cis <- function(n1_t, n1_c,
                       p1_t, p1_c,
                       y1m_t, y1m_c,
                       s1_t, s1_c,
                       minY, maxY, alpha){

  upper_bound_est <-
    gen_mean(y1m_t, p1_t, delta = 1, lower_bound = FALSE, minY = minY, maxY = maxY) -
    gen_mean(y1m_c, p1_c, delta = 1, lower_bound = TRUE, minY = minY, maxY = maxY)
  lower_bound_est <-
    gen_mean(y1m_t, p1_t, delta = 1, lower_bound = TRUE, minY = minY, maxY = maxY) -
    gen_mean(y1m_c, p1_c, delta = 1, lower_bound = FALSE, minY = minY, maxY = maxY)

  upper_bound_var_est <-
    gen_var(y1m_t, s1_t, p1_t, delta = 1, lower_bound = FALSE, minY = minY, maxY = maxY)/n1_t +
    gen_var(y1m_c, s1_c, p1_c, delta = 1, lower_bound = TRUE, minY = minY, maxY = maxY)/n1_c
  lower_bound_var_est <-
    gen_var(y1m_t, s1_t, p1_t, delta = 1, lower_bound = TRUE, minY = minY, maxY = maxY)/n1_t +
    gen_var(y1m_c, s1_c, p1_c, delta = 1, lower_bound = FALSE, minY = minY, maxY = maxY)/n1_c

  im_interval(lower_bound_est, upper_bound_est, lower_bound_var_est, upper_bound_var_est, alpha)
}

# Double-sampling bounds with the analytic variance of Coppock, Gerber, Green,
# and Kern (2017) and a joint Imbens-Manski interval. `m` is the list
# ds_moments() returns. delta is the share of the follow-up nonrespondents given
# worst-case treatment: at delta = 1 this is estimator_ds(), and at anything
# less it is estimator_ds_sens().
ds_manski_cis_2s <- function(m, minY, maxY, alpha, delta) {
  t <- m$t
  ctl <- m$c

  y2m_t_L <- gen_mean(t$y2m, t$p2, delta, lower_bound = TRUE, minY, maxY)
  y2m_t_U <- gen_mean(t$y2m, t$p2, delta, lower_bound = FALSE, minY, maxY)
  y2m_c_L <- gen_mean(ctl$y2m, ctl$p2, delta, lower_bound = TRUE, minY, maxY)
  y2m_c_U <- gen_mean(ctl$y2m, ctl$p2, delta, lower_bound = FALSE, minY, maxY)

  s2_t_L <- gen_var(t$y2m, t$s2, t$p2, delta, lower_bound = TRUE, minY, maxY)^0.5
  s2_t_U <- gen_var(t$y2m, t$s2, t$p2, delta, lower_bound = FALSE, minY, maxY)^0.5
  s2_c_L <- gen_var(ctl$y2m, ctl$s2, ctl$p2, delta, lower_bound = TRUE, minY, maxY)^0.5
  s2_c_U <- gen_var(ctl$y2m, ctl$s2, ctl$p2, delta, lower_bound = FALSE, minY, maxY)^0.5

  bounds <- construct_manski_bounds(t$p1, t$y1m,
                                    ctl$p1, ctl$y1m,
                                    y2m_t_L, y2m_t_U,
                                    y2m_c_L, y2m_c_U)

  lower_bound_var_est <-
    ds_var(t$n1, t$n2, t$p1, t$p2, t$s1, s2_t_L, t$y1m, y2m_t_L) +
    ds_var(ctl$n1, ctl$n2, ctl$p1, ctl$p2, ctl$s1, s2_c_U, ctl$y1m, y2m_c_U)
  upper_bound_var_est <-
    ds_var(t$n1, t$n2, t$p1, t$p2, t$s1, s2_t_U, t$y1m, y2m_t_U) +
    ds_var(ctl$n1, ctl$n2, ctl$p1, ctl$p2, ctl$s1, s2_c_L, ctl$y1m, y2m_c_L)

  im_interval(bounds[1], bounds[2], lower_bound_var_est, upper_bound_var_est, alpha)
}


# Nonparametric bootstrap for the trimming bounds. Units are resampled within
# treatment arm, holding the arm sizes fixed as complete random assignment does.
# Works for both the single-stage and the weighted double-sampling estimator,
# which is the reason it exists: Lee's analytic variance covers only the former.
#
# Replicates in which monotonicity fails contribute NA and are counted, since
# silently dropping them would understate the variance without saying so.
bootstrap_trim_variance <- function(compute_bounds, Z, sims) {
  idx_t <- which(Z == 1)
  idx_c <- which(Z == 0)

  reps <- vapply(seq_len(sims), function(i) {
    take <- c(sample(idx_t, length(idx_t), replace = TRUE),
              sample(idx_c, length(idx_c), replace = TRUE))
    out <- tryCatch(compute_bounds(take),
                    attrition_monotonicity_violation = function(e) c(NA_real_, NA_real_),
                    error = function(e) c(NA_real_, NA_real_))
    c(lower = unname(out[1]), upper = unname(out[2]))
  }, numeric(2))

  n_ok <- sum(stats::complete.cases(t(reps)))
  if (n_ok < 2) {
    stop("Bootstrap failed: fewer than two replicates produced usable bounds.")
  }
  if (n_ok < sims) {
    warning(sims - n_ok, " of ", sims,
            " bootstrap replicates did not yield bounds (monotonicity violated in the resample). ",
            "Standard errors are computed from the ", n_ok, " that did.",
            call. = FALSE)
  }

  c(lower_var = stats::var(reps["lower", ], na.rm = TRUE),
    upper_var = stats::var(reps["upper", ], na.rm = TRUE),
    n_boot = n_ok)
}

# Asymptotic variance of the Lee (2009) trimming bounds, Proposition 3 (eq. 7).
# Each bound estimate carries four independent contributions: the variance of
# the retained (trimmed) treated outcomes; the variance from estimating the
# trimming THRESHOLD, a quantile; the variance from estimating the trimming
# PROPORTION, which in Lee's own Job Corps example is the largest of the three;
# and the variance of the control-group respondent mean.
#
# The third term is written here in the algebraically equivalent form used by
# Tauchmann (2014, Stata Journal 14(4):884-894). Lee's published expression,
#   ((1 - p_c) - Q(1 - E[Z])) / (E[Z] p_c (1 - E[Z])),
# is identical but can look negative; this form is manifestly non-negative.
#
# Applies only to the single-stage, unweighted, monotonicity case. Lee's
# derivation assumes i.i.d. sampling and trimming of one group only.
lee_variance <- function(trim_out, n_treat, n_control) {
  Q   <- unname(trim_out["Q"])
  p_t <- unname(trim_out["pi_r_1"])
  p_c <- unname(trim_out["pi_r_0"])

  # Variance of the estimated trimming proportion Q = (p_t - p_c)/p_t, by the
  # delta method, scaled by the squared derivative of the trimmed mean in Q.
  trim_prop_term <- (1 - p_t)/(p_t * n_treat) + (1 - p_c)/(p_c * n_control)

  var_bound <- function(side) {
    y        <- unname(trim_out[if (side == "upper") "yU" else "yL"])
    mu       <- unname(trim_out[if (side == "upper") "Out1U_mono" else "Out1L_mono"])
    var_keep <- unname(trim_out[if (side == "upper") "var_keep_U" else "var_keep_L"])
    n_keep   <- unname(trim_out[if (side == "upper") "n_keep_U" else "n_keep_L"])
    var_keep/n_keep +                            # trimmed distribution
      (y - mu)^2 * Q/n_keep +                    # estimated threshold
      (y - mu)^2 * trim_prop_term +              # estimated trimming proportion
      unname(trim_out["var_control"])/unname(trim_out["control_group_N"])
  }

  c(lower_var = var_bound("lower"), upper_var = var_bound("upper"))
}

# Map a trimming result computed on relabelled arms back to the original ones.
# Reverse monotonicity is the forward estimator run with the arms swapped, so
# what comes back is the effect of control: each bound is negated and the pair
# swapped, and every quantity carrying a group label moves with them. The U/L
# suffix names the bound a quantity feeds, which is what it names in the
# forward case too, so Out0L_mono is the trimmed control mean behind the lower
# bound even though it is the higher of the two control means.
reverse_monotonicity_labels <- function(out) {
  c(upper_bound     = -unname(out["lower_bound"]),
    lower_bound     = -unname(out["upper_bound"]),
    Out1_mono       =  unname(out["Out0_mono"]),
    Out0L_mono      =  unname(out["Out1U_mono"]),
    Out0U_mono      =  unname(out["Out1L_mono"]),
    control_group_N =  unname(out["treat_group_N"]),
    treat_group_N   =  unname(out["control_group_N"]),
    Q               =  unname(out["Q"]),
    f1              =  unname(out["f0"]),
    f0              =  unname(out["f1"]),
    pi_r_1          =  unname(out["pi_r_0"]),
    pi_r_0          =  unname(out["pi_r_1"]),
    yU              =  unname(out["yL"]),
    yL              =  unname(out["yU"]))
}

trimming_bounds <-
  function(Out, Treat, Fail, Weight, monotonicity = FALSE) {

    dataf <- data.frame(Out, Treat, Fail, Weight)
    datafsort <- dataf[order(Out),]

    OutS0 <- datafsort[datafsort$Fail==0 & datafsort$Treat==0,]
    OutS1 <- datafsort[datafsort$Fail==0 & datafsort$Treat==1,]

    if(nrow(OutS0) == 0 | nrow(OutS1) == 0){
      stop("Trimming bounds require at least one observed outcome in each treatment group.")
    }

    OutS0$Weight <- OutS0$Weight/sum(OutS0$Weight)
    OutS1$Weight <- OutS1$Weight/sum(OutS1$Weight)

    OutS0.CDF <- cumsum(OutS0$Weight)
    OutS1.CDF <- cumsum(OutS1$Weight)

    f0 <- sum(Weight[Fail==1 & Treat==0])/sum(Weight[Treat==0])
    f1 <- sum(Weight[Fail==1 & Treat==1])/sum(Weight[Treat==1])

    if(monotonicity){
      Q <- ((1 - f1) - (1 - f0))/(1-f1)
      if(Q < 0){
        stop(structure(
          class = c("attrition_monotonicity_violation", "error", "condition"),
          list(message = "Monotonicity appears to be violated: the treatment group is more likely to be missing than the control group.",
               call = NULL)
        ))
      }

      Out0_mono <- weighted.mean(OutS0$Out, OutS0$Weight)

      if (Q == 0) {
        Out1_mean <- weighted.mean(OutS1$Out, OutS1$Weight)
        var_treat <- sum(OutS1$Weight * (OutS1$Out - Out1_mean)^2)
        return(c(upper_bound = Out1_mean - Out0_mono, lower_bound = Out1_mean - Out0_mono,
                 Out0_mono = Out0_mono, Out1L_mono = Out1_mean, Out1U_mono = Out1_mean,
                 control_group_N = nrow(OutS0), treat_group_N = nrow(OutS1),
                 Q = 0, f1 = f1, f0 = f0, pi_r_1 = 1 - f1, pi_r_0 = 1 - f0,
                 yU = min(OutS1$Out), yL = max(OutS1$Out),
                 var_keep_U = var_treat, var_keep_L = var_treat,
                 n_keep_U = nrow(OutS1), n_keep_L = nrow(OutS1),
                 var_control = sum(OutS0$Weight * (OutS0$Out - Out0_mono)^2)))
      }

      keep_U <- OutS1.CDF > Q
      keep_L <- OutS1.CDF < (1-Q)
      if (!any(keep_U) || !any(keep_L)) {
        stop("Trimming ", signif(Q, 3), " of the trimmed group leaves nothing behind on one ",
             "side. The group has too few distinct weighted observations to support a ",
             "trimming proportion this large.", call. = FALSE)
      }
      Out1U_mono <- weighted.mean(OutS1$Out[keep_U], OutS1$Weight[keep_U])
      Out1L_mono <- weighted.mean(OutS1$Out[keep_L], OutS1$Weight[keep_L])

      upper_bound <- Out1U_mono - Out0_mono
      lower_bound <- Out1L_mono - Out0_mono

      # Trimming thresholds actually used, and the variance of the retained
      # outcomes on each side. Lee (2009) Proposition 3 needs both.
      yU <- min(OutS1$Out[keep_U])
      yL <- max(OutS1$Out[keep_L])

      return(c(upper_bound = upper_bound, lower_bound = lower_bound,
               Out0_mono = Out0_mono, Out1L_mono=Out1L_mono, Out1U_mono = Out1U_mono,
               control_group_N = nrow(OutS0), treat_group_N = nrow(OutS1), Q = Q, f1 = f1, f0 = f0, pi_r_1 = 1 - f1, pi_r_0 = 1 - f0,
               yU = yU, yL = yL,
               var_keep_U = sum(OutS1$Weight[keep_U]/sum(OutS1$Weight[keep_U]) * (OutS1$Out[keep_U] - Out1U_mono)^2),
               var_keep_L = sum(OutS1$Weight[keep_L]/sum(OutS1$Weight[keep_L]) * (OutS1$Out[keep_L] - Out1L_mono)^2),
               n_keep_U = sum(keep_U), n_keep_L = sum(keep_L),
               var_control = sum(OutS0$Weight * (OutS0$Out - Out0_mono)^2)))

    }else{

      # Without monotonicity the always-reporter share is only bounded below, by the
      # Frechet-Hoeffding bound 1 - f0 - f1 (Imai 2008, Eq. 5). Each group is trimmed
      # by the largest share of its respondents that could fail to be always-reporters,
      # which is that bound subtracted from the group's own response rate. The bounds
      # exist only while the Frechet bound is positive.
      if (f0 + f1 >= 1) {
        stop("The two missingness rates sum to ", signif(f0 + f1, 4), ", which is not less ",
             "than one, so nothing bounds the always-reporter share away from zero and the ",
             "trimming bounds are undefined. Monotonicity, or a design that recovers some ",
             "of the missing outcomes, is what makes this case tractable.", call. = FALSE)
      }

      trim0 <- (f1)/(1-f0)
      trim1 <- (f0)/(1-f1)

      keep0_U <- OutS0.CDF > trim0
      keep0_L <- OutS0.CDF < (1-trim0)
      keep1_U <- OutS1.CDF > trim1
      keep1_L <- OutS1.CDF < (1-trim1)
      if (!any(keep0_U) || !any(keep0_L) || !any(keep1_U) || !any(keep1_L)) {
        stop("Trimming ", signif(trim0, 3), " of the control group and ", signif(trim1, 3),
             " of the treatment group leaves nothing behind on one side. The groups have too ",
             "few distinct weighted observations to support trimming proportions this large.",
             call. = FALSE)
      }

      Out0U <- weighted.mean(OutS0$Out[keep0_U], OutS0$Weight[keep0_U])
      Out0L <- weighted.mean(OutS0$Out[keep0_L], OutS0$Weight[keep0_L])

      Out1U <- weighted.mean(OutS1$Out[keep1_U], OutS1$Weight[keep1_U])
      Out1L <- weighted.mean(OutS1$Out[keep1_L], OutS1$Weight[keep1_L])

      upper_bound = Out1U - Out0L
      lower_bound = Out1L - Out0U

      return(c(upper_bound = upper_bound, lower_bound = lower_bound, Out0L=Out0L, Out0U=Out0U, Out1L=Out1L, Out1U=Out1U,
               control_group_N = nrow(OutS0), treat_group_N = nrow(OutS1), trim0 = trim0, trim1 = trim1,
               f1 = f1, f0 = f0, pi_r_1 = 1 - f1, pi_r_0 = 1 - f0))

    }
  }
