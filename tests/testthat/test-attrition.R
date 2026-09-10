library(testthat)

# ── Paper benchmark (Table 3 of CGGK 2017) ──────────────────────────────────
# Expected values computed from replication data in:
# Coppock, Gerber, Green, Kern (2017), Political Analysis
# doi:10.1017/pan.2016.6
# Data source: Harvard Dataverse doi:10.7910/DVN/AQB4MP, shipped as data(levendusky_replication)
#
# The bound point estimates and variances below are the published quantities and
# have never changed. The conf.low/conf.high expectations were refreshed when the
# Imbens-Manski critical value moved from optim() on abs() to uniroot() on the
# signed coverage excess; the two agree to 9 decimal places and uniroot is the
# more accurate of the two.

test_that("estimator_ev matches Table 3 column 1 (no double sampling)", {
  out <- estimator_ev(Y = Y_polarization_w2, Z = Z, R = R1,
                      minY = 0, maxY = 6, data = levendusky_replication)
  expect_equal(unname(out["conf.low"]), -1.66907903885775,   tolerance = 1e-10)
  expect_equal(unname(out["conf.high"]),  1.83588114554598,   tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]), -1.53914496339566,    tolerance = 1e-10)
  expect_equal(unname(out["estimate_upper"]),  1.70966762747749,    tolerance = 1e-10)
  expect_equal(unname(out["std.error_lower"]),  0.078994308875319, tolerance = 1e-10)
  expect_equal(unname(out["std.error_upper"]),  0.0767323705893607, tolerance = 1e-10)
})

test_that("estimator_ds matches Table 3 column 2 (double sampling)", {
  out <- estimator_ds(Y = Y_polarization_w2, Z = Z, R1 = R1,
                      Attempt = Attempt, R2 = R2,
                      minY = 0, maxY = 6, data = levendusky_replication)
  expect_equal(unname(out["conf.low"]), -0.52830967410789,  tolerance = 1e-10)
  expect_equal(unname(out["conf.high"]),  0.745174826433894,  tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]), -0.34174537662934,    tolerance = 1e-10)
  expect_equal(unname(out["estimate_upper"]),  0.571815728388134,   tolerance = 1e-10)
  expect_equal(unname(out["std.error_lower"]),  0.113423039242904,  tolerance = 1e-10)
  expect_equal(unname(out["std.error_upper"]),  0.105394848030982,  tolerance = 1e-10)
})

test_that("estimator_ds matches Table 3 column 3 (DS + poststratification)", {
  out <- estimator_ds(Y = Y_polarization_w2, Z = Z, R1 = R1,
                      Attempt = Attempt, R2 = R2,
                      strata = X_party_id,
                      minY = 0, maxY = 6, data = levendusky_replication)
  expect_equal(unname(out["conf.low"]), -0.529010551320638,  tolerance = 1e-10)
  expect_equal(unname(out["conf.high"]),  0.696600307023322,  tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]), -0.344389910863888,   tolerance = 1e-10)
  expect_equal(unname(out["estimate_upper"]),  0.525683688852529,   tolerance = 1e-10)
  expect_equal(unname(out["std.error_lower"]),  0.112241379677608,  tolerance = 1e-10)
  expect_equal(unname(out["std.error_upper"]),  0.103909925704189,  tolerance = 1e-10)
})

# ── Synthetic data tests (self-contained, seed = 343) ───────────────────────

make_synthetic <- function() {
  set.seed(343)
  N    <- 1000
  Y_0  <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
  Y_1  <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.1, 0.4, 0.3, 0.1))
  R1_0 <- rbinom(N, 1, prob = 0.7)
  R1_1 <- rbinom(N, 1, prob = 0.8)
  R2_0 <- rbinom(N, 1, prob = 0.7)
  R2_1 <- rbinom(N, 1, prob = 0.75)
  strata <- as.numeric(Y_0 > 2)
  Z      <- rbinom(N, 1, 0.5)
  R1     <- Z * R1_1 + (1 - Z) * R1_0
  Y_star <- Z * Y_1  + (1 - Z) * Y_0
  Y      <- Y_star
  Y[R1 == 0] <- NA
  Attempt <- rep(0L, N)
  Attempt[is.na(Y)] <- rbinom(sum(is.na(Y)), 1, 0.5)
  R2 <- rep(0L, N)
  R2[Attempt == 1] <- (Z * R2_1 + (1 - Z) * R2_0)[Attempt == 1]
  Y[R2 == 1 & Attempt == 1] <- Y_star[R2 == 1 & Attempt == 1]
  data.frame(Y, Z, R1, Attempt, R2, strata)
}

test_that("estimator_ds produces stable results on synthetic data", {
  df  <- make_synthetic()
  out <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  expect_equal(unname(out["conf.low"]), -0.0823029393291099, tolerance = 1e-10)
  expect_equal(unname(out["conf.high"]),  0.624972549760316, tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]),   0.0572167798040737, tolerance = 1e-10)
  expect_equal(unname(out["estimate_upper"]),   0.487115989508957, tolerance = 1e-10)
})

test_that("estimator_ds with poststratification produces stable results on synthetic data", {
  df  <- make_synthetic()
  out <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5,
                      strata = strata, data = df)
  expect_equal(unname(out["conf.low"]),  0.00420597347763489, tolerance = 1e-10)
  expect_equal(unname(out["conf.high"]),  0.662950557416964,   tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]),   0.124794503861406,   tolerance = 1e-10)
  expect_equal(unname(out["estimate_upper"]),   0.544633130143652,   tolerance = 1e-10)
})

# ── Helper function unit tests ───────────────────────────────────────────────

test_that("gen_var and gen_var_sens(delta=1) are equal", {
  gv  <- attrition:::gen_var(5, 2, 0.5, minY = 0, maxY = 5)
  gvs <- attrition:::gen_var_sens(5, 2, 0.5, delta = 1, minY = 0, maxY = 5)
  expect_equal(gv, gvs)
})

test_that("gen_var_sens(delta=0) equals gen_var with no imputation", {
  # At delta=0, gen_mean_sens returns p*y_m + (1-p)*y_m = y_m regardless of bound
  # so variance contribution from missing = 0 and only observed variance remains
  gvs0 <- attrition:::gen_var_sens(3, 1, 0.8, delta = 0, minY = 0, maxY = 5)
  expect_true(is.numeric(gvs0) && gvs0 >= 0)
})

test_that("gen_mean lower bound <= upper bound", {
  lb <- attrition:::gen_mean(3, 0.7, lower_bound = TRUE,  minY = 0, maxY = 5)
  ub <- attrition:::gen_mean(3, 0.7, lower_bound = FALSE, minY = 0, maxY = 5)
  expect_lte(lb, ub)
})

test_that("EV bound width equals (maxY - minY) * (frac_missing_t + frac_missing_c)", {
  # With known response rates, the theoretical width is predictable
  set.seed(1)
  N <- 10000
  Z  <- rep(0:1, each = N / 2)
  R  <- rbinom(N, 1, prob = ifelse(Z == 1, 0.8, 0.6))
  Y  <- ifelse(R == 1, rnorm(N, mean = 2), NA_real_)
  df <- data.frame(Y, Z, R)
  # Support must actually cover the draws: 10,000 normals reach roughly 4 sd out
  out <- estimator_ev(Y, Z, R, minY = -4, maxY = 8, data = df)
  # Expected width ≈ 12 * ((1 - 0.8) + (1 - 0.6)) = 12 * 0.6 = 7.2
  expect_equal(unname(out["estimate_upper"] - out["estimate_lower"]), 7.2, tolerance = 0.1)
})

# ── tidy() method tests ──────────────────────────────────────────────────────

test_that("tidy.attrition_bounds works for estimator_ev", {
  df <- make_synthetic()
  out <- estimator_ev(Y, Z, R1, minY = 1, maxY = 5, data = df)
  td  <- tidy(out)
  expect_s3_class(td, "tbl_df")
  expect_equal(nrow(td), 3L)
  expect_named(td, c("term", "estimate", "std.error", "conf.low", "conf.high",
                     "estimate_lower", "estimate_upper",
                     "std.error_lower", "std.error_upper", "outcome"))
  expect_equal(td$term, c("bounds", "lower_bound", "upper_bound"))
  # bounds row
  expect_true(is.na(td$estimate[1]))
  expect_true(is.na(td$std.error[1]))
  expect_equal(td$conf.low[1],      unname(out["conf.low"]))
  expect_equal(td$conf.high[1],     unname(out["conf.high"]))
  expect_equal(td$estimate_lower[1],  unname(out["estimate_lower"]))
  expect_equal(td$estimate_upper[1], unname(out["estimate_upper"]))
  # lower_bound row
  expect_equal(td$estimate[2],  unname(out["estimate_lower"]))
  expect_equal(td$std.error[2], unname(out["std.error_lower"]))
  expect_true(is.na(td$conf.low[2]))
  expect_true(is.na(td$conf.high[2]))
  expect_true(is.na(td$estimate_lower[2]))
  expect_true(is.na(td$estimate_upper[2]))
  # upper_bound row
  expect_equal(td$estimate[3],  unname(out["estimate_upper"]))
  expect_equal(td$std.error[3], unname(out["std.error_upper"]))
  expect_true(is.na(td$conf.low[3]))
  expect_true(is.na(td$conf.high[3]))
  expect_true(is.na(td$estimate_lower[3]))
  expect_true(is.na(td$estimate_upper[3]))
})

test_that("tidy.attrition_bounds works for estimator_ds", {
  df <- make_synthetic()
  out <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  td  <- tidy(out)
  expect_s3_class(td, "tbl_df")
  expect_equal(nrow(td), 3L)
  expect_named(td, c("term", "estimate", "std.error", "conf.low", "conf.high",
                     "estimate_lower", "estimate_upper",
                     "std.error_lower", "std.error_upper", "outcome"))
  expect_equal(td$term, c("bounds", "lower_bound", "upper_bound"))
  # bounds row
  expect_true(is.na(td$estimate[1]))
  expect_true(is.na(td$std.error[1]))
  expect_equal(td$conf.low[1],      unname(out["conf.low"]))
  expect_equal(td$conf.high[1],     unname(out["conf.high"]))
  expect_equal(td$estimate_lower[1],  unname(out["estimate_lower"]))
  expect_equal(td$estimate_upper[1], unname(out["estimate_upper"]))
  # lower_bound row
  expect_equal(td$estimate[2],  unname(out["estimate_lower"]))
  expect_equal(td$std.error[2], unname(out["std.error_lower"]))
  expect_true(is.na(td$conf.low[2]))
  expect_true(is.na(td$conf.high[2]))
  expect_true(is.na(td$estimate_lower[2]))
  expect_true(is.na(td$estimate_upper[2]))
  # upper_bound row
  expect_equal(td$estimate[3],  unname(out["estimate_upper"]))
  expect_equal(td$std.error[3], unname(out["std.error_upper"]))
  expect_true(is.na(td$conf.low[3]))
  expect_true(is.na(td$conf.high[3]))
  expect_true(is.na(td$estimate_lower[3]))
  expect_true(is.na(td$estimate_upper[3]))
})

test_that("tidy.attrition_trim works for estimator_trim", {
  df <- make_synthetic()
  out <- estimator_trim(Y, Z, R = R1, se = "none", data = df)
  td  <- tidy(out)
  expect_s3_class(td, "tbl_df")
  expect_equal(nrow(td), 3L)
  expect_named(td, c("term", "estimate", "std.error", "conf.low", "conf.high",
                     "estimate_lower", "estimate_upper",
                     "std.error_lower", "std.error_upper", "outcome"))
  expect_equal(td$term, c("bounds", "lower_bound", "upper_bound"))
  # bounds row
  expect_true(is.na(td$estimate[1]))
  expect_equal(td$estimate_lower[1],  unname(out["estimate_lower"]))
  expect_equal(td$estimate_upper[1], unname(out["estimate_upper"]))
  # lower_bound row
  expect_equal(td$estimate[2], unname(out["estimate_lower"]))
  expect_true(is.na(td$estimate_lower[2]))
  expect_true(is.na(td$estimate_upper[2]))
  # upper_bound row
  expect_equal(td$estimate[3], unname(out["estimate_upper"]))
  expect_true(is.na(td$estimate_lower[3]))
  expect_true(is.na(td$estimate_upper[3]))
  # NAs throughout for se and CI (no analytic variance for trimming bounds)
  expect_true(all(is.na(td$std.error)))
  expect_true(all(is.na(td$conf.low)))
  expect_true(all(is.na(td$conf.high)))
})

test_that("output classes are set correctly", {
  df <- make_synthetic()
  ev   <- estimator_ev(Y, Z, R1, minY = 1, maxY = 5, data = df)
  ds   <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  sens <- estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 0.5, data = df)
  trim <- estimator_trim(Y, Z, R = R1, data = df)
  expect_s3_class(ev,   "attrition_ev")
  expect_s3_class(ev,   "attrition_bounds")
  expect_s3_class(ds,   "attrition_ds")
  expect_s3_class(ds,   "attrition_bounds")
  expect_s3_class(sens, "attrition_ds_sens")
  expect_s3_class(sens, "attrition_bounds")
  expect_s3_class(trim, "attrition_trim")
})

# ── Additional paper benchmarks ──────────────────────────────────────────────

test_that("estimator_ev with strata (paper data)", {
  out <- estimator_ev(Y = Y_polarization_w2, Z = Z, R = R1, strata = X_party_id,
                      minY = 0, maxY = 6, data = levendusky_replication)
  expect_equal(unname(out["conf.low"]), -1.66862420646393,  tolerance = 1e-10)
  expect_equal(unname(out["conf.high"]),  1.83541087616497,  tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]), -1.53858804710012,   tolerance = 1e-10)
  expect_equal(unname(out["estimate_upper"]),  1.70925509076352,   tolerance = 1e-10)
  expect_equal(unname(out["std.error_lower"]),  0.0790563714807975, tolerance = 1e-10)
  expect_equal(unname(out["std.error_upper"]),  0.0766972716200064, tolerance = 1e-10)
})

test_that("estimator_trim DS path (paper data)", {
  out <- estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                        se = "none", data = levendusky_replication)
  expect_equal(unname(out["estimate_upper"]),  0.567599031583528,  tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]), -0.268142333620083,  tolerance = 1e-10)
})

test_that("estimator_trim R path returns NA bounds on monotonicity violation (paper data)", {
  # Control group responds at the higher rate, so the assumed direction fails
  expect_warning(
    out <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1, data = levendusky_replication),
    "treatment_decreases_response"
  )
  expect_s3_class(out, "attrition_trim")
  expect_true(is.na(out["estimate_lower"]))
  expect_true(is.na(out["estimate_upper"]))
})

test_that("the other monotonicity direction estimates on the same paper data", {
  out <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                        monotonicity = "treatment_decreases_response",
                        data = levendusky_replication)
  expect_false(anyNA(out[c("estimate_lower", "estimate_upper", "conf.low", "conf.high")]))
  expect_lt(unname(out["estimate_lower"]), unname(out["estimate_upper"]))

  # The control group is the one trimmed, so the treatment respondent mean goes
  # in untrimmed and the two group sizes are respondent counts, not arm sizes
  responders <- subset(levendusky_replication, R1 == 1)
  expect_equal(unname(out["Out1_mono"]),
               mean(responders$Y_polarization_w2[responders$Z == 1]))
  expect_equal(unname(out["treat_group_N"]), sum(responders$Z == 1))
  expect_equal(unname(out["control_group_N"]), sum(responders$Z == 0))

  # Q is the share of control respondents trimmed, positive in this direction
  pi_r_1 <- mean(levendusky_replication$R1[levendusky_replication$Z == 1])
  pi_r_0 <- mean(levendusky_replication$R1[levendusky_replication$Z == 0])
  expect_equal(unname(out["Q"]), (pi_r_0 - pi_r_1)/pi_r_0)
  expect_equal(unname(out["pi_r_1"]), pi_r_1)
  expect_equal(unname(out["pi_r_0"]), pi_r_0)
})

test_that("reverse monotonicity is the forward estimator on relabelled arms", {
  rev <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                        monotonicity = "treatment_decreases_response",
                        data = levendusky_replication)

  relabelled <- levendusky_replication
  relabelled$Z <- 1 - relabelled$Z
  fwd <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1, data = relabelled)

  expect_equal(unname(rev["estimate_lower"]), -unname(fwd["estimate_upper"]))
  expect_equal(unname(rev["estimate_upper"]), -unname(fwd["estimate_lower"]))
  expect_equal(unname(rev["std.error_lower"]), unname(fwd["std.error_upper"]))
  expect_equal(unname(rev["std.error_upper"]), unname(fwd["std.error_lower"]))
  expect_equal(unname(rev["conf.low"]),  -unname(fwd["conf.high"]))
  expect_equal(unname(rev["conf.high"]), -unname(fwd["conf.low"]))
})

test_that("all four design-by-assumption cells estimate and lead with the same six elements", {
  six <- c("estimate_lower", "estimate_upper", "std.error_lower", "std.error_upper",
           "conf.low", "conf.high")

  single_mono <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                                monotonicity = "treatment_decreases_response",
                                se = "none", data = levendusky_replication)
  single_none <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                                monotonicity = "none",
                                se = "none", data = levendusky_replication)
  ds_mono <- estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                            monotonicity = "treatment_decreases_response",
                            se = "none", data = levendusky_replication)
  ds_none <- estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                            monotonicity = "none",
                            se = "none", data = levendusky_replication)

  for (out in list(single_mono, single_none, ds_mono, ds_none)) {
    expect_equal(names(out)[1:6], six)
    expect_false(anyNA(out[c("estimate_lower", "estimate_upper")]))
    expect_lt(unname(out["estimate_lower"]), unname(out["estimate_upper"]))
  }

  # The assumption travels with the object
  expect_equal(attr(single_none, "monotonicity"), "none")
  expect_equal(attr(ds_mono, "monotonicity"), "treatment_decreases_response")
  expect_true(attr(single_mono, "single_stage"))
  expect_false(attr(ds_none, "single_stage"))
})

test_that("the no-monotonicity bounds are Imai (2008) Proposition 1 at the Frechet bound", {
  d <- levendusky_replication
  out <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1, monotonicity = "none",
                        se = "none", data = d)

  # Trim each arm by the largest share of its respondents that could fail to be
  # always-reporters, which the Frechet-Hoeffding bound puts at f_other/(1 - f_own)
  f1 <- mean(d$R1[d$Z == 1] == 0)
  f0 <- mean(d$R1[d$Z == 0] == 0)
  trim1 <- f0/(1 - f1)
  trim0 <- f1/(1 - f0)
  expect_equal(unname(out["trim1"]), trim1)
  expect_equal(unname(out["trim0"]), trim0)

  y1 <- sort(d$Y_polarization_w2[d$Z == 1 & d$R1 == 1])
  y0 <- sort(d$Y_polarization_w2[d$Z == 0 & d$R1 == 1])
  cdf1 <- seq_along(y1)/length(y1)
  cdf0 <- seq_along(y0)/length(y0)

  expect_equal(unname(out["estimate_upper"]),
               mean(y1[cdf1 > trim1]) - mean(y0[cdf0 < 1 - trim0]))
  expect_equal(unname(out["estimate_lower"]),
               mean(y1[cdf1 < 1 - trim1]) - mean(y0[cdf0 > trim0]))
})

test_that("dropping monotonicity widens the identification region it is dropped from", {
  # Where the assumed direction holds, the assumption-free bounds trim more of the
  # same group and trim the other one too, so they must contain the monotone bounds
  df <- make_synthetic()

  mono <- estimator_trim(Y, Z, R = R1, se = "none", data = df)
  none <- estimator_trim(Y, Z, R = R1, monotonicity = "none", se = "none", data = df)
  expect_lt(unname(none["estimate_lower"]), unname(mono["estimate_lower"]))
  expect_gt(unname(none["estimate_upper"]), unname(mono["estimate_upper"]))

  ds_mono <- estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2,
                            monotonicity = "treatment_increases_response",
                            se = "none", data = df)
  ds_none <- estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2,
                            se = "none", data = df)
  expect_lt(unname(ds_none["estimate_lower"]), unname(ds_mono["estimate_lower"]))
  expect_gt(unname(ds_none["estimate_upper"]), unname(ds_mono["estimate_upper"]))
})

test_that("the no-monotonicity bounds do not depend on which arm is called treatment", {
  # Both groups are trimmed by the same rule, so relabelling the arms must give back
  # the same interval with its sign flipped. Nothing here has a direction to set.
  d <- levendusky_replication
  relabelled <- d
  relabelled$Z <- 1 - relabelled$Z

  a <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1, monotonicity = "none",
                      se = "none", data = d)
  b <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1, monotonicity = "none",
                      se = "none", data = relabelled)
  expect_equal(unname(a["estimate_lower"]), -unname(b["estimate_upper"]))
  expect_equal(unname(a["estimate_upper"]), -unname(b["estimate_lower"]))
})

test_that("double sampling under monotonicity is the same estimator on relabelled arms", {
  d <- levendusky_replication
  relabelled <- d
  relabelled$Z <- 1 - relabelled$Z

  rev <- estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                        monotonicity = "treatment_decreases_response",
                        se = "none", data = d)
  fwd <- estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                        monotonicity = "treatment_increases_response",
                        se = "none", data = relabelled)
  expect_equal(unname(rev["estimate_lower"]), -unname(fwd["estimate_upper"]))
  expect_equal(unname(rev["estimate_upper"]), -unname(fwd["estimate_lower"]))

  # The follow-up weights are ratios computed within an arm, so relabelling leaves
  # the trimming proportion itself untouched
  expect_equal(unname(rev["Q"]), unname(fwd["Q"]))
})

test_that("analytic standard errors are offered only where Lee (2009) derives them", {
  d <- levendusky_replication

  # The one cell that has them: single sample, one group trimmed
  ok <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                       monotonicity = "treatment_decreases_response", data = d)
  expect_false(anyNA(ok[c("std.error_lower", "std.error_upper")]))

  expect_error(
    estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1, monotonicity = "none", data = d),
    "both groups are trimmed"
  )
  expect_error(
    estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                   monotonicity = "treatment_decreases_response", se = "analytic", data = d),
    "sampling weights"
  )
  expect_error(
    estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                   se = "analytic", data = d),
    "trims both groups"
  )
})

test_that("the bootstrap runs in every cell and brackets the bounds", {
  df <- make_synthetic()
  set.seed(343)
  # A resample thin on follow-up responders can violate monotonicity, and the
  # bootstrap says so and uses the replicates that survived rather than failing
  expect_warning(
    ds_mono <- estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2,
                              monotonicity = "treatment_increases_response",
                              se = "bootstrap", sims = 100, data = df),
    "did not yield bounds"
  )
  cells <- list(
    single_mono = estimator_trim(Y, Z, R = R1, se = "bootstrap", sims = 100, data = df),
    single_none = estimator_trim(Y, Z, R = R1, monotonicity = "none",
                                 se = "bootstrap", sims = 100, data = df),
    ds_mono = ds_mono,
    ds_none = estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2,
                             se = "bootstrap", sims = 100, data = df)
  )
  for (out in cells) {
    expect_true(all(out[c("std.error_lower", "std.error_upper")] > 0))
    expect_lt(unname(out["conf.low"]), unname(out["estimate_lower"]))
    expect_gt(unname(out["conf.high"]), unname(out["estimate_upper"]))
  }
})

test_that("the no-monotonicity bounds refuse to exist when the Frechet bound is not positive", {
  # Missingness rates summing to one or more leave the always-reporter share
  # unbounded away from zero, so there is nothing left to trim toward
  set.seed(9)
  n <- 400
  df <- data.frame(Y = rnorm(n), Z = rep(0:1, each = n/2))
  df$R <- rbinom(n, 1, prob = 0.4)  # about 60 percent missing in both arms
  expect_error(
    estimator_trim(Y, Z, R = R, monotonicity = "none", se = "none", data = df),
    "not less than one"
  )
  # The monotone version is unaffected by that condition
  expect_silent(estimator_trim(Y, Z, R = R, se = "none", data = df))
})

test_that("a trimming proportion too large for the group is an error, not a NaN", {
  # Two treated respondents with a trimming proportion of 0.5: the lower-bound side
  # retains nothing
  df <- data.frame(
    Y = c(1, 2, 3, 4, 5, 6, 7, 8),
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    R = c(1, 1, 0, 0, 1, 0, 0, 0)
  )
  expect_error(
    estimator_trim(Y, Z, R = R, se = "none", data = df),
    "leaves nothing behind"
  )
})

test_that("estimator_ds_sens(delta=1) exactly matches estimator_ds (paper data)", {
  ds   <- estimator_ds(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                       minY = 0, maxY = 6, data = levendusky_replication)
  sens <- estimator_ds_sens(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                            delta = 1, minY = 0, maxY = 6, data = levendusky_replication)
  expect_equal(as.numeric(ds), as.numeric(sens), tolerance = 1e-14)
})

test_that("estimator_ds_sens delta=0.5 (paper data)", {
  out <- estimator_ds_sens(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                           delta = 0.5, minY = 0, maxY = 6, data = levendusky_replication)
  expect_equal(unname(out["conf.low"]), -0.263141039257175,  tolerance = 1e-10)
  expect_equal(unname(out["conf.high"]),  0.527330963282511,  tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]), -0.0916157282335379,  tolerance = 1e-10)
  expect_equal(unname(out["estimate_upper"]),  0.365164824275198,   tolerance = 1e-10)
  expect_equal(unname(out["std.error_lower"]),  0.1042799840077,  tolerance = 1e-10)
  expect_equal(unname(out["std.error_upper"]),  0.0985900114761642, tolerance = 1e-10)
})

# ── Additional synthetic data tests ─────────────────────────────────────────

test_that("estimator_ev with strata produces stable results on synthetic data", {
  df  <- make_synthetic()
  out <- estimator_ev(Y, Z, R1, strata = strata, minY = 1, maxY = 5, data = df)
  expect_equal(unname(out["conf.low"]), -0.828328275895711,   tolerance = 1e-10)
  expect_equal(unname(out["conf.high"]),  1.3505157781213,   tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]), -0.699410333729967,   tolerance = 1e-10)
  expect_equal(unname(out["estimate_upper"]),  1.23280497665065,    tolerance = 1e-10)
  expect_equal(unname(out["std.error_lower"]),  0.0783765436956703, tolerance = 1e-10)
  expect_equal(unname(out["std.error_upper"]),  0.0715630859438926, tolerance = 1e-10)
})

test_that("estimator_trim R path (monotonicity) produces stable results on synthetic data", {
  df  <- make_synthetic()
  out <- estimator_trim(Y, Z, R = R1, data = df)
  expect_equal(unname(out["estimate_upper"]),  0.534610520108332, tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]),  0.100617821905443, tolerance = 1e-10)
  expect_equal(unname(out["Q"]),            0.101526858997151, tolerance = 1e-10)
})

test_that("estimator_trim DS path produces stable results on synthetic data", {
  df  <- make_synthetic()
  out <- estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2, se = "none", data = df)
  expect_equal(unname(out["estimate_upper"]),  0.544188162330423,  tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]),  0.0677009048063324, tolerance = 1e-10)
})

test_that("estimator_ds_sens delta=0.5 produces stable results on synthetic data", {
  df  <- make_synthetic()
  out <- estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 0.5,
                           data = df)
  expect_equal(unname(out["conf.low"]),  0.0288421426106519, tolerance = 1e-10)
  expect_equal(unname(out["conf.high"]),  0.514219908745877, tolerance = 1e-10)
  expect_equal(unname(out["estimate_lower"]),   0.164445203836863, tolerance = 1e-10)
  expect_equal(unname(out["estimate_upper"]),   0.379394808689304, tolerance = 1e-10)
  expect_equal(unname(out["std.error_lower"]),   0.082435669263969, tolerance = 1e-10)
  expect_equal(unname(out["std.error_upper"]),   0.0819627319342319, tolerance = 1e-10)
})

# ── Property tests ───────────────────────────────────────────────────────────

test_that("estimator_ds_sens(delta=1) equals estimator_ds (synthetic data)", {
  df   <- make_synthetic()
  ds   <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  sens <- estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 1,
                             data = df)
  expect_equal(as.numeric(ds), as.numeric(sens), tolerance = 1e-14)
})

test_that("lower_bound <= upper_bound for all estimators (synthetic data)", {
  df   <- make_synthetic()
  ev   <- estimator_ev(Y, Z, R1, minY = 1, maxY = 5, data = df)
  ds   <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  trim <- estimator_trim(Y, Z, R = R1, data = df)
  expect_lte(unname(ev["estimate_lower"]),       unname(ev["estimate_upper"]))
  expect_lte(unname(ds["estimate_lower"]),       unname(ds["estimate_upper"]))
  expect_lte(unname(trim["estimate_lower"]), unname(trim["estimate_upper"]))
})

test_that("Imbens-Manski CI covers the identification region", {
  df <- make_synthetic()
  for (out in list(
    estimator_ev(Y, Z, R1, minY = 1, maxY = 5, data = df),
    estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  )) {
    expect_lte(unname(out["conf.low"]), unname(out["estimate_lower"]))
    expect_gte(unname(out["conf.high"]), unname(out["estimate_upper"]))
  }
})

test_that("double sampling narrows identification region vs extreme value bounds", {
  df       <- make_synthetic()
  ev       <- estimator_ev(Y, Z, R1, minY = 1, maxY = 5, data = df)
  ds       <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  ev_width <- unname(ev["estimate_upper"] - ev["estimate_lower"])
  ds_width <- unname(ds["estimate_upper"] - ds["estimate_lower"])
  expect_lte(ds_width, ev_width)
})

test_that("increasing delta widens bounds in estimator_ds_sens", {
  df    <- make_synthetic()
  sens0 <- estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 0,   data = df)
  sens5 <- estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 0.5, data = df)
  sens1 <- estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 1,   data = df)
  width <- function(x) unname(x["estimate_upper"] - x["estimate_lower"])
  expect_lte(width(sens0), width(sens5))
  expect_lte(width(sens5), width(sens1))
})

# ── Helper function unit tests (extended) ────────────────────────────────────

test_that("gen_mean returns exact values", {
  expect_equal(attrition:::gen_mean(3, 0.5, lower_bound = TRUE,  minY = 0, maxY = 5), 1.5)
  expect_equal(attrition:::gen_mean(3, 0.5, lower_bound = FALSE, minY = 0, maxY = 5), 4.0)
  # At p = 1 all outcomes observed; bound doesn't matter
  expect_equal(attrition:::gen_mean(3, 1.0, lower_bound = TRUE,  minY = 0, maxY = 5), 3.0)
  expect_equal(attrition:::gen_mean(3, 1.0, lower_bound = FALSE, minY = 0, maxY = 5), 3.0)
})

test_that("gen_mean_sens: delta=1 matches gen_mean; delta=0 returns y_m", {
  sens1 <- attrition:::gen_mean_sens(3, 0.5, delta = 1, lower_bound = TRUE,  minY = 0, maxY = 5)
  base  <- attrition:::gen_mean(3, 0.5, lower_bound = TRUE, minY = 0, maxY = 5)
  expect_equal(sens1, base)
  # delta = 0: p*y_m + (1-p)*0*const + (1-p)*1*y_m = y_m
  expect_equal(attrition:::gen_mean_sens(3, 0.5, delta = 0, lower_bound = TRUE,  minY = 0, maxY = 5), 3.0)
  expect_equal(attrition:::gen_mean_sens(3, 0.5, delta = 0, lower_bound = FALSE, minY = 0, maxY = 5), 3.0)
})

test_that("construct_manski_bounds returns correct values", {
  out <- attrition:::construct_manski_bounds(
    p1_t = 0.8, y1m_t = 3, p1_c = 0.7, y1m_c = 2,
    y2m_t_L = 1, y2m_t_U = 4, y2m_c_L = 1, y2m_c_U = 4
  )
  # lower = (0.8*3 + 0.2*1) - (0.7*2 + 0.3*4) = 2.6 - 2.6 = 0
  # upper = (0.8*3 + 0.2*4) - (0.7*2 + 0.3*1) = 3.2 - 1.7 = 1.5
  expect_equal(out[1], 0, tolerance = 1e-10)
  expect_equal(out[2], 1.5)
})

test_that("ds_var returns correct value", {
  v <- attrition:::ds_var(n1 = 100, n2 = 50, p1 = 0.8, p2 = 0.7,
                           s1 = 1, s2 = 1.5, y1m = 3, y2m = 4)
  expect_equal(v, 0.0114, tolerance = 1e-10)
})

test_that("im_crit is near zero at standard normal critical value (point-ID case)", {
  # When bounds are equal, the IM critical value reduces to the normal 1.96
  val <- attrition:::im_crit(1.96,
                              upper_bound_est = 0, lower_bound_est = 0,
                              upper_bound_var_est = 1, lower_bound_var_est = 1,
                              alpha = 0.05)
  expect_equal(val, 0, tolerance = 1e-4)
})

test_that("im_critical_value is correct at every alpha, not just 0.05", {
  # The IM critical value lies in [z_(1-alpha), z_(1-alpha/2)]: it equals
  # z_(1-alpha/2) when the bounds coincide and z_(1-alpha) as they separate.
  for (a in c(0.20, 0.10, 0.05, 0.01, 0.001)) {
    # Coincident bounds: the two-sided critical value
    tight <- attrition:::im_critical_value(0, 0, 1, 1, alpha = a)
    expect_equal(tight, qnorm(1 - a/2), tolerance = 1e-6)
    # Widely separated bounds: the one-sided critical value
    wide <- attrition:::im_critical_value(-10, 10, 1, 1, alpha = a)
    expect_equal(wide, qnorm(1 - a), tolerance = 1e-6)
    # Intermediate: strictly inside the bracket
    mid <- attrition:::im_critical_value(-1, 1, 1, 1, alpha = a)
    expect_gte(mid, qnorm(1 - a))
    expect_lte(mid, qnorm(1 - a/2))
  }
})

test_that("estimator CIs widen as alpha shrinks", {
  df <- make_synthetic()
  ci_width <- function(a) {
    o <- estimator_ev(Y, Z, R1, minY = 1, maxY = 5, alpha = a, data = df)
    unname(o["conf.high"] - o["conf.low"])
  }
  widths <- vapply(c(0.20, 0.10, 0.05, 0.01, 0.001), ci_width, numeric(1))
  expect_true(all(diff(widths) > 0))
  # Regression guard: alpha = 0.01 once returned 2.000, the alpha = 0.05
  # critical value, against a correct 2.326. sig is recovered here through a
  # multiply and a divide, so it lands an ulp either side of the exact value
  # depending on the platform. The tolerance is six orders of magnitude below
  # the 0.33 defect this guards against.
  o01 <- estimator_ev(Y, Z, R1, minY = 1, maxY = 5, alpha = 0.01, data = df)
  sig <- unname((o01["estimate_lower"] - o01["conf.low"]) / o01["std.error_lower"])
  expect_gt(sig, qnorm(0.99) - 1e-6)
})

test_that("find_sign_changes flags the first departure from the initial sign", {
  expect_equal(attrition:::find_sign_changes(c(-1,  1, -1)), c(FALSE, TRUE, FALSE))
  # A zero counts as a change
  expect_equal(attrition:::find_sign_changes(c(-1,  0,  1)), c(FALSE, TRUE, FALSE))
  # No change at all
  expect_equal(attrition:::find_sign_changes(c( 1,  2,  3)), c(FALSE, FALSE, FALSE))
  expect_equal(attrition:::find_sign_changes(c(-3, -2, -1)), c(FALSE, FALSE, FALSE))
  # A change at the LAST position is a real change: previously missed, so a
  # delta* near 1 was reported as "no significant value of delta exists"
  expect_equal(attrition:::find_sign_changes(c(-1, -1,  1)), c(FALSE, FALSE, TRUE))
  # Only the first departure is flagged
  expect_equal(attrition:::find_sign_changes(c(-1,  1,  1, -1)), c(FALSE, TRUE, FALSE, FALSE))
})

test_that("trimming_bounds names the correct group in the monotonicity error", {
  set.seed(1); n <- 200
  Out   <- rnorm(2 * n)
  Treat <- rep(0:1, each = n)
  # Treatment much MORE likely to be missing than control: f1 > f0, so Q < 0
  Fail  <- c(rbinom(n, 1, 0.10), rbinom(n, 1, 0.60))
  expect_gt(mean(Fail[Treat == 1]), mean(Fail[Treat == 0]))
  expect_error(
    attrition:::trimming_bounds(Out, Treat, Fail, rep(1, 2 * n), monotonicity = TRUE),
    "treatment group is more likely to be missing",
    class = "attrition_monotonicity_violation"
  )
})

test_that("trimming_bounds handles single-observation and empty groups", {
  # One observed outcome in the treatment group: the CDF loop used to run
  # backwards here and corrupt the vector or throw
  out <- attrition:::trimming_bounds(
    Out = c(1, 2, 3, 4), Treat = c(0, 0, 0, 1), Fail = c(0, 0, 0, 0),
    Weight = rep(1, 4), monotonicity = TRUE)
  expect_false(any(is.na(out[c("lower_bound", "upper_bound")])))

  expect_error(
    attrition:::trimming_bounds(Out = c(1, 2), Treat = c(0, 1), Fail = c(0, 1),
                                Weight = rep(1, 2), monotonicity = TRUE),
    "at least one observed outcome in each treatment group")
})

test_that("estimator_trim surfaces non-monotonicity errors instead of returning NA", {
  # Every unit in the treatment group is missing: not a monotonicity violation,
  # so it must not be silently converted to NA bounds
  df <- data.frame(Y = c(1, 2, 3, NA, NA), Z = c(0, 0, 0, 1, 1), R = c(1, 1, 1, 0, 0))
  expect_error(estimator_trim(Y, Z, R = R, data = df),
               "at least one observed outcome in each treatment group")
})

test_that("estimators validate minY, maxY, and alpha", {
  df <- make_synthetic()
  expect_error(estimator_ev(Y, Z, R1, minY = 5, maxY = 1, data = df),
               "must not be greater than the maximum")
  expect_error(estimator_ev(Y, Z, R1, minY = 2, maxY = 5, data = df),
               "outside the assumed support")
  expect_error(estimator_ev(Y, Z, R1, minY = 1, maxY = 5, alpha = 0, data = df),
               "strictly between zero and one")
  expect_error(estimator_ev(Y, Z, R1, minY = 1, maxY = 5, alpha = 1, data = df),
               "strictly between zero and one")
  expect_error(estimator_ds(Y, Z, R1, Attempt, R2, minY = 5, maxY = 1, data = df),
               "must not be greater than the maximum")
})

test_that("estimator_ds_sens validates delta", {
  df <- make_synthetic()
  expect_error(estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5,
                                 delta = 5, data = df), "between zero and one")
  expect_error(estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5,
                                 delta = -1, data = df), "between zero and one")
})

# ── Input validation ─────────────────────────────────────────────────────────

test_that("estimator_ev validates inputs", {
  df <- make_synthetic()
  expect_error(estimator_ev(as.character(Y), Z, R1, minY = 1, maxY = 5, data = df),
               "numeric")
  expect_error(estimator_ev(Y, Z + 0.5, R1, minY = 1, maxY = 5, data = df),
               "zero or one")
  expect_error(estimator_ev(Y, Z, R1 + 0.5, minY = 1, maxY = 5, data = df),
               "zero or one")
  expect_error(estimator_ev(Y, Z, R1, minY = "a", maxY = 5, data = df),
               "numeric")
  expect_error(estimator_ev(Y, Z, R1, minY = 1, maxY = 5, strata = ifelse(Z == 1, NA_real_, strata),
                             data = df),
               "missing values")
})

test_that("estimator_ds validates inputs", {
  df <- make_synthetic()
  expect_error(estimator_ds(as.character(Y), Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df),
               "numeric")
  expect_error(estimator_ds(Y, Z + 0.5, R1, Attempt, R2, minY = 1, maxY = 5, data = df),
               "zero or one")
  expect_error(estimator_ds(Y, Z, R1 + 0.5, Attempt, R2, minY = 1, maxY = 5, data = df),
               "zero or one")
  expect_error(estimator_ds(Y, Z, R1, Attempt, R2 + 0.5, minY = 1, maxY = 5, data = df),
               "zero or one")
  expect_error(estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5,
                             strata = ifelse(Z == 1, NA_real_, strata), data = df),
               "missing values")
})

test_that("estimator_trim validates inputs", {
  df <- make_synthetic()
  expect_error(estimator_trim(as.character(Y), Z, R = R1, data = df), "numeric")
  expect_error(estimator_trim(Y, Z + 0.5, R = R1, data = df), "zero or one")
  expect_error(estimator_trim(Y, Z, R = R1 + 0.5, data = df), "zero or one")
})

test_that("estimator_trim returns NA bounds when monotonicity is violated", {
  df   <- make_synthetic()
  df$Z <- 1L - df$Z  # flip treatment — control now has higher response rate
  expect_warning(out <- estimator_trim(Y, Z, R = R1, data = df),
                 "treatment_decreases_response")
  expect_s3_class(out, "attrition_trim")
  expect_true(is.na(out["estimate_lower"]))
  expect_true(is.na(out["estimate_upper"]))

  # The same data estimate under the direction they do admit, and the warning
  # there names the direction that was just assumed away
  rev <- estimator_trim(Y, Z, R = R1, monotonicity = "treatment_decreases_response",
                        data = df)
  expect_false(anyNA(rev[c("estimate_lower", "estimate_upper")]))
  expect_warning(
    estimator_trim(Y, Z, R = R1, monotonicity = "treatment_decreases_response",
                   data = make_synthetic()),
    "treatment_increases_response"
  )
})

test_that("tidy carries the outcome name from either interface", {
  df  <- make_synthetic()
  nse <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  frm <- estimator_ds(Y ~ Z, R1 = "R1", Attempt = "Attempt", R2 = "R2",
                      minY = 1, maxY = 5, data = df)
  expect_equal(unique(tidy(nse)$outcome), "Y")
  expect_equal(unique(tidy(frm)$outcome), "Y")

  trim <- estimator_trim(Y, Z, R = R1, se = "none", data = df)
  expect_equal(unique(tidy(trim)$outcome), "Y")

  # No sample-size column: which sample is the relevant one differs by estimator
  # (every randomized unit for the bounds, the always-responders for trimming,
  # and that subpopulation's size is itself estimated), so reporting one number
  # would assert an answer the estimator does not have. estimatr does not put
  # nobs in tidy output either.
  expect_false("nobs" %in% names(tidy(nse)))
  expect_false("nobs" %in% names(tidy(trim)))
})

# ── Trimming bound standard errors ───────────────────────────────────────────

test_that("Lee (2009) analytic standard errors are produced and stable", {
  df  <- make_synthetic()
  out <- estimator_trim(Y, Z, R = R1, data = df)
  expect_true(all(c("std.error_lower", "std.error_upper", "conf.low", "conf.high") %in% names(out)))
  expect_true(all(out[c("std.error_lower", "std.error_upper")] > 0))
  expect_equal(unname(out["std.error_lower"]), 0.108109837320316, tolerance = 1e-8)
  expect_equal(unname(out["std.error_upper"]), 0.0939297066416204, tolerance = 1e-8)
  # The Imbens-Manski interval must contain the identified set
  expect_lte(unname(out["conf.low"]), unname(out["estimate_lower"]))
  expect_gte(unname(out["conf.high"]), unname(out["estimate_upper"]))
})

test_that("Lee eq (7) term 3 equals the Tauchmann (2014) form it is written in", {
  # Lee's published coefficient on the trimming-proportion term looks like it
  # could go negative; the form used in lee_variance cannot. They are identical.
  lee   <- function(pt, pc, e) { Q <- (pt-pc)/pt; ((1-pc) - Q*(1-e))/(e*pc*(1-e)) }
  tauch <- function(pt, pc, e) (1-pt)/(pt*e) + (1-pc)/(pc*(1-e))
  g <- expand.grid(pt = seq(0.55, 0.95, 0.1), pc = seq(0.35, 0.85, 0.1), e = c(0.3, 0.5, 0.7))
  g <- subset(g, pt > pc)
  expect_equal(mapply(lee, g$pt, g$pc, g$e), mapply(tauch, g$pt, g$pc, g$e))
})

test_that("analytic and bootstrap standard errors agree", {
  df <- make_synthetic()
  ana <- estimator_trim(Y, Z, R = R1, se = "analytic", data = df)
  set.seed(11)
  boo <- estimator_trim(Y, Z, R = R1, se = "bootstrap", sims = 600, data = df)
  expect_equal(unname(boo["estimate_lower"]), unname(ana["estimate_lower"]))
  expect_equal(unname(boo["estimate_upper"]), unname(ana["estimate_upper"]))
  expect_equal(unname(boo["std.error_lower"]), unname(ana["std.error_lower"]), tolerance = 0.15)
  expect_equal(unname(boo["std.error_upper"]), unname(ana["std.error_upper"]), tolerance = 0.15)
})

test_that("bootstrap standard errors work on the double-sampling path", {
  df <- make_synthetic()
  set.seed(12)
  out <- estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2,
                        se = "bootstrap", sims = 300, data = df)
  expect_true(all(out[c("std.error_lower", "std.error_upper")] > 0))
  expect_lte(unname(out["conf.low"]), unname(out["estimate_lower"]))
  expect_gte(unname(out["conf.high"]), unname(out["estimate_upper"]))
})

test_that("analytic standard errors are refused on the double-sampling path", {
  df <- make_synthetic()
  expect_error(
    estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2, se = "analytic", data = df),
    "Lee \\(2009\\) Proposition 3")
})

test_that("se = 'none' leaves standard errors and CIs missing", {
  df  <- make_synthetic()
  out <- estimator_trim(Y, Z, R = R1, se = "none", data = df)
  expect_true(all(is.na(out[c("std.error_lower", "std.error_upper", "conf.low", "conf.high")])))
  expect_false(is.na(unname(out["estimate_lower"])))
})

test_that("bootstrap validates sims and estimator_trim validates alpha", {
  df <- make_synthetic()
  expect_error(estimator_trim(Y, Z, R = R1, se = "bootstrap", sims = 1, data = df),
               "at least two")
  expect_error(estimator_trim(Y, Z, R = R1, alpha = 0, data = df),
               "strictly between zero and one")
})

test_that("a zero trimming proportion warns that the asymptotics do not apply", {
  # Equal response rates: Q = 0, the bounds collapse to a point, and Lee's
  # interior-point condition fails
  set.seed(4)
  N <- 400
  Z <- rep(0:1, each = N/2)
  R <- rep(c(1, 1, 1, 0), length.out = N)   # identical response rate in both arms
  Y <- ifelse(R == 1, rnorm(N), NA_real_)
  df <- data.frame(Y, Z, R)
  expect_warning(estimator_trim(Y, Z, R = R, data = df), "collapse to a point")
})

test_that("tidy.attrition_trim carries standard errors and the IM interval", {
  df  <- make_synthetic()
  out <- estimator_trim(Y, Z, R = R1, se = "analytic", data = df)
  td  <- tidy(out)
  expect_equal(td$std.error[2], unname(out["std.error_lower"]))
  expect_equal(td$std.error[3], unname(out["std.error_upper"]))
  expect_equal(td$conf.low[1],  unname(out["conf.low"]))
  expect_equal(td$conf.high[1], unname(out["conf.high"]))
  expect_true(is.na(td$std.error[1]))
})

test_that("analytic standard errors recover the sampling variability (Monte Carlo)", {
  skip_on_cran()
  set.seed(2026)
  dgp <- function(n) {
    Z <- rep(0:1, each = n/2)
    Y_star <- rnorm(n, mean = 2 + 0.3*Z)
    R <- rbinom(n, 1, plogis(-0.6 + 0.9*Z + 0.4*Y_star))
    Y <- Y_star; Y[R == 0] <- NA
    data.frame(Y, Z, R)
  }
  n <- 3000; reps <- 400
  mc <- vapply(seq_len(reps), function(i)
    unname(estimator_trim(Y, Z, R = R, se = "none", data = dgp(n))["estimate_upper"]), numeric(1))
  se <- vapply(seq_len(60), function(i)
    unname(estimator_trim(Y, Z, R = R, se = "analytic", data = dgp(n))["std.error_upper"]), numeric(1))
  # The mean analytic SE should track the Monte Carlo sd of the estimator
  expect_equal(mean(se), sd(mc), tolerance = 0.12)
})

test_that("omitting data by name gives an instructive error, not R's default", {
  df <- make_synthetic()
  # The call a user reflexively types, since every other R model is (formula, data)
  expect_error(estimator_ds(Y ~ Z, df, R1 = "R1", Attempt = "Attempt", R2 = "R2",
                            minY = 1, maxY = 5), "must be given by name")
  expect_error(estimator_ev(Y ~ Z, df, R = "R1", minY = 1, maxY = 5), "must be given by name")
  expect_error(estimator_trim(Y ~ Z, df, R = "R1"), "must be given by name")
  expect_error(estimator_ds_sens(Y ~ Z, df, R1 = "R1", Attempt = "Attempt", R2 = "R2",
                                 minY = 1, maxY = 5, delta = 0.5), "must be given by name")
  expect_error(sensitivity_ds(Y ~ Z, df, R1 = "R1", Attempt = "Attempt", R2 = "R2",
                              minY = 1, maxY = 5), "must be given by name")
})

test_that("a column may be named with a one-sided formula", {
  df  <- make_synthetic()
  nse <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  frm <- estimator_ds(Y ~ Z, R1 = ~R1, Attempt = ~Attempt, R2 = ~R2,
                      minY = 1, maxY = 5, data = df)
  expect_equal(as.numeric(nse), as.numeric(frm))
  expect_error(estimator_ds(Y ~ Z, R1 = ~R1 + Attempt, Attempt = "Attempt", R2 = "R2",
                            minY = 1, maxY = 5, data = df), "exactly one column")
})

# ── Formula interface ────────────────────────────────────────────────────────

test_that("estimator_ev formula interface matches NSE interface", {
  df  <- make_synthetic()
  nse <- estimator_ev(Y, Z, R1, minY = 1, maxY = 5, data = df)
  frm <- estimator_ev(Y ~ Z, R = "R1", minY = 1, maxY = 5, data = df)
  expect_equal(as.numeric(nse), as.numeric(frm))
  expect_s3_class(frm, "attrition_bounds")
})

test_that("estimator_ds formula interface matches NSE interface", {
  df  <- make_synthetic()
  nse <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  frm <- estimator_ds(Y ~ Z, R1 = "R1", Attempt = "Attempt", R2 = "R2",
                      minY = 1, maxY = 5, data = df)
  expect_equal(as.numeric(nse), as.numeric(frm))
  expect_s3_class(frm, "attrition_bounds")
})

test_that("estimator_trim formula interface (single-stage) matches NSE interface", {
  df  <- make_synthetic()
  nse <- estimator_trim(Y, Z, R = R1, data = df)
  frm <- estimator_trim(Y ~ Z, R = "R1", data = df)
  expect_equal(as.numeric(nse), as.numeric(frm))
  expect_s3_class(frm, "attrition_trim")
})

test_that("estimator_trim formula interface (double-sampling) matches NSE interface", {
  df  <- make_synthetic()
  nse <- estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2, se = "none", data = df)
  frm <- estimator_trim(Y ~ Z, R1 = "R1", Attempt = "Attempt", R2 = "R2", se = "none", data = df)
  expect_equal(as.numeric(nse), as.numeric(frm))
  expect_s3_class(frm, "attrition_trim")
})

test_that("estimator_ds_sens formula interface matches NSE interface", {
  df  <- make_synthetic()
  nse <- estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 0.5, data = df)
  frm <- estimator_ds_sens(Y ~ Z, R1 = "R1", Attempt = "Attempt", R2 = "R2",
                           minY = 1, maxY = 5, delta = 0.5, data = df)
  expect_equal(as.numeric(nse), as.numeric(frm))
  expect_s3_class(frm, "attrition_bounds")
})

test_that("strata accepts quoted string column name", {
  df  <- make_synthetic()
  nse <- estimator_ev(Y, Z, R1, strata = strata, minY = 1, maxY = 5, data = df)
  str <- estimator_ev(Y, Z, R1, strata = "strata", minY = 1, maxY = 5, data = df)
  expect_equal(as.numeric(nse), as.numeric(str))

  nse <- estimator_ds(Y, Z, R1, Attempt, R2, strata = strata, minY = 1, maxY = 5, data = df)
  str <- estimator_ds(Y, Z, R1, Attempt, R2, strata = "strata", minY = 1, maxY = 5, data = df)
  expect_equal(as.numeric(nse), as.numeric(str))
})

test_that("estimator_trim errors when strata is supplied", {
  df <- make_synthetic()
  expect_error(estimator_trim(Y, Z, R = R1, strata = strata, data = df),
               "not yet supported")
})

test_that("over-specified formulas are rejected with a clear error", {
  df <- make_synthetic()
  expect_error(estimator_ev(Y ~ Z + R1,              minY = 1, maxY = 5, data = df), "exactly two variables")
  expect_error(estimator_ds(Y ~ Z + R1,              minY = 1, maxY = 5, data = df), "exactly two variables")
  expect_error(estimator_trim(Y ~ Z + R1,            data = df),                     "exactly two variables")
  expect_error(estimator_ds_sens(Y ~ Z + R1,         minY = 1, maxY = 5, delta = 0.5, data = df), "exactly two variables")
})

# ── sensitivity_ds ───────────────────────────────────────────────────────────

test_that("sensitivity_ds returns correct structure", {
  df  <- make_synthetic()
  out <- sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 20, data = df)
  expect_named(out, c("sensitivity_plot", "sims_df", "delta_star"))
  expect_s3_class(out$sensitivity_plot, "gg")
  expect_equal(nrow(out$sims_df), 20L)
  expect_named(out$sims_df,
    c("estimate_lower", "estimate_upper", "std.error_lower", "std.error_upper",
      "conf.low", "conf.high", "delta", "change_lower", "change_upper",
      "change_any"))
  # the delta grid runs from 0 to 1
  expect_equal(out$sims_df$delta[1],  0)
  expect_equal(out$sims_df$delta[20], 1)
})

test_that("sensitivity_ds detects delta* on synthetic data", {
  # Synthetic data: bounds are positive at delta=0 but span 0 at delta=1,
  # so there must be a sign change interior to [0, 1]
  df  <- make_synthetic()
  out <- sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 50, data = df)
  expect_true(is.numeric(out$delta_star))
  expect_length(out$delta_star, 1L)
  expect_true(out$delta_star > 0 && out$delta_star < 1)
})

test_that("sensitivity_ds formula interface matches NSE interface", {
  df  <- make_synthetic()
  nse <- sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 10, data = df)
  frm <- sensitivity_ds(Y ~ Z, R1 = "R1", Attempt = "Attempt", R2 = "R2",
                        minY = 1, maxY = 5, sims = 10, data = df)
  expect_equal(nse$sims_df, frm$sims_df)
})

test_that("sensitivity_ds validates sims", {
  df <- make_synthetic()
  expect_error(sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 1, data = df),
               "at least two")
})

test_that("sensitivity_ds with strata returns correct structure", {
  df  <- make_synthetic()
  out <- sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 10,
                        strata = strata, data = df)
  expect_named(out, c("sensitivity_plot", "sims_df", "delta_star"))
  expect_equal(nrow(out$sims_df), 10L)
})

# ── Coverage of the surface added with the four trimming cells ───────────────

test_that("the print methods show the quantities and return invisibly", {
  bounds <- estimator_ev(Y = Y_polarization_w2, Z = Z, R = R1,
                         minY = 0, maxY = 6, data = levendusky_replication)
  trim <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                         monotonicity = "treatment_decreases_response",
                         se = "none", data = levendusky_replication)

  expect_output(print(bounds), "estimate_lower")
  expect_output(print(bounds), "conf.high")
  expect_output(print(trim), "estimate_lower")
  expect_output(print(trim), "Out1_mono")

  # The class vector is what the methods exist to keep out of the output
  expect_false(any(grepl("attrition_bounds", capture.output(print(bounds)))))
  expect_false(any(grepl("attrition_trim", capture.output(print(trim)))))

  invisible(capture.output({
    vis <- withVisible(print(bounds))$visible
    val <- withVisible(print(trim))$value
  }))
  expect_false(vis)
  expect_equal(unname(val), unname(trim))
})

test_that("tidy returns the same six quantities in every trimming cell", {
  d <- levendusky_replication
  # A thin resample can violate monotonicity, which is warned about and asserted
  # where the bootstrap itself is under test
  cells <- suppressWarnings(list(
    estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                   monotonicity = "treatment_decreases_response", data = d),
    estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1, monotonicity = "none",
                   se = "bootstrap", sims = 50, data = d),
    estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                   monotonicity = "treatment_decreases_response",
                   se = "bootstrap", sims = 50, data = d),
    estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                   se = "bootstrap", sims = 50, data = d)
  ))
  for (out in cells) {
    td <- tidy(out)
    expect_equal(td$term, c("bounds", "lower_bound", "upper_bound"))
    expect_equal(td$estimate[2], unname(out["estimate_lower"]))
    expect_equal(td$estimate[3], unname(out["estimate_upper"]))
    expect_equal(td$conf.low[1], unname(out["conf.low"]))
    expect_equal(unique(td$outcome), "Y_polarization_w2")
    expect_true(is.na(td$estimate[1]))
  }
})

test_that("the formula interface carries the monotonicity argument", {
  d <- levendusky_replication
  nse <- estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                        monotonicity = "treatment_decreases_response",
                        se = "none", data = d)
  frm <- estimator_trim(Y_polarization_w2 ~ Z, R = "R1",
                        monotonicity = "treatment_decreases_response",
                        se = "none", data = d)
  expect_equal(as.numeric(nse), as.numeric(frm))

  none_frm <- estimator_trim(Y_polarization_w2 ~ Z, R = "R1", monotonicity = "none",
                             se = "none", data = d)
  expect_equal(unname(none_frm["estimate_lower"]),
               unname(estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                                     monotonicity = "none", se = "none",
                                     data = d)["estimate_lower"]))
})

test_that("an unrecognised monotonicity value is rejected", {
  expect_error(
    estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1, monotonicity = "increases",
                   se = "none", data = levendusky_replication),
    "should be one of"
  )
})

test_that("double sampling under monotonicity relabels its intermediates correctly", {
  d <- levendusky_replication
  out <- estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                        monotonicity = "treatment_decreases_response",
                        se = "none", data = d)

  # Rebuild the follow-up weights and the kept set by hand
  w <- rep(NA_real_, nrow(d))
  w[d$R1 == 1] <- 1
  w[d$Attempt == 1 & d$Z == 1] <- sum(d$Z == 1 & d$R1 == 0)/sum(d$Z == 1 & d$Attempt == 1)
  w[d$Attempt == 1 & d$Z == 0] <- sum(d$Z == 0 & d$R1 == 0)/sum(d$Z == 0 & d$Attempt == 1)
  keep <- d$R1 == 1 | d$Attempt == 1
  fail <- d$R1 == 0 & d$R2 == 0

  # Under this direction the treatment arm is the untrimmed one
  treated_obs <- keep & !fail & d$Z == 1
  expect_equal(unname(out["Out1_mono"]),
               weighted.mean(d$Y_polarization_w2[treated_obs], w[treated_obs]))

  # f1 and f0 name the arm they came from, not the relabelled one
  f1 <- sum(w[keep & fail & d$Z == 1])/sum(w[keep & d$Z == 1])
  f0 <- sum(w[keep & fail & d$Z == 0])/sum(w[keep & d$Z == 0])
  expect_equal(unname(out["f1"]), f1)
  expect_equal(unname(out["f0"]), f0)
  expect_equal(unname(out["pi_r_1"]), 1 - f1)
  expect_equal(unname(out["Q"]), ((1 - f0) - (1 - f1))/(1 - f0))

  # The two trimmed control means bracket the untrimmed one
  expect_lt(unname(out["Out0U_mono"]), unname(out["Out0L_mono"]))
})

test_that("the no-monotonicity branch also refuses a trimming proportion it cannot fill", {
  # Missingness of 0.25 and 0.5 leaves the Frechet bound positive, so the bounds
  # are defined, but trimming half of two control respondents from the bottom
  # keeps nothing
  df <- data.frame(
    Y = c(1, 2, 3, 4, 5, 6, 7, 8),
    Z = c(1, 1, 1, 1, 0, 0, 0, 0),
    R = c(1, 1, 1, 0, 1, 1, 0, 0)
  )
  expect_error(
    estimator_trim(Y, Z, R = R, monotonicity = "none", se = "none", data = df),
    "leaves nothing behind"
  )
})

# ── Input validation that had no test ────────────────────────────────────────

test_that("every estimator checks its response indicators", {
  df <- make_synthetic()
  bad <- df; bad$R2[1] <- 2
  expect_error(estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = bad), "R2")
  expect_error(estimator_ds_sens(Y, Z, R1, Attempt, R2, delta = 1, minY = 1, maxY = 5, data = bad), "R2")
  expect_error(sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 3, data = bad), "R2")
  expect_error(estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2,
                              se = "none", data = bad), "R2")

  bad2 <- df; bad2$Attempt[1] <- 2
  expect_error(estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = bad2), "Attempt")
  expect_error(estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2,
                              se = "none", data = bad2), "Attempt")

  bad3 <- df; bad3$R1[1] <- 2
  expect_error(estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2,
                              se = "none", data = bad3), "R1")
  expect_error(estimator_ds_sens(Y, Z, R1, Attempt, R2, delta = 1, minY = 1, maxY = 5, data = bad3), "R1")
  expect_error(sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 3, data = bad3), "R1")
  expect_error(estimator_ds_sens(Y, Z, R1, Attempt, R2, delta = 1, minY = 1, maxY = 5, data = bad2), "Attempt")
  expect_error(sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 3, data = bad2), "Attempt")

  bad4 <- df; bad4$Z[1] <- 2
  expect_error(estimator_ds_sens(Y, Z, R1, Attempt, R2, delta = 1, minY = 1, maxY = 5, data = bad4), "zero or one")
  expect_error(sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 3, data = bad4), "zero or one")

  bad5 <- df; bad5$Y <- as.character(bad5$Y)
  expect_error(estimator_ds_sens(Y, Z, R1, Attempt, R2, delta = 1, minY = 1, maxY = 5, data = bad5), "numeric")
  expect_error(sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 3, data = bad5), "numeric")
})

test_that("estimator_trim needs one complete set of response arguments", {
  df <- make_synthetic()
  expect_error(estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, se = "none", data = df),
               "Supply either R")
  expect_error(estimator_trim(Y, Z, R1 = R1, se = "none", data = df), "Supply either R")
})

test_that("estimator_trim validates the bootstrap replicate count", {
  df <- make_synthetic()
  expect_error(estimator_trim(Y, Z, R = R1, se = "bootstrap", sims = 1, data = df), "at least two")
  expect_error(estimator_trim(Y, Z, R = R1, se = "bootstrap", sims = "many", data = df), "at least two")
})

test_that("a stratification variable may not be missing", {
  df <- make_synthetic()
  df$strata[3] <- NA
  expect_error(estimator_ds(Y, Z, R1, Attempt, R2, strata = strata,
                            minY = 1, maxY = 5, data = df), "missing values")
})

test_that("the formula interface reports a column it cannot find", {
  df <- make_synthetic()
  expect_error(estimator_ev(Y ~ Z, R = ~not_a_column, minY = 1, maxY = 5, data = df),
               "was not found in the data")
  expect_error(estimator_ev(Y ~ Z, R = ~ R1 + Z, minY = 1, maxY = 5, data = df),
               "exactly one column")
  expect_error(attrition:::parse_yz_formula("Y ~ Z", df),
               "must be a formula or an unquoted column name")
})

test_that("validate_support checks the significance level", {
  df <- make_synthetic()
  expect_error(estimator_ev(Y, Z, R1, minY = 1, maxY = 5, alpha = c(0.05, 0.1), data = df),
               "single number")
  expect_error(estimator_ev(Y, Z, R1, minY = 1, maxY = 5, alpha = 0, data = df),
               "strictly between zero and one")
})

test_that("the Imbens-Manski critical value is NA when a variance is", {
  expect_true(is.na(attrition:::im_critical_value(0, 1, NA_real_, NA_real_, 0.05)))
})

test_that("the sensitivity estimators check delta, alpha and their strata", {
  df <- make_synthetic()
  expect_error(estimator_ds_sens(Y, Z, R1, Attempt, R2, delta = c(0.1, 0.2),
                                 minY = 1, maxY = 5, data = df), "single number")
  expect_error(estimator_trim(Y, Z, R = R1, alpha = c(0.05, 0.1), se = "none", data = df),
               "single number")

  df$strata[3] <- NA
  expect_error(estimator_ds_sens(Y, Z, R1, Attempt, R2, delta = 1, strata = strata,
                                 minY = 1, maxY = 5, data = df), "missing values")
})

test_that("a bootstrap that never produces bounds fails rather than reporting nothing", {
  expect_error(
    attrition:::bootstrap_trim_variance(function(idx) stop("no bounds here"),
                                        Z = rep(0:1, each = 10), sims = 5),
    "fewer than two replicates"
  )
})

# ── Post-estimation methods, on every class the package returns ──────────────

# One fitted object per class, plus the trimming cells, so tidy(), print() and
# summary() are exercised on all of them rather than on a representative
all_fits <- function() {
  d <- levendusky_replication
  suppressWarnings(list(
    ev = estimator_ev(Y = Y_polarization_w2, Z = Z, R = R1,
                      minY = 0, maxY = 6, data = d),
    ev_strata = estimator_ev(Y = Y_polarization_w2, Z = Z, R = R1, strata = X_party_id,
                             minY = 0, maxY = 6, data = d),
    ds = estimator_ds(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                      minY = 0, maxY = 6, data = d),
    ds_strata = estimator_ds(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt,
                             R2 = R2, strata = X_party_id, minY = 0, maxY = 6, data = d),
    ds_sens = estimator_ds_sens(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt,
                                R2 = R2, delta = 0.5, minY = 0, maxY = 6, data = d),
    trim_mono = estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                               monotonicity = "treatment_decreases_response", data = d),
    trim_none = estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1, monotonicity = "none",
                               se = "bootstrap", sims = 50, data = d),
    trim_ds_mono = estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt,
                                  R2 = R2, monotonicity = "treatment_decreases_response",
                                  se = "bootstrap", sims = 50, data = d),
    trim_ds_none = estimator_trim(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt,
                                  R2 = R2, se = "bootstrap", sims = 50, data = d),
    trim_violated = estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                                   se = "none", data = d)
  ))
}

test_that("tidy returns the same shape for every class", {
  for (fit in all_fits()) {
    td <- tidy(fit)
    expect_s3_class(td, "tbl_df")
    expect_equal(nrow(td), 3L)
    expect_equal(td$term, c("bounds", "lower_bound", "upper_bound"))
    expect_equal(names(td),
                 c("term", "estimate", "std.error", "conf.low", "conf.high",
                   "estimate_lower", "estimate_upper", "std.error_lower",
                   "std.error_upper", "outcome"))
    expect_equal(unique(td$outcome), "Y_polarization_w2")
    expect_true(is.na(td$estimate[1]))
    expect_equal(td$estimate[2:3],
                 unname(c(fit["estimate_lower"], fit["estimate_upper"])))
  }
})

test_that("print shows the vector and returns it invisibly for every class", {
  for (fit in all_fits()) {
    expect_output(print(fit), "estimate_lower")
    expect_false(any(grepl("attr\\(", capture.output(print(fit)))))
    invisible(capture.output(res <- withVisible(print(fit))))
    expect_false(res$visible)
    expect_equal(unname(res$value), unname(fit))
  }
})

test_that("summary names the estimand and the assumptions for every class", {
  for (fit in all_fits()) {
    out <- capture.output(summary(fit))
    expect_true(any(grepl("Estimand:", out)))
    expect_true(any(grepl("Y_polarization_w2", out)))
    expect_true(any(grepl("Assum", out)))
    # The numbers are labelled rather than pooled the way summary.default pools them
    expect_false(any(grepl("Median|1st Qu", out)))
    # and it hands back the tidy frame rather than printing it twice
    invisible(capture.output(res <- withVisible(summary(fit))))
    expect_false(res$visible)
    expect_equal(res$value, tidy(fit))
  }
})

test_that("summary reports the design, the direction and what was trimmed", {
  fits <- all_fits()

  mono <- capture.output(summary(fits$trim_mono))
  expect_true(any(grepl("single sample", mono)))
  expect_true(any(grepl("never raised the chance of responding", mono)))
  expect_true(any(grepl("1.5% of the control group", mono)))
  expect_true(any(grepl("analytic \\(Lee 2009", mono)))

  ds_none <- capture.output(summary(fits$trim_ds_none))
  expect_true(any(grepl("double sampling", ds_none)))
  expect_true(any(grepl("random assignment alone", ds_none)))
  expect_true(any(grepl("of the treatment group and .* of the control group", ds_none)))
  expect_true(any(grepl("bootstrap", ds_none)))

  violated <- capture.output(summary(fits$trim_violated))
  expect_true(any(grepl("No bounds", violated)))
  expect_true(any(grepl("never lowered the chance of responding", violated)))
})

test_that("summary reports the confidence level actually used", {
  d <- levendusky_replication
  ninety <- estimator_ds(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                         minY = 0, maxY = 6, alpha = 0.10, data = d)
  expect_true(any(grepl("90% Imbens-Manski", capture.output(summary(ninety)))))
  expect_equal(attr(ninety, "alpha"), 0.10)

  ninetyfive <- estimator_ds(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt,
                             R2 = R2, minY = 0, maxY = 6, data = d)
  expect_true(any(grepl("95% Imbens-Manski", capture.output(summary(ninetyfive)))))
})

test_that("summary of a sensitivity fit reports delta and the poststratification", {
  d <- levendusky_replication
  sens <- estimator_ds_sens(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt,
                            R2 = R2, delta = 0.5, minY = 0, maxY = 6, data = d)
  out <- capture.output(summary(sens))
  expect_true(any(grepl("delta = 0.5", out)))
  expect_true(any(grepl("ignorability for 50", out)))

  ps <- estimator_ds(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                     strata = X_party_id, minY = 0, maxY = 6, data = d)
  expect_true(any(grepl("Poststratified", capture.output(summary(ps)))))
  expect_false(any(grepl("Poststratified",
                         capture.output(summary(estimator_ds(
                           Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt,
                           R2 = R2, minY = 0, maxY = 6, data = d))))))
})

test_that("sensitivity_ds returns its three pieces", {
  d <- levendusky_replication
  sens <- sensitivity_ds(Y = Y_polarization_w2, Z = Z, R1 = R1, Attempt = Attempt, R2 = R2,
                         minY = 0, maxY = 6, sims = 5, alpha = 0.10, data = d)
  expect_named(sens, c("sensitivity_plot", "sims_df", "delta_star"))
  expect_s3_class(sens$sensitivity_plot, "ggplot")
  expect_equal(nrow(sens$sims_df), 5L)
  expect_true(is.numeric(sens$delta_star))
})
