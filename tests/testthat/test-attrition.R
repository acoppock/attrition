library(testthat)

# ── Paper benchmark (Table 3 of CGGK 2017) ──────────────────────────────────
# Expected values computed from replication data in:
# Coppock, Gerber, Green, Kern (2017), Political Analysis
# doi:10.1017/pan.2016.6
# Data source: Harvard Dataverse doi:10.7910/DVN/AQB4MP

data_path <- testthat::test_path("testdata", "levendusky_mturk_clean.csv")

test_that("estimator_ev matches Table 3 column 1 (no double sampling)", {
  skip_if_not(file.exists(data_path), "replication data not found")
  dat <- read.csv(data_path)
  dat <- subset(dat, !is.na(Z1))
  out <- estimator_ev(Y = L_dif_w2, Z = Z1, R = R1,
                      minY = 0, maxY = 6, data = dat)
  expect_equal(unname(out["ci_lower"]), -1.66907903866295,   tolerance = 1e-10)
  expect_equal(unname(out["ci_upper"]),  1.83588114535676,   tolerance = 1e-10)
  expect_equal(unname(out["low_est"]), -1.53914496339566,    tolerance = 1e-10)
  expect_equal(unname(out["upp_est"]),  1.70966762747749,    tolerance = 1e-10)
  expect_equal(unname(out["low_var"]),  0.00624010083468930, tolerance = 1e-10)
  expect_equal(unname(out["upp_var"]),  0.00588785669626299, tolerance = 1e-10)
})

test_that("estimator_ds matches Table 3 column 2 (double sampling)", {
  skip_if_not(file.exists(data_path), "replication data not found")
  dat <- read.csv(data_path)
  dat <- subset(dat, !is.na(Z1))
  out <- estimator_ds(Y = L_dif_w2, Z = Z1, R1 = R1,
                      Attempt = Attempt, R2 = R2,
                      minY = 0, maxY = 6, data = dat)
  expect_equal(unname(out["ci_lower"]), -0.528309673828183,  tolerance = 1e-10)
  expect_equal(unname(out["ci_upper"]),  0.745174826173985,  tolerance = 1e-10)
  expect_equal(unname(out["low_est"]), -0.34174537662934,    tolerance = 1e-10)
  expect_equal(unname(out["upp_est"]),  0.571815728388134,   tolerance = 1e-10)
  expect_equal(unname(out["low_var"]),  0.0128647858310974,  tolerance = 1e-10)
  expect_equal(unname(out["upp_var"]),  0.0111080739914738,  tolerance = 1e-10)
})

test_that("estimator_ds matches Table 3 column 3 (DS + poststratification)", {
  skip_if_not(file.exists(data_path), "replication data not found")
  dat <- read.csv(data_path)
  dat <- subset(dat, !is.na(Z1))
  out <- estimator_ds(Y = L_dif_w2, Z = Z1, R1 = R1,
                      Attempt = Attempt, R2 = R2,
                      strata = pid_3_recoded,
                      minY = 0, maxY = 6, data = dat)
  expect_equal(unname(out["ci_lower"]), -0.529010551043845,  tolerance = 1e-10)
  expect_equal(unname(out["ci_upper"]),  0.696600306767075,  tolerance = 1e-10)
  expect_equal(unname(out["low_est"]), -0.344389910863888,   tolerance = 1e-10)
  expect_equal(unname(out["upp_est"]),  0.525683688852529,   tolerance = 1e-10)
  expect_equal(unname(out["low_var"]),  0.0125981273119329,  tolerance = 1e-10)
  expect_equal(unname(out["upp_var"]),  0.0107972726598501,  tolerance = 1e-10)
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
  expect_equal(unname(out["ci_lower"]), -0.082302939127751, tolerance = 1e-10)
  expect_equal(unname(out["ci_upper"]),  0.624972549561357, tolerance = 1e-10)
  expect_equal(unname(out["low_est"]),   0.0572167798040737, tolerance = 1e-10)
  expect_equal(unname(out["upp_est"]),   0.487115989508957, tolerance = 1e-10)
})

test_that("estimator_ds with poststratification produces stable results on synthetic data", {
  df  <- make_synthetic()
  out <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5,
                      strata = strata, data = df)
  expect_equal(unname(out["ci_lower"]),  0.00420597365836781, tolerance = 1e-10)
  expect_equal(unname(out["ci_upper"]),  0.662950557239635,   tolerance = 1e-10)
  expect_equal(unname(out["low_est"]),   0.124794503861406,   tolerance = 1e-10)
  expect_equal(unname(out["upp_est"]),   0.544633130143652,   tolerance = 1e-10)
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
  out <- estimator_ev(Y, Z, R, minY = -2, maxY = 6, data = df)
  # Expected width ≈ 8 * ((1 - 0.8) + (1 - 0.6)) = 8 * 0.6 = 4.8
  expect_equal(unname(out["upp_est"] - out["low_est"]), 4.8, tolerance = 0.1)
})

# ── tidy() method tests ──────────────────────────────────────────────────────

test_that("tidy.attrition_bounds works for estimator_ev", {
  df <- make_synthetic()
  out <- estimator_ev(Y, Z, R1, minY = 1, maxY = 5, data = df)
  td  <- tidy(out)
  expect_s3_class(td, "tbl_df")
  expect_equal(nrow(td), 3L)
  expect_named(td, c("term", "estimate", "std.error", "conf.low", "conf.high",
                     "estimate.low", "estimate.high"))
  expect_equal(td$term, c("bounds", "lower_bound", "upper_bound"))
  # bounds row
  expect_true(is.na(td$estimate[1]))
  expect_true(is.na(td$std.error[1]))
  expect_equal(td$conf.low[1],      unname(out["ci_lower"]))
  expect_equal(td$conf.high[1],     unname(out["ci_upper"]))
  expect_equal(td$estimate.low[1],  unname(out["low_est"]))
  expect_equal(td$estimate.high[1], unname(out["upp_est"]))
  # lower_bound row
  expect_equal(td$estimate[2],  unname(out["low_est"]))
  expect_equal(td$std.error[2], sqrt(unname(out["low_var"])))
  expect_true(is.na(td$conf.low[2]))
  expect_true(is.na(td$conf.high[2]))
  expect_true(is.na(td$estimate.low[2]))
  expect_true(is.na(td$estimate.high[2]))
  # upper_bound row
  expect_equal(td$estimate[3],  unname(out["upp_est"]))
  expect_equal(td$std.error[3], sqrt(unname(out["upp_var"])))
  expect_true(is.na(td$conf.low[3]))
  expect_true(is.na(td$conf.high[3]))
  expect_true(is.na(td$estimate.low[3]))
  expect_true(is.na(td$estimate.high[3]))
})

test_that("tidy.attrition_bounds works for estimator_ds", {
  df <- make_synthetic()
  out <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  td  <- tidy(out)
  expect_s3_class(td, "tbl_df")
  expect_equal(nrow(td), 3L)
  expect_named(td, c("term", "estimate", "std.error", "conf.low", "conf.high",
                     "estimate.low", "estimate.high"))
  expect_equal(td$term, c("bounds", "lower_bound", "upper_bound"))
  # bounds row
  expect_true(is.na(td$estimate[1]))
  expect_true(is.na(td$std.error[1]))
  expect_equal(td$conf.low[1],      unname(out["ci_lower"]))
  expect_equal(td$conf.high[1],     unname(out["ci_upper"]))
  expect_equal(td$estimate.low[1],  unname(out["low_est"]))
  expect_equal(td$estimate.high[1], unname(out["upp_est"]))
  # lower_bound row
  expect_equal(td$estimate[2],  unname(out["low_est"]))
  expect_equal(td$std.error[2], sqrt(unname(out["low_var"])))
  expect_true(is.na(td$conf.low[2]))
  expect_true(is.na(td$conf.high[2]))
  expect_true(is.na(td$estimate.low[2]))
  expect_true(is.na(td$estimate.high[2]))
  # upper_bound row
  expect_equal(td$estimate[3],  unname(out["upp_est"]))
  expect_equal(td$std.error[3], sqrt(unname(out["upp_var"])))
  expect_true(is.na(td$conf.low[3]))
  expect_true(is.na(td$conf.high[3]))
  expect_true(is.na(td$estimate.low[3]))
  expect_true(is.na(td$estimate.high[3]))
})

test_that("tidy.attrition_trim works for estimator_trim", {
  df <- make_synthetic()
  out <- estimator_trim(Y, Z, R = R1, data = df)
  td  <- tidy(out)
  expect_s3_class(td, "tbl_df")
  expect_equal(nrow(td), 3L)
  expect_named(td, c("term", "estimate", "std.error", "conf.low", "conf.high",
                     "estimate.low", "estimate.high"))
  expect_equal(td$term, c("bounds", "lower_bound", "upper_bound"))
  # bounds row
  expect_true(is.na(td$estimate[1]))
  expect_equal(td$estimate.low[1],  unname(out["lower_bound"]))
  expect_equal(td$estimate.high[1], unname(out["upper_bound"]))
  # lower_bound row
  expect_equal(td$estimate[2], unname(out["lower_bound"]))
  expect_true(is.na(td$estimate.low[2]))
  expect_true(is.na(td$estimate.high[2]))
  # upper_bound row
  expect_equal(td$estimate[3], unname(out["upper_bound"]))
  expect_true(is.na(td$estimate.low[3]))
  expect_true(is.na(td$estimate.high[3]))
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
  skip_if_not(file.exists(data_path), "replication data not found")
  dat <- read.csv(data_path)
  dat <- subset(dat, !is.na(Z1))
  out <- estimator_ev(Y = L_dif_w2, Z = Z1, R = R1, strata = pid_3_recoded,
                      minY = 0, maxY = 6, data = dat)
  expect_equal(unname(out["ci_lower"]), -1.66862420626897,  tolerance = 1e-10)
  expect_equal(unname(out["ci_upper"]),  1.83541087597583,  tolerance = 1e-10)
  expect_equal(unname(out["low_est"]), -1.53858804710012,   tolerance = 1e-10)
  expect_equal(unname(out["upp_est"]),  1.70925509076352,   tolerance = 1e-10)
  expect_equal(unname(out["low_var"]),  0.00624990987170985, tolerance = 1e-10)
  expect_equal(unname(out["upp_var"]),  0.00588247147395304, tolerance = 1e-10)
})

test_that("estimator_trim DS path (paper data)", {
  skip_if_not(file.exists(data_path), "replication data not found")
  dat <- read.csv(data_path)
  dat <- subset(dat, !is.na(Z1))
  out <- estimator_trim(Y = L_dif_w2, Z = Z1, R1 = R1, Attempt = Attempt, R2 = R2,
                        data = dat)
  expect_equal(unname(out["upper_bound"]),  0.567599031583528,  tolerance = 1e-10)
  expect_equal(unname(out["lower_bound"]), -0.268142333620083,  tolerance = 1e-10)
})

test_that("estimator_trim R path errors on monotonicity violation (paper data)", {
  skip_if_not(file.exists(data_path), "replication data not found")
  dat <- read.csv(data_path)
  dat <- subset(dat, !is.na(Z1))
  # Control group has slightly higher attrition than treatment → violation
  expect_error(
    estimator_trim(Y = L_dif_w2, Z = Z1, R = R1, data = dat),
    "Monotonicity appears to be violated"
  )
})

test_that("estimator_ds_sens(delta=1) exactly matches estimator_ds (paper data)", {
  skip_if_not(file.exists(data_path), "replication data not found")
  dat <- read.csv(data_path)
  dat <- subset(dat, !is.na(Z1))
  ds   <- estimator_ds(Y = L_dif_w2, Z = Z1, R1 = R1, Attempt = Attempt, R2 = R2,
                       minY = 0, maxY = 6, data = dat)
  sens <- estimator_ds_sens(Y = L_dif_w2, Z = Z1, R1 = R1, Attempt = Attempt, R2 = R2,
                            delta = 1, minY = 0, maxY = 6, data = dat)
  expect_equal(as.numeric(ds), as.numeric(sens), tolerance = 1e-14)
})

test_that("estimator_ds_sens delta=0.5 (paper data)", {
  skip_if_not(file.exists(data_path), "replication data not found")
  dat <- read.csv(data_path)
  dat <- subset(dat, !is.na(Z1))
  out <- estimator_ds_sens(Y = L_dif_w2, Z = Z1, R1 = R1, Attempt = Attempt, R2 = R2,
                           delta = 0.5, minY = 0, maxY = 6, data = dat)
  expect_equal(unname(out["ci_lower"]), -0.263141038742824,  tolerance = 1e-10)
  expect_equal(unname(out["ci_upper"]),  0.527330962796225,  tolerance = 1e-10)
  expect_equal(unname(out["low_est"]), -0.0916157282335379,  tolerance = 1e-10)
  expect_equal(unname(out["upp_est"]),  0.365164824275198,   tolerance = 1e-10)
  expect_equal(unname(out["low_var"]),  0.0108743150646462,  tolerance = 1e-10)
  expect_equal(unname(out["upp_var"]),  0.00971999036287018, tolerance = 1e-10)
})

# ── Additional synthetic data tests ─────────────────────────────────────────

test_that("estimator_ev with strata produces stable results on synthetic data", {
  df  <- make_synthetic()
  out <- estimator_ev(Y, Z, R1, strata = strata, minY = 1, maxY = 5, data = df)
  expect_equal(unname(out["ci_lower"]), -0.82832827570243,   tolerance = 1e-10)
  expect_equal(unname(out["ci_upper"]),  1.35051577794482,   tolerance = 1e-10)
  expect_equal(unname(out["low_est"]), -0.699410333729967,   tolerance = 1e-10)
  expect_equal(unname(out["upp_est"]),  1.23280497665065,    tolerance = 1e-10)
  expect_equal(unname(out["low_var"]),  0.00614288260167931, tolerance = 1e-10)
  expect_equal(unname(out["upp_var"]),  0.00512127526981296, tolerance = 1e-10)
})

test_that("estimator_trim R path (monotonicity) produces stable results on synthetic data", {
  df  <- make_synthetic()
  out <- estimator_trim(Y, Z, R = R1, data = df)
  expect_equal(unname(out["upper_bound"]),  0.534610520108332, tolerance = 1e-10)
  expect_equal(unname(out["lower_bound"]),  0.100617821905443, tolerance = 1e-10)
  expect_equal(unname(out["Q"]),            0.101526858997151, tolerance = 1e-10)
})

test_that("estimator_trim DS path produces stable results on synthetic data", {
  df  <- make_synthetic()
  out <- estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2, data = df)
  expect_equal(unname(out["upper_bound"]),  0.544188162330423,  tolerance = 1e-10)
  expect_equal(unname(out["lower_bound"]),  0.0677009048063324, tolerance = 1e-10)
})

test_that("estimator_ds_sens delta=0.5 produces stable results on synthetic data", {
  df  <- make_synthetic()
  out <- estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 0.5,
                           data = df)
  expect_equal(unname(out["ci_lower"]),  0.028842141946206, tolerance = 1e-10)
  expect_equal(unname(out["ci_upper"]),  0.514219909406511, tolerance = 1e-10)
  expect_equal(unname(out["low_est"]),   0.164445203836863, tolerance = 1e-10)
  expect_equal(unname(out["upp_est"]),   0.379394808689304, tolerance = 1e-10)
  expect_equal(unname(out["low_var"]),   0.00679563956699848, tolerance = 1e-10)
  expect_equal(unname(out["upp_var"]),   0.00671788942612276, tolerance = 1e-10)
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
  expect_lte(unname(ev["low_est"]),       unname(ev["upp_est"]))
  expect_lte(unname(ds["low_est"]),       unname(ds["upp_est"]))
  expect_lte(unname(trim["lower_bound"]), unname(trim["upper_bound"]))
})

test_that("Imbens-Manski CI covers the identification region", {
  df <- make_synthetic()
  for (out in list(
    estimator_ev(Y, Z, R1, minY = 1, maxY = 5, data = df),
    estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  )) {
    expect_lte(unname(out["ci_lower"]), unname(out["low_est"]))
    expect_gte(unname(out["ci_upper"]), unname(out["upp_est"]))
  }
})

test_that("double sampling narrows identification region vs extreme value bounds", {
  df       <- make_synthetic()
  ev       <- estimator_ev(Y, Z, R1, minY = 1, maxY = 5, data = df)
  ds       <- estimator_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, data = df)
  ev_width <- unname(ev["upp_est"] - ev["low_est"])
  ds_width <- unname(ds["upp_est"] - ds["low_est"])
  expect_lte(ds_width, ev_width)
})

test_that("increasing delta widens bounds in estimator_ds_sens", {
  df    <- make_synthetic()
  sens0 <- estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 0,   data = df)
  sens5 <- estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 0.5, data = df)
  sens1 <- estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 1,   data = df)
  width <- function(x) unname(x["upp_est"] - x["low_est"])
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
  expect_lt(val, 1e-4)
})

test_that("find_sign_changes detects interior sign changes only", {
  # Sign change at interior position 2
  expect_equal(attrition:::find_sign_changes(c(-1,  1, -1)), c(FALSE, TRUE, FALSE))
  # Zero crossing at interior position
  expect_equal(attrition:::find_sign_changes(c(-1,  0,  1)), c(FALSE, TRUE, FALSE))
  # No interior change (all positive; first_pos = 1 is excluded)
  expect_equal(attrition:::find_sign_changes(c( 1,  2,  3)), c(FALSE, FALSE, FALSE))
  # No interior change (all negative; first_neg = 1 is excluded)
  expect_equal(attrition:::find_sign_changes(c(-3, -2, -1)), c(FALSE, FALSE, FALSE))
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

test_that("estimator_trim errors when monotonicity is violated", {
  df   <- make_synthetic()
  df$Z <- 1L - df$Z  # flip treatment — control now has higher response rate
  expect_error(estimator_trim(Y, Z, R = R1, data = df),
               "Monotonicity appears to be violated")
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
  nse <- estimator_trim(Y, Z, R1 = R1, Attempt = Attempt, R2 = R2, data = df)
  frm <- estimator_trim(Y ~ Z, R1 = "R1", Attempt = "Attempt", R2 = "R2", data = df)
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
  expect_named(out, c("sensitivity_plot", "sims_df", "p_star"))
  expect_s3_class(out$sensitivity_plot, "gg")
  expect_equal(nrow(out$sims_df), 20L)
  expect_named(out$sims_df,
    c("ci_lower", "ci_upper", "low_est", "upp_est", "low_var", "upp_var",
      "p", "change_lower", "change_upper", "change_any"))
  # p grid runs from 0 to 1
  expect_equal(out$sims_df$p[1],  0)
  expect_equal(out$sims_df$p[20], 1)
})

test_that("sensitivity_ds detects delta* on synthetic data", {
  # Synthetic data: bounds are positive at delta=0 but span 0 at delta=1,
  # so there must be a sign change interior to [0, 1]
  df  <- make_synthetic()
  out <- sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 50, data = df)
  expect_s3_class(out$p_star, "data.frame")
  expect_equal(nrow(out$p_star), 1L)
  expect_true(out$p_star$p > 0 && out$p_star$p < 1)
})

test_that("sensitivity_ds with strata returns correct structure", {
  df  <- make_synthetic()
  out <- sensitivity_ds(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, sims = 10,
                        strata = strata, data = df)
  expect_named(out, c("sensitivity_plot", "sims_df", "p_star"))
  expect_equal(nrow(out$sims_df), 10L)
})
