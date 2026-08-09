# Changelog

## attrition 1.0.0

First release.

### New features

- Formula interface on all estimators:
  `estimator_ev(Y ~ Z, R = "R", ...)`,
  `estimator_ds(Y ~ Z, R1 = "R1", Attempt = "Attempt", R2 = "R2", ...)`,
  and likewise for
  [`estimator_trim()`](https://alexandercoppock.com/attrition/reference/estimator_trim.md),
  [`estimator_ds_sens()`](https://alexandercoppock.com/attrition/reference/estimator_ds_sens.md),
  and
  [`sensitivity_ds()`](https://alexandercoppock.com/attrition/reference/sensitivity_ds.md).
  Response, attempt, and stratification arguments accept either an
  unquoted column name or a quoted string. The formula form is what
  `DeclareDesign::declare_estimator(.method = ...)` expects.

- [`tidy()`](https://generics.r-lib.org/reference/tidy.html) methods for
  all estimator output classes (`attrition_bounds`, `attrition_trim`)
  via the `generics` package. Each method returns a three-row tibble
  with rows `"bounds"`, `"lower_bound"`, and `"upper_bound"`, suitable
  for use with DeclareDesign. The `"bounds"` row carries the joint
  Imbens-Manski confidence interval in `conf.low`/`conf.high` and the
  bound point estimates in `estimate.low`/`estimate.high`; the
  individual rows carry point estimates and standard errors selectable
  via the `term` argument.

- S3 classes on all estimator return values:
  [`estimator_ds()`](https://alexandercoppock.com/attrition/reference/estimator_ds.md)
  returns class `c("attrition_ds", "attrition_bounds")`,
  [`estimator_ev()`](https://alexandercoppock.com/attrition/reference/estimator_ev.md)
  returns `c("attrition_ev", "attrition_bounds")`,
  [`estimator_ds_sens()`](https://alexandercoppock.com/attrition/reference/estimator_ds_sens.md)
  returns `c("attrition_ds_sens", "attrition_bounds")`, and
  [`estimator_trim()`](https://alexandercoppock.com/attrition/reference/estimator_trim.md)
  returns `c("attrition_trim")`.

- Stratified estimators now index results by name rather than position,
  preventing silent errors if the output vector order ever changes.

- [`sensitivity_ds()`](https://alexandercoppock.com/attrition/reference/sensitivity_ds.md)
  no longer calls [`require()`](https://rdrr.io/r/base/library.html) at
  runtime; all dependencies are declared via `@importFrom`.

- `reshape2` dependency replaced by `tibble`; `generics` added for
  `tidy` re-export.

- [`estimator_trim()`](https://alexandercoppock.com/attrition/reference/estimator_trim.md)
  gains standard errors and confidence intervals, via a new `se`
  argument. `se = "analytic"` (the default) uses the closed-form
  asymptotic variance of Lee (2009), Proposition 3, and covers the
  single-stage `R` path. `se = "bootstrap"` resamples units within
  treatment arm and covers both that path and the weighted
  double-sampling path, for which no analytic variance exists in the
  literature. `se = "none"` returns bounds alone. The bound standard
  errors feed the existing Imbens-Manski machinery, so
  [`tidy.attrition_trim()`](https://alexandercoppock.com/attrition/reference/tidy.attrition_trim.md)
  now returns `std.error`, `conf.low` and `conf.high` instead of `NA`
  throughout.

  The analytic variance was validated three ways: Lee’s published term
  for the estimated trimming proportion agrees to machine precision with
  the algebraically distinct form used by Tauchmann’s Stata `leebounds`;
  the mean analytic standard error tracks the Monte Carlo standard
  deviation of the estimator across sample sizes from 1,000 to 64,000;
  and a conventional 95% interval around each bound endpoint covers the
  true population endpoint 95.4% and 94.9% of the time.

  Requesting `se = "analytic"` on the double-sampling path is an error
  rather than a silent substitution, since Lee’s derivation assumes
  i.i.d. sampling and trimming of one group only. A zero trimming
  proportion warns, because the bounds then collapse to a point on the
  boundary of the parameter space and Lee’s interior-point condition
  fails.

- Worked examples on every exported function.

### Bug fixes

- The Imbens-Manski critical value was found by minimizing an absolute
  value with `optim(..., method = "Brent", lower = 1, upper = 2)`. That
  interval only covers alpha near 0.05, so at other significance levels
  the optimizer returned a boundary value and the confidence intervals
  were silently wrong: at `alpha = 0.01` it returned 2.000 against a
  correct 2.326, and at `alpha = 0.001` it returned 2.000 against a
  correct 3.090. The search interval is now derived from alpha, and the
  root is found with [`uniroot()`](https://rdrr.io/r/stats/uniroot.html)
  on the signed coverage excess rather than by minimizing its absolute
  value. Bound point estimates and variances are unaffected.

- [`estimator_trim()`](https://alexandercoppock.com/attrition/reference/estimator_trim.md)
  reported every failure inside `trimming_bounds()` as a monotonicity
  violation, because it caught all errors and returned `NA` bounds. It
  now catches only a classed monotonicity condition.

- The monotonicity violation message named the wrong group. `Q < 0`
  means the treatment group is more likely to be missing than the
  control group.

- `trimming_bounds()` built its weighted CDFs with a loop that counted
  backwards when a treatment group had one observed outcome, silently
  appending an `NA`, and threw an opaque error when a group had none.
  Both are now [`cumsum()`](https://rdrr.io/r/base/cumsum.html), with an
  explicit check for empty groups.

- [`sensitivity_ds()`](https://alexandercoppock.com/attrition/reference/sensitivity_ds.md)
  could not detect a delta\* at the last point of the delta grid,
  reporting a genuine crossing near delta = 1 as no crossing at all.

- `minY`, `maxY`, `alpha`, `delta`, and `sims` are validated. Previously
  a reversed `minY`/`maxY` returned a lower bound above the upper bound,
  an assumed support that did not cover the observed outcomes returned
  bounds that are not bounds, and a `delta` outside `[0, 1]` returned
  `NaN` variances.

- [`estimator_ev()`](https://alexandercoppock.com/attrition/reference/estimator_ev.md):
  fixed a copy-paste error where `n1_c_s` was incorrectly referenced as
  `n1_c_c` in the unstratified path.

- [`estimator_ds()`](https://alexandercoppock.com/attrition/reference/estimator_ds.md)
  and
  [`estimator_ds_sens()`](https://alexandercoppock.com/attrition/reference/estimator_ds_sens.md):
  replaced a commented-out sentinel initialization (`-99`) with
  `NA_real_` for `c1a_t`, `c1r_t`, `c2a_t`, `c2r_t`. These arguments
  were never used and have since been removed from the internal
  estimators altogether.
