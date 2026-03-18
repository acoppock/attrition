# attrition 0.1.0

First release.

## New features

* `tidy()` methods for all estimator output classes (`attrition_bounds`,
  `attrition_trim`) via the `generics` package. Each method returns a
  three-row tibble with rows `"bounds"`, `"lower_bound"`, and `"upper_bound"`,
  suitable for use with DeclareDesign. The `"bounds"` row carries the joint
  Imbens-Manski confidence interval in `conf.low`/`conf.high` and the bound
  point estimates in `estimate.low`/`estimate.high`; the individual rows carry
  point estimates and standard errors selectable via the `term` argument.

* S3 classes on all estimator return values: `estimator_ds()` returns class
  `c("attrition_ds", "attrition_bounds")`, `estimator_ev()` returns
  `c("attrition_ev", "attrition_bounds")`, `estimator_ds_sens()` returns
  `c("attrition_ds_sens", "attrition_bounds")`, and `estimator_trim()` returns
  `c("attrition_trim")`.

* Stratified estimators now index results by name rather than position,
  preventing silent errors if the output vector order ever changes.

* `sensitivity_ds()` no longer calls `require()` at runtime; all dependencies
  are declared via `@importFrom`.

* `reshape2` dependency replaced by `tibble`; `generics` added for `tidy`
  re-export.

## Bug fixes

* `estimator_ev()`: fixed a copy-paste error where `n1_c_s` was incorrectly
  referenced as `n1_c_c` in the unstratified path.

* `estimator_ds()` and `estimator_ds_sens()`: replaced a commented-out
  sentinel initialisation (`-99`) with `NA_real_` for `c1a_t`, `c1r_t`,
  `c2a_t`, `c2r_t`.
