# Package index

## Bounds on the average treatment effect

Worst-case bounds from a single round of data collection, and the
double-sampling bounds that narrow them by chasing a random sample of
the nonrespondents. All three accept a stratification variable for
poststratification, which targets the same identification region and
estimates it more precisely.

- [`estimator_ev()`](https://alexandercoppock.com/attrition/reference/estimator_ev.md)
  : Extreme Value (Manski) Bounds
- [`estimator_ds()`](https://alexandercoppock.com/attrition/reference/estimator_ds.md)
  : Extreme Value Bounds with Double Sampling
- [`estimator_ds_sens()`](https://alexandercoppock.com/attrition/reference/estimator_ds_sens.md)
  : Extreme Value Bounds with Double Sampling with Sensitivity

## Sensitivity analysis

How far ignorability has to fail among the follow-up nonrespondents
before a finding stops holding.

- [`sensitivity_ds()`](https://alexandercoppock.com/attrition/reference/sensitivity_ds.md)
  : Sensitivity Analysis

## Trimming bounds

Lee (2009) bounds, which assume treatment moves response in one
direction only and do not require the outcome to be bounded.

- [`estimator_trim()`](https://alexandercoppock.com/attrition/reference/estimator_trim.md)
  : Trimming Bounds

## Working with results

Every estimator returns a named numeric vector. The tidy methods put the
same quantities in a data frame with the column names DeclareDesign
expects.

- [`tidy(`*`<attrition_bounds>`*`)`](https://alexandercoppock.com/attrition/reference/tidy.attrition_bounds.md)
  : Tidy an attrition bounds object
- [`tidy(`*`<attrition_trim>`*`)`](https://alexandercoppock.com/attrition/reference/tidy.attrition_trim.md)
  : Tidy a trimming bounds object
- [`print(`*`<attrition_bounds>`*`)`](https://alexandercoppock.com/attrition/reference/print.attrition_bounds.md)
  : Print bounds
- [`print(`*`<attrition_trim>`*`)`](https://alexandercoppock.com/attrition/reference/print.attrition_trim.md)
  : Print trimming bounds
- [`reexports`](https://alexandercoppock.com/attrition/reference/reexports.md)
  [`tidy`](https://alexandercoppock.com/attrition/reference/reexports.md)
  : Objects exported from other packages

## Data

- [`levendusky`](https://alexandercoppock.com/attrition/reference/levendusky.md)
  : Perceived polarization under double sampling

## Package overview

- [`attrition`](https://alexandercoppock.com/attrition/reference/attrition-package.md)
  [`attrition-package`](https://alexandercoppock.com/attrition/reference/attrition-package.md)
  : attrition: bounds for experiments with missing outcomes
