# attrition 1.0.0

First release. The package implements the estimators of Coppock, Gerber, Green, and Kern (2017) for randomized experiments with nonignorable missing outcomes.

## Estimators

* `estimator_ev()`: worst-case (Manski 1990) bounds from a single round of data collection, with a joint Imbens-Manski (2004) confidence interval.

* `estimator_ds()`: the double-sampling bounds of the paper, with its analytic variance. A random sample of the initial nonrespondents is pursued a second time, and only the subjects who refuse twice get worst-case treatment.

* `estimator_ds_sens()` and `sensitivity_ds()`: the interpolation between worst-case bounds and ignorability, indexed by `delta`, the share of follow-up nonrespondents left unmodeled, and a sweep over `delta` for delta*, the smallest value at which the confidence interval reaches zero.

* `estimator_trim()`: trimming bounds on the effect among always-reporters. The design and the selection assumption are separate choices. Which response arguments are supplied picks the design, single sample (`R`) or double sampling (`R1`, `Attempt`, `R2`); `monotonicity` picks the assumption, `"treatment_increases_response"` (Lee 2009), `"treatment_decreases_response"`, or `"none"`, which gives the sharp bounds of Imai (2008, Proposition 1) under random assignment alone. All four combinations estimate. Standard errors come from Lee (2009, Proposition 3) where that derivation applies, which is the single-sample, one-group-trimmed case, and from a bootstrap resampled within treatment arm everywhere else; asking for analytic standard errors where they do not apply is an error rather than a silent substitution.

* `estimator_ev()`, `estimator_ds()`, `estimator_ds_sens()` and `sensitivity_ds()` take a `strata` argument for poststratification on a discrete covariate, which targets the same identification region and estimates it more precisely.

## Output

* Every estimator returns the same six named elements in the same order, under broom's names: `estimate_lower` and `estimate_upper`, the two ends of the identification region; `std.error_lower` and `std.error_upper`, their standard errors; and `conf.low` and `conf.high`, the joint Imbens-Manski interval. `estimator_trim()` follows them with the intermediate quantities of the path taken.

* `print()`, `summary()` and `tidy()` methods for every result. `summary()` names the estimand and the assumptions that produced the numbers; for trimming bounds that includes the design, the direction assumed, how much of which group was trimmed, and where the standard errors came from. `tidy()` returns a three-row tibble, a `bounds` row carrying the whole vector and a row per endpoint, so `DeclareDesign::declare_estimator()` can select an endpoint with `term`. `sensitivity_ds()` returns a classed list whose `print()` reports delta* and whose `tidy()` returns the bounds at every value of `delta`.

* A formula interface on every estimator, `estimator_ds(Y ~ Z, R1 = "R1", Attempt = "Attempt", R2 = "R2", ...)`, which is what `declare_estimator(.method = ...)` expects. Response, attempt and stratification arguments accept an unquoted column name, a quoted string, or a one-sided formula.

## Guards

* The assumed support of the outcome must cover the observed outcomes, `alpha` and `delta` must lie in their ranges, and every treatment group needs at least one respondent, one follow-up attempt and one follow-up respondent (within every stratum, when strata are supplied). Each of these is an error naming the problem rather than a `NaN` or a bound that is not a bound.

* A monotonicity violation in the direction assumed returns `NA` bounds with a warning naming the direction the response rates do admit. A zero trimming proportion warns that Lee's interior-point condition fails. Missingness rates summing to one or more, or a trimming proportion that leaves one side of a group empty, are errors, and the bootstrap drops the replicates in which they occur and reports how many survived.

## Data

* `levendusky_replication`: the replication study reported in the paper, a two-wave survey experiment in which 536 of 1,980 subjects did not answer the second wave and 100 of them were followed up at a higher incentive. The tibble holds the polarized-versus-moderate contrast the paper analyzes, with columns named for the roles they play. `data-raw/levendusky_replication.R` rebuilds it from the Harvard Dataverse deposit.

* `vignette("attrition")` works through the design and all five estimators on those data, reproduces Table 3 of the paper, draws the imputation the worst-case bounds average over, and closes with every estimator on one axis, grouped by the population each one describes.

## Changes from the development version on GitHub

For anyone who installed the pre-release package from GitHub:

* Output names changed. The bounding estimators returned `low_est`, `upp_est`, `low_var`, `upp_var`, `ci_lower` and `ci_upper`, and `estimator_trim()` returned `lower_bound`, `upper_bound`, `lower_se` and `upper_se`. The third pair of the old names held variances; the new `std.error_lower` and `std.error_upper` are standard errors.

* The dataset was `levendusky` and held all three experimental conditions with the archive's column names. It is now `levendusky_replication`, the analyzed contrast alone.

* `sensitivity_ds()` returns `delta_star`, a single number or `NA`, in place of `p_star`, and its `sims_df` has a `delta` column in place of `p`.

* The Imbens-Manski critical value was found by a bounded optimizer whose search interval only covered `alpha` near 0.05, so confidence intervals at other levels were wrong: at `alpha = 0.01` it returned 2.000 against a correct 2.326. The root is now found with `uniroot()` on a bracket derived from `alpha`. Bound estimates and variances are unaffected.

* `estimator_trim()` with `R1`, `Attempt` and `R2` previously trimmed both groups without a way to assume monotonicity, and with `R` assumed monotonicity in one fixed direction. Both defaults are unchanged, and the `monotonicity` argument now reaches the other combinations.
