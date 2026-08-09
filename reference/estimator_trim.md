# Trimming Bounds

Trimming Bounds

## Usage

``` r
estimator_trim(
  Y,
  Z,
  R = NULL,
  R1 = NULL,
  Attempt = NULL,
  R2 = NULL,
  strata = NULL,
  alpha = 0.05,
  se = c("analytic", "bootstrap", "none"),
  sims = 1000,
  data
)
```

## Arguments

- Y:

  The (unquoted) outcome variable, or a formula `outcome ~ treatment`
  for use with `declare_estimator(.method = estimator_trim)`. Must be
  numeric.

- Z:

  The (unquoted) assignment indicator variable. Must be numeric and take
  values 0 or 1. Ignored when `Y` is a formula.

- R:

  The single-stage response indicator: unquoted column name, or a quoted
  string column name when using the formula interface. Must be numeric
  and take values 0 or 1. Supply either `R` (single-stage) or
  `R1`/`Attempt`/`R2` (double-sampling).

- R1:

  The initial sample response indicator. Unquoted or quoted string
  column name. Must be numeric and take values 0 or 1.

- Attempt:

  The follow-up attempt indicator. Unquoted or quoted string column
  name. Must be numeric and take values 0 or 1.

- R2:

  The follow-up response indicator. Unquoted or quoted string column
  name. Must be numeric and take values 0 or 1.

- strata:

  Not supported; supplying any value raises an error.

- alpha:

  The desired significance level. 0.05 by default.

- se:

  How to obtain standard errors. `"analytic"` (the default) uses the
  closed-form asymptotic variance of Lee (2009), Proposition 3, and is
  available for the single-stage `R` path only. `"bootstrap"` resamples
  units within treatment arm and works for both paths. `"none"` returns
  bounds alone.

- sims:

  Number of bootstrap replicates when `se = "bootstrap"`. 1000 by
  default.

- data:

  A dataframe

## Value

A named numeric vector containing `lower_bound` and `upper_bound`, the
trimming bound estimates; `lower_se` and `upper_se`, their standard
errors; `ci_lower` and `ci_upper`, the joint Imbens-Manski confidence
interval; and the intermediate quantities used to build them. All
elements are `NA` when monotonicity is violated. Pass to
[`tidy()`](https://alexandercoppock.com/attrition/reference/tidy.attrition_trim.md)
for a data frame.

## Details

The analytic variance has four contributions: the variance of the
retained (trimmed) outcomes, the variance from estimating the trimming
threshold, the variance from estimating the trimming proportion, and the
variance of the control-group respondent mean. The third of these is
often the largest, so treating the trimming proportion as known would
understate the uncertainty substantially.

Lee's derivation assumes the bounds are interior points, which fails
when the two response rates are equal: the trimming proportion is then
zero, the bounds collapse to a point, and the standard errors are not
trustworthy. This case warns.

## References

Lee, David S. (2009). Training, Wages, and Sample Selection: Estimating
Sharp Bounds on Treatment Effects. *Review of Economic Studies*
76(3):1071-1102.

Tauchmann, Harald (2014). Lee (2009) Treatment-Effect Bounds for
Nonrandom Sample Selection. *Stata Journal* 14(4):884-894.

## Examples

``` r
set.seed(343)
N <- 1000
Y_0 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
Y_1 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.1, 0.4, 0.3, 0.1))
Z <- rbinom(N, 1, 0.5)
Y_star <- Z * Y_1 + (1 - Z) * Y_0

# Treated units respond at a higher rate, so the missingness is nonignorable
R <- rbinom(N, 1, prob = 0.7 + 0.1 * Z)
Y <- Y_star
Y[R == 0] <- NA
df <- data.frame(Y, Z, R)

# Single-stage: trimming bounds under monotonicity, with Lee (2009) standard errors
estimator_trim(Y, Z, R = R, data = df)
#>     upper_bound     lower_bound       Out0_mono      Out1L_mono      Out1U_mono 
#>      0.64239649      0.01761468      2.93353474      2.95114943      3.57593123 
#> control_group_N   treat_group_N               Q              f1              f0 
#>    331.00000000    417.00000000      0.16385742      0.18713450      0.32032854 
#>          pi_r_1          pi_r_0              yU              yL        lower_se 
#>      0.81286550      0.67967146      2.00000000      4.00000000      0.09007908 
#>        upper_se        ci_lower        ci_upper 
#>      0.09972250     -0.13055222      0.80642541 

# Bootstrap standard errors instead
estimator_trim(Y, Z, R = R, se = "bootstrap", sims = 200, data = df)
#>     upper_bound     lower_bound       Out0_mono      Out1L_mono      Out1U_mono 
#>      0.64239649      0.01761468      2.93353474      2.95114943      3.57593123 
#> control_group_N   treat_group_N               Q              f1              f0 
#>    331.00000000    417.00000000      0.16385742      0.18713450      0.32032854 
#>          pi_r_1          pi_r_0              yU              yL        lower_se 
#>      0.81286550      0.67967146      2.00000000      4.00000000      0.08467663 
#>        upper_se        ci_lower        ci_upper 
#>      0.09880524     -0.12166598      0.80491665 
```
