# Extreme Value Bounds with Double Sampling with Sensitivity

This function yields extreme value bounds under the assumption that the
outcomes of 1-delta of the missing second-round units are ignorable,
that is, that they are drawn from an unknown distribution with mean and
variance equal to the observed second-round groups.

## Usage

``` r
estimator_ds_sens(
  Y,
  Z,
  R1,
  Attempt,
  R2,
  minY,
  maxY,
  delta,
  strata = NULL,
  alpha = 0.05,
  data
)
```

## Arguments

- Y:

  The (unquoted) outcome variable, or a formula `outcome ~ treatment`
  for use with `declare_estimator(.method = estimator_ds_sens)`. Must be
  numeric.

- Z:

  The (unquoted) assignment indicator variable. Must be numeric and take
  values 0 or 1. Ignored when `Y` is a formula.

- R1:

  The initial sample response indicator: unquoted column name, or a
  quoted string column name when using the formula interface. Must be
  numeric and take values 0 or 1.

- Attempt:

  The follow-up attempt indicator: unquoted column name, or quoted
  string. Must be numeric and take values 0 or 1.

- R2:

  The follow-up response indicator: unquoted column name, or quoted
  string. Must be numeric and take values 0 or 1.

- minY:

  The minimum possible value of the outcome (Y) variable.

- maxY:

  The maximum possible value of the outcome (Y) variable.

- delta:

  Sensitivity parameter in \[0, 1\]. At delta = 1 (default) worst-case
  bounds apply; at delta = 0 ignorability holds for all follow-up
  non-responders.

- strata:

  Stratification variable: unquoted column name or a quoted string
  column name.

- alpha:

  The desired significance level. 0.05 by default.

- data:

  A dataframe

## Value

A named numeric vector with elements `ci_lower` and `ci_upper`, the
joint Imbens-Manski confidence interval; `low_est` and `upp_est`, the
bound point estimates; and `low_var` and `upp_var`, their variances.
Pass to
[`tidy()`](https://alexandercoppock.com/attrition/reference/tidy.attrition_bounds.md)
for a data frame.

## Examples

``` r
set.seed(343)
N <- 1000
Y_0 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
Y_1 <- sample(1:5, N, replace = TRUE, prob = c(0.1, 0.1, 0.4, 0.3, 0.1))
Z <- rbinom(N, 1, 0.5)
Y_star <- Z * Y_1 + (1 - Z) * Y_0
R1 <- rbinom(N, 1, prob = 0.7 + 0.1 * Z)
Y <- Y_star
Y[R1 == 0] <- NA

# Follow up intensively with a random half of the initial non-responders
Attempt <- rep(0, N)
Attempt[R1 == 0] <- rbinom(sum(R1 == 0), 1, 0.5)
R2 <- rep(0, N)
R2[Attempt == 1] <- rbinom(sum(Attempt == 1), 1, 0.9)
Y[Attempt == 1 & R2 == 1] <- Y_star[Attempt == 1 & R2 == 1]
df <- data.frame(Y, Z, R1, Attempt, R2)

# delta = 1 reproduces the worst-case double-sampling bounds
estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 1, data = df)
#>      ci_lower      ci_upper       low_est       upp_est       low_var 
#> -0.0003261488  0.5435113299  0.1356982015  0.4100265693  0.0068387546 
#>       upp_var 
#>  0.0065857778 

# delta = 0 assumes ignorability among follow-up non-responders
estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 0, data = df)
#>    ci_lower    ci_upper     low_est     upp_est     low_var     upp_var 
#> 0.128264350 0.433193885 0.280729118 0.280729118 0.006051218 0.006051218 
```
