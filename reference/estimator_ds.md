# Extreme Value Bounds with Double Sampling

Extreme Value Bounds with Double Sampling

## Usage

``` r
estimator_ds(
  Y,
  Z,
  R1,
  Attempt,
  R2,
  minY,
  maxY,
  strata = NULL,
  alpha = 0.05,
  data
)
```

## Arguments

- Y:

  The (unquoted) outcome variable, or a formula `outcome ~ treatment`
  for use with `declare_estimator(.method = estimator_ds)`. Must be
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

- strata:

  Stratification variable: unquoted column name or a quoted string
  column name.

- alpha:

  The desired significance level. 0.05 by default.

- data:

  A dataframe. Must be given by name: `data` is the last argument, so
  passing it positionally assigns it to another argument.

## Value

A named numeric vector with elements `ci_lower` and `ci_upper`, the
joint Imbens-Manski confidence interval; `low_est` and `upp_est`, the
bound point estimates; and `low_var` and `upp_var`, their variances.
Pass to
[`tidy()`](https://alexandercoppock.com/attrition/reference/tidy.attrition_bounds.md)
for a data frame.

## Examples

``` r
set.seed(343) # For reproducibility
N <- 1000

# Potential Outcomes
Y_0 <- sample(1:5, N, replace=TRUE, prob = c(0.1, 0.3, 0.3, 0.2, 0.1))
Y_1 <- sample(1:5, N, replace=TRUE, prob = c(0.1, 0.1, 0.4, 0.3, 0.1))

R1_0 <- rbinom(N, 1, prob = 0.7)
R1_1 <- rbinom(N, 1, prob = 0.8)

R2_0 <- rbinom(N, 1, prob = 0.9)
R2_1 <- rbinom(N, 1, prob = 0.95)

# Covariate
strata <- as.numeric(Y_0 > 2)

# Random Assignment
Z <- rbinom(N, 1, .5)

# Reveal Initial Sample Outcomes
R1 <- Z*R1_1 + (1-Z)*R1_0 # Initial sample response
Y_star <- Z*Y_1 + (1-Z)*Y_0 # True outcomes
Y <- Y_star
Y[R1==0] <- NA # Mask outcome of non-responders

# Conduct Double Sampling
Attempt <- rep(0, N)
Attempt[is.na(Y)] <- rbinom(sum(is.na(Y)), 1, .5)

R2 <- rep(0, N)
R2[Attempt==1] <-  (Z*R2_1 + (1-Z)*R2_0)[Attempt==1]

Y[R2==1 & Attempt==1] <- Y_star[R2==1 & Attempt==1]

df <- data.frame(Y, Z, R1, Attempt, R2, strata)

# Without post-stratification
estimator_ds(Y, Z, R1, Attempt, R2, minY=1, maxY=5, data=df)
#>    ci_lower    ci_upper     low_est     upp_est     low_var     upp_var 
#> 0.048650402 0.437209789 0.181568350 0.304099481 0.006471977 0.006490723 

# With post-stratification
estimator_ds(Y, Z, R1, Attempt, R2, minY=1, maxY=5, strata=strata, data=df)
#>    ci_lower    ci_upper     low_est     upp_est     low_var     upp_var 
#> 0.134213533 0.482968974 0.246987534 0.375801359 0.004689110 0.004234475 
```
