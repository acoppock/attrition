# Extreme Value (Manski) Bounds

Extreme Value (Manski) Bounds

## Usage

``` r
estimator_ev(Y, Z, R, minY, maxY, strata = NULL, alpha = 0.05, data)
```

## Arguments

- Y:

  The (unquoted) outcome variable, or a formula `outcome ~ treatment`
  for use with `declare_estimator(.method = estimator_ev)`. Must be
  numeric.

- Z:

  The (unquoted) assignment indicator variable. Must be numeric and take
  values 0 or 1. Ignored when `Y` is a formula.

- R:

  The response indicator variable: unquoted column name, or a quoted
  string column name when using the formula interface. Must be numeric
  and take values 0 or 1.

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

estimator_ev(Y, Z, R, minY = 1, maxY = 5, data = df)
#>     ci_lower     ci_upper      low_est      upp_est      low_var      upp_var 
#> -0.918861985  1.375461565 -0.782617049  1.247235131  0.006860981  0.006077159 

# Equivalently, via the formula interface
estimator_ev(Y ~ Z, R = "R", minY = 1, maxY = 5, data = df)
#>     ci_lower     ci_upper      low_est      upp_est      low_var      upp_var 
#> -0.918861985  1.375461565 -0.782617049  1.247235131  0.006860981  0.006077159 
```
