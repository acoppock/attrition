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

A named numeric vector with elements `estimate_lower` and
`estimate_upper`, the two ends of the identification region;
`std.error_lower` and `std.error_upper`, their standard errors; and
`conf.low` and `conf.high`, the joint Imbens-Manski confidence interval.
The names are those of the `bounds` row of
[`tidy()`](https://alexandercoppock.com/attrition/reference/tidy.attrition_bounds.md),
which returns the same quantities as a data frame.

## References

Coppock, Alexander, Alan S. Gerber, Donald P. Green, and Holger L. Kern
(2017). Combining Double Sampling and Bounds to Address Nonignorable
Missing Outcomes in Randomized Experiments. *Political Analysis*
25(2):188-206.
[doi:10.1017/pan.2016.6](https://doi.org/10.1017/pan.2016.6)

Imbens, Guido W., and Charles F. Manski (2004). Confidence Intervals for
Partially Identified Parameters. *Econometrica* 72(6):1845-1857.
[doi:10.1111/j.1468-0262.2004.00555.x](https://doi.org/10.1111/j.1468-0262.2004.00555.x)

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
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>    0.1356982015    0.4100265693    0.0826967632    0.0811528054   -0.0003261488 
#>       conf.high 
#>    0.5435113299 

# delta = 0 assumes ignorability among follow-up non-responders
estimator_ds_sens(Y, Z, R1, Attempt, R2, minY = 1, maxY = 5, delta = 0, data = df)
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>      0.28072912      0.28072912      0.07778958      0.07778958      0.12826435 
#>       conf.high 
#>      0.43319388 
```
