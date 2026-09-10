# Extreme Value (Manski) Bounds

Bounds the average treatment effect when some outcomes are missing and
nothing is assumed about why. Filling every missing outcome in the
treatment group with the lowest value the outcome can take and every
missing outcome in the control group with the highest gives the smallest
average effect the data can support; reversing the fills gives the
largest. Reports the resulting identification region with a joint
Imbens-Manski confidence interval.

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

A named numeric vector with elements `estimate_lower` and
`estimate_upper`, the two ends of the identification region;
`std.error_lower` and `std.error_upper`, their standard errors; and
`conf.low` and `conf.high`, the joint Imbens-Manski confidence interval.
The names are those of the `bounds` row of
[`tidy()`](https://alexandercoppock.com/attrition/reference/tidy.attrition_bounds.md),
which returns the same quantities as a data frame.

## References

Manski, Charles F. (1990). Nonparametric Bounds on Treatment Effects.
*American Economic Review Papers and Proceedings* 80(2):319-323.

Imbens, Guido W., and Charles F. Manski (2004). Confidence Intervals for
Partially Identified Parameters. *Econometrica* 72(6):1845-1857.
[doi:10.1111/j.1468-0262.2004.00555.x](https://doi.org/10.1111/j.1468-0262.2004.00555.x)

Miratrix, Luke W., Jasjeet S. Sekhon, and Bin Yu (2013). Adjusting
Treatment Effect Estimates by Post-Stratification in Randomized
Experiments. *Journal of the Royal Statistical Society, Series B*
75(2):369-396.
[doi:10.1111/j.1467-9868.2012.01048.x](https://doi.org/10.1111/j.1467-9868.2012.01048.x)

Coppock, Alexander, Alan S. Gerber, Donald P. Green, and Holger L. Kern
(2017). Combining Double Sampling and Bounds to Address Nonignorable
Missing Outcomes in Randomized Experiments. *Political Analysis*
25(2):188-206.
[doi:10.1017/pan.2016.6](https://doi.org/10.1017/pan.2016.6)

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
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>     -0.78261705      1.24723513      0.08283104      0.07795614     -0.91886198 
#>       conf.high 
#>      1.37546156 

# Equivalently, via the formula interface
estimator_ev(Y ~ Z, R = "R", minY = 1, maxY = 5, data = df)
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>     -0.78261705      1.24723513      0.08283104      0.07795614     -0.91886198 
#>       conf.high 
#>      1.37546156 
```
