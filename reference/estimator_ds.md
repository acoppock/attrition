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

Manski, Charles F. (1990). Nonparametric Bounds on Treatment Effects.
*American Economic Review Papers and Proceedings* 80(2):319-323.

Neyman, Jerzy (1938). Contribution to the Theory of Sampling Human
Populations. *Journal of the American Statistical Association*
33(201):101-116.
[doi:10.1080/01621459.1938.10503378](https://doi.org/10.1080/01621459.1938.10503378)

Hansen, Morris H., and William N. Hurwitz (1946). The Problem of
Non-Response in Sample Surveys. *Journal of the American Statistical
Association* 41(236):517-529.
[doi:10.1080/01621459.1946.10501894](https://doi.org/10.1080/01621459.1946.10501894)

Miratrix, Luke W., Jasjeet S. Sekhon, and Bin Yu (2013). Adjusting
Treatment Effect Estimates by Post-Stratification in Randomized
Experiments. *Journal of the Royal Statistical Society, Series B*
75(2):369-396.
[doi:10.1111/j.1467-9868.2012.01048.x](https://doi.org/10.1111/j.1467-9868.2012.01048.x)

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
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>      0.18156835      0.30409948      0.08044860      0.08056503      0.04865040 
#>       conf.high 
#>      0.43720979 

# With post-stratification
estimator_ds(Y, Z, R1, Attempt, R2, minY=1, maxY=5, strata=strata, data=df)
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>      0.24698753      0.37580136      0.06847708      0.06507285      0.13421353 
#>       conf.high 
#>      0.48296897 
```
