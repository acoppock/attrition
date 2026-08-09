# attrition

Bounds on average treatment effects for randomized experiments with
missing outcomes, implementing Coppock, Gerber, Green, and Kern (2017),
[*Political Analysis*
25(2):188-206](https://doi.org/10.1017/pan.2016.6).

Attrition is the problem covariate adjustment cannot solve. If subjects
leave your experiment for reasons related to the outcomes they would
have reported, the respondents you are left with are not a random sample
of the subjects you started with, and no regression on the covariates
you happened to measure will put that right. The standard response is to
assume the problem away: missingness is ignorable, conditional on
treatment and a few covariates. The assumption is untestable, and when
it fails it fails silently.

The package takes the other route. It reports the full range of
treatment effects consistent with the data, and it implements a research
design, double sampling, that makes that range narrow enough to be worth
reporting.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("acoppock/attrition")
```

## What double sampling buys

The package ships the replication study from the paper: a two-wave
survey experiment in which 1,980 subjects were asked about perceived
polarization, and 536 of them did not answer the second wave.

Refusing any assumption about the missing outcomes gives worst-case
bounds. Because the outcome runs from 0 to 6, filling every missing
value in one group with 0 and the other with 6 gives the lowest effect
the data can support, and reversing the fills gives the highest.

``` r

library(attrition)
dat <- subset(levendusky, !is.na(Z1))

estimator_ev(L_dif_w2, Z1, R1, minY = 0, maxY = 6, data = dat)
#>  ci_lower  ci_upper   low_est   upp_est   low_var   upp_var 
#> -1.669079  1.835881 -1.539145  1.709668  0.006240  0.005888
```

The effect lies somewhere between -1.54 and 1.71, which is honest and
nearly useless. Note how little of that width is sampling error: the
confidence interval is barely wider than the bounds. A larger sample
would not have helped, because the problem is 536 unknown outcomes
rather than noise.

Double sampling attacks the unknown outcomes directly. After the first
round of data collection, draw a random sample of the nonrespondents and
pursue them harder. Here, 50 nonrespondents per condition were offered
\$4.00 rather than the original \$1.00, and 72 of those 100 answered.
Because they are a random sample of the nonrespondents, their outcomes
stand in for all 536, and only the 28 who refused twice still need
worst-case treatment.

``` r

estimator_ds(L_dif_w2, Z1, R1, Attempt, R2, minY = 0, maxY = 6, data = dat)
#> ci_lower ci_upper  low_est  upp_est  low_var  upp_var 
#> -0.52831  0.74517 -0.34175  0.57182  0.01286  0.01111
```

The identification region shrinks by a factor of 3.6, from 3.25 points
wide to 0.91. Chasing 100 subjects bought all of it.

## The estimators

| Function | What it assumes |
|----|----|
| [`estimator_ev()`](https://alexandercoppock.com/attrition/reference/estimator_ev.md) | The outcome has known lower and upper limits. |
| [`estimator_ds()`](https://alexandercoppock.com/attrition/reference/estimator_ds.md) | The same, plus a random follow-up sample of nonrespondents. |
| [`estimator_ds_sens()`](https://alexandercoppock.com/attrition/reference/estimator_ds_sens.md) | The same, with ignorability allowed to fail for a fraction delta of the follow-up nonrespondents. |
| [`sensitivity_ds()`](https://alexandercoppock.com/attrition/reference/sensitivity_ds.md) | A search over delta for the point where the interval starts to include zero. |
| [`estimator_trim()`](https://alexandercoppock.com/attrition/reference/estimator_trim.md) | Treatment moves response in one direction only ([Lee 2009](https://doi.org/10.1111/j.1467-937X.2009.00536.x)). The outcome need not be bounded. |

The first three take a `strata` argument for poststratification on a
discrete covariate, which targets the same identification region but
estimates it more precisely. All five have
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) methods and a
formula interface for use with
[DeclareDesign](https://declaredesign.org).

``` r

tidy(estimator_ds(L_dif_w2 ~ Z1, R1 = "R1", Attempt = "Attempt", R2 = "R2",
                  minY = 0, maxY = 6, data = dat))
#> # A tibble: 3 × 7
#>   term        estimate std.error conf.low conf.high estimate.low estimate.high
#>   <chr>          <dbl>     <dbl>    <dbl>     <dbl>        <dbl>         <dbl>
#> 1 bounds        NA        NA       -0.528     0.745       -0.342         0.572
#> 2 lower_bound   -0.342     0.113   NA        NA           NA            NA    
#> 3 upper_bound    0.572     0.105   NA        NA           NA            NA
```

Bounds have no single point estimate, so `estimate` is `NA` on the
`bounds` row: the identification region sits in `estimate.low` and
`estimate.high`, and the joint Imbens-Manski interval in `conf.low` and
`conf.high`.

## Learning more

[`vignette("attrition")`](https://alexandercoppock.com/attrition/articles/attrition.md)
works through the design and all five estimators on the shipped data,
reproducing the paper’s Table 3 along the way.

## References

Coppock, Alexander, Alan S. Gerber, Donald P. Green, and Holger L. Kern
(2017). Combining Double Sampling and Bounds to Address Nonignorable
Missing Outcomes in Randomized Experiments. *Political Analysis*
25(2):188-206.

Imbens, Guido W., and Charles F. Manski (2004). Confidence Intervals for
Partially Identified Parameters. *Econometrica* 72(6):1845-1857.

Lee, David S. (2009). Training, Wages, and Sample Selection: Estimating
Sharp Bounds on Treatment Effects. *Review of Economic Studies*
76(3):1071-1102.

Manski, Charles F. (1990). Nonparametric Bounds on Treatment Effects.
*American Economic Review Papers and Proceedings* 80(2):319-323.
