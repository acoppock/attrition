
<!-- README.md is generated from README.Rmd. Please edit that file -->

# attrition

Addressing nonignorable attrition with double sampling and bounds: the
estimators of Coppock, Gerber, Green, and Kern (2017), [*Political
Analysis* 25(2):188-206](https://doi.org/10.1017/pan.2016.6), for
randomized experiments in which some outcomes are missing.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("acoppock/attrition")
```

## When these estimators apply

A randomized experiment has been run and some subjects have no outcome
recorded. A two-wave survey where part of the sample never comes back
for the second wave. A field experiment whose endline cannot reach
everyone. A study where the outcome exists only for subjects who cleared
an earlier hurdle, as wages do for the subset who found work.

Each estimator reports the range of average treatment effects consistent
with the data. Which one fits depends on what researchers are able to
assume about the subjects they did not observe:

- The outcome has a known minimum and maximum, such as a 0 to 6 scale, a
  binary indicator, or a bounded index: `estimator_ev()`.
- The same, and a random sample of the nonrespondents was pursued in a
  second round: `estimator_ds()`, which is what the paper is about.
- No minimum and maximum can be fixed, but treatment can be assumed to
  have only raised a subject’s chance of responding, never lowered it:
  `estimator_trim()`.
- The double-sampling design was run, and the question is what fraction
  of the follow-up nonrespondents would have to violate ignorability
  before the finding stops holding: `estimator_ds_sens()` and
  `sensitivity_ds()`.

## The estimators at a glance

| Function | What it assumes |
|----|----|
| `estimator_ev()` | The outcome has known lower and upper limits (Manski 1990). |
| `estimator_ds()` | The same, plus a random follow-up sample of nonrespondents. |
| `estimator_ds_sens()` | The same, with ignorability allowed to fail for a fraction `delta` of the follow-up nonrespondents. |
| `sensitivity_ds()` | A search over `delta` for the point where the interval starts to include zero. |
| `estimator_trim()` | Treatment moves response in one direction only (Lee 2009). The outcome need not be bounded. |

The first four take a `strata` argument for poststratification on a
discrete covariate, which targets the same identification region and
estimates it more precisely; `estimator_trim()` does not. All five have
`tidy()` methods and a formula interface for use with
[DeclareDesign](https://declaredesign.org).

## When a double-sampling design is worth running

The decision belongs before fielding, because the design is the part
that does the work. Double sampling is worth planning whenever
substantial attrition is expected and the nonrespondents could be
reached by spending more on them than the first round spent: a larger
incentive, more callbacks, an in-person visit, a switch from web to
phone.

The procedure is to close the first round, draw a random sample of
whoever did not respond, and pursue that sample hard. Drawing at random
is what makes the recovered outcomes stand in for every nonrespondent
rather than for the subset who happen to be easy to reach, so only the
subjects who refuse twice remain unknown.

The budget for it belongs to design rather than to analysis. A larger
initial sample buys precision and does nothing about attrition; a
follow-up sample reduces the attrition bias directly.

## What double sampling buys

The package ships the replication study from the paper as
`levendusky_replication`: a two-wave survey experiment in which 1,980
subjects were asked about perceived polarization, and 536 of them did
not answer the second wave. A third condition read nothing on the topic
and is not analyzed, which is what the subset below removes: `Z1` is
defined only for the two conditions being compared, so filtering on its
missingness drops the unused group and nothing else.

Refusing any assumption about the missing outcomes gives worst-case
bounds. Because the outcome runs from 0 to 6, filling every missing
value in one group with 0 and the other with 6 gives the lowest effect
the data can support, and reversing the fills gives the highest.

``` r
library(attrition)
dat <- droplevels(subset(levendusky_replication, !is.na(Z1)))

estimator_ev(L_dif_w2, Z1, R1, minY = 0, maxY = 6, data = dat)
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>        -1.53914         1.70967         0.07899         0.07673        -1.66908 
#>       conf.high 
#>         1.83588
```

The effect lies somewhere between -1.54 and 1.71, which is honest and
nearly useless. Little of that width is sampling error: the confidence
interval is barely wider than the bounds. A larger sample would not have
helped, because the problem is 536 unknown outcomes rather than noise.

Double sampling addresses the unknown outcomes directly. After the first
round of data collection, a random sample of the nonrespondents is drawn
and pursued with more effort than the first attempt. Here, 50
nonrespondents per condition were offered \$4.00 rather than the
original \$1.00, and 72 of those 100 answered. Because they are a random
sample of the nonrespondents, their outcomes stand in for all 536, and
only the 28 who refused twice still need worst-case treatment.

``` r
estimator_ds(L_dif_w2, Z1, R1, Attempt, R2, minY = 0, maxY = 6, data = dat)
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>         -0.3417          0.5718          0.1134          0.1054         -0.5283 
#>       conf.high 
#>          0.7452
```

The identification region shrinks by a factor of 3.6, from 3.25 points
wide to 0.91.

## Reading the output

Every estimator returns the same six named elements: `estimate_lower`
and `estimate_upper`, the two ends of the identification region;
`std.error_lower` and `std.error_upper`, their standard errors; and
`conf.low` and `conf.high`, the joint Imbens-Manski interval. Those are
broom’s names, and `tidy()` returns the same six quantities as a data
frame under the same names.

``` r
tidy(estimator_ds(L_dif_w2 ~ Z1, R1 = "R1", Attempt = "Attempt", R2 = "R2",
                  minY = 0, maxY = 6, data = dat))
#> # A tibble: 3 × 10
#>   term       estimate std.error conf.low conf.high estimate_lower estimate_upper
#>   <chr>         <dbl>     <dbl>    <dbl>     <dbl>          <dbl>          <dbl>
#> 1 bounds       NA        NA       -0.528     0.745         -0.342          0.572
#> 2 lower_bou…   -0.342     0.113   NA        NA             NA             NA    
#> 3 upper_bou…    0.572     0.105   NA        NA             NA             NA    
#> # ℹ 3 more variables: std.error_lower <dbl>, std.error_upper <dbl>,
#> #   outcome <chr>
```

Bounds have no single point estimate, so `estimate` is `NA` on the
`bounds` row, which carries the whole vector across its columns. The two
rows below split it, one endpoint each, so `estimate` and `std.error`
mean there what broom means by them and `declare_estimator()` can select
an endpoint with `term`.

## Learning more

`vignette("attrition")` works through the design and all five estimators
on the shipped data, reproducing the paper’s Table 3 along the way.
`vignette("drawing-the-bounds")` draws the imputation that
`estimator_ev()` averages over and checks the picture against the
estimates.

## AI statement

attrition 1.0.0 was prepared by Alexander Coppock working with Claude
(Anthropic), across the rewrite of the estimators, the test suite and
the documentation. The method and the original implementation come from
Coppock, Gerber, Green, and Kern (2017). The 1.0.0 release rewrote the
internals, added a formula interface and `tidy()` methods, corrected
several defects in the earlier code, and wrote both vignettes.

The code base has been reviewed but not written line by line, so the
guarantee offered is not that every line has been vouched for. It is
that the package reproduces the published analysis it implements.
`vignette("attrition")` reproduces Table 3 of the paper, and the test
suite holds each estimator to the published quantities, which have not
changed. The analytic variance added to `estimator_trim()` was checked
three ways: against the algebraically distinct form in Tauchmann’s Stata
`leebounds`, against the Monte Carlo standard deviation of the estimator
across sample sizes from 1,000 to 64,000, and by the coverage of a
conventional interval around each bound endpoint.

## References

Coppock, Alexander, Alan S. Gerber, Donald P. Green, and Holger L. Kern
(2017). Combining Double Sampling and Bounds to Address Nonignorable
Missing Outcomes in Randomized Experiments. *Political Analysis*
25(2):188-206. <https://doi.org/10.1017/pan.2016.6>

Imbens, Guido W., and Charles F. Manski (2004). Confidence Intervals for
Partially Identified Parameters. *Econometrica* 72(6):1845-1857.
<https://doi.org/10.1111/j.1468-0262.2004.00555.x>

Lee, David S. (2009). Training, Wages, and Sample Selection: Estimating
Sharp Bounds on Treatment Effects. *Review of Economic Studies*
76(3):1071-1102. <https://doi.org/10.1111/j.1467-937X.2009.00536.x>

Manski, Charles F. (1990). Nonparametric Bounds on Treatment Effects.
*American Economic Review Papers and Proceedings* 80(2):319-323.
