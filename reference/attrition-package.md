# attrition: bounds for experiments with missing outcomes

When subjects go missing from an experiment and the reason they went
missing is related to what their outcome would have been, no amount of
covariate adjustment will fix the problem. The package takes the other
route: rather than assume the missingness away, it reports the range of
average treatment effects consistent with the data, and it offers a
research design that makes that range small enough to be useful.

## The design

Worst-case bounds fill in every missing outcome with the smallest and
largest values the outcome could take. The resulting interval is honest,
and it is usually far too wide to settle anything. Double sampling
narrows it. After the first round of data collection, draw a random
sample of the nonrespondents and pursue them harder: pay more, call
again, send an interviewer. Because those subjects are a random sample
of the nonrespondents, their recovered outcomes stand in for all of
them, and only the residual group who refuse twice needs worst-case
treatment. In the application shipped with the package, chasing 100 of
536 nonrespondents cut the width of the 95 percent confidence interval
from 3.50 to 1.23.

## The estimators

- [`estimator_ev`](https://alexandercoppock.com/attrition/reference/estimator_ev.md):

  Worst-case (Manski) bounds from a single round of data collection.

- [`estimator_ds`](https://alexandercoppock.com/attrition/reference/estimator_ds.md):

  Double-sampling bounds, with analytic variances and Imbens-Manski
  confidence intervals. The estimator of Coppock, Gerber, Green, and
  Kern (2017).

- [`estimator_ds_sens`](https://alexandercoppock.com/attrition/reference/estimator_ds_sens.md):

  Double-sampling bounds at a chosen value of delta, the fraction of
  follow-up nonrespondents for whom ignorability is allowed to fail.

- [`sensitivity_ds`](https://alexandercoppock.com/attrition/reference/sensitivity_ds.md):

  A search over delta for the point at which the confidence interval
  starts to include zero.

- [`estimator_trim`](https://alexandercoppock.com/attrition/reference/estimator_trim.md):

  Lee (2009) trimming bounds, which assume monotone selection instead of
  a bounded outcome.

`estimator_ev`, `estimator_ds`, and `estimator_ds_sens` accept a
`strata` argument for poststratification on a discrete covariate. The
identified set is the same either way; poststratification estimates it
more precisely, and by the law of total variance the asymptotic variance
is no larger. Every estimator has a
[`tidy()`](https://alexandercoppock.com/attrition/reference/tidy.attrition_bounds.md)
method and a formula interface for use with DeclareDesign.

## Where to start

[`vignette("attrition")`](https://alexandercoppock.com/attrition/articles/attrition.md)
walks through the design and all five estimators on the replication data
in
[`levendusky_replication`](https://alexandercoppock.com/attrition/reference/levendusky_replication.md),
reproducing the published table as it goes and drawing the imputation
the worst-case bounds average over.

## References

Coppock, Alexander, Alan S. Gerber, Donald P. Green, and Holger L. Kern
(2017). Combining Double Sampling and Bounds to Address Nonignorable
Missing Outcomes in Randomized Experiments. *Political Analysis*
25(2):188-206.
[doi:10.1017/pan.2016.6](https://doi.org/10.1017/pan.2016.6)

Hansen, Morris H., and William N. Hurwitz (1946). The Problem of
Non-Response in Sample Surveys. *Journal of the American Statistical
Association* 41(236):517-529.
[doi:10.1080/01621459.1946.10501894](https://doi.org/10.1080/01621459.1946.10501894)

Horowitz, Joel L., and Charles F. Manski (1995). Identification and
Robustness with Contaminated and Corrupted Data. *Econometrica*
63(2):281-302. [doi:10.2307/2951627](https://doi.org/10.2307/2951627)

Imai, Kosuke (2008). Sharp Bounds on the Causal Effects in Randomized
Experiments with "Truncation-by-Death". *Statistics & Probability
Letters* 78(2):144-149.
[doi:10.1016/j.spl.2007.05.015](https://doi.org/10.1016/j.spl.2007.05.015)

Imbens, Guido W., and Charles F. Manski (2004). Confidence Intervals for
Partially Identified Parameters. *Econometrica* 72(6):1845-1857.
[doi:10.1111/j.1468-0262.2004.00555.x](https://doi.org/10.1111/j.1468-0262.2004.00555.x)

Lee, David S. (2009). Training, Wages, and Sample Selection: Estimating
Sharp Bounds on Treatment Effects. *Review of Economic Studies*
76(3):1071-1102.
[doi:10.1111/j.1467-937X.2009.00536.x](https://doi.org/10.1111/j.1467-937X.2009.00536.x)

Manski, Charles F. (1990). Nonparametric Bounds on Treatment Effects.
*American Economic Review Papers and Proceedings* 80(2):319-323.

Miratrix, Luke W., Jasjeet S. Sekhon, and Bin Yu (2013). Adjusting
Treatment Effect Estimates by Post-Stratification in Randomized
Experiments. *Journal of the Royal Statistical Society, Series B*
75(2):369-396.
[doi:10.1111/j.1467-9868.2012.01048.x](https://doi.org/10.1111/j.1467-9868.2012.01048.x)

Neyman, Jerzy (1938). Contribution to the Theory of Sampling Human
Populations. *Journal of the American Statistical Association*
33(201):101-116.
[doi:10.1080/01621459.1938.10503378](https://doi.org/10.1080/01621459.1938.10503378)

Zhang, Junni L., and Donald B. Rubin (2003). Estimation of Causal
Effects via Principal Stratification When Some Outcomes are Truncated by
"Death". *Journal of Educational and Behavioral Statistics*
28(4):353-368.
[doi:10.3102/10769986028004353](https://doi.org/10.3102/10769986028004353)

## See also

Useful links:

- <https://alexandercoppock.com/attrition/>

- <https://github.com/acoppock/attrition>

- Report bugs at <https://github.com/acoppock/attrition/issues>

## Author

**Maintainer**: Alexander Coppock <acoppock@gmail.com>

Authors:

- Alexander Coppock <acoppock@gmail.com>

- Alan S. Gerber

- Donald P. Green

- Holger L. Kern
