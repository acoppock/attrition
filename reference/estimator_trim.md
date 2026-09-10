# Trimming Bounds

Bounds the average treatment effect among the subjects who would report
an outcome under either assignment, the always-reporters, by trimming a
tail of the arm with respondents to spare. The outcome need not be
bounded, which is what distinguishes this from
[`estimator_ev`](https://alexandercoppock.com/attrition/reference/estimator_ev.md).
Two separate choices shape the estimate: which response arguments are
supplied picks the design, and `monotonicity` picks the selection
assumption, which decides which arm is trimmed and by how much.

## Usage

``` r
estimator_trim(
  Y,
  Z,
  R = NULL,
  R1 = NULL,
  Attempt = NULL,
  R2 = NULL,
  monotonicity = c("treatment_increases_response", "treatment_decreases_response",
    "none"),
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

- monotonicity:

  The selection assumption, which is separate from the choice of design
  and is available on both paths. `"treatment_increases_response"`
  assumes \\R_i(1) \ge R_i(0)\\, which makes the control respondents the
  always-reporters and trims the treatment group;
  `"treatment_decreases_response"` assumes the reverse and trims the
  control group; `"none"` assumes neither and trims both groups, by the
  largest share of each that could fail to be always-reporters. The
  default is `"treatment_increases_response"` on the single-stage `R`
  path, which is Lee (2009), and `"none"` on the double-sampling path,
  which is what the follow-up makes affordable; either default gives way
  to an explicit value. The assumption is the researcher's to make
  rather than the data's to choose, so nothing here picks a direction
  from the observed response rates.

- strata:

  Not supported; supplying any value raises an error.

- alpha:

  The desired significance level. 0.05 by default.

- se:

  How to obtain standard errors. `"analytic"` (the default) uses the
  closed-form asymptotic variance of Lee (2009), Proposition 3, which
  covers the single-stage, unweighted, one-group-trimmed case and so is
  available only for the single-stage path under an assumed direction.
  `"bootstrap"` resamples units within treatment arm and works
  everywhere. `"none"` returns bounds alone. Asking for analytic
  standard errors where they do not apply is an error rather than a
  silent substitution.

- sims:

  Number of bootstrap replicates when `se = "bootstrap"`. 1000 by
  default.

- data:

  A dataframe. Must be given by name: `data` is the last argument, so
  passing it positionally assigns it to another argument.

## Value

A named numeric vector leading with the same six elements as the
bounding estimators, in the same order whichever path was taken:
`estimate_lower` and `estimate_upper`, the two trimming bounds;
`std.error_lower` and `std.error_upper`, their standard errors; and
`conf.low` and `conf.high`, the joint Imbens-Manski confidence interval.
The intermediate quantities used to build them follow, and differ
between the two paths. All six lead elements are `NA`, with a warning
naming the other direction, when monotonicity is violated in the
direction assumed. Pass to
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

The design and the assumption are separate choices, and `estimator_trim`
takes them separately. Which response arguments are supplied picks the
design: `R` is the single-stage estimator, `R1`, `Attempt` and `R2` the
double-sampling one, which recovers outcomes from a random sample of the
nonrespondents and so has many fewer subjects left to trim for.
`monotonicity` picks the assumption, in either design.

Monotonicity has a direction, and the two directions are not two ways of
writing the same assumption: each names a different group as the
always-reporters and trims the other, so they give different bounds on
the same data, and only one of them is consistent with any given pair of
response rates. The estimator under the reverse direction is the forward
estimator run with the arms relabelled, so the bounds it returns are
negated and swapped back onto the original contrast.

`monotonicity = "none"` assumes only random assignment. The share of
always-reporters is then bounded below by the Frechet-Hoeffding bound
\\1 - f_0 - f_1\\, where \\f_z\\ is the missingness rate in arm \\z\\,
and each arm is trimmed by the largest share of its respondents that
could fail to be always-reporters: \\f_0/(1 - f_1)\\ of the treatment
group and \\f_1/(1 - f_0)\\ of the control group. These are the sharp
bounds of Imai (2008), Proposition 1, building on Zhang and Rubin (2003)
and Horowitz and Manski (1995). They exist only while \\f_0 + f_1 \<
1\\; beyond that nothing keeps the always-reporter share away from zero,
and the estimator says so rather than returning a number. Because both
arms are trimmed by the same rule, this case has no direction to set and
is invariant to which arm is called treatment.

## References

Horowitz, Joel L., and Charles F. Manski (1995). Identification and
Robustness with Contaminated and Corrupted Data. *Econometrica*
63(2):281-302.

Imai, Kosuke (2008). Sharp Bounds on the Causal Effects in Randomized
Experiments with "Truncation-by-Death". *Statistics & Probability
Letters* 78(2):144-149.

Lee, David S. (2009). Training, Wages, and Sample Selection: Estimating
Sharp Bounds on Treatment Effects. *Review of Economic Studies*
76(3):1071-1102.

Tauchmann, Harald (2014). Lee (2009) Treatment-Effect Bounds for
Nonrandom Sample Selection. *Stata Journal* 14(4):884-894.
[doi:10.1177/1536867X1401400411](https://doi.org/10.1177/1536867X1401400411)

Zhang, Junni L., and Donald B. Rubin (2003). Estimation of Causal
Effects via Principal Stratification When Some Outcomes are Truncated by
"Death". *Journal of Educational and Behavioral Statistics*
28(4):353-368.

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
estimator_trim(Y = Y, Z = Z, R = R, data = df)
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>      0.01761468      0.64239649      0.09007908      0.09972250     -0.13055222 
#>       conf.high       Out0_mono      Out1L_mono      Out1U_mono control_group_N 
#>      0.80642541      2.93353474      2.95114943      3.57593123    331.00000000 
#>   treat_group_N               Q              f1              f0          pi_r_1 
#>    417.00000000      0.16385742      0.18713450      0.32032854      0.81286550 
#>          pi_r_0              yU              yL 
#>      0.67967146      2.00000000      4.00000000 

# Bootstrap standard errors instead
estimator_trim(Y = Y, Z = Z, R = R, se = "bootstrap", sims = 200, data = df)
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>      0.01761468      0.64239649      0.08467663      0.09880524     -0.12166598 
#>       conf.high       Out0_mono      Out1L_mono      Out1U_mono control_group_N 
#>      0.80491665      2.93353474      2.95114943      3.57593123    331.00000000 
#>   treat_group_N               Q              f1              f0          pi_r_1 
#>    417.00000000      0.16385742      0.18713450      0.32032854      0.81286550 
#>          pi_r_0              yU              yL 
#>      0.67967146      2.00000000      4.00000000 

# The other direction, on data built the other way round: here treatment
# lowers response, so the control group holds the extra respondents and is
# the group trimmed
R_rev <- rbinom(N, 1, prob = 0.8 - 0.1 * Z)
Y_rev <- Y_star
Y_rev[R_rev == 0] <- NA
df_rev <- data.frame(Y = Y_rev, Z = Z, R = R_rev)

estimator_trim(Y = Y, Z = Z, R = R,
               monotonicity = "treatment_decreases_response", data = df_rev)
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>      0.10825411      0.64774281      0.09089730      0.09435331     -0.04125865 
#>       conf.high       Out1_mono      Out0L_mono      Out0U_mono control_group_N 
#>      0.80294018      3.30056180      3.19230769      2.65281899    397.00000000 
#>   treat_group_N               Q              f1              f0          pi_r_1 
#>    356.00000000      0.14872263      0.30604288      0.18480493      0.69395712 
#>          pi_r_0              yU              yL 
#>      0.81519507      4.00000000      2.00000000 

# No direction at all: both groups trimmed, randomization the only assumption.
# Wider, and available in either design.
estimator_trim(Y = Y, Z = Z, R = R, monotonicity = "none",
               se = "bootstrap", sims = 200, data = df)
#>  estimate_lower  estimate_upper std.error_lower std.error_upper        conf.low 
#>      -0.8567460       1.4592918       0.1331254       0.1180387      -1.0757179 
#>       conf.high           Out0L           Out0U           Out1L           Out1U 
#>       1.6534482       2.4142259       3.4083333       2.5515873       3.8735178 
#> control_group_N   treat_group_N           trim0           trim1              f1 
#>     331.0000000     417.0000000       0.2753308       0.3940732       0.1871345 
#>              f0          pi_r_1          pi_r_0 
#>       0.3203285       0.8128655       0.6796715 
```
