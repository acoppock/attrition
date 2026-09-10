# Summarize trimming bounds

Reports the identification region, its standard errors and the joint
Imbens-Manski interval, together with the design and the selection
assumption that produced them. The two are separate choices in
\[estimator_trim()\] and both change what the numbers mean, so both are
named here.

## Usage

``` r
# S3 method for class 'attrition_trim'
summary(object, ...)
```

## Arguments

- object:

  An object of class \`"attrition_trim"\`, produced by
  \[estimator_trim()\].

- ...:

  Unused; included for S3 compatibility.

## Value

The \[tidy()\] data frame, invisibly. The report is printed.

## Examples

``` r
summary(estimator_trim(Y = Y_polarization_w2, Z = Z, R = R1,
                       monotonicity = "treatment_decreases_response",
                       data = levendusky_replication))
#> Trimming bounds on Y_polarization_w2
#> Estimand: the average treatment effect among always-reporters
#> 
#>   Identification region [0.079, 0.163]
#>   Standard errors       0.102, 0.095
#>   95% Imbens-Manski CI  [-0.095, 0.326]
#> 
#> Design:     single sample
#> Assumption: treatment never raised the chance of responding, so the treated
#>             respondents are the always-reporters and the control group is
#>             trimmed
#> Trimmed:    1.5% of the control group
#> Std errors: analytic (Lee 2009, Proposition 3)
```
