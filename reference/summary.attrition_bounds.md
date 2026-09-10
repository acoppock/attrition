# Summarize bounds

Reports the identification region, its standard errors and the joint
Imbens-Manski interval with a line naming the estimand and the
assumptions that produced them. Where \[print()\] shows the returned
vector, \`summary()\` says what the numbers are of.

## Usage

``` r
# S3 method for class 'attrition_bounds'
summary(object, ...)
```

## Arguments

- object:

  An object of class \`"attrition_bounds"\`, produced by
  \[estimator_ev()\], \[estimator_ds()\], or \[estimator_ds_sens()\].

- ...:

  Unused; included for S3 compatibility.

## Value

The \[tidy()\] data frame, invisibly. The report is printed.

## Examples

``` r
summary(estimator_ev(Y = Y_polarization_w2, Z = Z, R = R1,
                     minY = 0, maxY = 6, data = levendusky_replication))
#> Extreme value (Manski) bounds on Y_polarization_w2
#> Estimand: the average treatment effect among all subjects
#> 
#>   Identification region [-1.539, 1.710]
#>   Standard errors       0.079, 0.077
#>   95% Imbens-Manski CI  [-1.669, 1.836]
#> 
#> Assuming the outcome lies within the stated minimum and maximum.
```
