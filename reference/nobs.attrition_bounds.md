# Number of observations

The number of units the estimator was given, counting those whose
outcome is missing. Every variance in the package is an asymptotic
result indexed to the full randomized sample rather than to the
respondents: the worst-case bounds divide by the number assigned to each
arm, and Lee (2009) Proposition 3 is a root-n result in which the
response and trimming rates appear as proportions of that same n. Counts
of respondents, of retained observations after trimming, and of the
groups behind them are returned as named elements of the estimator
output.

## Usage

``` r
# S3 method for class 'attrition_bounds'
nobs(object, ...)

# S3 method for class 'attrition_trim'
nobs(object, ...)
```

## Arguments

- object:

  An object of class \`"attrition_bounds"\` or \`"attrition_trim"\`.

- ...:

  Unused; included for S3 compatibility.

## Value

An integer.
