# Tidy a trimming bounds object

Returns a three-row tibble with the same columns as
\[tidy.attrition_bounds()\]. Standard errors and the joint Imbens-Manski
confidence interval are filled in when \[estimator_trim()\] was called
with \`se = "analytic"\` or \`se = "bootstrap"\`, and are \`NA\` when it
was called with \`se = "none"\` or when monotonicity failed.

## Usage

``` r
# S3 method for class 'attrition_trim'
tidy(x, ...)
```

## Arguments

- x:

  An object of class \`"attrition_trim"\` (produced by
  \[estimator_trim()\]).

- ...:

  Unused; included for S3 compatibility.

## Value

A \[tibble::tibble()\] with columns \`term\`, \`estimate\`,
\`std.error\`, \`conf.low\`, \`conf.high\`, \`estimate_lower\`,
\`estimate_upper\`, \`std.error_lower\`, \`std.error_upper\`,
\`outcome\`.
