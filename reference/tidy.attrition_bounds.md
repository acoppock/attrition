# Tidy an attrition bounds object

Returns a one-row tibble representing the identified interval. Designed
for use with DeclareDesign: \`estimate\` and \`std.error\` are \`NA\`
because bounds do not yield a single point estimate, while
\`conf.low\`/\`conf.high\` hold the Imbens-Manski joint confidence
interval and \`estimate.low\`/\`estimate.high\` hold the bound point
estimates.

## Usage

``` r
# S3 method for class 'attrition_bounds'
tidy(x, ...)
```

## Arguments

- x:

  An object of class \`"attrition_bounds"\` (produced by
  \[estimator_ev()\] or \[estimator_ds()\]).

- ...:

  Unused; included for S3 compatibility.

## Value

A \[tibble::tibble()\] with columns \`term\`, \`estimate\`,
\`std.error\`, \`conf.low\`, \`conf.high\`, \`estimate.low\`,
\`estimate.high\`.
