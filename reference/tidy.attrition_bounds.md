# Tidy an attrition bounds object

Returns a three-row tibble. The \`bounds\` row is the estimator's own
named vector transcribed, so \`estimate_lower\`, \`estimate_upper\`,
\`std.error_lower\`, \`std.error_upper\`, \`conf.low\` and \`conf.high\`
are the same names in both places. \`estimate\` and \`std.error\` are
\`NA\` on that row, because bounds do not yield a single point estimate;
the \`lower_bound\` and \`upper_bound\` rows carry one endpoint each in
broom's long form, which is what DeclareDesign selects on with \`term\`.

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
\`std.error\`, \`conf.low\`, \`conf.high\`, \`estimate_lower\`,
\`estimate_upper\`, \`std.error_lower\`, \`std.error_upper\`,
\`outcome\`.
