# Tidy a sensitivity analysis

Returns the bounds and joint Imbens-Manski interval at every value of
the sensitivity parameter, one row per \`delta\`, under the same names
the estimators return.

## Usage

``` r
# S3 method for class 'attrition_sensitivity'
tidy(x, ...)
```

## Arguments

- x:

  An object of class \`"attrition_sensitivity"\` (produced by
  \[sensitivity_ds()\]).

- ...:

  Unused; included for S3 compatibility.

## Value

A \[tibble::tibble()\] with columns \`delta\`, \`estimate_lower\`,
\`estimate_upper\`, \`std.error_lower\`, \`std.error_upper\`,
\`conf.low\`, \`conf.high\`, \`outcome\`.
