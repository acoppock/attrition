# Print a sensitivity analysis

Reports delta\*, the smallest value of the sensitivity parameter at
which the confidence interval reaches zero, and names the components of
the object. Printing the list itself would draw the plot, which is not
what a glance at the result should do.

## Usage

``` r
# S3 method for class 'attrition_sensitivity'
print(x, ...)
```

## Arguments

- x:

  An object of class \`"attrition_sensitivity"\`, produced by
  \[sensitivity_ds()\].

- ...:

  Unused; included for S3 compatibility.

## Value

\`x\`, invisibly.
