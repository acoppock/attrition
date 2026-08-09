# Print bounds

Prints the named vector of results on its own. \`print.default()\` would
otherwise append the class vector to every result, which is noise.

## Usage

``` r
# S3 method for class 'attrition_bounds'
print(x, ...)
```

## Arguments

- x:

  An object of class \`"attrition_bounds"\`, produced by
  \[estimator_ev()\], \[estimator_ds()\], or \[estimator_ds_sens()\].

- ...:

  Passed to \`print.default()\`.

## Value

\`x\`, invisibly.
