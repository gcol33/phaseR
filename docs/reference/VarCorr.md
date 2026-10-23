# Extract variance components

Posterior summaries of the SD of every random-effect coefficient and the
posterior mean of the correlation matrix of each term (the identity for
intercept and uncorrelated terms).

## Usage

``` r
VarCorr(object)
```

## Arguments

- object:

  A `phase_fit` object

## Value

A list with one element per predictor and term, each with `sd` (a data
frame with one row per coefficient) and `cor` (a matrix)
