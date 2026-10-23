# Extract variance components

Returns variance-covariance information for random effects. For random
intercepts, returns the SD. For random slopes, returns the full
covariance matrix (correlated) or diagonal variances (uncorrelated).

## Usage

``` r
VarCorr(object)
```

## Arguments

- object:

  A `phase_fit` object

## Value

A list with variance component estimates by term
