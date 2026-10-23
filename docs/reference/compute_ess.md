# Compute effective sample size

Multi-chain effective sample size of Geyer's initial monotone sequence
estimator, as in Stan: autocovariances are averaged over chains,
combined with the between-chain variance, and summed in pairs of lags
until a pair turns negative.

## Usage

``` r
compute_ess(draws, chain_id)
```

## Arguments

- draws:

  Matrix of draws, one column per parameter

- chain_id:

  Chain of each row of `draws`

## Value

Numeric vector of effective sample sizes, one per column
