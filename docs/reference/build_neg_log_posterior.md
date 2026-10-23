# Build negative log posterior function for optimization

Creates a function that computes -log p(theta \| y) for use with
optim(), from the same model description and prior the NUTS backend
samples.

## Usage

``` r
build_neg_log_posterior(stan_data)
```

## Arguments

- stan_data:

  Prepared data list

## Value

Function of parameter vector returning negative log posterior
