# Build negative log posterior function for optimization

Creates -log p(theta \| y) and its exact gradient for use with optim(),
from the same model description and prior the NUTS backend samples.

## Usage

``` r
build_neg_log_posterior(stan_data)
```

## Arguments

- stan_data:

  Prepared data list

## Value

List of two functions of the parameter vector: `fn`, the negative log
posterior, and `gr`, its gradient
