# Laplace approximation backend for phase models

Fast approximate inference using the Laplace approximation. Finds the
posterior mode and approximates the posterior as Gaussian with
covariance equal to the inverse Hessian.

## Usage

``` r
fit_laplace(stan_data, n_samples = 1000, seed, ...)
```

## Arguments

- stan_data:

  Prepared data list

- n_samples:

  Number of samples to draw from Gaussian approximation

- seed:

  Integer seed for the draws from the Gaussian approximation

- ...:

  Additional arguments (ignored)

## Value

Raw fit object with same structure as HMC backend
