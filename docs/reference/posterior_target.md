# Describe the log posterior of a phase model

Both backends evaluate the log posterior through this description: the
NUTS sampler in `phaseR_nuts_sampler()` and the Laplace approximation
through `phaseR_log_posterior()`.

## Usage

``` r
posterior_target(stan_data)
```

## Arguments

- stan_data:

  Prepared data list

## Value

List with `kind` (which C++ log posterior to use), `n_params`, `data`,
`args` (the dimensions and random-effect structures that log posterior
takes) and `prior` (the prior of each parameter class, from
[`prior_to_cpp()`](https://gillescolling.com/phaseR/reference/prior_to_cpp.md))
