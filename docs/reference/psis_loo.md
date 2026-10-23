# Pareto-Smoothed Importance Sampling for LOO

Core PSIS-LOO algorithm. Computes smoothed importance weights and LOO
estimates from a log-likelihood matrix.

## Usage

``` r
psis_loo(ll_matrix)
```

## Arguments

- ll_matrix:

  Matrix of log-likelihoods (draws x observations)

## Value

List with elpd_loo, p_loo, and pareto_k vectors
