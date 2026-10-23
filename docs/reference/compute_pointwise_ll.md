# Compute pointwise log-likelihood

The log-likelihood of each observation given the earlier observations of
its unit, marginalised over the latent phases, for every posterior draw
(random effects at their drawn values). The forward recursion is the one
both backends fit with, so every model family is covered.

## Usage

``` r
compute_pointwise_ll(fit)
```

## Arguments

- fit:

  A `phase_fit` object

## Value

Matrix with one row per draw and one column per row of `fit$data`
