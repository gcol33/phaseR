# Simulate data from a k-phase model

Simulate data from a k-phase model

## Usage

``` r
sim_phaseR_k(
  n_units = 100,
  n_times = 8,
  k = 3,
  beta_trans = NULL,
  beta_dyn = NULL,
  sigma = NULL,
  seed = NULL
)
```

## Arguments

- n_units:

  Number of units

- n_times:

  Number of time points per unit

- k:

  Number of phases (must be \>= 2)

- beta_trans:

  List of transition coefficient vectors (length k-1)

- beta_dyn:

  List of dynamics coefficient vectors (length k)

- sigma:

  Vector of residual SDs (length k)

- seed:

  Random seed

## Value

A data frame suitable for
[`fit_phaseR()`](https://gillescolling.com/phaseR/reference/fit_phaseR.md)

## Examples

``` r
# 3-phase model
dat <- sim_phaseR_k(n_units = 50, n_times = 8, k = 3)
table(dat$true_phase)
```
