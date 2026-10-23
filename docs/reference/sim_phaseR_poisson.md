# Simulate Poisson data from a phase model

Simulate Poisson data from a phase model

## Usage

``` r
sim_phaseR_poisson(
  n_units = 100,
  n_times = 5,
  beta_trans = c(-1, 0.5),
  beta_0 = c(2, 0),
  beta_1 = c(1.5, -0.3),
  seed = NULL
)
```

## Arguments

- n_units:

  Number of units

- n_times:

  Number of time points per unit

- beta_trans:

  Transition coefficients

- beta_0:

  Phase 0 dynamics coefficients (log scale)

- beta_1:

  Phase 1 dynamics coefficients (log scale)

- seed:

  Random seed

## Value

A data frame suitable for
[`fit_phaseR()`](https://gillescolling.com/phaseR/reference/fit_phaseR.md)
with `family = poisson()`

## Examples

``` r
dat <- sim_phaseR_poisson(n_units = 50, n_times = 5)
head(dat)
```
