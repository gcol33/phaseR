# Simulate data from a phase model

Simulate data from a phase model

## Usage

``` r
sim_phaseR(
  n_units = 100,
  n_times = 5,
  beta_trans = c(-1, 0.5),
  beta_0 = c(10, 0),
  beta_1 = c(8, -0.5),
  sigma_0 = 1,
  sigma_1 = 1.5,
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

  Phase 0 dynamics coefficients

- beta_1:

  Phase 1 dynamics coefficients

- sigma_0:

  Phase 0 residual SD

- sigma_1:

  Phase 1 residual SD

- seed:

  Random seed

## Value

A data frame suitable for
[`fit_phaseR()`](https://gillescolling.com/phaseR/reference/fit_phaseR.md)

## Examples

``` r
dat <- sim_phaseR(n_units = 100, n_times = 5)
head(dat)
```
