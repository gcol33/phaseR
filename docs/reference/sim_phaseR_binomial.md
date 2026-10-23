# Simulate Binomial data from a phase model

Simulate Binomial data from a phase model

## Usage

``` r
sim_phaseR_binomial(
  n_units = 100,
  n_times = 5,
  n_trials = 10,
  beta_trans = c(-1, 0.5),
  beta_0 = c(0.5, 0),
  beta_1 = c(-0.5, 0.3),
  seed = NULL
)
```

## Arguments

- n_units:

  Number of units

- n_times:

  Number of time points per unit

- n_trials:

  Number of trials per observation

- beta_trans:

  Transition coefficients

- beta_0:

  Phase 0 dynamics coefficients (logit scale)

- beta_1:

  Phase 1 dynamics coefficients (logit scale)

- seed:

  Random seed

## Value

A data frame suitable for
[`fit_phaseR()`](https://gillescolling.com/phaseR/reference/fit_phaseR.md)
with `family = binomial()`

## Examples

``` r
dat <- sim_phaseR_binomial(n_units = 50, n_times = 5, n_trials = 10)
head(dat)
```
