# Fit a phase model

Fit a phase model

## Usage

``` r
fit_phaseR(
  model,
  data,
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  backend = c("hmc", "laplace"),
  cores = 1,
  prior = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- model:

  A `phase_model` object

- data:

  A data frame with columns: id, time, outcome(s), covariates

- chains:

  Number of MCMC chains (default 4)

- iter:

  Total iterations per chain (default 2000)

- warmup:

  Warmup iterations per chain (default iter/2)

- backend:

  Estimation backend: "hmc" (default), the No-U-Turn sampler with
  step-size and diagonal-metric adaptation during warmup, or "laplace",
  a Gaussian approximation at the posterior mode

- cores:

  Number of cores for parallel chains (default 1)

- prior:

  Prior specification from
  [`prior()`](https://gillescolling.com/phaseR/reference/prior.md)
  (default: weakly informative); both backends evaluate it for every
  model family

- seed:

  Integer seed. Two fits with the same `seed`, data and settings return
  identical draws. With `seed = NULL` the seed is drawn from R's random
  number generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  makes the fit reproducible; the seed used is stored in the returned
  object as `seed`.

- ...:

  Backend settings. For `backend = "hmc"`: `adapt_delta`, the target
  acceptance statistic of step-size adaptation (default 0.8; raise it
  towards 1 when there are divergent transitions), and `max_treedepth`
  (default 10). For `backend = "laplace"`: `n_samples` (default 1000).

## Value

A `phase_fit` object

## Examples

``` r
if (FALSE) { # \dontrun{
model <- phase_model(
  phases("inactive", "active"),
  transition("inactive", "active", ~ x),
  dynamics("inactive", y ~ 1),
  dynamics("active", y ~ x)
)
dat <- sim_phaseR(n_units = 50, n_times = 5)
fit <- fit_phaseR(model, data = dat)
} # }
```
