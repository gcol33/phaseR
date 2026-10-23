# HMC backend for phase models

Runs the No-U-Turn sampler on the log posterior of the model, one chain
after another or in parallel.

## Usage

``` r
fit_hmc(
  stan_data,
  chains,
  iter,
  warmup,
  cores,
  seed,
  adapt_delta = 0.8,
  max_treedepth = 10L,
  ...
)
```

## Arguments

- stan_data:

  Prepared data list

- chains:

  Number of chains

- iter:

  Iterations per chain

- warmup:

  Warmup iterations

- cores:

  Number of cores

- seed:

  Integer seed; chain `c` samples from the stream seeded by `(seed, c)`,
  which also perturbs its initial values. The initial values are the
  data-informed starting points of
  [`initial_values()`](https://gillescolling.com/phaseR/reference/initial_values.md),
  drawn with `seed`.

- adapt_delta:

  Target average acceptance statistic of the step-size adaptation, in
  (0, 1)

- max_treedepth:

  Maximum tree depth; a transition takes at most `2^max_treedepth - 1`
  leapfrog steps

- ...:

  Unused

## Value

Raw fit object
