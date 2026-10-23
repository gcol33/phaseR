# Mode of the fixed effects with the random effects at zero

Maximises the log posterior over the coefficients and residual SDs with
every random effect at zero and its covariance parameters at their
starting values. Each step costs one evaluation of the log posterior and
its gradient instead of an inner mode search, and it brings a starting
point of the Laplace search close to the mode of the phase structure it
leads to.

## Usage

``` r
fixed_effects_mode(start, target)
```

## Arguments

- start:

  Parameter vector from
  [`initial_values()`](https://gillescolling.com/phaseR/reference/initial_values.md)

- target:

  Model description from
  [`posterior_target()`](https://gillescolling.com/phaseR/reference/posterior_target.md)

## Value

`start` with the coefficients and log residual SDs replaced by their
mode
