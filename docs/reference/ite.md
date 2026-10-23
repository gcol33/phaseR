# Individual treatment effects

The effect of the treatment on the expected outcome, as in
[`ate()`](https://gillescolling.com/phaseR/reference/ate.md), averaged
over the observations of each unit.

## Usage

``` r
ite(object, treatment_var, treat_value = 1, control_value = 0, n_sims = NULL)
```

## Arguments

- object:

  A `phase_fit` object

- treatment_var:

  Name of the treatment variable

- treat_value:

  Value representing treatment (default 1)

- control_value:

  Value representing control (default 0)

- n_sims:

  Number of posterior draws (default: all)

## Value

A data frame with one row per unit: `id` and the posterior mean, SD and
95% interval of its effect
