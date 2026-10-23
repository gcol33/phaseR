# Individual treatment effects

Compute treatment effects for each unit, integrating over phase
uncertainty.

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

  Number of posterior draws

## Value

A data frame with unit-level treatment effects
