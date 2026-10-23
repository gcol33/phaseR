# Effects of a treatment on the expected outcome

Effects of a treatment on the expected outcome

## Usage

``` r
treatment_effects(object, treatment_var, treat_value, control_value, n_sims)
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

Matrix with one row per posterior draw and one column per row of the
data: expected outcome under treatment minus under control
