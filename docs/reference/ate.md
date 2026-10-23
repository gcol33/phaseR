# Average treatment effect

The effect of setting `treatment_var` to `treat_value` against
`control_value` for every observation, on the expected outcome given the
covariates: for each posterior draw, the expected outcome integrates
over the phase each observation is in, so a treatment that changes the
transition probabilities acts through the phase as well as through the
dynamics. The effect is averaged over the observations of the data.

## Usage

``` r
ate(object, treatment_var, treat_value = 1, control_value = 0, n_sims = NULL)
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

A list with:

- `ate`: Posterior mean of ATE

- `ate_sd`: Posterior SD of ATE

- `ate_ci`: 95% credible interval

- `samples`: Vector of ATE samples, one per draw

## Examples

``` r
if (FALSE) { # \dontrun{
effect <- ate(fit, treatment_var = "treatment")
} # }
```
