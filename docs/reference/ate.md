# Average treatment effect

Compute the average treatment effect by comparing counterfactual
outcomes under treatment vs control conditions.

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

  Number of posterior draws

## Value

A list with:

- `ate`: Posterior mean of ATE

- `ate_sd`: Posterior SD of ATE

- `ate_ci`: 95% credible interval

- `samples`: Vector of ATE samples

## Examples

``` r
if (FALSE) { # \dontrun{
effect <- ate(fit, treatment_var = "treatment")
} # }
```
