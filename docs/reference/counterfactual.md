# Counterfactual prediction

Predict outcomes under hypothetical covariate values, integrating over
phase uncertainty.

## Usage

``` r
counterfactual(object, newdata, n_sims = NULL)
```

## Arguments

- object:

  A `phase_fit` object

- newdata:

  Data frame with counterfactual covariate values

- n_sims:

  Number of posterior draws to use

## Value

A list with:

- `mean`: Mean predicted outcome per observation

- `sd`: SD of predicted outcome per observation

- `samples`: Matrix of samples (rows = obs, cols = draws)

## Examples

``` r
if (FALSE) { # \dontrun{
# Predict under treatment = 1 for all
cf_data <- data
cf_data$treatment <- 1
cf <- counterfactual(fit, newdata = cf_data)
} # }
```
