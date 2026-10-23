# Counterfactual prediction

Simulates outcomes under hypothetical covariate values with
[`posterior_predict()`](https://gillescolling.com/phaseR/reference/posterior_predict.md):
every unit moves through the phases with the transition probabilities of
`newdata` and draws its outcomes from the dynamics of its current phase.

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

- `lower`, `upper`: 2.5% and 97.5% quantiles per observation

- `samples`: Matrix of samples (rows = rows of `newdata`, cols = draws)

## Examples

``` r
if (FALSE) { # \dontrun{
# Predict under treatment = 1 for all
cf_data <- data
cf_data$treatment <- 1
cf <- counterfactual(fit, newdata = cf_data)
} # }
```
