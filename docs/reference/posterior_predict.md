# Posterior predictive simulation

Simulate outcomes from the posterior predictive distribution.

## Usage

``` r
posterior_predict(
  object,
  newdata = NULL,
  n_sims = NULL,
  type = c("response", "phase")
)
```

## Arguments

- object:

  A `phase_fit` object

- newdata:

  Optional new data frame. If NULL, uses original data.

- n_sims:

  Number of posterior draws to use (default: all)

- type:

  "response" for outcome predictions, "phase" for phase probabilities

## Value

A matrix of predictions (rows = observations, cols = posterior samples)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_phaseR(model, data)
preds <- posterior_predict(fit)
} # }
```
