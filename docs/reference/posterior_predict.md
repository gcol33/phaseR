# Posterior predictive simulation

Simulates outcomes from the posterior predictive distribution given the
covariates of `newdata`, or returns phase probabilities given its
outcomes.

## Usage

``` r
posterior_predict(
  object,
  newdata = NULL,
  n_sims = NULL,
  type = c("response", "phase"),
  phase = NULL
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

- phase:

  For `type = "phase"`, the name of the phase whose probability is
  returned (default: the last phase)

## Value

A matrix with one row per row of `newdata` (in its order) and one column
per posterior draw

## Details

For `type = "response"`, each posterior draw simulates every unit
forward: the unit starts in the first phase, moves to the next phase
with the transition probability its covariates imply at each later
observation, and draws each outcome from the dynamics of its current
phase (Gaussian, Poisson or binomial, random effects included; levels
the fit has not seen get no effect). The observed outcomes are not used,
so changing a covariate in `newdata` changes both who transitions and
the outcome in each phase.

For `type = "phase"`, the probability that each observation is in
`phase` given the unit's outcomes up to and including that observation
(the forward filter), which needs the outcome in `newdata`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_phaseR(model, data)
preds <- posterior_predict(fit)
} # }
```
