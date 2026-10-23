# Predict from a phase model

Predictions given the covariates of `newdata` alone, averaged over the
posterior draws: each unit starts in the first phase and moves through
the phases with the transition probabilities its covariates imply, and
each phase's outcome follows its dynamics, random effects included.
Random-effect levels that the fit has not seen contribute no effect.

## Usage

``` r
# S3 method for class 'phase_fit'
predict(
  object,
  newdata = NULL,
  type = c("response", "phase", "transition"),
  summary = TRUE,
  ...
)
```

## Arguments

- object:

  A `phase_fit` object

- newdata:

  Optional new data for prediction, with columns `id`, `time` and the
  covariates of the model; defaults to the data of the fit

- type:

  `"response"` (default): the expected outcome, marginalised over the
  phase; `"phase"`: the probability of each phase; `"transition"`: the
  probability of each move between consecutive phases at each
  observation (zero at a unit's first observation)

- summary:

  Return summary statistics (default `TRUE`) or the draws (one row per
  draw, one column per row of `newdata`)

- ...:

  Additional arguments (unused)

## Value

For `"response"`, a data frame of summaries (or the matrix of draws).
For `"phase"`, a data frame with one column per phase holding the
posterior mean probability. For `"transition"`, the summaries (or draws)
of the single move of a two-phase model, or a named list of them with
one element per move of a model with more phases.
