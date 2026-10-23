# Phase summaries of a fit on new data

Evaluates the shared C++ model at every draw on the rows of `newdata`.

## Usage

``` r
phase_summaries(fit, newdata, draws)
```

## Arguments

- fit:

  A `phase_fit` object

- newdata:

  Data frame with `id`, `time` and the model's covariates; the outcome
  is needed only for the outcome-conditional summaries

- draws:

  Matrix of parameter draws (one row per draw)

## Value

List of arrays indexed `[draw, row of newdata, phase]`: `mean` (the
outcome mean of each phase, on the response scale), `marginal` (the
phase probability given the covariates), `predictive` (given the unit's
earlier outcomes), `filtered` (given the unit's outcomes up to and
including the row) and `move` (the probability of moving from each phase
to the next at the row), plus `sigma` (`[draw, phase]`, the residual
SDs), `family` (the family code of each phase) and `has_response`
