# Model description of a fit evaluated on new data

The parameter layout is the fit's; the outcome, the units and the design
matrices come from `newdata`. Random-effect levels are matched to the
fit's by name; a level the fit has not seen gets no effect.

## Usage

``` r
prediction_target(fit, newdata)
```

## Arguments

- fit:

  A `phase_fit` object

- newdata:

  Data frame with `id`, `time` and the model's covariates

## Value

List with `target` (as from
[`posterior_target()`](https://gillescolling.com/phaseR/reference/posterior_target.md)),
`order` (the rows of `newdata` in the order of the target) and
`has_response`
