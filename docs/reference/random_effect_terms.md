# Random-effect terms of a fit with their posterior draws

Random-effect terms of a fit with their posterior draws

## Usage

``` r
random_effect_terms(object)
```

## Arguments

- object:

  A `phase_fit` object

## Value

List with one element per term, in the order of the parameter vector:
`predictor` (e.g. "transition", "dynamics_active"), `component`,
`phase`, `name` (grouping factor), `label`, `coefficients` (draws
`[draw, level, coefficient]`), `sd` (`[draw, coefficient]`) and `cor`
(`[draw, coefficient, coefficient]`)
