# Evaluate the log posterior of a model description

Evaluate the log posterior of a model description

## Usage

``` r
log_posterior(target, params, gradient = FALSE)
```

## Arguments

- target:

  Model description from
  [`posterior_target()`](https://gillescolling.com/phaseR/reference/posterior_target.md)

- params:

  Parameter vector

- gradient:

  If `TRUE`, also return the exact gradient

## Value

The log posterior at `params`; with `gradient = TRUE` a list with
`value` and `gradient`
