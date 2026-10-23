# Extract fixed effects

Returns posterior summaries of fixed effect coefficients.

## Usage

``` r
fixef(object, component = c("all", "transition", "dynamics"), phase = NULL)
```

## Arguments

- object:

  A `phase_fit` object

- component:

  Which component: "all" (default), "transition", or "dynamics"

- phase:

  For dynamics, which phase (NULL for all)

## Value

A data frame with fixed effect estimates
