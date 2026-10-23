# Extract residuals from phase_fit

One-step prediction residuals: the observed outcome minus
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), divided for
`type = "pearson"` by the SD of the one-step predictive distribution
(phase mixture of the family variances, plus the posterior variance of
the predictive mean).

## Usage

``` r
# S3 method for class 'phase_fit'
residuals(object, type = c("response", "pearson"), ...)
```

## Arguments

- object:

  A `phase_fit` object

- type:

  Type of residuals: "response" (default) or "pearson"

- ...:

  Additional arguments (unused)

## Value

Numeric vector of residuals
