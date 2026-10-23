# Extract residuals from phase_fit

Computes residuals (observed - fitted) marginalized over phase
uncertainty.

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
