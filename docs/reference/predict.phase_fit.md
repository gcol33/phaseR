# Predict from a phase model

Predict from a phase model

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

  Optional new data for prediction

- type:

  Type of prediction: "response" (default), "phase", or "transition"

- summary:

  Return summary statistics (default TRUE) or full draws

- ...:

  Additional arguments (unused)

## Value

A data frame with predictions and credible intervals
