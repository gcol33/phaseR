# Confidence intervals for phase_fit parameters

Returns credible intervals from the posterior distribution.

## Usage

``` r
# S3 method for class 'phase_fit'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `phase_fit` object

- parm:

  Character vector of parameter names (NULL for all)

- level:

  Credible interval level (default 0.95)

- ...:

  Additional arguments (unused)

## Value

Matrix with lower and upper bounds
