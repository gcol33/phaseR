# Extract fitted values from phase_fit

The posterior mean of the one-step predictive mean of each observation,
E(y_t \| y_0, ..., y\_(t-1)): the outcome mean of each phase weighted by
the probability of the phase given the unit's earlier outcomes.

## Usage

``` r
# S3 method for class 'phase_fit'
fitted(object, ...)
```

## Arguments

- object:

  A `phase_fit` object

- ...:

  Additional arguments (unused)

## Value

Numeric vector of fitted values, one per row of the data
