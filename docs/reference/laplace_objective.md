# Negative Laplace log marginal posterior for optimization

Negative Laplace log marginal posterior for optimization

## Usage

``` r
laplace_objective(laplace)
```

## Arguments

- laplace:

  External pointer from `phaseR_laplace()`

## Value

List of two functions of the parameters outside the random effects:
`fn`, the negative log marginal posterior with the random effects
integrated out by the Laplace approximation, and `gr`, its exact
gradient. The C++ object caches the last evaluation, so `fn` and `gr` at
the same point cost one evaluation.
