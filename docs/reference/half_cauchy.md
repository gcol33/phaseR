# Half-Cauchy prior distribution

A Cauchy(0, scale) distribution truncated to positive values. Has
heavier tails than half-normal.

## Usage

``` r
half_cauchy(scale = 1)
```

## Arguments

- scale:

  Scale parameter (must be positive)

## Value

A prior specification object

## Examples

``` r
prior(sigma ~ half_cauchy(1))
```
