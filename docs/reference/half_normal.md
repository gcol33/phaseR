# Half-Normal prior distribution

A Normal(0, sd) distribution truncated to positive values. Suitable for
scale parameters like sigma.

## Usage

``` r
half_normal(sd = 1)
```

## Arguments

- sd:

  Standard deviation of the underlying normal (must be positive)

## Value

A prior specification object

## Examples

``` r
prior(sigma ~ half_normal(1))
```
