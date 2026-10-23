# Cauchy prior distribution

Cauchy prior distribution

## Usage

``` r
cauchy(location = 0, scale = 1)
```

## Arguments

- location:

  Location parameter

- scale:

  Scale parameter (must be positive)

## Value

A prior specification object

## Examples

``` r
prior(beta ~ cauchy(0, 2.5))
```
