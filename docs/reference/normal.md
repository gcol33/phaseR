# Normal prior distribution

Normal prior distribution

## Usage

``` r
normal(mean = 0, sd = 1)
```

## Arguments

- mean:

  Mean of the distribution

- sd:

  Standard deviation (must be positive)

## Value

A prior specification object

## Examples

``` r
prior(beta ~ normal(0, 2.5))
```
