# Student-t prior distribution

Student-t prior distribution

## Usage

``` r
student_t(df = 3, mean = 0, scale = 2.5)
```

## Arguments

- df:

  Degrees of freedom (must be positive)

- mean:

  Location parameter

- scale:

  Scale parameter (must be positive)

## Value

A prior specification object

## Examples

``` r
# Robust prior for coefficients
prior(beta ~ student_t(3, 0, 2.5))
```
