# Generalized inverse (Moore-Penrose pseudo-inverse)

Simple implementation without MASS dependency.

## Usage

``` r
MASS_ginv(X, tol = sqrt(.Machine$double.eps))
```

## Arguments

- X:

  Matrix to invert

- tol:

  Tolerance for singular values

## Value

Pseudo-inverse of X
