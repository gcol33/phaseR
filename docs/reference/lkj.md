# LKJ prior for correlation matrices

The LKJ(eta) prior places more mass on the identity matrix as eta
increases. eta = 1 is uniform over correlation matrices.

## Usage

``` r
lkj(eta = 2)
```

## Arguments

- eta:

  Shape parameter (must be positive). eta = 1 is uniform, larger values
  favor weaker correlations.

## Value

A prior specification object

## Examples

``` r
prior(cor ~ lkj(2))
```
