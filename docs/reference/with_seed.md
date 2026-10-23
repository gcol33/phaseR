# Evaluate an expression with a seeded random number generator

Seeds R's generator with `seed` for `expr` and restores the caller's
generator state afterwards.

## Usage

``` r
with_seed(seed, expr)
```

## Arguments

- seed:

  Integer seed

- expr:

  Expression to evaluate

## Value

The value of `expr`
