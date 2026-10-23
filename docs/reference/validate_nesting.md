# Validate nested structure

Checks that nested groups are properly nested (each level of b is within
one level of a).

## Usage

``` r
validate_nesting(data, outer, inner)
```

## Arguments

- data:

  Data frame

- outer:

  Outer grouping variable name

- inner:

  Inner grouping variable name

## Value

TRUE if valid, otherwise throws an error
