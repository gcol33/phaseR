# Expand nested RE syntax

Expands `(1|a/b)` to `(1|a) + (1|a:b)` and `(1|a/b/c)` to
`(1|a) + (1|a:b) + (1|a:b:c)`.

## Usage

``` r
expand_nested_re(formula_str)
```

## Arguments

- formula_str:

  Formula string

## Value

List with expanded formula_str, has_nested flag, and interactions
