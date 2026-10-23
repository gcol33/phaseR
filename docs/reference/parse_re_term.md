# Parse a single RE term string

Parse a single RE term string

## Usage

``` r
parse_re_term(term_str)
```

## Arguments

- term_str:

  String like "(1 + x \| group)" or "(1 + x \|\| group)"

## Value

List with group, type, vars, correlated
