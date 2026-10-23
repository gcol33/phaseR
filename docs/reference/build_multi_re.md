# Build multiple random effect structures

Build multiple random effect structures

## Usage

``` r
build_multi_re(data, re_terms)
```

## Arguments

- data:

  Data frame

- re_terms:

  List of RE term specifications from parse_formula_re()

## Value

List with:

- `n_re_terms`: Number of RE terms

- `re_list`: List of RE structures (one per term)

- `total_re_params`: Total number of RE parameters

- `has_slopes`: Whether any term has slopes
