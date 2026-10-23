# Parse formula for random effects

Parses formulas with random intercept and/or slope terms. Supports:

- `y ~ x + (1|site)` - random intercept

- `y ~ x + (1 + x|site)` - correlated random intercept and slope

- `y ~ x + (1 + x||site)` - uncorrelated random intercept and slope

- `y ~ x + (0 + x|site)` or `(x|site)` - random slope only

- `y ~ x + (1|site/plot)` - nested random effects (expands to
  (1\|site) + (1\|site:plot))

## Usage

``` r
parse_formula_re(formula)
```

## Arguments

- formula:

  A formula potentially containing random effect terms

## Value

A list with:

- `fixed_formula`: Formula without RE terms

- `has_re`: Logical, TRUE if any RE terms found

- `re_terms`: List of RE specifications (one per term)

- `n_re_terms`: Number of RE terms

- `re_groups`: Character vector of group names (for backward compat)

- `has_slopes`: Logical, TRUE if any RE term has slopes

- `has_nested`: Logical, TRUE if nested RE terms were expanded

- `nested_interactions`: Character vector of interaction terms to create
