# Compare models using WAIC

Compare multiple fitted phase models using WAIC.

## Usage

``` r
compare_waic(..., model_names = NULL)
```

## Arguments

- ...:

  `phase_fit` objects to compare

- model_names:

  Optional names for the models

## Value

A data frame with WAIC comparisons

## Examples

``` r
if (FALSE) { # \dontrun{
compare_waic(fit1, fit2, model_names = c("Simple", "Complex"))
} # }
```
