# Compare models using LOO-CV

Compare multiple fitted phase models using leave-one-out
cross-validation.

## Usage

``` r
loo_compare(..., model_names = NULL)
```

## Arguments

- ...:

  `phase_fit` objects to compare, or `loo_phaseR` objects

- model_names:

  Optional names for the models

## Value

A `loo_compare` object (data frame with comparison)

## Details

Models are ranked by expected log pointwise predictive density
(elpd_loo), with higher values indicating better predictive performance.
The difference in elpd and its standard error are computed relative to
the best model.

## Examples

``` r
if (FALSE) { # \dontrun{
fit1 <- fit_phaseR(model1, data)
fit2 <- fit_phaseR(model2, data)
loo_compare(fit1, fit2, model_names = c("Simple", "Complex"))
} # }
```
