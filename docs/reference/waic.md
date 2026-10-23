# Compute WAIC for model comparison

Computes the Widely Applicable Information Criterion (WAIC) for a fitted
phase model. This provides an estimate of out-of-sample predictive
accuracy.

## Usage

``` r
waic(object, ...)
```

## Arguments

- object:

  A `phase_fit` object

- ...:

  Additional arguments (unused)

## Value

A list with components:

- `waic`: The WAIC value (lower is better)

- `lppd`: Log pointwise predictive density

- `p_waic`: Effective number of parameters

- `se`: Standard error of WAIC

## Details

WAIC is computed as: \$\$WAIC = -2 \* (lppd - p\_{WAIC})\$\$

where lppd is the log pointwise predictive density and p_WAIC is the
effective number of parameters (variance of the log-likelihood).

For phase models, the log-likelihood of each observation is that given
the earlier observations of its unit, marginalised over the latent
phases, at every posterior draw.

## References

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
Journal of Machine Learning Research 11, 3571-3594.

## Examples

``` r
if (FALSE) { # \dontrun{
fit1 <- fit_phaseR(model1, data)
fit2 <- fit_phaseR(model2, data)
waic(fit1)$waic
waic(fit2)$waic
} # }
```
