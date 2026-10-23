# Leave-One-Out Cross-Validation via PSIS

Computes approximate leave-one-out cross-validation using
Pareto-smoothed importance sampling (PSIS-LOO). This provides an
estimate of out-of-sample predictive accuracy that is more robust than
WAIC.

## Usage

``` r
loo(object, ...)
```

## Arguments

- object:

  A `phase_fit` object

- ...:

  Additional arguments (unused)

## Value

A `loo_phaseR` object with components:

- `elpd_loo`: Expected log pointwise predictive density for LOO

- `p_loo`: Effective number of parameters

- `looic`: LOO information criterion (-2 \* elpd_loo)

- `se_elpd_loo`: Standard error of elpd_loo

- `pointwise`: Data frame with per-observation values

- `diagnostics`: List with k-hat values and convergence info

## Details

PSIS-LOO approximates leave-one-out cross-validation without refitting
the model. It uses importance sampling with Pareto-smoothed weights to
estimate the predictive density for each held-out observation.

The k-hat diagnostic indicates the reliability of the estimate:

- k \< 0.5: Good, estimates are reliable

- 0.5 \<= k \< 0.7: Okay, estimates may have some error

- k \>= 0.7: Bad, estimates are unreliable for these observations

If many observations have k \> 0.7, consider using K-fold
cross-validation instead.

## References

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. Statistics and
Computing 27, 1413-1432.

## See also

[`waic()`](https://gillescolling.com/phaseR/reference/waic.md),
[`loo_compare()`](https://gillescolling.com/phaseR/reference/loo_compare.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_phaseR(model, data)
loo_result <- loo(fit)
print(loo_result)

# Check diagnostics
plot(loo_result$diagnostics$pareto_k)
} # }
```
