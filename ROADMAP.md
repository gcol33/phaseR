# phaseR Roadmap to v1.1.0

## Current State (v1.0.0)

phaseR provides phase-structured inference for longitudinal data with:
- 2-phase and k-phase models
- GLM families (Gaussian, Poisson, Binomial)
- Single random intercept per component: `(1 | group)`
- Native HMC/NUTS sampler (no external dependencies)
- Counterfactual inference (ATE, ITE)
- WAIC for model comparison

---

## Feature Gap (vs ratiod)

### Random Effects

| Feature | ratiod | phaseR | Version |
|---------|--------|--------|---------|
| Single random intercept | ✓ | ✓ | v1.0.0 |
| Multiple crossed RE | ✓ | ✓ | v1.0.1 |
| Random slopes (correlated) | ✓ | ✓ | v1.0.2 |
| Random slopes (uncorrelated) | ✓ | ✓ | v1.0.2 |
| Nested RE | ✓ | ✓ | v1.0.3 |

### Inference Backends

| Feature | ratiod | phaseR | Version |
|---------|--------|--------|---------|
| HMC/NUTS | ✓ | ✓ | v1.0.0 |
| Laplace approximation | ✓ | ✓ | v1.0.4 |

### Priors & Model Comparison

| Feature | ratiod | phaseR | Version |
|---------|--------|--------|---------|
| Default weakly informative | ✓ | ✓ | v1.0.0 |
| Custom prior specification | ✓ | ✓ | v1.0.5 |
| LOO-CV | ✓ | ✗ | v1.0.6 |

### Future (v2.0.0+)

| Feature | ratiod | phaseR | Notes |
|---------|--------|--------|-------|
| Spatial effects (ICAR/CAR) | ✓ | ✗ | If needed |
| Temporal effects (AR1/RW1) | ✓ | ✗ | If needed |
| Zero-inflation | ✓ | ✗ | If needed |
| Polya-Gamma sampler | ✓ | ✗ | If needed |

---

## Implementation Plan

### v1.0.1: Multiple Crossed Random Effects

**Goal**: Support `y ~ x + (1|site) + (1|year)` syntax

**R code**:
- `R/formula.R`: Parse multiple `(1|group)` terms, return list
- `R/fit.R`: Build multiple RE index vectors, track n_groups per term
- `R/ranef.R`: Return named list by grouping factor

**C++ code**:
- `src/phase_likelihood.cpp`: Sum RE contributions across groups
- `src/hmc_sampler.cpp`: Variable parameter vector length

**Tests**:
- Crossed RE recovery with simulated data
- Comparison with single RE model

---

### v1.0.2: Random Slopes

**Goal**: Support `(1 + x | group)` and `(1 + x || group)`

**R code**:
- `R/formula.R`: Parse slope terms, detect `||` syntax
- `R/fit.R`: Build Z matrix (RE design matrix)

**C++ code**:
- Likelihood with Z matrix multiplication
- Cholesky factor for correlated slopes (LKJ prior)
- Diagonal covariance for uncorrelated slopes

**Output**:
- `VarCorr()` returns full covariance matrix
- Correlation estimates in summary

---

### v1.0.3: Nested Random Effects

**Goal**: Support `(1|site/plot)` syntax

**R code**:
- `R/formula.R`: Expand `(1|a/b)` to `(1|a) + (1|a:b)`
- Validation for proper nesting structure

**Implementation**: Reuses v1.0.1 crossed RE machinery

---

### v1.0.4: Laplace Approximation Backend

**Goal**: Fast approximate inference for large datasets

**R code**:
- `R/backend_laplace.R`: Full implementation
- `backend = "laplace"` option in `fit_phaseR()`

**C++ code**:
- Marginal likelihood (integrate out RE)
- Hessian computation at RE mode
- Sparse matrix support

**Output**:
- Consistent structure with HMC backend
- Standard errors via observed information

---

### v1.0.5: Custom Priors

**Goal**: User-specified priors on all parameters

**R code**:
- `R/prior.R`: Extend DSL
  ```r
  prior(
    beta ~ normal(0, 2.5),
    sigma ~ half_cauchy(0, 1),
    cor ~ lkj(2)
  )
  ```
- Prior classes: Normal, Student-t, Cauchy, Half-normal, Half-Cauchy, LKJ

**C++ code**:
- Pass prior parameters from R
- Dispatch to log-prior functions

---

### v1.0.6: LOO-CV

**Goal**: Leave-one-out cross-validation via PSIS

**R code**:
- `R/loo.R`: `loo()` and `loo_compare()` functions
- PSIS-LOO with diagnostic k-hat values

**Dependencies**: Builds on existing `compute_pointwise_ll()`

---

## Version Summary

| Version | Feature | Complexity |
|---------|---------|------------|
| v1.0.0 | Initial release | Done |
| v1.0.1 | Multiple crossed RE | Done |
| v1.0.2 | Random slopes | Done |
| v1.0.3 | Nested RE | Done |
| v1.0.4 | Laplace backend | Done |
| v1.0.5 | Custom priors | Done |
| v1.0.6 | LOO-CV | Medium |
| **v1.1.0** | **Feature parity with ratiod** | - |

---

## References

- Gelman et al. (2013). Bayesian Data Analysis, 3rd ed.
- Vehtari et al. (2017). Practical Bayesian model evaluation using LOO-CV and WAIC
- Stan Development Team. Stan User's Guide (random effects chapters)
