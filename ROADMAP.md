# phaseR Development Roadmap

## Current State (v1.0.0)

phaseR provides phase-structured inference for longitudinal data with:
- 2-phase and k-phase models
- GLM families (Gaussian, Poisson, Binomial)
- Single random intercept per component: `(1 | group)`
- Native HMC/NUTS sampler (no external dependencies)
- Counterfactual inference (ATE, ITE)
- WAIC for model comparison

---

## Feature Gap Analysis (vs ratiod)

### Random Effects

| Feature | ratiod | phaseR | Priority |
|---------|--------|--------|----------|
| Single random intercept | ✓ | ✓ | Done |
| Multiple crossed RE | ✓ | ✗ | High |
| Random slopes (correlated) | ✓ | ✗ | High |
| Random slopes (uncorrelated) | ✓ | ✗ | Medium |
| Nested RE | ✓ | ✗ | Medium |

### Inference Backends

| Feature | ratiod | phaseR | Priority |
|---------|--------|--------|----------|
| HMC/NUTS | ✓ | ✓ | Done |
| Laplace approximation | ✓ | ✗ | High |
| Polya-Gamma (binomial) | ✓ | ✗ | Low |

### Priors

| Feature | ratiod | phaseR | Priority |
|---------|--------|--------|----------|
| Default weakly informative | ✓ | ✓ | Done |
| Custom prior specification | ✓ | ✗ | Medium |
| Prior predictive checks | ✓ | ✗ | Low |

### Model Extensions

| Feature | ratiod | phaseR | Priority |
|---------|--------|--------|----------|
| Spatial effects (ICAR/CAR) | ✓ | ✗ | Low |
| Temporal effects (AR1/RW1) | ✓ | ✗ | Low |
| Zero-inflation | ✓ | ✗ | Low |

### Diagnostics & Output

| Feature | ratiod | phaseR | Priority |
|---------|--------|--------|----------|
| Rhat, ESS | ✓ | ✓ | Done |
| Trace plots | ✓ | ✓ | Done |
| LOO-CV | ✓ | ✗ | Medium |
| Posterior predictive checks | ✓ | partial | Medium |

---

## Implementation Plan

### Phase 1: Multiple Crossed Random Effects (v1.1.0)

**Goal**: Support `y ~ x + (1|site) + (1|year)` syntax

**Tasks**:

1. **Update formula parser** (`R/formula.R`)
   - Parse multiple `(1|group)` terms
   - Return list of RE specifications
   - Handle different RE terms for transition vs dynamics

2. **Update data preparation** (`R/fit.R`)
   - Build multiple RE index vectors
   - Track n_groups for each RE term
   - Update parameter naming convention

3. **Update C++ likelihood** (`src/phase_likelihood.cpp`)
   - Accept variable number of RE terms
   - Sum contributions from each RE group
   - Update gradient computation

4. **Update C++ sampler** (`src/hmc_sampler.cpp`)
   - Handle variable parameter vector length
   - Proper indexing for multiple RE blocks

5. **Update extractors** (`R/ranef.R`)
   - Return named list of RE estimates by grouping factor

**Estimated complexity**: Medium-High

---

### Phase 2: Random Slopes (v1.2.0)

**Goal**: Support `y ~ x + (1 + x | group)` and `y ~ x + (1 + x || group)`

**Tasks**:

1. **Extend formula parser**
   - Detect slope terms in RE specification
   - Parse `||` for uncorrelated slopes
   - Extract slope variable names

2. **Design matrix for RE**
   - Build Z matrix (RE design matrix) per group
   - Handle intercept + slope columns

3. **Covariance structure**
   - For correlated: estimate Cholesky factor of 2x2 (or larger) covariance
   - For uncorrelated: diagonal covariance (separate sigmas)

4. **Update C++ code**
   - Likelihood with Z matrix multiplication
   - Prior on Cholesky factors (LKJ or similar)
   - Gradient computation for correlation parameters

5. **Update output**
   - Report correlation estimates
   - `VarCorr()` returns full covariance matrix

**Estimated complexity**: High

---

### Phase 3: Laplace Approximation Backend (v1.3.0)

**Goal**: Fast approximate inference for large datasets

**Tasks**:

1. **Implement marginal likelihood**
   - Integrate out RE using Laplace approximation
   - Compute Hessian at RE mode

2. **Optimization**
   - Inner optimization: find RE mode given fixed effects
   - Outer optimization: maximize marginal likelihood

3. **Standard errors**
   - Compute observed information matrix
   - Delta method for transformed parameters

4. **C++ implementation**
   - Efficient Hessian computation
   - Sparse matrix support for large RE

5. **R interface**
   - `backend = "laplace"` option
   - Consistent output structure with HMC

**Estimated complexity**: High

---

### Phase 4: Nested Random Effects (v1.4.0)

**Goal**: Support `y ~ x + (1|site/plot)` syntax

**Tasks**:

1. **Parse nested syntax**
   - Expand `(1|site/plot)` to `(1|site) + (1|site:plot)`
   - Create interaction grouping factor

2. **Validation**
   - Check nesting structure is valid
   - Warn if groups are not properly nested

3. **Update internals**
   - Reuse crossed RE machinery
   - Proper labeling of nested levels

**Estimated complexity**: Low (builds on Phase 1)

---

### Phase 5: Custom Priors (v1.5.0)

**Goal**: User-specified priors on all parameters

**Tasks**:

1. **Prior specification DSL**
   ```r
   prior(
     beta ~ normal(0, 2.5),
     sigma ~ half_cauchy(0, 1),
     sigma_re ~ half_normal(0, 0.5)
   )
   ```

2. **Prior classes**
   - Normal, Student-t, Cauchy
   - Half-normal, Half-Cauchy, Exponential
   - LKJ for correlations

3. **Update C++ prior computation**
   - Pass prior parameters from R
   - Dispatch to appropriate log-prior function

4. **Prior predictive simulation**
   - Sample from prior
   - Generate fake data for model checking

**Estimated complexity**: Medium

---

### Phase 6: LOO-CV (v1.6.0)

**Goal**: Leave-one-out cross-validation via PSIS

**Tasks**:

1. **Pointwise log-likelihood storage**
   - Already have `compute_pointwise_ll()`
   - Ensure it works for all model types

2. **PSIS-LOO implementation**
   - Pareto smoothed importance sampling
   - Diagnostic k-hat values

3. **Model comparison**
   - `loo()` function
   - `loo_compare()` for multiple models

**Estimated complexity**: Medium

---

### Future Considerations (v2.0.0+)

These features are lower priority but may be added based on user demand:

- **Spatial random effects** (ICAR/CAR/BYM2)
  - Requires adjacency matrix input
  - Sparse precision matrix computation
  - Most useful for areal phase data

- **Temporal autocorrelation** (AR1/RW1)
  - Within-phase temporal structure
  - Requires careful handling with phase transitions

- **Zero-inflated models**
  - Mixture of point mass at zero and count distribution
  - Additional parameters for zero-inflation probability

- **Negative binomial family**
  - Overdispersed counts
  - Estimate dispersion parameter

- **Multivariate outcomes**
  - Multiple response variables
  - Shared latent phases across outcomes

---

## Version Timeline

| Version | Features | Target |
|---------|----------|--------|
| v1.1.0 | Multiple crossed RE | TBD |
| v1.2.0 | Random slopes | TBD |
| v1.3.0 | Laplace backend | TBD |
| v1.4.0 | Nested RE | TBD |
| v1.5.0 | Custom priors | TBD |
| v1.6.0 | LOO-CV | TBD |
| v2.0.0 | Spatial/temporal (if needed) | TBD |

---

## Contributing

Contributions welcome. Priority areas:
1. C++ optimization for HMC
2. Additional GLM families
3. Vignettes and documentation
4. Test coverage for edge cases

## References

- Gelman et al. (2013). Bayesian Data Analysis, 3rd ed.
- Vehtari et al. (2017). Practical Bayesian model evaluation using LOO-CV and WAIC
- Stan Development Team. Stan User's Guide (random effects chapters)
