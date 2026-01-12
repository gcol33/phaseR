# phaseR 0.8.0

Experimental release with feature parity to ratiod.

## Features

* **Phase-structured models** for longitudinal data
  - k-phase support (arbitrary number of sequential phases)
  - GLM families: Gaussian, Poisson, Binomial
  - Counterfactual inference: `ate()`, `ite()`, `counterfactual()`

* **Random effects**
  - Random intercepts: `(1 | group)`
  - Multiple crossed RE: `(1|site) + (1|year)`
  - Random slopes: `(1 + x | group)` and `(1 + x || group)`
  - Nested RE: `(1|site/plot)`

* **Inference backends**
  - Native HMC/NUTS sampler (no external dependencies)
  - Laplace approximation: `backend = "laplace"`

* **Custom priors**
  - Formula syntax: `prior(beta ~ normal(0, 2.5), sigma ~ half_cauchy(1))`
  - Distributions: `normal()`, `student_t()`, `cauchy()`, `half_normal()`, `half_cauchy()`, `lkj()`

* **Model comparison**
  - WAIC: `waic()`, `compare_waic()`
  - LOO-CV via PSIS: `loo()`, `loo_compare()`

## Core Functions

* `phases()`: Declare ordered phases
* `transition()`: Specify transition models
* `dynamics()`: Specify phase-specific outcome models
* `phase_model()`: Compose complete model
* `fit_phaseR()`: Fit model
* `posterior_predict()`: Posterior predictive simulation
* `fixef()`, `ranef()`, `VarCorr()`: Extract estimates
* `residuals()`, `confint()`: Diagnostics

---

# phaseR 0.5.0

Initial release.
