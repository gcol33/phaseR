# phaseR

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Phase-structured inference for longitudinal change.

## Installation

```r
# From GitHub
devtools::install_github("gcol33/phaseR")
```

## The Problem

A common analytical pattern is to compare outcomes "before" vs "after" a change
(disease onset, treatment adoption, regime shift). This approach is **fundamentally
flawed** when the probability of transitioning depends on covariates that also
affect outcomes.

**Inference conditional on realized change is not a valid estimand.**

## The Solution

phaseR models the complete generative process:

1. **Transition model**: How units move between phases
2. **Dynamics model**: How outcomes behave within each phase

By jointly modeling both components, phaseR produces valid causal estimates.

## Quick Start

```r
library(phaseR)

# Simulate data
dat <- sim_phaseR(n_units = 100, n_times = 5, seed = 42)

# Define model
model <- phase_model(
  phases("inactive", "active"),
  transition("inactive", "active", ~ x),
  dynamics("inactive", y ~ x),
  dynamics("active", y ~ x)
)

# Fit
fit <- fit_phaseR(model, data = dat, chains = 2, iter = 1000)

# Summarize
summary(fit)

# Counterfactual: effect of setting x = 1 vs x = 0
cf_high <- counterfactual(fit, newdata = transform(dat, x = 1))
cf_low <- counterfactual(fit, newdata = transform(dat, x = 0))
```

## Features

- **k-phase models**: Support for arbitrary number of sequential phases
- **GLM families**: Gaussian, Poisson, Binomial outcomes
- **Random effects**: Hierarchical models with `(1 | group)` syntax
- **Counterfactual inference**: `ate()`, `ite()`, `counterfactual()`
- **Native HMC**: No external dependencies (Stan, JAGS, etc.)

## Model Specification

```r
# Phases (must be sequential)
phases("susceptible", "infected", "recovered")

# Transitions (each non-absorbing phase -> next)
transition("susceptible", "infected", ~ exposure + age)
transition("infected", "recovered", ~ treatment)

# Dynamics (outcome model per phase)
dynamics("susceptible", y ~ 1)
dynamics("infected", y ~ time_since_infection)
dynamics("recovered", y ~ 1)

# Compose
model <- phase_model(
  phases(...),
  transition(...),
  transition(...),
  dynamics(...),
  dynamics(...),
  dynamics(...)
)
```

## Design Principles

1. **Change (delta) is never modelled directly** - we model outcomes, not differences
2. **Phases are regimes, not time points** - before/after are states, not calendar time
3. **Phase membership is latent** - we don't condition on observed transitions
4. **Transitions are modelled** - the selection mechanism is part of the model

## Applications

phaseR is domain-agnostic. Example applications:

- **Ecology**: Species invasions, ecosystem regime shifts
- **Epidemiology**: Disease onset, recovery, relapse
- **Economics**: Technology adoption, market entry
- **Psychology**: Behavioral change, skill acquisition
- **Medicine**: Treatment response, disease progression

## Citation

If you use phaseR in your research, please cite:

```
@software{phaseR,
  title = {phaseR: Phase-Structured Inference for Longitudinal Change},
  author = {Colling, Gilles},
  year = {2025},
  url = {https://github.com/gcol33/phaseR}
}
```

## License

MIT
