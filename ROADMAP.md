# phaseR Development Roadmap

## Vision

**phaseR** is a phase-structured framework for analysing change in longitudinal data that treats "before" and "after" as regimes, not time points. It separates regime activation from regime-specific dynamics, preventing outcome-dependent stratification and enabling valid counterfactual inference.

**Core methodological claim**: Inference conditional on realised change is not an estimand.

**Application domains**: Ecology (invasion, colonization, recovery), epidemiology (disease onset, treatment initiation), economics (regime switches, policy adoption), psychology (behavioral state changes), medicine (disease progression), and any field where longitudinal data involves latent phase transitions.

## Version Overview

### Pre-Release (v0.X.0)

| Version | Title | Description |
|---------|-------|-------------|
| [v0.1.0](docs/roadmap/v0.1.0-skeleton.md) | Package Skeleton | Core S3 classes, directory structure, CI setup |
| [v0.2.0](docs/roadmap/v0.2.0-two-phase.md) | Two-Phase Model | Minimal working model: inactive → active |
| [v0.3.0](docs/roadmap/v0.3.0-random-intercepts.md) | Random Intercepts | Site-level RE in transitions and dynamics |
| [v0.4.0](docs/roadmap/v0.4.0-k-phase.md) | k-Phase Support | Ordered phases, absorbing states, explicit skips |
| [v0.5.0](docs/roadmap/v0.5.0-counterfactual.md) | Counterfactual | "Never initiated" counterfactual inference |
| [v0.6.0](docs/roadmap/v0.6.0-families.md) | GLM Families | Gaussian, Poisson, Binomial outcomes |
| [v0.7.0](docs/roadmap/v0.7.0-estimand-vignette.md) | Estimand Failure Demo | Null simulation vignette |
| [v0.8.0](docs/roadmap/v0.8.0-documentation.md) | Documentation | Full docs, pkgdown site, vignettes |
| [v0.9.0](docs/roadmap/v0.9.0-cran-prep.md) | CRAN Preparation | R CMD check clean, final polish |
| [v1.0.0](docs/roadmap/v1.0.0-release.md) | CRAN Release | First stable release |

### Post-Release (v1.X.0)

| Version | Title | Description |
|---------|-------|-------------|
| [v1.1.0](docs/roadmap/v1.1.0-diagnostics.md) | Model Diagnostics | LOO, WAIC, posterior predictive checks |
| [v1.2.0](docs/roadmap/v1.2.0-random-slopes.md) | Random Slopes | Correlated random effects |
| [v1.2.1](docs/roadmap/v1.2.1-ratiod-parity.md) | ratiod Parity | Nested RE, crossed RE, interactions, offsets |
| [v1.3.0](docs/roadmap/v1.3.0-extended-families.md) | Extended Families | Negative binomial, zero-inflated, beta |
| [v1.4.0](docs/roadmap/v1.4.0-time-varying-hazards.md) | Time-Varying Hazards | Duration dependence, smooth time effects |
| [v1.5.0](docs/roadmap/v1.5.0-spatial-temporal.md) | Spatial & Temporal | AR(1), CAR, BYM2, Gaussian process |
| [v1.6.0](docs/roadmap/v1.6.0-multivariate.md) | Multivariate Outcomes | Multiple correlated responses |
| [v1.7.0](docs/roadmap/v1.7.0-timing-counterfactuals.md) | Timing Counterfactuals | "Initiate at time t" interventions |
| [v1.8.0](docs/roadmap/v1.8.0-path-counterfactuals.md) | Path Counterfactuals | Fixed phase path interventions |
| [v1.9.0](docs/roadmap/v1.9.0-performance.md) | Performance | Variational inference, parallelization |

## Dependency Graph

```
v0.1.0 (skeleton)
    │
    v
v0.2.0 (2-phase) ──────────────────────────────────────┐
    │                                                   │
    v                                                   │
v0.3.0 (random intercepts)                              │
    │                                                   │
    v                                                   │
v0.4.0 (k-phase)                                        │
    │                                                   │
    v                                                   │
v0.5.0 (counterfactual)                                 │
    │                                                   │
    v                                                   │
v0.6.0 (families)                                       │
    │                                                   │
    v                                                   │
v0.7.0 (vignette) ◄─────────────────────────────────────┘
    │                     (uses all prior features)
    v
v0.8.0 (documentation)
    │
    v
v0.9.0 (CRAN prep)
    │
    v
v1.0.0 (CRAN RELEASE)
    │
    ├── v1.1.0 (diagnostics)
    │
    ├── v1.2.0 (random slopes)
    │       │
    │       └── v1.5.0 (spatial/temporal)
    │
    ├── v1.3.0 (extended families)
    │       │
    │       └── v1.6.0 (multivariate)
    │
    ├── v1.4.0 (time-varying hazards)
    │       │
    │       └── v1.7.0 (timing counterfactuals)
    │               │
    │               └── v1.8.0 (path counterfactuals)
    │
    └── v1.9.0 (performance)

v2.0.0 (reserved for breaking changes)
```

## Design Principles

### Non-Negotiable Constraints

1. **Change (Δ) is never modelled directly** — Δ is always a derived quantity
2. **"Before" and "after" are phases (regimes), not time points**
3. **Phases are declared a priori** — never defined from realised outcomes
4. **Movement between phases is modelled via explicit transition processes**
5. **Phase membership is latent or probabilistic**
6. **Group comparisons based on realised transitions are structurally impossible**

### API Guardrails

| Invalid Operation | How Blocked |
|-------------------|-------------|
| Δ as response | Formula parser rejects difference terms |
| Phases from outcome | Phase membership always latent |
| Condition on realised transitions | No function returns transition groups |
| Post-hoc transition slopes | Extractors return phase-conditional only |

### Architecture Patterns (following ratiod)

- **S3 classes** for all objects (not S4)
- **Native C++/Rcpp** with custom HMC sampler (not TMB)
- **Draws matrix** as primary output format
- **Formula-as-config** pattern with structured objects
- **Backend-agnostic API** with late dispatch
- **OpenMP** for parallelization

## File Organization

```
phaseR/
├── R/
│   ├── phaseR-package.R      # Package docs, useDynLib
│   ├── phases.R              # Phase declaration
│   ├── transition.R          # Transition specification
│   ├── dynamics.R            # Phase-specific dynamics
│   ├── formula.R             # Formula parsing
│   ├── fit.R                 # Main fitting function
│   ├── backend_hmc.R         # HMC backend
│   ├── backend_laplace.R     # Laplace approximation
│   ├── methods.R             # print, summary, plot
│   ├── predict.R             # Predictions
│   ├── counterfactual.R      # Counterfactual inference
│   └── validate.R            # Input validation
├── src/
│   ├── hmc_core.cpp          # HMC gradient & trajectory
│   ├── hmc_sampler.cpp       # NUTS sampler
│   ├── laplace_core.cpp      # Laplace approximation
│   ├── phase_likelihood.cpp  # Phase-switching likelihood
│   └── Makevars              # OpenMP config
├── tests/testthat/
├── vignettes/
├── docs/roadmap/             # This roadmap
├── man/
├── DESCRIPTION
├── NAMESPACE
└── _pkgdown.yml
```

## Explicitly Not Planned

| Feature | Reason |
|---------|--------|
| Arbitrary transition graphs | Violates ordered-phase discipline |
| Outcome-defined phases | Core methodological prohibition |
| Δ as response | Core methodological prohibition |
| Bayesian model averaging over phase counts | Complicates interpretation |
| Real-time/streaming inference | Different use case |

## Contributing

Each version milestone has a detailed specification in `docs/roadmap/`. Before starting work on a version:

1. Read the version spec completely
2. Ensure all prerequisites are met
3. Follow the validation criteria
4. Update the spec if implementation diverges

## References

- ratiod package: architectural patterns
- Methodology paper: theoretical grounding (in progress)
