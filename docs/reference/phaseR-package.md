# phaseR: Phase-Structured Inference for Longitudinal Change

A framework for analysing change in longitudinal data that treats
"before" and "after" as regimes (phases), not time points. The key
insight is that **inference conditional on realized change is not a
valid estimand**.

## The Problem

Conventional analyses often:

1.  Observe when a unit "changes" (e.g., disease onset, treatment
    adoption)

2.  Compare outcomes before vs. after the change

3.  Attribute the difference to the "effect" of the change

This approach conflates the true phase effect with selection effects
from the transition mechanism.

## The Solution

phaseR models the complete generative process:

- **Transition model**: How units move between phases

- **Dynamics model**: How outcomes behave within each phase

By jointly modeling both components and integrating over phase
uncertainty, phaseR produces valid causal estimates.

## Core Functions

**Model specification:**

- [`phases`](https://gillescolling.com/phaseR/reference/phases.md):
  Declare ordered phases

- [`transition`](https://gillescolling.com/phaseR/reference/transition.md):
  Specify transition model

- [`dynamics`](https://gillescolling.com/phaseR/reference/dynamics.md):
  Specify phase-specific dynamics

- [`phase_model`](https://gillescolling.com/phaseR/reference/phase_model.md):
  Compose a complete model

**Fitting:**

- [`fit_phaseR`](https://gillescolling.com/phaseR/reference/fit_phaseR.md):
  Fit the model with the No-U-Turn sampler or a Laplace approximation

**Inference:**

- [`posterior_predict`](https://gillescolling.com/phaseR/reference/posterior_predict.md):
  Posterior predictive simulation

- [`counterfactual`](https://gillescolling.com/phaseR/reference/counterfactual.md):
  Counterfactual prediction

- [`ate`](https://gillescolling.com/phaseR/reference/ate.md): Average
  treatment effects

- [`ite`](https://gillescolling.com/phaseR/reference/ite.md): Individual
  treatment effects

**Simulation:**

- [`sim_phaseR`](https://gillescolling.com/phaseR/reference/sim_phaseR.md):
  Simulate 2-phase Gaussian data

- [`sim_phaseR_k`](https://gillescolling.com/phaseR/reference/sim_phaseR_k.md):
  Simulate k-phase data

- [`sim_phaseR_poisson`](https://gillescolling.com/phaseR/reference/sim_phaseR_poisson.md):
  Simulate Poisson data

## Design Principles

- Change (delta) is never modelled directly

- Phases are regimes, not time points

- Phase membership is latent or probabilistic

- Transitions are modelled, not conditioned on

## Author

**Maintainer**: Gilles Colling <gilles.colling051@gmail.com>
([ORCID](https://orcid.org/0000-0003-3070-6066))

Authors:

- Gilles Colling <gilles.colling051@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-3070-6066))
