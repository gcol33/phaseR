# Story: phaseR

## Claim
Bayesian models of longitudinal change that treat "before" and "after" as latent regimes, fitting the transition and the within-regime dynamics jointly.

## Audience
Applied researchers (ecology, epidemiology, economics, psychology, medicine) who compare outcomes before and after a change that units select into. They know regression formulas and mixed-model syntax.

## Relation
- to: none
- home: none

## Arcs

### A1 Before/after splits are biased
- question: why not compare before and after the observed change?
- answer: conditioning on the realized transition mixes the regime effect with selection
- evidence: vignettes/estimand-failure.Rmd (known truth on simulated data)
- home: vignettes/estimand-failure.Rmd
- also: README, DESCRIPTION

### A2 One joint model
- question: how do I fit it?
- answer: `phase_model()` from `phases()`, `transition()`, `dynamics()`; `fit_phaseR()` with the HMC or Laplace backend
- home: README Usage; ?fit_phaseR

### A3 Counterfactual effects
- question: what would outcomes be under other covariate values?
- answer: `ate()`, `ite()`, `counterfactual()`
- home: README; ?ate

### A4 Model extensions
- question: does it handle my design?
- answer: k regimes, GLM outcomes, random effects, custom priors, WAIC and PSIS-LOO
- home: help pages
- also: README

## Documents
- DESCRIPTION: claim, A1
- README: claim, A1 to A4, first example, routing
- vignette estimand-failure: A1

## Open
- The HMC samplers in src/hmc_sampler*.cpp are named `phaseR_nuts_sampler*` but run fixed-length leapfrog trajectories ("simplified, not full NUTS"); the README says Hamiltonian Monte Carlo.
- The samplers seed from `std::random_device`, so `seed` does not make HMC fits reproducible.
- No pkgdown site; the vignette is linked by its repository path.
