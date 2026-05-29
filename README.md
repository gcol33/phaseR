# phaseR

*the before/after split keeps lying to you*

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/gcol33/phaseR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/phaseR/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Joint Bayesian inference for longitudinal change, fit by a native NUTS sampler in C++.**

`phaseR` treats "before" and "after" as latent regimes and models two things at once: the transition that moves a unit between regimes, and the outcome dynamics within each regime. The sampler marginalizes over which regime each observation belongs to, so estimates never condition on the realized transition. Splitting your data at the observed change point and comparing the two halves does condition on it, which biases the comparison whenever the chance of transitioning depends on the same covariates that drive the outcome.

```r
library(phaseR)

dat <- sim_phaseR(n_units = 100, n_times = 5, seed = 42)

model <- phase_model(
  phases("inactive", "active"),
  transition("inactive", "active", ~ x),   # what moves a unit between regimes
  dynamics("inactive", y ~ x),             # outcome model in each regime
  dynamics("active",   y ~ x)
)

fit <- fit_phaseR(model, data = dat, chains = 2, iter = 1000, seed = 42)
summary(fit)
```

## Why before/after comparison is biased

Compare the outcomes of units that transitioned, before versus after, and you have stratified on a realized outcome. When units with certain covariate values are both more likely to transition and have different outcomes, the difference you measure mixes the regime effect with that selection. The bias survives any amount of data.

`phaseR` writes down the full generative process instead: a logistic transition model for regime membership and a per-regime outcome model, fit jointly. The likelihood integrates over the latent regime of every observation (a log-sum-exp over phase assignments), so the selection mechanism is part of the model rather than a thing you condition away. The [`estimand-failure`](vignettes/estimand-failure.Rmd) vignette runs the biased before/after estimate and the `phaseR` estimate against a known truth on simulated data.

## Counterfactuals once you have a fit

Because the transition and dynamics are modeled together, you can ask what outcomes would have looked like under covariate values you never observed, integrating over phase uncertainty:

```r
ate(fit, treatment_var = "x")                          # average treatment effect, with a credible interval

counterfactual(fit, newdata = transform(dat, x = 1))   # predicted outcomes under x = 1
ite(fit, treatment_var = "x")                          # per-unit effects
```

## What the sampler handles

- **k regimes**: arbitrary sequential phases, declared with `phases("susceptible", "infected", "recovered")` and one `transition()` per consecutive pair.
- **GLM outcomes**: Gaussian, Poisson, and Binomial dynamics, set per phase.
- **Random effects**: random intercepts `(1 | group)`, crossed `(1|site) + (1|year)`, random slopes `(1 + x | group)`, and nested `(1|site/plot)`.
- **Custom priors**: `prior(beta ~ normal(0, 2.5), sigma ~ half_cauchy(1))`, with `normal()`, `student_t()`, `cauchy()`, `half_normal()`, `half_cauchy()`, and `lkj()`.
- **Model comparison**: `waic()` / `compare_waic()`, and LOO-CV via PSIS with `loo()` / `loo_compare()`.

The NUTS sampler is written from scratch in C++ via Rcpp. There is no Stan, JAGS, or other external inference engine to install.

## Specifying a model

```r
phases("susceptible", "infected", "recovered")        # ordered regimes

transition("susceptible", "infected", ~ exposure + age)
transition("infected", "recovered", ~ treatment)      # one per consecutive pair

dynamics("susceptible", y ~ 1)
dynamics("infected",    y ~ time_since_infection)      # outcome model per regime
dynamics("recovered",   y ~ 1)

model <- phase_model(
  phases("susceptible", "infected", "recovered"),
  transition("susceptible", "infected", ~ exposure + age),
  transition("infected", "recovered", ~ treatment),
  dynamics("susceptible", y ~ 1),
  dynamics("infected", y ~ time_since_infection),
  dynamics("recovered", y ~ 1)
)
```

## Where it applies

Any longitudinal setting with a latent before/after where the chance of crossing over depends on covariates that also shape the outcome: disease onset and recovery, treatment adoption, technology uptake, behavioral change, ecological regime shifts.

## Installation

```r
# install.packages("pak")
pak::pak("gcol33/phaseR")
```

## Documentation

- [Why conditioning on realized change fails](vignettes/estimand-failure.Rmd)

## License

MIT (see the LICENSE file)

## Citation

```bibtex
@software{phaseR,
  title  = {phaseR: Phase-Structured Inference for Longitudinal Change},
  author = {Colling, Gilles},
  year   = {2025},
  url    = {https://github.com/gcol33/phaseR}
}
```
