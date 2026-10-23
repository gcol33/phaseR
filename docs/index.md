# phaseR

*longitudinal change as a switch between regimes*

phaseR fits Bayesian models of change in longitudinal data, such as
disease onset and recovery, treatment adoption or an ecological regime
shift, in which each unit moves from a “before” regime to an “after”
regime.

The package treats “before” and “after” as latent regimes and models two
things jointly: the transition that moves a unit between regimes, and
the outcome dynamics within each regime. The likelihood sums over the
regime each observation may belong to, so the estimates do not condition
on the realized transition. Splitting the data at the observed change
point and comparing the two halves does condition on it, which biases
the comparison whenever the chance of transitioning depends on the same
covariates that drive the outcome. From a fit you can estimate regime
effects and counterfactual outcomes, and decide which comparison answers
your question. phaseR supports:

- **More than two regimes.** Ordered phases such as
  `phases("susceptible", "infected", "recovered")`, with one
  [`transition()`](https://gillescolling.com/phaseR/reference/transition.md)
  per consecutive pair.
- **Non-Gaussian outcomes.** Gaussian, Poisson and binomial dynamics,
  set per phase.
- **Random effects.** Random intercepts `(1 | group)`, crossed
  `(1 | site) + (1 | year)`, random slopes `(1 + x | group)` and nested
  `(1 | site/plot)`.
- **Custom priors.** [`prior()`](#usage) with
  [`normal()`](https://gillescolling.com/phaseR/reference/normal.md),
  [`student_t()`](https://gillescolling.com/phaseR/reference/student_t.md),
  [`cauchy()`](https://gillescolling.com/phaseR/reference/cauchy.md),
  [`half_normal()`](https://gillescolling.com/phaseR/reference/half_normal.md),
  [`half_cauchy()`](https://gillescolling.com/phaseR/reference/half_cauchy.md)
  and [`lkj()`](https://gillescolling.com/phaseR/reference/lkj.md).
- **Counterfactual effects.**
  [`ate()`](https://gillescolling.com/phaseR/reference/ate.md),
  [`ite()`](https://gillescolling.com/phaseR/reference/ite.md) and
  [`counterfactual()`](https://gillescolling.com/phaseR/reference/counterfactual.md)
  predict outcomes under covariate values that were never observed,
  integrating over phase uncertainty.
- **Model comparison.** WAIC with
  [`waic()`](https://gillescolling.com/phaseR/reference/waic.md) and
  [`compare_waic()`](https://gillescolling.com/phaseR/reference/compare_waic.md),
  and leave-one-out cross-validation by Pareto-smoothed importance
  sampling with
  [`loo()`](https://gillescolling.com/phaseR/reference/loo.md) and
  [`loo_compare()`](https://gillescolling.com/phaseR/reference/loo_compare.md).

The default backend is the No-U-Turn sampler (NUTS), an adaptive
Hamiltonian Monte Carlo method, written in C++ with Rcpp, so no Stan,
JAGS or other inference engine has to be installed. [Why Conditioning on
Realized Change
Fails](https://gillescolling.com/phaseR/vignettes/estimand-failure.Rmd)
runs the before/after estimate and the phaseR estimate against a known
truth on simulated data.

## Installation

The package compiles C++ code. Install the build tools for your
operating system first:

| Operating system | Build tools                                              |
|------------------|----------------------------------------------------------|
| Windows          | [Rtools](https://cran.r-project.org/bin/windows/Rtools/) |
| macOS            | Run `xcode-select --install` in Terminal                 |
| Debian / Ubuntu  | Install `r-base-dev`                                     |

Check that R can find the compiler:

``` r

pkgbuild::has_build_tools(debug = TRUE)
```

Then install phaseR:

``` r

# install.packages("pak")
pak::pak("gcol33/phaseR")
```

The C++ code is compiled with Rcpp to the C++17 standard and needs no
system libraries.

## Usage

1.  Declare the regimes with
    [`phases()`](https://gillescolling.com/phaseR/reference/phases.md),
    the transition model with
    [`transition()`](https://gillescolling.com/phaseR/reference/transition.md)
    and the outcome model of each regime with
    [`dynamics()`](https://gillescolling.com/phaseR/reference/dynamics.md),
    and combine them with
    [`phase_model()`](https://gillescolling.com/phaseR/reference/phase_model.md).
2.  Fit the model to long-format data (one row per unit and time) with
    [`fit_phaseR()`](https://gillescolling.com/phaseR/reference/fit_phaseR.md).
3.  Read the coefficients with
    [`summary()`](https://rdrr.io/r/base/summary.html), and estimate
    effects with
    [`ate()`](https://gillescolling.com/phaseR/reference/ate.md),
    [`ite()`](https://gillescolling.com/phaseR/reference/ite.md) or
    [`counterfactual()`](https://gillescolling.com/phaseR/reference/counterfactual.md).

[`sim_phaseR()`](https://gillescolling.com/phaseR/reference/sim_phaseR.md)
simulates data with a unit identifier `id`, a time index `time`, a
covariate `x`, an outcome `y` and the true phase. The transition depends
on `x`, and each regime has its own regression of `y` on `x`:

``` r

library(phaseR)

dat <- sim_phaseR(n_units = 100, n_times = 5, seed = 42)

model <- phase_model(
  phases("inactive", "active"),
  transition("inactive", "active", ~ x),
  dynamics("inactive", y ~ x),
  dynamics("active", y ~ x)
)

fit <- fit_phaseR(model, data = dat, backend = "laplace", seed = 42)
summary(fit)
#> Phase Model Summary
#> Phases: inactive -> active
#> 
#> Transition coefficients:
#>                         mean   sd  2.5% median 97.5% Rhat    ESS
#> beta_trans_(Intercept) -0.92 0.23 -1.39  -0.92 -0.49   NA 878.85
#> beta_trans_x            0.42 0.25 -0.05   0.42  0.93   NA 818.19
#> 
#> Dynamics coefficients (inactive):
#>                           mean   sd  2.5% median 97.5% Rhat     ESS
#> beta_inactive_(Intercept) 9.88 0.08  9.72   9.88 10.03   NA  986.47
#> beta_inactive_x           0.06 0.07 -0.06   0.06  0.19   NA 1053.52
#> 
#> Dynamics coefficients (active):
#>                          mean   sd  2.5% median 97.5% Rhat    ESS
#> beta_active_(Intercept)  7.95 0.15  7.65   7.96  8.24   NA 966.92
#> beta_active_x           -0.47 0.13 -0.71  -0.47 -0.22   NA 986.75
#> 
#> Variance parameters:
#>                    mean   sd  2.5% median 97.5% Rhat     ESS
#> log_sigma_inactive 0.04 0.05 -0.05   0.04  0.14   NA  865.83
#> log_sigma_active   0.41 0.06  0.30   0.40  0.52   NA 1025.50
```

`backend = "laplace"` approximates the posterior by a Gaussian at its
mode, with any random effects integrated out first, and returns within a
second. The default backend, `"hmc"`, samples the posterior with the
No-U-Turn sampler, and `chains` and `iter` set the number of chains and
the iterations per chain:

``` r

fit <- fit_phaseR(model, data = dat, chains = 2, iter = 1000, seed = 42)
```

## Why splitting at the change point is biased

Comparing the outcomes of units that transitioned, before against after,
stratifies on a realized outcome. When units with certain covariate
values are both more likely to transition and have different outcomes,
the measured difference mixes the regime effect with that selection, and
the bias does not shrink with more data.

phaseR writes down the full generative process instead: a logistic
transition model for regime membership and an outcome model per regime,
fitted jointly. The likelihood sums over the latent regime of every
observation with a log-sum-exp over phase assignments, so the selection
mechanism is part of the model. [Why Conditioning on Realized Change
Fails](https://gillescolling.com/phaseR/vignettes/estimand-failure.Rmd)
shows the bias on simulated data with a known effect.

## Counterfactual effects

Because the transition and the dynamics are modelled together, a fit
predicts outcomes under covariate values that were never observed,
integrating over phase uncertainty.
[`ate()`](https://gillescolling.com/phaseR/reference/ate.md) returns the
average effect of setting `x` to 1 against 0, with a 95% credible
interval:

``` r

effect <- ate(fit, treatment_var = "x")
effect$ate
#> [1] -0.3790779
effect$ate_ci
#>       2.5%      97.5% 
#> -0.7154333 -0.0678609
```

[`counterfactual()`](https://gillescolling.com/phaseR/reference/counterfactual.md)
predicts outcomes for new covariate values, and
[`ite()`](https://gillescolling.com/phaseR/reference/ite.md) returns the
effect per unit:

``` r

counterfactual(fit, newdata = transform(dat, x = 1))
ite(fit, treatment_var = "x")
```

## Models with more regimes

A model with three regimes declares the phases in order, one transition
per consecutive pair and one outcome model per regime:

``` r

model <- phase_model(
  phases("susceptible", "infected", "recovered"),
  transition("susceptible", "infected", ~ exposure + age),
  transition("infected", "recovered", ~ treatment),
  dynamics("susceptible", y ~ 1),
  dynamics("infected", y ~ time_since_infection),
  dynamics("recovered", y ~ 1)
)
```

## Documentation

- [Why Conditioning on Realized Change
  Fails](https://gillescolling.com/phaseR/articles/estimand-failure.html):
  the before/after estimate and the phaseR estimate against a known
  effect on simulated data.
- The help pages document each function, starting with
  [`?phase_model`](https://gillescolling.com/phaseR/reference/phase_model.md)
  and
  [`?fit_phaseR`](https://gillescolling.com/phaseR/reference/fit_phaseR.md).

## Getting help

Questions and bug reports go to the [GitHub
issues](https://github.com/gcol33/phaseR/issues). A bug report is
easiest to act on with a small reproducible example, for instance made
with [reprex](https://reprex.tidyverse.org/).

## Support

> “Software is like sex: it’s better when it’s free.” — Linus Torvalds

I’m a PhD student who builds R packages in my free time because I
believe good tools should be free and open. I started these projects for
my own work and figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to
say thanks. It helps with my coffee addiction.

[![Buy Me A
Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

phaseR is licensed under MIT.

## Citation

``` text
Colling G (2026). _phaseR: Phase-Structured Inference for Longitudinal
Change_. R package version 0.8.0, <https://github.com/gcol33/phaseR>.
```

BibTeX:

``` bibtex
@Manual{colling2026phaser,
  title = {phaseR: Phase-Structured Inference for Longitudinal Change},
  author = {Gilles Colling},
  year = {2026},
  note = {R package version 0.8.0},
  url = {https://github.com/gcol33/phaseR},
}
```

Download the citation as
[BibTeX](https://gillescolling.com/phaseR/phaseR.bib) or
[RIS](https://gillescolling.com/phaseR/phaseR.ris).
