# Prepare data for C++ likelihood

Builds the design of a phase model with `K` phases: one transition
linear predictor per consecutive pair of phases (from the
[`transition()`](https://gillescolling.com/phaseR/reference/transition.md)
between them, intercept-only if there is none) and one dynamics linear
predictor per phase (from its
[`dynamics()`](https://gillescolling.com/phaseR/reference/dynamics.md),
an intercept-only Gaussian model of the response if there is none). Each
may contain random-effect terms, and each phase has its own outcome
family.

## Usage

``` r
prepare_stan_data(model, data, prior = NULL)
```

## Arguments

- model:

  A `phase_model`

- data:

  Data frame

- prior:

  Prior specification from
  [`prior()`](https://gillescolling.com/phaseR/reference/prior.md)

## Value

List with the outcome, the unit ranges, the design matrices
(`X_trans_list`, `X_dyn_list`), the random-effect structures
(`re_trans`, `re_dyn`, from
[`build_multi_re()`](https://gillescolling.com/phaseR/reference/build_multi_re.md)),
the family code of each phase, the parameter count and names, and the
prior in C++ form

## Details

The parameter vector holds the transition coefficients, the dynamics
coefficients and the log residual SD of each phase, followed by the
random-effect parameters of each predictor in the same order.
