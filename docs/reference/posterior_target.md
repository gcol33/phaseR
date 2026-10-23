# Describe the log posterior of a phase model

Both backends evaluate the log posterior through this description: the
NUTS sampler in `phaseR_nuts_sampler()` and the Laplace approximation. A
model with `K` phases has `K - 1` transition linear predictors (the
logit of moving from phase `k` to `k + 1`) and `K` dynamics linear
predictors (the mean of each phase on the scale of its link). Each is a
fixed-effect design matrix plus random-effect terms, whose parameters
sit at the offsets of the parameter vector that
[`prepare_stan_data()`](https://gillescolling.com/phaseR/reference/prepare_stan_data.md)
names.

## Usage

``` r
posterior_target(stan_data)
```

## Arguments

- stan_data:

  Prepared data list

## Value

List with `kind = "phase"`, `n_params`, `n_phases`, the outcome `y`, the
0-based observation range of each unit (`unit_start`, `unit_end`), the
outcome family code and number of trials of each phase, the 0-based
offset of the log residual SDs (`sigma_offset`), the linear predictors
(`transitions`, `dynamics`; each a list with design matrix `X`, 0-based
coefficient `offset` and random-effect terms `re`) and the prior of each
parameter class (`prior`, from
[`prior_to_cpp()`](https://gillescolling.com/phaseR/reference/prior_to_cpp.md))
