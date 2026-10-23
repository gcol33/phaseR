# Laplace approximation backend for phase models

Approximates the posterior by a Gaussian at its mode. In a model with
random effects the effects u are integrated out first: for the other
parameters theta (coefficients, log SDs, Cholesky factors) the marginal
posterior is approximated by the Laplace formula
`p(theta | y) ~ p(theta, u*(theta) | y) (2 pi)^(n_u / 2) det(A)^(-1/2)`,
with `u*(theta)` the mode of the effects given theta and `A` the
negative Hessian of the log posterior in u there (Tierney and Kadane
1986; the approach of INLA and TMB). The joint mode of (theta, u) does
not exist: as a random-effect SD shrinks with its effects at zero, the
Normal density of the effects grows without bound, faster than any
proper SD prior falls. The marginal posterior of theta is bounded, so
its mode is well defined. theta is drawn from the Gaussian at that mode,
with covariance the inverse negative Hessian of the log marginal
(central differences of its exact gradient), and u given theta from
`u*(theta_hat) + du*/dtheta (theta - theta_hat) + N(0, A^-1)`. In a
model without random effects this is the Gaussian approximation at the
joint mode.

## Usage

``` r
fit_laplace(stan_data, n_samples = 1000, seed, n_starts = 5L, ...)
```

## Arguments

- stan_data:

  Prepared data list

- n_samples:

  Number of samples to draw from Gaussian approximation

- seed:

  Integer seed for the random starting points and the draws from the
  Gaussian approximation

- n_starts:

  Number of starting points of the mode search

- ...:

  Additional arguments (ignored)

## Value

Raw fit object with same structure as HMC backend

## Details

Given theta the log posterior splits over blocks of effects that share
no unit (one block per unit for unit-level effects, one block per level
of a grouping factor nested in the units), so `A` is block diagonal and
the mode of the effects is found by Newton's method from `u = 0` in each
block separately. When the effects of a block have several modes given
theta, for example unit-level effects in phases whose timing is
uncertain, the approximation can become degenerate where two modes
merge; the fit then stops or, if several modes remain at the optimum,
warns, and the HMC backend is the reliable choice.

The log posterior of a phase model is multimodal: a mode in which a
phase is never entered leaves that phase's parameters at their prior.
The mode is therefore searched from `n_starts` starting points built
from provisional phase assignments
([`initial_values()`](https://gillescolling.com/phaseR/reference/initial_values.md)),
and the start reaching the highest log marginal posterior is kept. With
random effects, each start is first moved to the mode of the fixed
effects with the effects at zero
([`fixed_effects_mode()`](https://gillescolling.com/phaseR/reference/fixed_effects_mode.md)),
and starts that arrive at the same point share one search.
