# Data-informed starting points

Starting points of the Laplace mode search and initial values of the
NUTS chains. A phase model's posterior can have a mode in which a phase
is never entered and its parameters stay at their prior; starting all
parameters near zero, far from the data, can end there. Each starting
point is built from a provisional phase assignment: every unit follows a
path that starts in the first phase and moves at most one phase forward
per observation. Given the paths, each phase's dynamics are fitted by a
weighted ridge regression (Gaussian, Poisson or logistic) of the outcome
on that phase's design matrix, using the observations assigned to the
phase (weight 1) and all others (weight `1e-3`, which keeps the fit
defined for a phase no unit enters), and each transition by a ridge
logistic regression of "moved" on the transition design matrix over the
observations at risk of that move. Random effects start at zero, with
their SDs at the scale of the random-effect prior and their correlations
at zero.

## Usage

``` r
initial_values(target, n_starts)
```

## Arguments

- target:

  Model description from
  [`posterior_target()`](https://gillescolling.com/phaseR/reference/posterior_target.md)

- n_starts:

  Number of starting points

## Value

List of `n_starts` parameter vectors

## Details

The first path splits each unit's observations into `K` consecutive runs
of equal length, so every phase is entered; the others draw the number
and the times of each unit's moves at random from R's random number
generator.
