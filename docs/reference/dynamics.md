# Specify phase-specific dynamics

Specify phase-specific dynamics

## Usage

``` r
dynamics(phase, formula, family = gaussian())
```

## Arguments

- phase:

  Character name of the phase

- formula:

  A two-sided formula specifying the outcome model (e.g.,
  `y ~ x1 + x2`).

- family:

  A family object (e.g.,
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`poisson()`](https://rdrr.io/r/stats/family.html)). Default is
  [`gaussian()`](https://rdrr.io/r/stats/family.html).

## Value

A `phase_dynamics` object

## Examples

``` r
# Gaussian outcome
dynamics("active", y ~ treatment + time)

# Poisson outcome
dynamics("infected", cases ~ age + region, family = poisson())
```
