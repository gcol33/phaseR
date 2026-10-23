# Specify a phase transition

Specify a phase transition

## Usage

``` r
transition(from, to, formula = ~1)
```

## Arguments

- from:

  Character name of the origin phase

- to:

  Character name of the destination phase

- formula:

  A one-sided formula specifying covariates for the transition
  probability (e.g., `~ x1 + x2`). Default is `~ 1` (intercept only).

## Value

A `phase_transition` object

## Examples

``` r
# Transition with covariates
transition("inactive", "active", ~ treatment + age)

# Intercept-only transition
transition("stage1", "stage2")
```
