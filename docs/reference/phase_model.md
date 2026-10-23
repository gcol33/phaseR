# Compose a phase model

Compose a phase model

## Usage

``` r
phase_model(phases, ...)
```

## Arguments

- phases:

  A `phase_spec` object from
  [`phases()`](https://gillescolling.com/phaseR/reference/phases.md)

- ...:

  `phase_transition` and `phase_dynamics` objects

## Value

A `phase_model` object

## Examples

``` r
model <- phase_model(
  phases("inactive", "active"),
  transition("inactive", "active", ~ treatment),
  dynamics("inactive", y ~ 1),
  dynamics("active", y ~ treatment)
)
```
