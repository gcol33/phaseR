# Declare ordered phases

Declare ordered phases

## Usage

``` r
phases(..., absorbing = NULL)
```

## Arguments

- ...:

  Character names of phases in order (e.g., "inactive", "active")

- absorbing:

  Character vector of phase names that are absorbing states. Defaults to
  the last declared phase.

## Value

A `phase_spec` object

## Examples

``` r
# Two-phase model
phases("inactive", "active")

# Multi-phase with explicit absorbing state
phases("susceptible", "infected", "recovered", absorbing = "recovered")
```
