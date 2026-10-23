# Provisional phase of every observation

Provisional phase of every observation

## Usage

``` r
phase_paths(target, random)
```

## Arguments

- target:

  Model description from
  [`posterior_target()`](https://gillescolling.com/phaseR/reference/posterior_target.md)

- random:

  If `FALSE`, split each unit into `K` runs of equal length; if `TRUE`,
  draw the number of moves uniformly from those possible and their times
  uniformly without replacement

## Value

Integer vector of 0-based phases, one per observation
