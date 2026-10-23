# Simulate outcomes forward through the phases

Simulate outcomes forward through the phases

## Usage

``` r
simulate_outcomes(s, data)
```

## Arguments

- s:

  Phase summaries from
  [`phase_summaries()`](https://gillescolling.com/phaseR/reference/phase_summaries.md)

- data:

  The data the summaries were computed on

## Value

Matrix of simulated outcomes, one row per draw and one column per row of
`data`
