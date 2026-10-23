# Combine diagnostics from multiple chains

Combine diagnostics from multiple chains

## Usage

``` r
combine_diagnostics(chain_list, max_treedepth)
```

## Arguments

- chain_list:

  List of per-chain results of `phaseR_nuts_sampler()`

- max_treedepth:

  Maximum tree depth the chains ran with

## Value

List with the post-warmup number of divergent transitions, the number of
transitions that hit `max_treedepth`, the mean acceptance statistic, and
per chain the step size, mean number of leapfrog steps, E-BFMI and
adapted inverse metric
