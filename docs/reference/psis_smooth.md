# Apply Pareto smoothing to importance ratios

Fits a generalized Pareto distribution to the tail of the importance
ratios and replaces extreme values with expected order statistics.

## Usage

``` r
psis_smooth(log_ratios)
```

## Arguments

- log_ratios:

  Log importance ratios

## Value

List with smoothed log_weights and Pareto k estimate
