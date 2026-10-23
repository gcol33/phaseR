# Extract random effects

Posterior summaries of the random-effect coefficients of every level:
the effects themselves for intercept and uncorrelated terms, and
`diag(sigma) L z` for correlated terms.

## Usage

``` r
ranef(
  object,
  component = c("all", "transition", "dynamics"),
  phase = NULL,
  summary = TRUE
)
```

## Arguments

- object:

  A `phase_fit` object

- component:

  Which component: "transition", "dynamics", or "all" (default)

- phase:

  For dynamics RE, which phase (ignored for transition)

- summary:

  If TRUE (default), return posterior summary; if FALSE, return draws

## Value

With `summary = TRUE`, a data frame with one row per predictor, grouping
term, level and coefficient and the posterior mean, SD and quantiles.
With `summary = FALSE`, a named list of arrays of draws indexed
`[draw, level, coefficient]`, one per predictor and term.

## Examples

``` r
if (FALSE) { # \dontrun{
ranef(fit)
ranef(fit, component = "transition")
ranef(fit, component = "dynamics", phase = "active")
} # }
```
