# Extract random effects

Extract random effects

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

A data frame with random effect estimates

## Examples

``` r
if (FALSE) { # \dontrun{
ranef(fit)
ranef(fit, component = "transition")
ranef(fit, component = "dynamics", phase = "active")
} # }
```
