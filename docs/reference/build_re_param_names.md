# Build RE parameter names

Generates parameter names for RE terms including slopes

## Usage

``` r
build_re_param_names(
  component,
  re_name,
  re_info,
  transition = component == "trans"
)
```

## Arguments

- component:

  Label of the linear predictor (e.g. "trans", a phase name)

- re_name:

  Name of grouping variable

- re_info:

  RE info structure from build_multi_re

- transition:

  Whether the term belongs to a transition predictor (effects named
  `u_...`) or a dynamics predictor (`v_...`)

## Value

Character vector of parameter names
