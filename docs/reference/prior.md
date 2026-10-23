# Specify priors for phase model

Creates a prior specification object with either simple SD parameters or
full distribution specifications using formula syntax.

## Usage

``` r
prior(..., beta_sd = 2.5, sigma_sd = 1, re_sd = 0.5)
```

## Arguments

- ...:

  Prior specifications using formula syntax (e.g.,
  `beta ~ normal(0, 2.5)`) or named SD parameters for backward
  compatibility (`beta_sd = 2.5`).

- beta_sd:

  Prior SD for regression coefficients (default 2.5, ignored if formula
  syntax is used)

- sigma_sd:

  Prior SD for log(sigma) parameters (default 1, ignored if formula
  syntax is used)

- re_sd:

  Prior SD for random effect hyperparameter (default 0.5, ignored if
  formula syntax is used)

## Value

A `phase_prior` object

## Details

The default priors are weakly informative:

- Regression coefficients: Normal(0, 2.5)

- Residual SD (sigma): Half-Normal(0, 1)

- Random effect SD: Half-Normal(0, 0.5)

- Correlation (LKJ): LKJ(2)

### Supported distributions

- `normal(mean, sd)`:

  Normal distribution

- `student_t(df, mean, scale)`:

  Student-t distribution

- `cauchy(location, scale)`:

  Cauchy distribution

- `half_normal(sd)`:

  Half-Normal (positive values only)

- `half_cauchy(scale)`:

  Half-Cauchy (positive values only)

- `lkj(eta)`:

  LKJ prior for correlation matrices

### Parameter classes

- `beta`:

  All regression coefficients (fixed effects)

- `sigma`:

  Residual standard deviations

- `re`:

  Random effect standard deviations

- `cor`:

  Correlation matrices (for correlated random slopes)

## Examples

``` r
# Default priors
prior()

# Simple SD specification (backward compatible)
prior(beta_sd = 1)

# Full distribution specification
prior(
  beta ~ normal(0, 2.5),
  sigma ~ half_cauchy(1),
  re ~ half_normal(0.5),
  cor ~ lkj(2)
)

# Student-t for robust regression coefficients
prior(beta ~ student_t(3, 0, 2.5))
```
