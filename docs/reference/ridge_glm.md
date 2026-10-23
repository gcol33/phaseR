# Weighted ridge regression with a canonical link

Maximises `sum(w * loglik(y | X beta)) - lambda / 2 * sum(beta^2)` by
Newton steps with step halving, for the Gaussian (identity link),
Poisson (log link) and binomial (logit link, `y` a proportion of one
trial) families.

## Usage

``` r
ridge_glm(X, y, w, family, lambda = 0.001)
```

## Arguments

- X:

  Design matrix

- y:

  Outcome

- w:

  Non-negative observation weights

- family:

  `"gaussian"`, `"poisson"` or `"binomial"`

- lambda:

  Ridge penalty

## Value

Coefficient vector
