reference_data <- function() {
  dat <- sim_phaseR(n_units = 40, n_times = 4, seed = 21)
  set.seed(22)
  dat$site <- factor((dat$id - 1) %% 8)
  dat$region <- factor((dat$id - 1) %% 8 %% 3)
  dat$z <- sin(dat$id + dat$time)
  dat$count <- rpois(nrow(dat), 3)
  dat$bin <- rbinom(nrow(dat), 1, 0.4)
  dat
}

reference_two_phase <- function(trans, dyn_0, dyn_1, family_0 = gaussian(),
                                family_1 = family_0) {
  phase_model(
    phases("inactive", "active"),
    transition("inactive", "active", trans),
    dynamics("inactive", dyn_0, family = family_0),
    dynamics("active", dyn_1, family = family_1)
  )
}

# Models of every posterior target kind, each with groupings that are
# constant within a unit, and the prior each is evaluated under.
reference_cases <- function() {
  informative <- prior(beta ~ student_t(4, 0.5, 2), sigma ~ half_cauchy(1.5),
                       re ~ half_cauchy(0.7), cor ~ lkj(2))
  default <- prior()
  list(
    two_phase = list(
      model = reference_two_phase(~ x, y ~ x, y ~ x), prior = default),
    two_phase_designs = list(
      model = reference_two_phase(~ x + z, y ~ 1, y ~ x + z),
      prior = informative),
    poisson = list(
      model = reference_two_phase(~ x, count ~ x, count ~ 1,
                                  family_0 = poisson()),
      prior = default),
    binomial = list(
      model = reference_two_phase(~ 1, bin ~ x, bin ~ x + z,
                                  family_0 = binomial()),
      prior = informative),
    mixed_family = list(
      model = reference_two_phase(~ x, count ~ x, count ~ x,
                                  family_0 = poisson(),
                                  family_1 = gaussian()),
      prior = default),
    re = list(
      model = reference_two_phase(~ x + (1 | site), y ~ x + (1 | site),
                                  y ~ x + (1 | site)),
      prior = informative),
    multi_re = list(
      model = reference_two_phase(~ x + (1 | region),
                                  y ~ x + (1 | site) + (1 | region), y ~ x),
      prior = default),
    slopes_correlated = list(
      model = reference_two_phase(~ x, y ~ x + (1 + x | site), y ~ x),
      prior = informative),
    slopes_uncorrelated = list(
      model = reference_two_phase(~ x, y ~ x + (1 + x || site),
                                  y ~ x + (1 | region)),
      prior = default),
    slopes_transition = list(
      model = reference_two_phase(~ x + (1 + x | region), y ~ x + (1 | site),
                                  y ~ 1),
      prior = default),
    k_phase_3 = list(
      model = phase_model(
        phases("A", "B", "C"),
        transition("A", "B", ~ 1),
        transition("B", "C", ~ x),
        dynamics("A", y ~ 1),
        dynamics("B", y ~ x),
        dynamics("C", y ~ 1)
      ),
      prior = default),
    k_phase_4 = list(
      model = phase_model(
        phases("A", "B", "C", "D"),
        transition("A", "B", ~ x),
        transition("B", "C", ~ z),
        transition("C", "D", ~ 1),
        dynamics("A", y ~ x),
        dynamics("B", y ~ 1),
        dynamics("C", y ~ x + z),
        dynamics("D", y ~ 1)
      ),
      prior = informative)
  )
}

reference_target <- function(case, dat = reference_data()) {
  stan_data <- phaseR:::prepare_stan_data(case$model, dat, case$prior)
  phaseR:::posterior_target(stan_data)
}

# Central differences with one Richardson extrapolation step: the error is
# O(h^4) in the step and O(eps |f| / h) in rounding.
richardson_gradient <- function(f, x, h = 1e-4) {
  vapply(seq_along(x), function(i) {
    step <- h * max(1, abs(x[i]))
    central <- function(s) {
      e <- replace(numeric(length(x)), i, s)
      (f(x + e) - f(x - e)) / (2 * s)
    }
    (4 * central(step / 2) - central(step)) / 3
  }, numeric(1))
}

gradient_cases <- function() {
  dat <- reference_data()
  dat$year <- factor(dat$time)
  cases <- reference_cases()
  cases$crossed_year <- list(
    model = reference_two_phase(~ x + (1 | year), y ~ x + (1 | site),
                                y ~ x + (1 | year)),
    prior = prior(re ~ half_cauchy(1)))
  cases$slope_only <- list(
    model = reference_two_phase(~ x, y ~ 1 + (0 + x | site), y ~ x),
    prior = prior())
  cases$correlated_three <- list(
    model = reference_two_phase(~ x, y ~ x + z + (1 + x + z | region), y ~ x),
    prior = prior(cor ~ lkj(3), re ~ half_normal(1)))
  cases$binomial_re <- list(
    model = reference_two_phase(~ x, bin ~ x + (1 | site), bin ~ 1,
                                family_0 = binomial()),
    prior = prior())
  list(dat = dat, cases = cases)
}
