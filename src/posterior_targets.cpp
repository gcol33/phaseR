#include <Rcpp.h>
#include <memory>
#include <string>
#include "log_density.h"
#include "nuts.h"
#include "phase_model.h"
#include "rng.h"

using namespace Rcpp;

namespace {

// Log posterior of a phase model.
class PhaseDensity : public phaser::LogDensity {
public:
  explicit PhaseDensity(const List& target)
    : phaser::LogDensity(as<int>(target["n_params"])), model_(target) {}

protected:
  double evaluate(const double* q, double* grad) const {
    return model_.log_posterior(q, grad);
  }

private:
  phaser::PhaseModel model_;
};

// Multivariate normal with mean `mu` and precision matrix `precision`; a
// target with known moments for validating the sampler.
class MvNormalDensity : public phaser::LogDensity {
public:
  explicit MvNormalDensity(const List& data)
    : phaser::LogDensity(as<NumericVector>(data["mu"]).size()),
      mu_(as<std::vector<double> >(data["mu"])),
      precision_(as<NumericMatrix>(data["precision"])) {}

protected:
  double evaluate(const double* q, double* grad) const {
    const int n = dim();
    double quad = 0.0;
    for (int i = 0; i < n; i++) {
      double row = 0.0;
      for (int j = 0; j < n; j++) row += precision_(i, j) * (q[j] - mu_[j]);
      quad += (q[i] - mu_[i]) * row;
      if (grad != nullptr) grad[i] = -row;
    }
    return -0.5 * quad;
  }

private:
  std::vector<double> mu_;
  NumericMatrix precision_;
};

// Log density described by `target`: a phase model ("phase", the list built
// by posterior_target() in R, which both the NUTS and the Laplace backend
// evaluate) or a multivariate normal with known moments ("mvnormal").
std::unique_ptr<phaser::LogDensity> make_density(const List& target) {
  const std::string kind = as<std::string>(target["kind"]);
  if (kind == "phase") {
    return std::unique_ptr<phaser::LogDensity>(new PhaseDensity(target));
  }
  if (kind == "mvnormal") {
    return std::unique_ptr<phaser::LogDensity>(
      new MvNormalDensity(as<List>(target["data"])));
  }
  stop("Unknown posterior target kind '%s'.", kind);
}

void check_length(const phaser::LogDensity& density,
                  const NumericVector& params) {
  if (params.size() != density.dim()) {
    stop("`params` has length %d; the model has %d parameters.",
         static_cast<int>(params.size()), density.dim());
  }
}

}  // namespace


// [[Rcpp::export]]
SEXP phaseR_density(List target) {
  return XPtr<phaser::LogDensity>(make_density(target).release(), true);
}


// [[Rcpp::export]]
double phaseR_log_density(SEXP density, NumericVector params) {
  XPtr<phaser::LogDensity> d(density);
  check_length(*d, params);
  return d->log_density(params.begin());
}


// [[Rcpp::export]]
List phaseR_log_density_gradient(SEXP density, NumericVector params) {
  XPtr<phaser::LogDensity> d(density);
  check_length(*d, params);
  NumericVector grad(d->dim());
  double value = d->log_density_gradient(params.begin(), grad.begin());
  return List::create(Named("value") = value, Named("gradient") = grad);
}


// [[Rcpp::export]]
List phaseR_nuts_sampler(List target, int n_iter, int n_warmup, int seed,
                         int chain, double adapt_delta, int max_treedepth) {
  if (n_warmup < 0 || n_warmup >= n_iter) {
    stop("`warmup` must be non-negative and smaller than `iter`.");
  }
  if (!(adapt_delta > 0.0 && adapt_delta < 1.0)) {
    stop("`adapt_delta` must lie strictly between 0 and 1.");
  }
  if (max_treedepth < 1) {
    stop("`max_treedepth` must be at least 1.");
  }

  std::unique_ptr<phaser::LogDensity> density = make_density(target);
  phaser::Rng rng(seed, chain);
  phaser::NutsSettings settings = {n_iter, n_warmup, max_treedepth,
                                   adapt_delta};
  return phaser::run_nuts(*density, rng, settings);
}
