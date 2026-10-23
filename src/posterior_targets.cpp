#include <Rcpp.h>
#include <functional>
#include <memory>
#include <string>
#include "nuts.h"
#include "phase_model.h"
#include "rng.h"

using namespace Rcpp;

namespace {

typedef std::function<double(const NumericVector&)> LogPosterior;

class FunctionLogDensity : public phaser::LogDensity {
public:
  FunctionLogDensity(int dim, LogPosterior f)
    : phaser::LogDensity(dim), f_(f) {}

protected:
  double evaluate(const NumericVector& q) { return f_(q); }

private:
  LogPosterior f_;
};

// Log density described by `target`: a phase model ("phase", the list built
// by posterior_target() in R, which both the NUTS and the Laplace backend
// evaluate) or a multivariate normal with known moments ("mvnormal").
LogPosterior make_log_posterior(const List& target) {
  const std::string kind = as<std::string>(target["kind"]);
  if (kind == "mvnormal") {
    // Multivariate normal with mean `mu` and precision matrix `precision`;
    // a target with known moments for validating the sampler.
    const List data = target["data"];
    NumericVector mu = data["mu"];
    NumericMatrix precision = data["precision"];
    int n = mu.size();
    return [=](const NumericVector& q) {
      double quad = 0.0;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          quad += (q[i] - mu[i]) * precision(i, j) * (q[j] - mu[j]);
        }
      }
      return -0.5 * quad;
    };
  }

  if (kind == "phase") {
    std::shared_ptr<const phaser::PhaseModel> model =
      std::make_shared<const phaser::PhaseModel>(target);
    return [model](const NumericVector& q) {
      return model->log_posterior(q.begin());
    };
  }
  stop("Unknown posterior target kind '%s'.", kind);
}

}  // namespace


// [[Rcpp::export]]
double phaseR_log_posterior(List target, NumericVector params) {
  if (params.size() != as<int>(target["n_params"])) {
    stop("`params` has length %d; the model has %d parameters.",
         static_cast<int>(params.size()), as<int>(target["n_params"]));
  }
  return make_log_posterior(target)(params);
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

  FunctionLogDensity density(as<int>(target["n_params"]),
                             make_log_posterior(target));
  phaser::Rng rng(seed, chain);
  phaser::NutsSettings settings = {n_iter, n_warmup, max_treedepth,
                                   adapt_delta};
  return phaser::run_nuts(density, rng, settings);
}
