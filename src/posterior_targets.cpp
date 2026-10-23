#include <Rcpp.h>
#include <functional>
#include <string>
#include "nuts.h"
#include "prior_functions.h"
#include "rng.h"

using namespace Rcpp;

double phase_log_posterior(NumericVector params, List data,
                           int n_trans_coef, int n_dyn_coef_0, int n_dyn_coef_1,
                           const prior::PriorSet& pr);

double phase_log_posterior_re(NumericVector params, List data,
                              int n_trans_coef, int n_dyn_coef_0, int n_dyn_coef_1,
                              int n_trans_re, int n_dyn_re_0, int n_dyn_re_1,
                              bool has_trans_re, bool has_dyn_re_0, bool has_dyn_re_1,
                              const prior::PriorSet& pr);

double phase_log_posterior_glm(NumericVector params, List data,
                               int n_trans_coef, int n_dyn_coef_0, int n_dyn_coef_1,
                               int family_0, int family_1,
                               int n_trials_0, int n_trials_1,
                               const prior::PriorSet& pr);

double phase_log_posterior_k(NumericVector params, List data,
                             IntegerVector n_dyn_coef, IntegerVector n_trans_coef,
                             int k_phases, const prior::PriorSet& pr);

double phase_log_posterior_multi_re(NumericVector params, List data,
                                    int n_trans_coef, int n_dyn_coef_0, int n_dyn_coef_1,
                                    List trans_re_info, List dyn_re_info_0,
                                    List dyn_re_info_1, const prior::PriorSet& pr);

double phase_log_posterior_slopes(NumericVector params, List data,
                                  int n_trans_coef, int n_dyn_coef_0, int n_dyn_coef_1,
                                  List trans_re_info, List dyn_re_info_0,
                                  List dyn_re_info_1, const prior::PriorSet& pr);

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

// Log posterior of the model described by `target`, the list built by
// posterior_target() in R: `kind` (the model family), `n_params`, `data`,
// `args` (the model dimensions and random-effect structures) and `prior`
// (the prior of each parameter class, from prior_to_cpp()). Both the NUTS
// and the Laplace backend evaluate this function.
LogPosterior make_log_posterior(const List& target) {
  const std::string kind = as<std::string>(target["kind"]);
  const List data = target["data"];
  const List args = target["args"];

  if (kind == "mvnormal") {
    // Multivariate normal with mean `mu` and precision matrix `precision`;
    // a target with known moments for validating the sampler.
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

  const prior::PriorSet pr(as<List>(target["prior"]));

  if (kind == "two_phase") {
    int nt = args["n_trans_coef"], n0 = args["n_dyn_coef_0"],
        n1 = args["n_dyn_coef_1"];
    return [=](const NumericVector& q) {
      return phase_log_posterior(q, data, nt, n0, n1, pr);
    };
  }
  if (kind == "re") {
    int nt = args["n_trans_coef"], n0 = args["n_dyn_coef_0"],
        n1 = args["n_dyn_coef_1"];
    int rt = args["n_trans_re"], r0 = args["n_dyn_re_0"],
        r1 = args["n_dyn_re_1"];
    bool ht = args["has_trans_re"], h0 = args["has_dyn_re_0"],
         h1 = args["has_dyn_re_1"];
    return [=](const NumericVector& q) {
      return phase_log_posterior_re(q, data, nt, n0, n1, rt, r0, r1,
                                    ht, h0, h1, pr);
    };
  }
  if (kind == "glm") {
    int nt = args["n_trans_coef"], n0 = args["n_dyn_coef_0"],
        n1 = args["n_dyn_coef_1"];
    int f0 = args["family_0"], f1 = args["family_1"];
    int t0 = args["n_trials_0"], t1 = args["n_trials_1"];
    return [=](const NumericVector& q) {
      return phase_log_posterior_glm(q, data, nt, n0, n1, f0, f1, t0, t1, pr);
    };
  }
  if (kind == "k_phase") {
    IntegerVector nd = args["n_dyn_coef"], nt = args["n_trans_coef"];
    int k = args["k_phases"];
    return [=](const NumericVector& q) {
      return phase_log_posterior_k(q, data, nd, nt, k, pr);
    };
  }
  if (kind == "multi_re" || kind == "slopes") {
    int nt = args["n_trans_coef"], n0 = args["n_dyn_coef_0"],
        n1 = args["n_dyn_coef_1"];
    List rt = args["trans_re_info"], r0 = args["dyn_re_info_0"],
         r1 = args["dyn_re_info_1"];
    if (kind == "multi_re") {
      return [=](const NumericVector& q) {
        return phase_log_posterior_multi_re(q, data, nt, n0, n1, rt, r0, r1, pr);
      };
    }
    return [=](const NumericVector& q) {
      return phase_log_posterior_slopes(q, data, nt, n0, n1, rt, r0, r1, pr);
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
