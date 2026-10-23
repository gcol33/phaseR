#ifndef PHASER_PHASE_DATA_H
#define PHASER_PHASE_DATA_H

#include <Rcpp.h>
#include <vector>

namespace phaser {

// Observations and design matrices of a two-phase model. Each phase has its
// own dynamics design matrix, built from its own dynamics() formula, so the
// two phases may use different covariates and different numbers of columns.
struct TwoPhaseData {
  std::vector<double> y;
  Rcpp::NumericMatrix X_trans;
  Rcpp::NumericMatrix X_dyn_0;
  Rcpp::NumericMatrix X_dyn_1;
  int n_units;
  std::vector<int> unit_start;
  std::vector<int> unit_end;

  explicit TwoPhaseData(const Rcpp::List& data)
    : y(Rcpp::as<std::vector<double> >(data["y"])),
      X_trans(Rcpp::as<Rcpp::NumericMatrix>(data["X_trans"])),
      X_dyn_0(Rcpp::as<Rcpp::NumericMatrix>(data["X_dyn_0"])),
      X_dyn_1(Rcpp::as<Rcpp::NumericMatrix>(data["X_dyn_1"])),
      n_units(Rcpp::as<int>(data["n_units"])),
      unit_start(Rcpp::as<std::vector<int> >(data["unit_start"])),
      unit_end(Rcpp::as<std::vector<int> >(data["unit_end"])) {
    check_rows(X_trans, "X_trans");
    check_rows(X_dyn_0, "X_dyn_0");
    check_rows(X_dyn_1, "X_dyn_1");
  }

private:
  void check_rows(const Rcpp::NumericMatrix& X, const char* name) const {
    if (X.nrow() != static_cast<int>(y.size())) {
      Rcpp::stop("Design matrix %s has %d rows for %d observations.",
                 name, X.nrow(), static_cast<int>(y.size()));
    }
  }
};

// Linear predictor x_row' beta.
inline double linear_predictor(const Rcpp::NumericMatrix& X, int row,
                               const Rcpp::NumericVector& beta) {
  if (X.ncol() != beta.size()) {
    Rcpp::stop("Design matrix has %d columns for %d coefficients.",
               X.ncol(), static_cast<int>(beta.size()));
  }
  double eta = 0.0;
  for (int j = 0; j < beta.size(); j++) eta += beta[j] * X(row, j);
  return eta;
}

}  // namespace phaser

#endif
