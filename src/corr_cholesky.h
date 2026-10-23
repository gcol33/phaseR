#ifndef PHASER_CORR_CHOLESKY_H
#define PHASER_CORR_CHOLESKY_H

#include <cmath>
#include <vector>
#include "dual.h"
#include "scalar_math.h"

namespace phaser {

// Cholesky factor L (n x n, lower triangular, row-major) of a correlation
// matrix from n (n - 1) / 2 unconstrained values y, in the transform of
// Lewandowski, Kurowicka and Joe (2009) that Stan uses: canonical partial
// correlations z_ij = tanh(y_ij) for i > j, taken row by row, and
//
//   L_00 = 1,  L_i0 = z_i0,  L_ij = z_ij sqrt(1 - sum_{k<j} L_ik^2),
//   L_ii = sqrt(1 - sum_{k<i} L_ik^2).
//
// Returns the log density of y when L L^T ~ LKJ(eta): the LKJ density of the
// Cholesky factor, sum_{i>=1} (n - i - 1 + 2 (eta - 1)) log L_ii, plus the log
// Jacobian of y -> L, sum_{i>j} log(1 - z_ij^2) +
// sum_{i>j>0} log(1 - sum_{k<j} L_ik^2) / 2.
template <class T>
T corr_cholesky(const T* y, int n, double eta, T* L) {
  using std::log;
  using std::log1p;
  using std::sqrt;
  for (int i = 0; i < n * n; i++) L[i] = T(0.0);
  L[0] = T(1.0);
  T log_density(0.0);
  int k = 0;
  for (int i = 1; i < n; i++) {
    T sum_sq(0.0);
    for (int j = 0; j < i; j++) {
      const T z = tanh_stable(y[k++]);
      log_density = log_density + log1p(-(z * z));
      if (j == 0) {
        L[i * n] = z;
      } else {
        log_density = log_density + 0.5 * log1p(-sum_sq);
        L[i * n + j] = z * sqrt(1.0 - sum_sq);
      }
      sum_sq = sum_sq + L[i * n + j] * L[i * n + j];
    }
    L[i * n + i] = sqrt(1.0 - sum_sq);
    log_density = log_density + (n - i - 1 + 2.0 * (eta - 1.0)) * log(L[i * n + i]);
  }
  return log_density;
}

// corr_cholesky() and, unless derivatives is false, its derivatives with
// respect to y by forward-mode dual numbers: dL[m * n * n + i * n + j] is
// dL_ij / dy_m and d_log_density[m] the derivative of the log density.
template <class T>
struct CorrFactor {
  std::vector<T> L;
  std::vector<T> dL;
  std::vector<T> d_log_density;
  T log_density;
};

template <class T>
void corr_factor(const T* y, int n, double eta, bool derivatives,
                 CorrFactor<T>& out) {
  const int n_y = n * (n - 1) / 2;
  out.L.assign(n * n, T(0.0));
  out.log_density = corr_cholesky(y, n, eta, out.L.data());
  if (!derivatives) return;
  out.dL.assign(static_cast<std::size_t>(n_y) * n * n, T(0.0));
  out.d_log_density.assign(n_y, T(0.0));
  std::vector<Dual<T> > yd(n_y), Ld(n * n);
  for (int m = 0; m < n_y; m++) {
    for (int k = 0; k < n_y; k++) yd[k] = Dual<T>(y[k], T(k == m ? 1.0 : 0.0));
    Dual<T> lp = corr_cholesky(yd.data(), n, eta, Ld.data());
    for (int i = 0; i < n * n; i++) out.dL[m * n * n + i] = Ld[i].d;
    out.d_log_density[m] = lp.d;
  }
}

}  // namespace phaser

#endif
