#ifndef PHASER_DENSE_CHOLESKY_H
#define PHASER_DENSE_CHOLESKY_H

#include <cmath>
#include <vector>

namespace phaser {

// Cholesky factor L (lower triangular, A = L L^T) of a symmetric positive
// definite n x n matrix stored row-major. Overwrites the lower triangle of A
// with L and zeroes the upper triangle. Returns false if A is not positive
// definite.
inline bool cholesky(std::vector<double>& A, int n) {
  for (int j = 0; j < n; j++) {
    double diag = A[j * n + j];
    for (int k = 0; k < j; k++) diag -= A[j * n + k] * A[j * n + k];
    if (!(diag > 0.0)) return false;
    const double L_jj = std::sqrt(diag);
    A[j * n + j] = L_jj;
    for (int i = j + 1; i < n; i++) {
      double s = A[i * n + j];
      for (int k = 0; k < j; k++) s -= A[i * n + k] * A[j * n + k];
      A[i * n + j] = s / L_jj;
    }
    for (int i = 0; i < j; i++) A[i * n + j] = 0.0;
  }
  return true;
}

// Solves L L^T x = b in place, L from cholesky().
inline void cholesky_solve(const std::vector<double>& L, int n, double* b) {
  for (int i = 0; i < n; i++) {
    double s = b[i];
    for (int k = 0; k < i; k++) s -= L[i * n + k] * b[k];
    b[i] = s / L[i * n + i];
  }
  for (int i = n - 1; i >= 0; i--) {
    double s = b[i];
    for (int k = i + 1; k < n; k++) s -= L[k * n + i] * b[k];
    b[i] = s / L[i * n + i];
  }
}

// log det(L L^T).
inline double cholesky_log_det(const std::vector<double>& L, int n) {
  double s = 0.0;
  for (int i = 0; i < n; i++) s += std::log(L[i * n + i]);
  return 2.0 * s;
}

// (L L^T)^{-1}, row-major.
inline std::vector<double> cholesky_inverse(const std::vector<double>& L,
                                            int n) {
  std::vector<double> inv(static_cast<std::size_t>(n) * n, 0.0);
  std::vector<double> col(n);
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) col[i] = (i == j) ? 1.0 : 0.0;
    cholesky_solve(L, n, col.data());
    for (int i = 0; i < n; i++) inv[i * n + j] = col[i];
  }
  return inv;
}

}  // namespace phaser

#endif
