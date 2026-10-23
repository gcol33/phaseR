#ifndef PHASER_BLOCK_PRECISION_H
#define PHASER_BLOCK_PRECISION_H

#include <algorithm>
#include <cmath>
#include <vector>
#include "dense_cholesky.h"

namespace phaser {

// Symmetric matrix over n indices that is block diagonal once the indices
// are grouped into blocks (the connected components of its sparsity
// pattern), held and factorised as one dense matrix per block. For random
// effects nested in units, such as one intercept per unit, the blocks are a
// few effects each and the work is linear in the number of effects.
class BlockPrecision {
public:
  BlockPrecision() : n_(0) {}

  // Blocks from the sparsity pattern: rows[a] lists the indices coupled
  // with index a (symmetric).
  explicit BlockPrecision(const std::vector<std::vector<int> >& rows)
    : n_(static_cast<int>(rows.size())), block_of_(n_, -1), local_(n_, -1) {
    for (int start = 0; start < n_; start++) {
      if (block_of_[start] >= 0) continue;
      const int b = static_cast<int>(members_.size());
      members_.push_back(std::vector<int>());
      std::vector<int> stack(1, start);
      block_of_[start] = b;
      while (!stack.empty()) {
        const int a = stack.back();
        stack.pop_back();
        members_[b].push_back(a);
        for (int r : rows[a]) {
          if (block_of_[r] < 0) {
            block_of_[r] = b;
            stack.push_back(r);
          }
        }
      }
      std::sort(members_[b].begin(), members_[b].end());
      for (std::size_t i = 0; i < members_[b].size(); i++) {
        local_[members_[b][i]] = static_cast<int>(i);
      }
    }
    A_.resize(members_.size());
    for (std::size_t b = 0; b < members_.size(); b++) {
      A_[b].assign(members_[b].size() * members_[b].size(), 0.0);
    }
    L_ = A_;
  }

  int size() const { return n_; }
  int n_blocks() const { return static_cast<int>(members_.size()); }
  int block_of(int a) const { return block_of_[a]; }
  const std::vector<int>& members(int b) const { return members_[b]; }

  double max_abs_diagonal(int b) const {
    const int m = static_cast<int>(members_[b].size());
    double out = 0.0;
    for (int i = 0; i < m; i++) out = std::max(out, std::fabs(A_[b][i * m + i]));
    return out;
  }

  // Cholesky factor of block b of A + damping I; false if it is not
  // positive definite.
  bool factor_block(int b, double damping) {
    const int m = static_cast<int>(members_[b].size());
    L_[b] = A_[b];
    for (int i = 0; i < m; i++) L_[b][i * m + i] += damping;
    inverse_.clear();
    return cholesky(L_[b], m);
  }

  // x <- (block b of A + damping I)^{-1} x on the indices of block b, after
  // factor_block().
  void solve_block(int b, double* x) const {
    const std::vector<int>& idx = members_[b];
    std::vector<double> local(idx.size());
    for (std::size_t i = 0; i < idx.size(); i++) local[i] = x[idx[i]];
    cholesky_solve(L_[b], static_cast<int>(idx.size()), local.data());
    for (std::size_t i = 0; i < idx.size(); i++) x[idx[i]] = local[i];
  }

  void set_zero() {
    for (std::vector<double>& block : A_) std::fill(block.begin(), block.end(), 0.0);
  }

  // Entry (r, a); r and a must lie in one block.
  double& at(int r, int a) {
    const int b = block_of_[a];
    const int m = static_cast<int>(members_[b].size());
    return A_[b][local_[r] * m + local_[a]];
  }
  double at(int r, int a) const {
    const int b = block_of_[a];
    const int m = static_cast<int>(members_[b].size());
    return A_[b][local_[r] * m + local_[a]];
  }

  // Cholesky factors of the blocks of A + damping I; false if one is not
  // positive definite.
  bool factor(double damping) {
    for (std::size_t b = 0; b < A_.size(); b++) {
      const int m = static_cast<int>(members_[b].size());
      L_[b] = A_[b];
      for (int i = 0; i < m; i++) L_[b][i * m + i] += damping;
      if (!cholesky(L_[b], m)) return false;
    }
    inverse_.clear();
    return true;
  }

  // x <- (A + damping I)^{-1} x, after factor().
  void solve(double* x) const {
    std::vector<double> local;
    for (std::size_t b = 0; b < members_.size(); b++) {
      const std::vector<int>& idx = members_[b];
      local.resize(idx.size());
      for (std::size_t i = 0; i < idx.size(); i++) local[i] = x[idx[i]];
      cholesky_solve(L_[b], static_cast<int>(idx.size()), local.data());
      for (std::size_t i = 0; i < idx.size(); i++) x[idx[i]] = local[i];
    }
  }

  // x <- L^{-T} x with L L^T the factorised matrix: standard normal x
  // becomes a draw with covariance (L L^T)^{-1}.
  void solve_factor_transpose(double* x) const {
    std::vector<double> local;
    for (std::size_t b = 0; b < members_.size(); b++) {
      const std::vector<int>& idx = members_[b];
      const int m = static_cast<int>(idx.size());
      local.resize(m);
      for (int i = 0; i < m; i++) local[i] = x[idx[i]];
      for (int i = m - 1; i >= 0; i--) {
        double s = local[i];
        for (int k = i + 1; k < m; k++) s -= L_[b][k * m + i] * local[k];
        local[i] = s / L_[b][i * m + i];
      }
      for (int i = 0; i < m; i++) x[idx[i]] = local[i];
    }
  }

  double log_det() const {
    double out = 0.0;
    for (std::size_t b = 0; b < L_.size(); b++) {
      out += cholesky_log_det(L_[b], static_cast<int>(members_[b].size()));
    }
    return out;
  }

  // Entry (r, a) of the inverse, r and a in one block, after factor().
  double inverse_at(int r, int a) {
    if (inverse_.empty()) {
      inverse_.resize(L_.size());
      for (std::size_t b = 0; b < L_.size(); b++) {
        inverse_[b] = cholesky_inverse(L_[b], static_cast<int>(members_[b].size()));
      }
    }
    const int b = block_of_[a];
    const int m = static_cast<int>(members_[b].size());
    return inverse_[b][local_[r] * m + local_[a]];
  }

private:
  int n_;
  std::vector<int> block_of_;
  std::vector<int> local_;
  std::vector<std::vector<int> > members_;
  std::vector<std::vector<double> > A_;
  std::vector<std::vector<double> > L_;
  std::vector<std::vector<double> > inverse_;
};

}  // namespace phaser

#endif
