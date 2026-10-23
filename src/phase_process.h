#ifndef PHASER_PHASE_PROCESS_H
#define PHASER_PHASE_PROCESS_H

#include "scalar_math.h"

namespace phaser {

// Latent phase process of one unit observed at n times, with K ordered
// phases: the unit is in phase 0 at its first observation; between
// consecutive observations it either stays in its phase or moves to the
// next one, and the last phase is absorbing. The arrays are indexed by
// observation s = 0, ..., n - 1 and phase k as [s * K + k]:
//
//   log_emit  log density of observation s given phase k;
//   log_stay  log probability of staying in phase k between observations
//             s - 1 and s (s >= 1, k < K - 1);
//   log_move  log probability of moving from phase k to phase k + 1 between
//             observations s - 1 and s (s >= 1, k < K - 1).
//
// The emission and transition terms are computed by the model from its
// families and linear predictors; the recursion only combines them.

// Phases reachable at observation s.
inline int n_reachable(int s, int K) { return (s + 1 < K) ? s + 1 : K; }

// Forward recursion: alpha[s * K + k] = log p(y_0, ..., y_s, phase_s = k)
// for the reachable phases k < n_reachable(s, K). Returns the log
// likelihood of the unit's observations, log sum_k exp(alpha[(n-1) * K + k]).
template <class T>
T forward(int n, int K, const T* log_emit, const T* log_stay,
          const T* log_move, T* alpha) {
  alpha[0] = log_emit[0];
  for (int s = 1; s < n; s++) {
    const T* prev = alpha + (s - 1) * K;
    T* cur = alpha + s * K;
    const int reach_prev = n_reachable(s - 1, K);
    const int reach = n_reachable(s, K);
    for (int k = 0; k < reach; k++) {
      const bool from_same = k < reach_prev;
      const bool from_previous = k > 0;
      T in_same, in_previous;
      if (from_same) {
        in_same = (k < K - 1) ? T(prev[k] + log_stay[s * K + k]) : prev[k];
      }
      if (from_previous) {
        in_previous = prev[k - 1] + log_move[s * K + k - 1];
      }
      T log_in;
      if (from_same && from_previous) {
        log_in = log_sum_exp(in_same, in_previous);
      } else if (from_same) {
        log_in = in_same;
      } else {
        log_in = in_previous;
      }
      cur[k] = log_in + log_emit[s * K + k];
    }
  }
  return log_sum_exp(alpha + (n - 1) * K, n_reachable(n - 1, K));
}

// Derivatives of the log likelihood returned by forward() with respect to
// its inputs, from the backward recursion
// beta[s * K + k] = log p(y_{s+1}, ..., y_{n-1} | phase_s = k).
// The log likelihood is the log of a sum over phase paths of the exponent
// of the summed terms along each path, so its derivative with respect to a
// term is the posterior probability that the path uses it:
//
//   d_emit[s * K + k]  p(phase_s = k | y);
//   d_stay[s * K + k]  p(phase_{s-1} = k, phase_s = k | y);
//   d_move[s * K + k]  p(phase_{s-1} = k, phase_s = k + 1 | y).
//
// alpha and log_lik are the outputs of forward(); beta is workspace of the
// size of alpha. Entries of unreachable phases are set to zero.
template <class T>
void backward(int n, int K, const T* log_emit, const T* log_stay,
              const T* log_move, const T* alpha, const T& log_lik, T* beta,
              T* d_emit, T* d_stay, T* d_move) {
  using std::exp;
  for (int k = 0; k < K; k++) beta[(n - 1) * K + k] = T(0.0);
  for (int s = n - 1; s >= 1; s--) {
    const T* next = beta + s * K;
    const T* e = log_emit + s * K;
    const int reach_prev = n_reachable(s - 1, K);
    for (int k = 0; k < reach_prev; k++) {
      if (k < K - 1) {
        beta[(s - 1) * K + k] =
          log_sum_exp(T(log_stay[s * K + k] + e[k] + next[k]),
                      T(log_move[s * K + k] + e[k + 1] + next[k + 1]));
      } else {
        beta[(s - 1) * K + k] = e[k] + next[k];
      }
    }
  }

  for (int s = 0; s < n; s++) {
    const int reach = n_reachable(s, K);
    for (int k = 0; k < K; k++) {
      const int i = s * K + k;
      d_emit[i] = (k < reach) ? T(exp(alpha[i] + beta[i] - log_lik)) : T(0.0);
      d_stay[i] = T(0.0);
      d_move[i] = T(0.0);
    }
    if (s == 0) continue;
    const T* prev = alpha + (s - 1) * K;
    const T* e = log_emit + s * K;
    const T* next = beta + s * K;
    const int reach_prev = n_reachable(s - 1, K);
    for (int k = 0; k < reach_prev && k < K - 1; k++) {
      d_stay[s * K + k] =
        exp(prev[k] + log_stay[s * K + k] + e[k] + next[k] - log_lik);
      d_move[s * K + k] =
        exp(prev[k] + log_move[s * K + k] + e[k + 1] + next[k + 1] - log_lik);
    }
  }
}

}  // namespace phaser

#endif
