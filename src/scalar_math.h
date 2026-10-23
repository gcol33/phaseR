#ifndef PHASER_SCALAR_MATH_H
#define PHASER_SCALAR_MATH_H

#include <cmath>

// Elementary functions shared by the log posterior code. Every function is a
// template on the scalar type, so the same code evaluates the log posterior
// in double precision and, through the overloads of exp(), log() and log1p()
// for other scalar types found by argument-dependent lookup, its derivatives.

namespace phaser {

const double kLog2Pi = 1.8378770664093454836;  // log(2 pi)
const double kLogPi = 1.1447298858494001741;   // log(pi)
const double kLog2 = 0.69314718055994530942;   // log(2)

// The double value of a scalar, used for branching (which operand of a
// log-sum-exp is larger, the sign of an argument).
inline double value(double x) { return x; }

// log(exp(a) + exp(b)).
template <class T>
T log_sum_exp(const T& a, const T& b) {
  using std::exp;
  using std::log1p;
  if (value(a) > value(b)) return a + log1p(exp(b - a));
  return b + log1p(exp(a - b));
}

// log(sum_i exp(x[i])) for n >= 1 terms.
template <class T>
T log_sum_exp(const T* x, int n) {
  using std::exp;
  using std::log;
  int arg_max = 0;
  for (int i = 1; i < n; i++) {
    if (value(x[i]) > value(x[arg_max])) arg_max = i;
  }
  T sum(1.0);
  for (int i = 0; i < n; i++) {
    if (i != arg_max) sum = sum + exp(x[i] - x[arg_max]);
  }
  return x[arg_max] + log(sum);
}

// log(1 / (1 + exp(-x))), the log of the inverse logit, without overflow.
template <class T>
T log_inv_logit(const T& x) {
  using std::exp;
  using std::log1p;
  if (value(x) < 0.0) return x - log1p(exp(x));
  return -log1p(exp(-x));
}

// log(1 - 1 / (1 + exp(-x))) = log_inv_logit(-x).
template <class T>
T log1m_inv_logit(const T& x) {
  return log_inv_logit(T(-x));
}

// tanh(x) through exp(-2 |x|), without overflow.
template <class T>
T tanh_stable(const T& x) {
  using std::exp;
  if (value(x) >= 0.0) {
    T e = exp(-2.0 * x);
    return (1.0 - e) / (1.0 + e);
  }
  T e = exp(2.0 * x);
  return (e - 1.0) / (e + 1.0);
}

}  // namespace phaser

#endif
