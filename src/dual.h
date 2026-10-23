#ifndef PHASER_DUAL_H
#define PHASER_DUAL_H

#include <cmath>
#include "scalar_math.h"

namespace phaser {

// Forward-mode dual number v + d * epsilon with epsilon^2 = 0. Evaluating a
// function on Dual<double> arguments seeded with d = direction gives its
// directional derivative in d; nesting, Dual<Dual<double> >, gives second
// directional derivatives. Applied to the exact gradient of the log
// posterior, Dual<double> gives Hessian-vector products and
// Dual<Dual<double> > third derivatives contracted with two directions.
template <class T>
struct Dual {
  T v;
  T d;

  Dual() : v(0.0), d(0.0) {}
  Dual(double value) : v(value), d(0.0) {}
  Dual(const T& value, const T& derivative) : v(value), d(derivative) {}
};

template <class T>
inline double value(const Dual<T>& x) { return value(x.v); }

template <class T>
inline Dual<T> operator-(const Dual<T>& a) { return Dual<T>(-a.v, -a.d); }

template <class T>
inline Dual<T> operator+(const Dual<T>& a, const Dual<T>& b) {
  return Dual<T>(a.v + b.v, a.d + b.d);
}
template <class T>
inline Dual<T> operator+(const Dual<T>& a, double b) {
  return Dual<T>(a.v + b, a.d);
}
template <class T>
inline Dual<T> operator+(double a, const Dual<T>& b) {
  return Dual<T>(a + b.v, b.d);
}

template <class T>
inline Dual<T> operator-(const Dual<T>& a, const Dual<T>& b) {
  return Dual<T>(a.v - b.v, a.d - b.d);
}
template <class T>
inline Dual<T> operator-(const Dual<T>& a, double b) {
  return Dual<T>(a.v - b, a.d);
}
template <class T>
inline Dual<T> operator-(double a, const Dual<T>& b) {
  return Dual<T>(a - b.v, -b.d);
}

template <class T>
inline Dual<T> operator*(const Dual<T>& a, const Dual<T>& b) {
  return Dual<T>(a.v * b.v, a.d * b.v + a.v * b.d);
}
template <class T>
inline Dual<T> operator*(const Dual<T>& a, double b) {
  return Dual<T>(a.v * b, a.d * b);
}
template <class T>
inline Dual<T> operator*(double a, const Dual<T>& b) {
  return Dual<T>(a * b.v, a * b.d);
}

template <class T>
inline Dual<T> operator/(const Dual<T>& a, const Dual<T>& b) {
  T q = a.v / b.v;
  return Dual<T>(q, (a.d - q * b.d) / b.v);
}
template <class T>
inline Dual<T> operator/(const Dual<T>& a, double b) {
  return Dual<T>(a.v / b, a.d / b);
}
template <class T>
inline Dual<T> operator/(double a, const Dual<T>& b) {
  T q = a / b.v;
  return Dual<T>(q, -q * b.d / b.v);
}

template <class T>
inline Dual<T> exp(const Dual<T>& a) {
  using std::exp;
  T e = exp(a.v);
  return Dual<T>(e, e * a.d);
}

template <class T>
inline Dual<T> log(const Dual<T>& a) {
  using std::log;
  return Dual<T>(log(a.v), a.d / a.v);
}

template <class T>
inline Dual<T> log1p(const Dual<T>& a) {
  using std::log1p;
  return Dual<T>(log1p(a.v), a.d / (1.0 + a.v));
}

template <class T>
inline Dual<T> sqrt(const Dual<T>& a) {
  using std::sqrt;
  T r = sqrt(a.v);
  return Dual<T>(r, a.d / (2.0 * r));
}

}  // namespace phaser

#endif
