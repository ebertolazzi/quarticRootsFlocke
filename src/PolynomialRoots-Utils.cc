/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-function"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wvla-extension"
#pragma clang diagnostic ignored "-Wvla"
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

#include "PolynomialRoots.hh"
#include "PolynomialRoots-Utils.hh"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace PolynomialRoots {

  // static real_type const machepsi = std::numeric_limits<real_type>::epsilon();

  using std::abs;
  using std::pow;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of polinomial
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // (op[0] + op[1]/x.... + op[n]/x^n)*x^n
  //
  real_type
  evalPoly(
    real_type const op[],
    integer   const Degree,
    real_type const x
  ) {
    if ( std::abs(x) > 1 ) {
      real_type res(op[Degree]);
      real_type xn(1);
      for ( integer i{1}; i <= Degree; ++i ) {
        res /= x;
        res += op[Degree-i];
        xn  *= x;
      }
      res *= xn;
      return res;
    }
    real_type res(op[0]);
    for ( integer i{1}; i <= Degree; ++i ) {
      res *= x;
      res += op[i];
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of Newton step
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // (op[0] + op[1]/x.... + op[n]/x^n)*x^n
  //
  bool
  NewtonStep(
    real_type const op[],
    integer   const Degree,
    real_type     & x
  ) {
    real_type p, dp;
    if ( std::abs(x) > 1 ) {
      real_type xn{1};
      p = op[Degree];
      for ( integer i{1}; i <= Degree; ++i ) {
        p  /= x;
        p  += op[Degree-i];
        xn *= x;
      }
      p *= xn;

      dp = op[Degree];
      xn = 1;
      for ( integer i = 2; i <= Degree; ++i ) {
        dp /= x;
        dp += (Degree-i)*op[Degree-i];
        xn *= x;
      }
      dp *= xn;
    } else {
      dp = Degree*op[0];
      p  = op[0];
      for ( integer i = 1; i < Degree; ++i ) {
        p  *= x;
        p  += op[i];
        dp *= x;
        dp += (Degree-i)*op[i];
      }
      p *= x;
      p += op[Degree];
    }
    x -= p/dp;
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of polinomial and its derivative
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // (op[0] + op[1]/x.... + op[n]/x^n)*x^n
  //
  void
  evalPolyDPoly(
    real_type const op[],
    integer   const Degree,
    real_type const x,
    real_type     & p,
    real_type     & dp
  ) {
    if ( std::abs(x) > 1 ) {
      real_type const y{ real_type(1)/x };
      real_type const * pc{ op+Degree };

      p  = *pc--;
      dp = 0;
      for ( integer i{1}; i <= Degree; ++i ) {
        dp = dp*y + p;
        p  = p*y + *pc--;
      }

      real_type xn{1};
      for ( integer i{0}; i < Degree; ++i ) xn *= x;

      dp = static_cast<real_type>(Degree)*p*(xn/x) - dp*(xn/(x*x));
      p *= xn;
    } else {
      p  = op[0];
      dp = op[0];
      p  = p*x + op[1];
      for ( integer i = 2; i <= Degree; ++i ) {
        dp = dp*x + p;   // c_i = c_{i-1}*x + b_{i-1}
        p  = p*x + op[i];
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  complex_type
  evalPolyC(
    real_type    const op[],
    integer      const Degree,
    complex_type const x
  ) {
    if ( std::abs(x) > 1 ) {
      complex_type res(op[Degree]);
      complex_type xn(1,0);
      for ( integer i{1}; i <= Degree; ++i ) {
        res /= x;
        res += op[Degree-i];
        xn  *= x;
      }
      res *= xn;
      return res;
    }
    complex_type res(op[0]);
    for ( integer i{1}; i <= Degree; ++i ) {
      res *= x;
      res += op[i];
    }
    return res;
  }

}

#endif

// EOF: PolynomialRoots-Utils.cc
