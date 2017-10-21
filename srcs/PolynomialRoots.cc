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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "PolynomialRoots.hh"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

namespace PolynomialRoots {

  using std::pair ;
  using std::abs ;
  using std::pow ;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of polinomial
  // p0 + p1*x + p2*x^2 + ... + pn*x^n
  //
  // (p0/x^n + p1/x^(n-1) + p2/(x^(n-2) + ... + pn)*x^n
  //
  valueType
  evalPoly( valueType const op[],
            indexType       Degree,
            valueType       x,
            bool            reverse ) {
    valueType res(0) ;
    valueType xabs = std::abs(x) ;
    if ( xabs > 1 ) { x = valueType(1)/x ; reverse = !reverse ; }
    if ( reverse ) {
      for ( indexType i = 0 ; i <= Degree ; ++i ) res = res*x + op[i] ;
    } else {
      for ( indexType i = 0 ; i <= Degree ; ++i ) res = res*x + op[Degree-i] ;
    }
    return res*pow(x,Degree) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  std::complex<valueType>
  evalPolyC( valueType const         op[],
             indexType               Degree,
             std::complex<valueType> x,
             bool                    reverse ) {
    std::complex<valueType> res(0,0) ;
    valueType xabs = std::norm(x) ;
    if ( xabs > 1 ) { x = valueType(1)/x ; reverse = !reverse ; }
    if ( reverse ) {
      for ( indexType i = 0 ; i <= Degree ; ++i ) res = res*x + op[i] ;
    } else {
      for ( indexType i = 0 ; i <= Degree ; ++i ) res = res*x + op[Degree-i] ;
    }
    return res*pow(x,Degree) ;
  }

  //============================================================================

  /*
  ..  scale roots by pair.second, after scaling the polinomial has the form
  ..
  ..  p[0] + p[1]*x + ... + p[n]*x^n
  ..
  ..  pair.first is the index such that p[pair.first] = +-1
  ..
  */

  static
  pair<indexType,valueType>
  scalePolynomial( indexType       n, // degree
                   valueType const p[],
                   valueType       ps[] ) {
    indexType i_max = n ;
    valueType an    = p[n] ;
    valueType scale = -1 ;
    indexType i = n ;
    while ( --i >= 0 ) {
      ps[i] = p[i]/an ;
      valueType scale1 = pow( abs(ps[i]), 1.0/(n-i) ) ;
      if ( scale1 > scale ) { scale = scale1 ; i_max = i ; }
    }
    // scale coeffs
    pair<indexType,valueType> res(i_max,scale) ;
    valueType s = scale ;
    for ( i = 1 ; i <= n ; ++i, s *= scale ) ps[n-i] /= s ;
    ps[n] = 1 ;
    return res ;
  }

  //============================================================================

  /*
  .. divide a(x)  by (x-r) so that a(x) = (x-r) * b(x)
  */

  static
  void
  deflatePolynomial( indexType       n, // degree
                     valueType const a[],
                     valueType       r,
                     valueType       b[] ) {
    // crossover index for forward/backward deflation
    // G. Peters and J. H. Wilkinson.
    // Practical problems arising in the solution of polynomial equations.
    // J. Inst. Math. Appl. 8 (1971), 16–35.
    indexType i_cross = 0 ;
    valueType v_cross = abs(a[0]) ;
    valueType r1 = r ;
    for ( indexType i = 1 ; i < n ; ++i, r1 *= r ) {
      valueType v_cross1 = abs(a[i]*r1) ;
      if ( v_cross1 > v_cross )
        { v_cross = v_cross1 ; i_cross = i ; }
    }
    b[n-1] = 1 ;
    if ( a[n] == 1 ) {
      if ( i_cross > 0 ) {
        b[0] = -a[0] / r ;
        for ( indexType j = 1 ; j < i_cross ; ++j )
          b[j] = (a[j]-b[j-1]) / r ;
      }
      for ( indexType j = n-2 ; j >= i_cross ; --j )
        b[j] = a[j+1]+r*b[j+1] ;
    } else {
      valueType an = a[n] ;
      if ( i_cross > 0 ) {
        b[0] = -(a[0]/an) / r ;
        for ( indexType j = 1 ; j < i_cross ; ++j )
          b[j] = (a[j]/an-b[j-1]) / r ;
      }
      for ( indexType j = n-2 ; j >= i_cross ; --j )
        b[j] = a[j+1]/an+r*b[j+1] ;
    }
  }

}

/*
extern "C"
void
quartic_solver_( double op[], int & degree, double zeror[], double zeroi[] ) {
  double r1, r2, r3, r4 ;
  int nr, nc ;
  PolynomialRoots::solveQuartic( op[0], op[1], op[2], op[3], op[4], r1, r2, r3, r4, nr, nc ) ;
  switch ( nr ) {
  case 4:
    zeror[0] = r1 ; zeroi[0] = 0 ;
    zeror[1] = r2 ; zeroi[1] = 0 ;
    zeror[2] = r3 ; zeroi[2] = 0 ;
    zeror[3] = r4 ; zeroi[3] = 0 ;
    break ;
  case 2:
    zeror[0] = r1 ; zeroi[0] = r2 ;
    zeror[1] = r1 ; zeroi[1] = -r2 ;
    zeror[2] = r3 ; zeroi[2] = r4 ;
    zeror[3] = r3 ; zeroi[3] = -r4 ;
    break ;
  case 0:
    zeror[0] = r1 ; zeroi[0] = 0 ;
    zeror[1] = r2 ; zeroi[1] = 0 ;
    zeror[2] = r3 ; zeroi[2] = r4 ;
    zeror[3] = r3 ; zeroi[3] = -r4 ;
    break ;
  }
}
*/

// EOF: PolynomialRoots.cc
