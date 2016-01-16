/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2013                                                      |
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

  static valueType const epsi = 10*std::numeric_limits<valueType>::epsilon() ;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*
  ..  scale roots by pair.second, after scaling the polinomial has the form
  ..
  ..  x^n + p[1]*x^{n-1} + ... + p[n-1]*x + p[n]
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
    ps[n] = 1 ;
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
    return res ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
    i_cross = 0 ;
    b[n-1] = 1 ;
    if ( a[n] == 1 ) {
      b[0] = -a[0] / r ;
      for ( indexType j = 1 ; j < i_cross ; ++j )
        b[j] = (a[j]-b[j-1]) / r ;
      for ( indexType j = n-2 ; j >= i_cross ; --j )
        b[j] = a[j+1]+r*b[j+1] ;
    } else {
      valueType an = a[n] ;
      b[0] = -(a[0]/an) / r ;
      for ( indexType j = 1 ; j < i_cross ; ++j )
        b[j] = (a[j]/an-b[j-1]) / r ;
      for ( indexType j = n-2 ; j >= i_cross ; --j )
        b[j] = a[j+1]/an+r*b[j+1] ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  int
  zeroCubicByNewton( valueType const a[], valueType & x ) {
    indexType iter = 0 ;
    do {
      valueType f   = a[0]+x*(a[1]+x*(a[2]+x*a[3])) ;
      valueType df  = a[1]+x*(2*a[2]+x*(3*a[3])) ;
      valueType ddf = 2*a[2]+x*(6*a[3]) ;
      valueType dx = f/df ;
      //dx = f/(df-0.5*ddf*dx) ;
      dx /= 1-0.5*dx*ddf/df ;
      x -= dx ;
      if ( abs(dx) < epsi*abs(x) || abs(f) < epsi ) return iter ;
    } while ( ++iter < 20 ) ;
    return -1 ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*\
   *  Calculate the zeros of the quadratic a*z^2 + b*z + c.
   *  The quadratic formula, modified to avoid overflow, is used
   *  to find the larger zero if the zeros are real and both
   *  are complex. The smaller real zero is found directly from
   *  the product of the zeros c/a.
  \*/
  bool
  solveQuadratic( valueType   a,
                  valueType   b,
                  valueType   c,
                  valueType & r1,
                  valueType & r2 ) {
    bool real_root = true ;
    r1 = r2 = 0 ;
    if ( a == 0 ) { // less than two roots
      r1 = b == 0 ? 0 : -c/b ;
    } else if ( c == 0 ) { // one real root, one zero root
      r1 = -b/a;
    } else {
      // Compute discriminant avoiding overflow.
      b /= 2 ; // b now b/2
      valueType abs_b = abs(b) ;
      valueType abs_c = abs(c) ;
      valueType e, d ;
      if ( abs_b < abs_c ) {
        e = c < 0 ? -a : a ;
        e = b*(b/abs_c) - e ;
        d = sqrt(abs(e))*sqrt(abs_c);
      } else {
        e = 1 - (a/b)*(c/b);
        d = sqrt(abs(e))*abs_b ;
      }
      real_root = e >= 0 ;
      if ( real_root ) { // complex conjugate zeros
        if ( b >= 0 ) d = -d ; // real zeros
        r1 = (d-b)/a;
        if ( r1 != 0 ) r2 = (c/r1)/a;
      } else {
        r1 = -b/a ; // real part
        r2 = std::abs(d/a) ; // immaginary part
      }
    }
    return real_root ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  static
  inline
  valueType
  guess1( valueType const a[4] ) {
    valueType const p =  1.09574 ;
    valueType const q = -3.239E-1 ;
    valueType const r = -3.239E-1 ;
    valueType const s =  9.57439E-2 ;
    return p+q*a[1]+r*a[2]+s*a[1]*a[2] ;
  }

  static
  inline
  valueType
  guess2( valueType const a[4] ) {
    valueType const p = -1.09574 ;
    valueType const q =  3.239E-1 ;
    valueType const r = -3.239E-1 ;
    valueType const s =  9.57439E-2 ;
    return p+q*a[1]+r*a[2]+s*a[1]*a[2] ;
  }

  static
  inline
  valueType
  guess3( valueType const a[4] ) {
    valueType const p =  1.14413    ;
    valueType const q = -2.75509E-1 ;
    valueType const r = -4.45578E-1 ;
    valueType const s = -2.59342E-2 ;
    valueType t = a[2]/3 ;
    if ( a[0] < t*(2*t*t-1) ) return  p+q*a[0]+r*a[2]+s*a[0]*a[2] ;
    else                      return -p+q*a[0]+r*a[2]-s*a[0]*a[2] ;
  }

  static
  inline
  valueType
  guess4( valueType const a[4] ) {
    valueType const q = -7.71845E-1 ;
    valueType const s = -2.28155E-1 ;
    if ( a[0] > 0 ) return (q+s*a[2])*a[0] ;
    else            return (q-s*a[2])*a[0] ;
  }

  static
  inline
  valueType
  guess5( valueType const a[4], bool & use_shifted ) {
    valueType p, q, r, s ;
    use_shifted = abs(a[1]-1./3.) <= 0.01 && abs(a[0]+1./27.) ;
    if ( a[1] <= 1./3. ) {
      if ( a[0] < -a[1]/3+2./27. ) {
        p =  8.78558E-1 ;
        q = -5.71888E-1 ;
        r = -7.11154E-1 ;
        s = -3.22313E-1 ;
      } else {
        p = -1.92823E-1 ;
        q = -5.66324E-1 ;
        r = +5.05734E-1 ;
        s = -2.64881E-1 ;
      }
    } else {
      if ( a[0] < -a[1]/3+2./27. ) {
        p = 1.19748 ;
        q = -2.83772E-1 ;
        r = -8.37476E-1 ;
        s = -3.56228E-1 ;
      } else {
        p = -3.45219E-1 ;
        q = -4.01231E-1 ;
        r =  2.07216E-1 ;
        s = -4.45532E-3 ;
      }
    }
    return p+q*a[0]+r*a[1]+s*a[0]*a[1] ;
  }

  static
  inline
  valueType guess6( valueType const a[4], bool & use_shifted ) {
    valueType p, q, r, s ;
    use_shifted = abs(a[1]-1./3.) <= 0.01 && abs(a[0]-1./27.) ;
    if ( a[1] <= 1./3. ) {
      if ( a[0] > a[1]/3-2./27. ) {
        p = -8.78558E-1 ;
        q = -5.71888E-1 ;
        r =  7.11154E-1 ;
        s = -3.22313E-1 ;
      } else {
        p =  1.92823E-1 ;
        q = -5.66324E-1 ;
        r = -5.05734E-1 ;
        s = -2.64881E-1 ;
      }
    } else {
      if ( a[0] > a[1]/3-2./27. ) {
        p = -1.19748 ;
        q = -2.83772E-1 ;
        r =  8.37476E-1 ;
        s = -3.56228E-1 ;
      } else {
        p =  3.45219E-1 ;
        q = -4.01231E-1 ;
        r = -2.07216E-1 ;
        s = -4.45532E-3 ;
      }
    }
    return p+q*a[0]+r*a[1]+s*a[0]*a[1] ;
  }

  /*\
  ... Calculate the zeros of the cubic a*z^3 + b*z^2 + c*z + d.
  ...
  ... N. FLOCKE, Flash Center for Computational Science, University of Chicago
  ... Algorithm 954: An Accurate and Efficient Cubic and Quartic Equation Solver
  ... for Physical Applications
  ... ACM Transactions on Mathematical Software, Vol. 41, No. 4, 2015.
  ... DOI: http://dx.doi.org/10.1145/2699468
  \*/
  bool
  solveCubic( valueType   A,
              valueType   B,
              valueType   C,
              valueType   D,
              valueType & r1,
              valueType & r2,
              valueType & r3 ) {

    valueType a[4], p_unscaled[4] = { D, C, B, A } ;

    pair<indexType,valueType> icase = scalePolynomial( 3, p_unscaled, a ) ;

    // Class1: a[0] = −1, −1 <= a[1],a[2] <= +1
    // Class2: a[0] = +1, −1 <= a[1],a[2] <= +1
    // Class3: a[1] = −1, −1 <= a[0],a[2] <= +1
    // Class4: a[1] = +1, −1 <= a[0],a[2] <= +1
    // Class5: a[2] = −1, −1 <= a[0],a[1] <= +1
    // Class6: a[2] = +1, −1 <= a[0],a[1] <= +1
    indexType iclass = -1 ;
    switch ( icase.first ) {
      case 0: iclass = a[0] > 0 ? 2 : 1 ; break ;
      case 1: iclass = a[1] > 0 ? 4 : 3 ; break ;
      case 2: iclass = a[2] > 0 ? 6 : 5 ; break ;
    }
    bool use_shifted = false ;
    switch ( iclass ) {
      case 1: r1 = guess1(a) ; break ;
      case 2: r1 = guess2(a) ; break ;
      case 3: r1 = guess3(a) ; break ;
      case 4: r1 = guess4(a) ; break ;
      case 5: r1 = guess5(a,use_shifted) ; break ;
      case 6: r1 = guess6(a,use_shifted) ; break ;
    }
    indexType iter ;
    if ( use_shifted ) {
      valueType const A = a[1]-1./3. ;
      if ( iclass == 5 ) {
        // y^3 + A * y + (B+A/3), y = x-1/3
        // B = a[0]+1./27. ;
        valueType const p[4] = { a[0]+a[1]/3-2./27, A, 0, 1 } ;
        r1 -= 1./3. ;
        iter = zeroCubicByNewton( p, r1 ) ;
        r1 += 1./3. ;
      } else {
        // y^3 + A * y + (B-A/3), y = x+1/3
        // B = a[0]-1./27. ;
        valueType const p[4] = { a[0]-a[1]/3+2./27, A, 0, 1 } ;
        r1 += 1./3. ;
        iter = zeroCubicByNewton( p, r1 ) ;
        r1 -= 1./3. ;
      }
    } else {
      iter = zeroCubicByNewton( a, r1 ) ;
    }
    // scale
    r1 *= icase.second ;
    
    // una extra correzione con Newton dopo riscalatura
    valueType f  = p_unscaled[0]+r1*(p_unscaled[1]+r1*(p_unscaled[2]+r1*p_unscaled[3])) ;
    valueType df = p_unscaled[1]+r1*(2*p_unscaled[2]+r1*(3*p_unscaled[3])) ;
    r1 -= f/df ;

    // deflate
    valueType b[3] ;
    deflatePolynomial( 3, p_unscaled, r1, b ) ;
    bool real_roots = solveQuadratic( b[2], b[1], b[0], r2, r3 ) ;
    if ( real_roots ) { // if real roots sort it!
      if ( r1 > r2 ) std::swap(r1,r2) ;
      if ( r2 > r3 ) std::swap(r2,r3) ;
      if ( r1 > r2 ) std::swap(r1,r2) ;
    }

    return real_roots ;
  }
}

// EOF: PolynomialRoots.cc
