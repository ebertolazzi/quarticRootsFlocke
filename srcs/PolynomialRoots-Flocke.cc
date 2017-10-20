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

  static valueType const machepsi      = std::numeric_limits<valueType>::epsilon() ;
  static valueType const third         = 1./3. ;
  static valueType const one27th       = 1./27. ;
  static valueType const two27th       = 2./27. ;
  static int       const bitsValueType = std::numeric_limits<valueType>::digits ;
  static valueType const splitFactor   = (long(1)<<(bitsValueType-1))+1 ;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // a + b = x + err
  static
  inline
  void
  TwoSum( valueType   a,
          valueType   b,
          valueType & x,
          valueType & y ) {
    x = a+b ;
    valueType z = x-a ;
    y = (a-(x-z))+(b-z) ;
    //if ( abs(a) < abs(b) ) std::swap(a,b) ;
    //x = a+b ;
    //valueType z = x - a;
    //y = b-z ;
  }

  // a = x + y
  static
  inline
  void
  Split( valueType a, valueType & x, valueType & y ) {
    valueType c = splitFactor*a ;
    x = c-(c-a) ;
    y = a-x ;
  }

  // a * b = x + err
  static
  inline
  void
  TwoProduct( valueType   a,
              valueType   b,
              valueType & x,
              valueType & y ) {
    valueType a1, a2, b1, b2 ;
    Split( a, a1, a2 ) ;
    Split( b, b1, b2 ) ;
    x = a*b ;
    y = a2*b2-(((x-a1*b1)-a2*b1)-a1*b2) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of polinomial
  // p0 + p1*x + p2*x^2 + ... + pn*x^n
  //
  // (p0/x^n + p1/x^(n-1) + p2/(x^(n-2) + ... + pn)*x^n
  //
  valueType
  CompHorner( valueType const p[],
              indexType       Degree,
              valueType       x,
              bool            reverse ) {

    valueType xabs = std::abs(x) ;
    if ( xabs > 1 ) { x = valueType(1)/x ; reverse = !reverse ; }
    indexType ii0 = reverse ? 0 : Degree ;
    valueType res(p[ii0]) ;
    valueType c = 0 ;
    for ( indexType i = 1 ; i <= Degree ; ++i ) {
      indexType ii = reverse ? i : Degree-i ;
      valueType tmp, pi, sigma ;
      TwoProduct( res, x, tmp, pi ) ;
      TwoSum( tmp, p[ii], res, sigma ) ;
      c = c * x + (pi+sigma) ;
    }
    res += c ;
    if ( xabs > 1 ) res *= pow(x,Degree) ;
    return res ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // x^3 + A x^2 + B x + C
  static
  inline
  void
  scaleCubicMonicPolynomial( valueType   A,
                             valueType   B,
                             valueType   C,
                             valueType & AS,
                             valueType & BS,
                             valueType & CS,
                             indexType & i_case,
                             valueType & scale ) {

    valueType a = abs(A) ;
    valueType b = sqrt(abs(B)) ;
    valueType c = cbrt(abs(C)) ;

    if ( a < b ) {
      if ( b < c ) i_case = 0 ; // a < b < c --> c MAX
      else         i_case = 1 ; // a < b and c <= b --> b MAX
    } else {
      if ( a < c ) i_case = 0 ; // b <= a < c --> c MAX
      else         i_case = 2 ; // b <= a  and c <= a --> a MAX
    }

    switch ( i_case ) {
      case 0:
        scale = c ;
        AS    = A/c ;
        BS    = (B/c)/c ;
        CS    = C > 0 ? 1 : -1 ;
      break ;
      case 1:
        scale = b ;
        AS    = A/b ;
        BS    = B > 0 ? 1 : -1 ;
        CS    = ((C/b)/b)/b ;
      break ;
      case 2:
        scale = a ;
        AS    = A > 0 ? 1 : -1 ;
        BS    = (B/a)/a ;
        CS    = ((C/a)/a)/a ;
      break ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // x^4 + A x^3 + B x^2 + C x + D
  static
  inline
  void
  scaleQuarticMonicPolynomial( valueType   A,
                               valueType   B,
                               valueType   C,
                               valueType   D,
                               valueType & AS,
                               valueType & BS,
                               valueType & CS,
                               valueType & DS,
                               indexType & i_case,
                               valueType & scale ) {

    valueType a = abs(A) ;
    valueType b = sqrt(abs(B)) ;
    valueType c = cbrt(abs(C)) ;
    valueType d = sqrt(sqrt(abs(D))) ;

    if ( a < b ) {
      if ( b < c ) {
        if ( c < d ) i_case = 0 ; // a < b < c < d        --> d MAX
        else         i_case = 1 ; // a < b < c and d <= c --> c MAX
      } else {
        if ( b < d ) i_case = 0 ; // a < b and c <= b < d --> d MAX
        else         i_case = 2 ; // a < b and c <= b and d <= b --> b MAX
      }
    } else {
      if ( a < c ) {
        if ( c < d ) i_case = 0 ; // b <= a < c < d        --> d MAX
        else         i_case = 1 ; // b <= a < c and d <= c --> c MAX
      } else {
        if ( a < d ) i_case = 0 ; // b <= a and c <= a and a < d --> d MAX
        else         i_case = 3 ; // b <= a and c <= a and d <= a --> a MAX
      }
    }

    switch ( i_case ) {
      case 0:
        scale = d ;
        AS    = A/d ;
        BS    = (B/d)/d ;
        CS    = ((C/d)/d)/d ;
        DS    = D > 0 ? 1 : -1 ;
      break ;
      case 1:
        scale = c ;
        AS    = A/c ;
        BS    = (B/c)/c ;
        CS    = C > 0 ? 1 : -1 ;
        DS    = (((D/c)/c)/c)/c ;
      break ;
      case 2:
        scale = b ;
        AS    = A/b ;
        BS    = B > 0 ? 1 : -1 ;
        CS    = ((C/b)/b)/b ;
        DS    = (((D/b)/b)/b)/b ;
      break ;
      case 3:
        scale = a ;
        AS    = A > 0 ? 1 : -1 ;
        BS    = (B/a)/a ;
        CS    = ((C/a)/a)/a ;
        DS    = (((D/a)/a)/a)/a ;
      break ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // a3*x^3 + a2*x^2 + a1*x + a0 = (x-r)*(a3*x^2+b1*x+b0)
  static
  void
  deflateCubicPolynomial( valueType   a3,
                          valueType   a2,
                          valueType   a1,
                          valueType   a0,
                          valueType   r,
                          valueType & b1,
                          valueType & b0 ) {
    indexType i_cross  = 0 ;
    valueType r2       = r*r ;
    valueType v_cross  = abs(a0) ;
    valueType v_cross1 = abs(a1*r) ;
    if ( v_cross1 > v_cross ) { v_cross = v_cross1 ; i_cross = 1 ; }
    v_cross1 = abs(a2*r2) ;
    if ( v_cross1 > v_cross ) { v_cross = v_cross1 ; i_cross = 2 ; }
    v_cross1 = abs(a3*r*r2) ;
    if ( v_cross1 > v_cross ) i_cross = 3 ;
    switch ( i_cross ) {
      case 0: b1 = a2+a3*r ; b0 = a1+r*b1 ; break;
      case 1: b1 = a2+a3*r ; b0 = -a0/r   ; break;
      case 2:
      case 3: b0 = -a0/r ; b1 = (b0-a1)/r ; break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = (x-r)*(a4*x^3+b2*x^2+b1*x+b0)
  static
  void
  deflateQuarticPolynomial( valueType   a4,
                            valueType   a3,
                            valueType   a2,
                            valueType   a1,
                            valueType   a0,
                            valueType   r,
                            valueType & b2,
                            valueType & b1,
                            valueType & b0 ) {
    indexType i_cross  = 0 ;
    valueType r2       = r*r ;
    valueType v_cross  = abs(a0) ;
    valueType v_cross1 = abs(a1*r) ;
    if ( v_cross1 > v_cross ) { v_cross = v_cross1 ; i_cross = 1 ; }
    v_cross1 = abs(a2*r2) ;
    if ( v_cross1 > v_cross ) { v_cross = v_cross1 ; i_cross = 2 ; }
    v_cross1 = abs(a3*r*r2) ;
    if ( v_cross1 > v_cross ) { v_cross = v_cross1 ; i_cross = 3 ; }
    v_cross1 = abs(a4*r2*r2) ;
    if ( v_cross1 > v_cross ) i_cross = 4 ;
    switch ( i_cross ) {
      case 0: b2 = a3+a4*r ; b1 = a2+r*b2 ; b0 = a1+r*b1   ; break;
      case 1: b2 = a3+a4*r ; b1 = a2+r*b2 ; b0 = -a0/r     ; break;
      case 2: b2 = a3+a4*r ; b0 = -a0/r   ; b1 = (b0-a1)/r ; break;
      case 3:
      case 4: b0 = -a0/r ; b1 = (b0-a1)/r ; b2 = (b1-a2)/r ; break;
    }
  }

  // x^3 + a*x^2 + b*x + c
  static
  inline
  valueType
  evalMonicCubic( valueType x,
                  valueType a,
                  valueType b,
                  valueType c ) {
    return ((x+a)*x+b)*x+c ;
  }

  static
  inline
  void
  evalMonicCubic( valueType   x,
                  valueType   a,
                  valueType   b,
                  valueType   c,
                  valueType & p,
                  valueType & dp ) {
    p  = x + a ;
    dp = x + p ;
    p  = p  * x + b ;
    dp = dp * x + p ;
    p  = p  * x + c ;
  }

  // x^4 + a*x^3 + b*x^2 + c*x + d
  static
  inline
  valueType
  evalMonicQuartic( valueType x,
                    valueType a,
                    valueType b,
                    valueType c,
                    valueType d ) {
    return (((x+a)*x+b)*x+c)*x+d ;
  }

  static
  inline
  void
  evalMonicQuartic( valueType   x,
                    valueType   a,
                    valueType   b,
                    valueType   c,
                    valueType   d,
                    valueType & p,
                    valueType & dp ) {
    p  = x + a ;
    dp = x + p ;
    p  = p  * x + b ;
    dp = dp * x + p ;
    p  = p  * x + c ;
    dp = dp * x + p ;
    p  = p  * x + d ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static
  inline
  valueType
  evalHexic( valueType x,
             valueType q3,
             valueType q2,
             valueType q1,
             valueType q0 ) {
    valueType t1 = x + q3 ;
    valueType t2 = x + t1 ;
    valueType t3 = x + t2 ;
    t1 = t1 * x + q2 ;
    t2 = t2 * x + t1  ;
    t1 = t1 * x + q1 ;
    valueType Q    = t1 * x + q0 ; // Q(x)
    valueType dQ   = t2 * x + t1 ; // Q'(x)
    valueType ddQ  = t3 * x + t2 ; // Q''(x) / 2
    valueType dddQ = x + t3 ;      // Q'''(x) / 6
    return  Q * dddQ * dddQ - dQ * ddQ * dddQ + dQ * dQ ; // H(x), usually < 0
  }

  static
  inline
  void
  evalHexic( valueType   x,
             valueType   q3,
             valueType   q2,
             valueType   q1,
             valueType   q0,
             valueType & p,
             valueType & dp ) {
    valueType t1 = x + q3 ;
    valueType t2 = x + t1 ;
    valueType t3 = x + t2 ;
    t1 = t1 * x + q2 ;
    t2 = t2 * x + t1  ;
    t1 = t1 * x + q1 ;
    valueType Q    = t1 * x + q0 ; // Q(x)
    valueType dQ   = t2 * x + t1 ; // Q'(x)
    valueType ddQ  = t3 * x + t2 ; // Q''(x) / 2
    valueType dddQ = x + t3 ;      // Q'''(x) / 6
    p  = Q * dddQ * dddQ - dQ * ddQ * dddQ + dQ * dQ ; // H(x), usually < 0
    dp = 2 * dddQ * (4 * Q - dQ * dddQ - ddQ * ddQ) ; // H'(x)
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Translate to C from Polynomial234RootSolvers
  static
  indexType
  zeroCubicByNewtonBisection( valueType const a,
                              valueType const b,
                              valueType const c,
                              valueType     & x ) {

    valueType p, dp ;
    evalMonicCubic( x, a, b, c, p, dp ) ;
    valueType t = p ; // save p(x) for sign comparison
    x -= p/dp ; // 1st improved root

    indexType iter      = 1 ;
    indexType oscillate = 0 ;
    bool      bisection = false ;
    bool      converged = false ;
    valueType s(0), u(0) ; // to mute warning
    while ( ! (converged||bisection) ) {
      ++iter ;
      evalMonicCubic( x, a, b, c, p, dp ) ;
      if ( p*t < 0 ) { // does Newton start oscillating ?
        if ( p < 0 ) {
          ++oscillate ; // increment oscillation counter
          s = x ;       // save lower bisection bound
        } else {
          u = x ; // save upper bisection bound
        }
        t = p ; // save current p(x)
      }
      dp = p/dp ; // Newton correction
      x -= dp ;   // new Newton root
      bisection = oscillate > 2 ; // activate bisection
      converged = abs(dp) <= abs(x) * machepsi ; // Newton convergence indicator
    }
    if ( bisection ) {
      t = u - s ; // initial bisection interval
      while ( abs(t) > abs(x) * machepsi ) { // bisection iterates
        ++iter ;
        p = evalMonicCubic( x, a, b, c ) ;
        if ( p < 0 ) s = x ;
        else         u = x ; // keep bracket on root
        t = (u-s)/2 ; // new bisection interval
        x = s + t ;   // new bisection root
      }
    }
    return iter ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Translate to C from Polynomial234RootSolvers
  static
  indexType
  zeroQuarticByNewtonBisection( valueType const a,
                                valueType const b,
                                valueType const c,
                                valueType const d,
                                valueType     & x ) {

    valueType p, dp ;
    evalMonicQuartic( x, a, b, c, d, p, dp ) ;
    valueType t = p ; // save p(x) for sign comparison
    x -= p/dp ; // 1st improved root

    indexType iter      = 1 ;
    indexType oscillate = 0 ;
    bool      bisection = false ;
    bool      converged = false ;
    valueType s(0), u(0) ; // to mute warning
    while ( ! (converged||bisection) ) {
      ++iter ;
      evalMonicQuartic( x, a, b, c, d, p, dp ) ;
      if ( p*t < 0 ) { // does Newton start oscillating ?
        if ( p < 0 ) {
          ++oscillate ; // increment oscillation counter
          s = x ;       // save lower bisection bound
        } else {
          u = x ; // save upper bisection bound
        }
        t = p ; // save current p(x)
      }
      dp = p/dp ; // Newton correction
      x -= dp ; // new Newton root
      bisection = oscillate > 2 ; // activate bisection
      converged = abs(dp) <= abs(x) * machepsi ; // Newton convergence indicator
    }
    if ( bisection ) {
      t = u - s ; // initial bisection interval
      while ( abs(t) > abs(x) * machepsi ) { // bisection iterates
        ++iter ;
        p = evalMonicQuartic( x, a, b, c, d ) ;
        if ( p < 0 ) s = x ;
        else         u = x ; // keep bracket on root
        t = (u-s)/2 ; // new bisection interval
        x = s + t ;   // new bisection root
      }
    }
    return iter ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Translate to C from Polynomial234RootSolvers
  static
  indexType
  zeroHexicByNewtonBisection( valueType   q3,
                              valueType   q2,
                              valueType   q1,
                              valueType   q0,
                              valueType & x ) {
    valueType p, dp ;
    evalHexic( x, q3, q2, q1, q0, p, dp ) ;
    valueType t = p ; // save p(x) for sign comparison
    x -= p/dp ; // 1st improved root

    indexType iter      = 1 ;
    indexType oscillate = 0 ;
    bool      bisection = false ;
    bool      converged = false ;
    valueType s(0), u(0) ; // to mute warning
    while ( ! (converged||bisection) ) {
      ++iter ;
      evalHexic( x, q3, q2, q1, q0, p, dp ) ;
      if ( p*t < 0 ) { // does Newton start oscillating ?
        if ( p < 0 ) {
          ++oscillate ; // increment oscillation counter
          s = x ;       // save lower bisection bound
        } else {
          u = x ; // save upper bisection bound
        }
        t = p ; // save current p(x)
      }
      dp = p/dp ; // Newton correction
      x -= dp ; // new Newton root
      bisection = oscillate > 2 ; // activate bisection
      converged = abs(dp) <= abs(x) * machepsi ; // Newton convergence indicator
    }
    if ( bisection ) {
      t = u - s ; // initial bisection interval
      while ( abs(t) > abs(x) * machepsi ) { // bisection iterates
        ++iter ;
        p = evalHexic( x, q3, q2, q1, q0 ) ;
        if ( p < 0 ) s = x ;
        else         u = x ; // keep bracket on root
        t = (u-s)/2 ; // new bisection interval
        x = s + t ;   // new bisection root
      }
    }
    return iter ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if 1
  /*\
   *  Calculate the zeros of the quadratic a*z^2 + b*z + c.
   *  The quadratic formula, modified to avoid overflow, is used
   *  to find the larger zero if the zeros are real and both
   *  are complex. The smaller real zero is found directly from
   *  the product of the zeros c/a.
  \*/
  void
  solveQuadratic( valueType   a,
                  valueType   b,
                  valueType   c,
                  valueType & r1,
                  valueType & r2,
                  indexType & nr,
                  indexType & nc ) {
    r1 = r2 = 0 ;
    nr = nc = 0 ;
    if ( a == 0 ) { // less than two roots b*z + c = 0
      if ( b != 0 ) { nr = 1 ; r1 = -c/b ; }
    } else if ( c == 0 ) { // a*z^2 + b*z  = 0
      nr = 2 ;
      r1 = -b/a ;
      if ( r1 > 0 ) std::swap(r1,r2) ;
    } else { // Compute discriminant avoiding overflow.
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
      bool real_root = e >= 0 ;
      if ( real_root ) {       // complex conjugate zeros
        if ( b >= 0 ) d = -d ; // real zeros
        r1 = (d-b)/a;
        if ( r1 != 0 ) {
          r2 = (c/r1)/a ;
          if ( r1 > r2 ) std::swap(r1,r2) ; // order roots
        }
        nr = 2 ;
      } else {
        r1 = -b/a ;          // real part
        r2 = std::abs(d/a) ; // immaginary part
        nc = 2 ;
      }
    }
  }
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  static
  inline
  valueType
  guess1( valueType const a[3] ) {
    valueType const p =  1.09574 ;
    valueType const q = -3.239E-1 ;
    valueType const r = -3.239E-1 ;
    valueType const s =  9.57439E-2 ;
    return p+q*a[1]+r*a[2]+s*a[1]*a[2] ;
  }

  static
  inline
  valueType
  guess2( valueType const a[3] ) {
    valueType const p = -1.09574 ;
    valueType const q =  3.239E-1 ;
    valueType const r = -3.239E-1 ;
    valueType const s =  9.57439E-2 ;
    return p+q*a[1]+r*a[2]+s*a[1]*a[2] ;
  }

  static
  inline
  valueType
  guess3( valueType const a[3] ) {
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
  guess4( valueType const a[3] ) {
    valueType const q = -7.71845E-1 ;
    valueType const s = -2.28155E-1 ;
    if ( a[0] > 0 ) return (q+s*a[2])*a[0] ;
    else            return (q-s*a[2])*a[0] ;
  }

  static
  inline
  valueType
  guess5( valueType const a[3] ) {
    valueType p, q, r, s ;
    valueType tmp = two27th-a[1]/3 ;
    if ( a[1] <= third ) {
      if ( a[0] < tmp ) {
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
      if ( a[0] < tmp ) {
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
  valueType guess6( valueType const a[3] ) {
    valueType p, q, r, s ;
    valueType tmp = a[1]/3-two27th ;
    if ( a[1] <= third ) {
      if ( a[0] > tmp ) {
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
      if ( a[0] > tmp ) {
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
  ... Calculate the zeros of the cubic A*z^3 + B*z^2 + C*z + D.
  ...
  ... N. FLOCKE, Flash Center for Computational Science, University of Chicago
  ... Algorithm 954: An Accurate and Efficient Cubic and Quartic Equation Solver
  ... for Physical Applications
  ... ACM Transactions on Mathematical Software, Vol. 41, No. 4, 2015.
  ... DOI: http://dx.doi.org/10.1145/2699468
  \*/

  indexType
  solveCubic( valueType   A,
              valueType   B,
              valueType   C,
              valueType   D,
              valueType & r1,
              valueType & r2,
              valueType & r3,
              indexType & nr,
              indexType & nc ) {

    // special cases
    if ( A == 0 ) {
      solveQuadratic( B, C, D, r1, r2, nr, nc ) ;
      return 0 ;
    }
    if ( D == 0 ) {
      r1 = 0 ;
      solveQuadratic( A, B, C, r2, r3, nr, nc ) ;
      if ( nr == 1 ) { // caso degenere
        if ( r1 > r2 ) std::swap(r1,r2) ;
      } else if ( nr == 2 ) {
        if ( r1 > r2 ) std::swap(r1,r2) ;
        if ( r2 > r3 ) std::swap(r2,r3) ;
      }
      ++nr ;
      return 0 ;
    }

    valueType scale, a[3] ;
    indexType i_case ;
    scaleCubicMonicPolynomial( B/A, C/A, D/A, a[2], a[1], a[0], i_case, scale ) ;

    // Class1: a[0] = −1, −1 <= a[1],a[2] <= +1
    // Class2: a[0] = +1, −1 <= a[1],a[2] <= +1
    // Class3: a[1] = −1, −1 <= a[0],a[2] <= +1
    // Class4: a[1] = +1, −1 <= a[0],a[2] <= +1
    // Class5: a[2] = −1, −1 <= a[0],a[1] <= +1
    // Class6: a[2] = +1, −1 <= a[0],a[1] <= +1
    indexType iclass = -1 ;
    switch ( i_case ) {
      case 0: iclass = a[0] > 0 ? 2 : 1 ; break ;
      case 1: iclass = a[1] > 0 ? 4 : 3 ; break ;
      case 2: iclass = a[2] > 0 ? 6 : 5 ; break ;
    }
    bool use_shifted = false ;
    bool triple_root = false ;
    switch ( iclass ) {
      case 1: r1 = guess1(a) ; break ;
      case 2: r1 = guess2(a) ; break ;
      case 3: r1 = guess3(a) ; break ;
      case 4: r1 = guess4(a) ; break ;
      case 5:
        r2 = a[1]-third ;
        r3 = a[0]+one27th ;
        use_shifted = abs(r2) <= 0.01 && abs(r3) <= 0.01 ;
        triple_root = abs(r2) <= machepsi && abs(r3) <= machepsi ;
        r1 = guess5(a) ;
        break ;
      case 6:
        r2 = a[1]-third ;
        r3 = a[0]-one27th ;
        use_shifted = abs(r2) <= 0.01 && abs(r3) <= 0.01 ;
        triple_root = abs(r2) <= machepsi && abs(r3) <= machepsi ;
        r1 = guess6(a) ;
        break ;
    }
    indexType iter = 0 ;
    if ( triple_root ) {
      nr = 3 ;
      if ( iclass == 5 ) r1 = r2 = r3 = -third * scale ;
      else               r1 = r2 = r3 =  third * scale ;
      return iter ;
    } else if ( use_shifted ) {
      if ( iclass == 5 ) {
        // y^3 + A * y + (B+A/3), y = x-1/3
        // B = a[0]+1./27. ;
        r1 -= third ;
        r3 += third * r2 ;
        //if ( abs(r3) < machepsi ) r3 = 0 ;
        iter = zeroCubicByNewtonBisection( 0, r2, r3, r1 ) ;
        r1 += third ;
      } else {
        // y^3 + A * y + (B-A/3), y = x+1/3
        // B = a[0]-1./27. ;
        r1 += third ;
        r3 -= third * r2 ;
        //if ( abs(r3) < machepsi ) r3 = 0 ;
        iter = zeroCubicByNewtonBisection( 0, r2, r3, r1 ) ;
        r1 -= third ;
      }
    } else {
      iter = zeroCubicByNewtonBisection( a[2], a[1], a[0], r1 ) ;
    }
    // scale
    r1 *= scale ;
    
    valueType p  = ((A*r1+B)*r1+C)*r1+D ;
    valueType dp = ((4*A*r1+3*B)*r1+2*C)*r1 ;
    
    r1 -= p/dp ;
/*
    valueType const pp[]  = { A, B, C, D } ;
    valueType pH = CompHorner( pp, 3, r1, true ) ;
    std::cout << "pH = " << pH << "\n" ;

    // una extra correzione con Newton dopo riscalatura
    valueType const dpp[] = { 3*A, 2*B, C } ;
    for ( int k = 0 ; k < 10 ; ++k ) {
      valueType pH  = CompHorner( pp, 3, r1, true ) ;
      valueType dpH = (3*A*r1+2*B)*r1+C ;
      r1 -= pH/dpH ;
    }

    pH = CompHorner( pp, 3, r1, true ) ;
    std::cout << "pH = " << pH << "\n" ;
*/
    // deflate
    valueType b0, b1 ;
    deflateCubicPolynomial( A, B, C, D, r1, b1, b0 ) ;
    solveQuadratic( A, b1, b0, r2, r3, nr, nc ) ;
    if ( nr == 2 ) { // if real roots sort it!
      if ( r1 > r2 ) std::swap(r1,r2) ;
      if ( r2 > r3 ) std::swap(r2,r3) ;
    }
    ++nr ; // one more real root
    return iter ;
  }
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  indexType
  solveQuartic( valueType   A,
                valueType   B,
                valueType   C,
                valueType   D,
                valueType   E,
                valueType & r1,
                valueType & r2,
                valueType & r3,
                valueType & r4,
                indexType & nr,
                indexType & nc ) {
    
    // special cases
    if ( A == 0 ) {
      return solveCubic( B, C, D, E, r1, r2, r3, nr, nc ) ;
    }
    if ( E == 0 ) {
      r1 = 0 ;
      indexType iter = solveCubic( A, B, C, D, r2, r3, r4, nr, nc ) ;
      if ( nr == 1 ) { // caso degenere o regolare
        if ( r1 > r2 ) std::swap(r1,r2) ;
      } else if ( nr == 2 ) { // caso degenere
        if ( r1 > r2 ) std::swap(r1,r2) ;
        if ( r2 > r3 ) std::swap(r2,r3) ;
      } else if ( nr == 3 ) { // caso regolare, 3 radici reali
        if ( r1 > r2 ) std::swap(r1,r2) ;
        if ( r2 > r3 ) std::swap(r2,r3) ;
        if ( r3 > r4 ) std::swap(r3,r4) ;
      }
      ++nr ;
      return iter ;
    }
    if ( B == 0 && D == 0 ) { // biquadratic case
      // A x^4 + C x^2 + E
      valueType x, y ;
      solveQuadratic( A, C, E, x, y, nr, nc ) ;
      if ( nr == 2 ) { // real roots of quadratic are ordered x <= y
        if ( x >= 0 ) { // y >= 0
          nr = 4 ; nc = 0 ;
          x = sqrt(x) ; y = sqrt(y) ;
          r1 = -y ; r2 = -x ; r3 =  x ; r4 =  y ;
        } else if ( y >= 0 ) { // x < 0 && y >= 0
          nr = nc = 2 ;
          x = sqrt(-x) ; y = sqrt(y) ;
          r1 = -y ; r2 = y ; r3 = 0 ; r4 = x ; // (real,imaginary)
        } else { // x < 0 && y < 0
          nr = 0 ; nc = 4 ;
          x = sqrt(-x) ; y = sqrt(-y) ;
          r1 = 0 ; r2 = x ; r3 = 0 ; r4 = y ; // 2 x (real,imaginary)
        }
      } else { // complex conjugate pair biquadratic roots x +/- iy.
        nr = 0 ; nc = 4 ;
        x /= 2 ; // re
        y /= 2 ; // im
        valueType z = hypot(x,y) ;
        y = sqrt(z - x) ;
        x = sqrt(z + x) ;
        r1 = -x ;
        r2 = y ;
        r3 = x ;
        r4 = y ;
      }
      return 0 ;
    }
 
    /*
    .. The general case. Rescale quartic polynomial, such that largest absolute
    .. coefficient is (exactly!) equal to 1. Honor the presence of a special
    .. quartic case that might have been obtained during the rescaling process
    .. (due to underflow in the coefficients).
    */

    valueType const A3 = B/A ;
    valueType const A2 = C/A ;
    valueType const A1 = D/A ;
    valueType const A0 = E/A ;
    valueType q0, q1, q2, q3, scale ;
    indexType i_case ;
    scaleQuarticMonicPolynomial( A3, A2, A1, A0, q3, q2, q1, q0, i_case, scale ) ;

    /*
    ..  The general quartic case. Search for stationary points. Set the first
    ..  derivative polynomial (cubic) equal to zero and find its roots.
    ..  Check, if any minimum point of Q(x) is below zero, in which case we
    ..  must have real roots for Q(x). Hunt down only the real root, which
    ..  will potentially converge fastest during Newton iterates. The remaining
    ..  roots will be determined by deflation Q(x) -> cubic.
    ..
    ..  The best roots for the Newton iterations are the two on the opposite
    ..  ends, i.e. those closest to the +2 and -2. Which of these two roots
    ..  to take, depends on the location of the Q(x) minima x = s and x = u,
    ..  with s > u. There are three cases:
    ..
    ..     1) both Q(s) and Q(u) < 0
    ..        ----------------------
    ..
    ..        The best root is the one that corresponds to the lowest of
    ..        these minima. If Q(s) is lowest -> start Newton from +2
    ..        downwards (or zero, if s < 0 and a0 > 0). If Q(u) is lowest
    ..        -> start Newton from -2 upwards (or zero, if u > 0 and a0 > 0).
    ..
    ..     2) only Q(s) < 0
    ..        -------------
    ..
    ..        With both sides +2 and -2 possible as a Newton starting point,
    ..        we have to avoid the area in the Q(x) graph, where inflection
    ..        points are present. Solving Q''(x) = 0, leads to solutions
    ..        x = -a3/4 +/- discriminant, i.e. they are centered around -a3/4.
    ..        Since both inflection points must be either on the r.h.s or l.h.s.
    ..        from x = s, a simple test where s is in relation to -a3/4 allows
    ..        us to avoid the inflection point area.
    ..
    ..     3) only Q(u) < 0
    ..        -------------
    ..
    ..        Same of what has been said under 2) but with x = u.
    */

    // Q(x) = a0 + a1 * x + a2 * x^2 + a3 * x^3 + x^4 ;
    // Q'(x) = a1 + 2*a2 * x + 3*a3 * x^2 + 4*x ;
    valueType c2 = 3 * (q3/4) ;
    valueType c1 = q2/2 ;
    valueType c0 = q1/4 ;
    
    valueType u, t, s ; // root according to paper
    indexType iter = solveCubic( 1, c2, c1, c0, u, t, s, nr, nc ) ;

    valueType Qs, Qu ;
    if ( nr > 1 ) {
      Qs = evalMonicQuartic( s, q3, q2, q1, q0 ) ;
      Qu = evalMonicQuartic( u, q3, q2, q1, q0 ) ;
    } else {
      Qs = evalMonicQuartic( u, q3, q2, q1, q0 ) ;
      Qu = 1 ;
    }

    valueType tmp = q0 >= 0 ? 0 : 2 ;
    nr = 1 ;
    if ( Qs < 0 && Qu < 0 ) {
      if ( Qs < Qu ) r1 = s < 0 ?  tmp :  2 ;
      else           r1 = u > 0 ? -tmp : -2 ;
    } else if ( Qs < 0 ) {
      if ( 4*s < -q3 ) r1 = s > 0 ? tmp : -2 ;
      else             r1 = s < 0 ? tmp :  2 ;
    } else if ( Qu < 0 ) {
      if ( 4*u < -q3 ) r1 = u > 0 ? tmp : -2 ;
      else             r1 = u < 0 ? tmp :  2 ;
    } else {
      nr = 0 ;
    }
    /*
    ..  Do all necessary Newton iterations. In case we have more than 2 oscillations,
    ..  exit the Newton iterations and switch to bisection. Note, that from the
    ..  definition of the Newton starting point, we always have Q(x) > 0 and Q'(x)
    ..  starts (-ve/+ve) for the (-2/+2) starting points and (increase/decrease) smoothly
    ..  and staying (< 0 / > 0). In practice, for extremely shallow Q(x) curves near the
    ..  root, the Newton procedure can overshoot slightly due to rounding errors when
    ..  approaching the root. The result are tiny oscillations around the root. If such
    ..  a situation happens, the Newton iterations are abandoned after 3 oscillations
    ..  and further location of the root is done using bisection starting with the
    ..  oscillation brackets.
    */
    if ( nr > 0 ) {
      iter += zeroQuarticByNewtonBisection( q3, q2, q1, q0, r1 ) ;
      r1 *= scale ;

      /*
      ..  Find remaining roots -> reduce to cubic. The reduction to a cubic polynomial
      ..  is done using composite deflation to minimize rounding errors. Also, while
      ..  the composite deflation analysis is done on the reduced quartic, the actual
      ..  deflation is being performed on the original quartic again to avoid enhanced
      ..  propagation of root errors.
      */
      deflateQuarticPolynomial( A, B, C, D, E, r1, q2, q1, q0 );
      iter += solveCubic( A, q2, q1, q0, r2, r3, r4, nr, nc ) ;
      if ( nr > 0 && r1 > r2 ) std::swap( r1, r2 ) ;
      if ( nr > 1 && r2 > r3 ) std::swap( r2, r3 ) ;
      if ( nr > 2 && r3 > r4 ) std::swap( r3, r4 ) ;
      ++nr ;
    } else {
      /*
      .. If no real roots have been found by now, only complex roots are
      .. possible. Find real parts of roots first, followed by imaginary
      .. components.
      */
      s = q3/2 ;
      t = s * s - q2 ;
      u = s * t + q1 ; // value of Q'(-q3/4) at stationary point -q3/4
      bool notZero = (abs(u) >= machepsi) ; // H(-a3/4) is considered > 0 at stationary point
      bool minimum ;
      if ( q3 != 0 ) {
        s = q1 / q3 ;
        minimum = q0 > s * s ; // H''(-q3/4) > 0 -> minimum
      } else {
        minimum = 4 * q0 > q2 * q2 ; // H''(-q3/4) > 0 -> minimum
      }

      bool iterate = notZero || (!notZero && minimum) ;

      valueType a, b, c, d ;
      if ( iterate ) {
        valueType x = q3 >= 0 ? 2 : -2 ; // initial root -> target = smaller mag root
        iter += zeroHexicByNewtonBisection( q3, q2, q1, q0, x ) ;

        a = x*scale ;   // 1st real component -> a
        b = -A3/2 - a ; // 2nd real component -> b

        x = 4 * a + A3 ; // Q'''(a)
        valueType y = x + 2*A3 ;
        y = y * a + 2*A2 ; // Q'(a)
        y = y * a + A1 ;
        y /= x ;
        if ( y < 0 ) y = 0 ; // ensure Q'(a) / Q'''(a) >= 0

        x = 4 * b + A3 ; // Q'''(b)
        valueType z = x + 2*A3 ;
        z = z * b + 2*A2 ; // Q'(b)
        z = z * b + A1 ;
        z /= x ;
        if ( z < 0 ) z = 0 ; // ensure Q'(b) / Q'''(b) >= 0

        c = a * a ;              // store a^2 for later
        d = b * b ;              // store b^2 for later
        valueType s = c + y ;    // magnitude^2 of (a + iy) root
        valueType t = d + z ;    // magnitude^2 of (b + iz) root

        if ( s > t ) {           // minimize imaginary error
          c = sqrt(y) ;          // 1st imaginary component -> c
          d = sqrt(A0 / s - d) ; // 2nd imaginary component -> d
        } else {
          c = sqrt(A0 / t - c) ; // 1st imaginary component -> c
          d = sqrt(z) ;          // 2nd imaginary component -> d
        }

      } else { // no bisection -> real components equal

        a = -A3/4 ; // 1st real component -> a
        b = a ;     // 2nd real component -> b = a

        valueType x = a + A3 ;
        x = x * a + A2 ; // Q(a)
        x = x * a + A1 ;
        x = x * a + A0 ;
        valueType y = A2/2 - 3*(a*a) ; // Q''(a) / 2
        valueType z = y * y - x ;
        z = z > 0 ? sqrt(z) : 0 ;     // force discriminant to be >= 0
                                      // square root of discriminant
        y = y > 0 ? y + z : y - z ;   // larger magnitude root
        x /= y  ;                     // smaller magnitude root
        c = y < 0 ? 0 : sqrt(y) ;     // ensure root of biquadratic > 0
        d = x < 0 ? 0 : sqrt(x) ;     // large magnitude imaginary component
      }

      nr = 0 ;
      nc = 4 ;
      if      (a > b) { r1 = a ; r2 = c ; r3 = b ; r4 = d ; }
      else if (a < b) { r1 = b ; r2 = d ; r3 = a ; r4 = c ; }
      else            { r1 = a ; r2 = c ; r3 = a ; r4 = d ; }
    }
    return iter ;
  }

}

// EOF: PolynomialRoots.cc
