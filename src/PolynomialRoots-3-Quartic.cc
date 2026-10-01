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

#include "PolynomialRoots.hh"

#include <ostream>

#define MAX_ITER_SAFETY 50

//
// Safety factor on the machine epsilon used to decide that the quartic touches
// the real axis at a stationary point, that is that it has a real root of even
// multiplicity there.  Calibrated on (x-r)^2 (x-s)^2 and (x-r)^4 families.
//
#define QRF_TOUCHING_ROOT_FACTOR 64

//
// Largest relative residual accepted for an isolated real root, eps^(3/4).
//
#define QRF_ACCEPT_RESIDUAL ( sqrt( machine_eps<T_real>() ) * sqrt( sqrt( machine_eps<T_real>() ) ) )

namespace PolynomialRoots
{

  using std::abs;
  using std::cbrt;
  using std::hypot;
  using std::sqrt;

  //!
  //! Collect the computed real roots.
  //!
  //! The storage layout depends on `m_ncplx`:
  //!   - 0 : m_r0..m_r3 hold up to `m_nreal` real roots;
  //!   - 2 : (m_r0,m_r1) is a conjugate pair, m_r2/m_r3 hold up to `m_nreal`
  //!         real roots;
  //!   - 4 : two conjugate pairs, no real root.
  //!
  //! `m_nreal` is smaller than 4 only for degenerate quartics (A == 0), hence
  //! the explicit bound: without it a degenerate cubic would report a fourth,
  //! never computed, "root".
  //!
  template <typename T_real, typename T_complex> integer QuarticT<T_real, T_complex>::get_real_roots( T_real r[] ) const
  {
    integer nr = 0;
    switch ( m_ncplx )
    {
      case 0:
        if ( m_nreal > 0 ) r[nr++] = m_r0;
        if ( m_nreal > 1 ) r[nr++] = m_r1;
        if ( m_nreal > 2 ) r[nr++] = m_r2;
        if ( m_nreal > 3 ) r[nr++] = m_r3;
        break;
      case 2:
        if ( m_nreal > 0 ) r[nr++] = m_r2;
        if ( m_nreal > 1 ) r[nr++] = m_r3;
        break;
      default: break;  // 4 complex roots: nothing real to report
    }
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer QuarticT<T_real, T_complex>::get_positive_roots( T_real r[] ) const
  {
    T_real        tmp[4];
    integer const n  = get_real_roots( tmp );
    integer       nr = 0;
    for ( integer i{ 0 }; i < n; ++i )
      if ( tmp[i] > 0 ) r[nr++] = tmp[i];
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer QuarticT<T_real, T_complex>::get_negative_roots( T_real r[] ) const
  {
    T_real        tmp[4];
    integer const n  = get_real_roots( tmp );
    integer       nr = 0;
    for ( integer i{ 0 }; i < n; ++i )
      if ( tmp[i] < 0 ) r[nr++] = tmp[i];
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer QuarticT<T_real, T_complex>::get_roots_in_range( T_real const & a, T_real const & b, T_real r[] ) const
  {
    T_real        tmp[4];
    integer const n  = get_real_roots( tmp );
    integer       nr = 0;
    for ( integer i{ 0 }; i < n; ++i )
      if ( tmp[i] >= a && tmp[i] <= b ) r[nr++] = tmp[i];
    return nr;
  }

  template <typename T_real, typename T_complex>
  integer QuarticT<T_real, T_complex>::get_roots_in_open_range( T_real const & a, T_real const & b, T_real r[] ) const
  {
    T_real        tmp[4];
    integer const n  = get_real_roots( tmp );
    integer       nr = 0;
    for ( integer i{ 0 }; i < n; ++i )
      if ( tmp[i] > a && tmp[i] < b ) r[nr++] = tmp[i];
    return nr;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // Compute (p/q) * 2^shift without any intermediate overflow or underflow:
  // the mantissas are divided first, the binary exponents are added afterwards.
  //
  template <typename T_real> static T_real scaled_ratio( T_real const & p, T_real const & q, int const shift )
  {
    using std::frexp;
    using std::ldexp;
    if ( p == 0 ) return T_real( 0 );
    int          ep{ 0 }, eq{ 0 };
    T_real const mp{ frexp( p, &ep ) };
    T_real const mq{ frexp( q, &eq ) };
    return ldexp( mp / mq, ep - eq + shift );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // x^4 + A x^3 + B x^2 + C x + D
  template <typename T_real> static void scaleQuarticMonicPolynomial(
    T_real const & A,
    T_real const & B,
    T_real const & C,
    T_real const & D,
    T_real &       AS,
    T_real &       BS,
    T_real &       CS,
    T_real &       DS,
    integer &      i_case,
    T_real &       scale )
  {
    T_real const a = abs( A );
    T_real const b = sqrt( abs( B ) );
    T_real const c = cbrt( abs( C ) );
    T_real const d = sqrt( sqrt( abs( D ) ) );

    if ( a < b )
    {
      if ( b < c )
      {
        if ( c < d )
          i_case = 0;  // a < b < c < d        --> d MAX
        else
          i_case = 1;  // a < b < c and d <= c --> c MAX
      }
      else
      {
        if ( b < d )
          i_case = 0;  // a < b and c <= b < d --> d MAX
        else
          i_case = 2;  // a < b and c <= b and d <= b --> b MAX
      }
    }
    else
    {
      if ( a < c )
      {
        if ( c < d )
          i_case = 0;  // b <= a < c < d        --> d MAX
        else
          i_case = 1;  // b <= a < c and d <= c --> c MAX
      }
      else
      {
        if ( a < d )
          i_case = 0;  // b <= a and c <= a and a < d --> d MAX
        else
          i_case = 3;  // b <= a and c <= a and d <= a --> a MAX
      }
    }

    switch ( i_case )
    {
      case 0:
        scale = d;
        AS    = A / d;
        BS    = ( B / d ) / d;
        CS    = ( ( C / d ) / d ) / d;
        DS    = D > 0 ? 1 : -1;
        break;
      case 1:
        scale = c;
        AS    = A / c;
        BS    = ( B / c ) / c;
        CS    = C > 0 ? 1 : -1;
        DS    = ( ( ( D / c ) / c ) / c ) / c;
        break;
      case 2:
        scale = b;
        AS    = A / b;
        BS    = B > 0 ? 1 : -1;
        CS    = ( ( C / b ) / b ) / b;
        DS    = ( ( ( D / b ) / b ) / b ) / b;
        break;
      case 3:
        scale = a;
        AS    = A > 0 ? 1 : -1;
        BS    = ( B / a ) / a;
        CS    = ( ( C / a ) / a ) / a;
        DS    = ( ( ( D / a ) / a ) / a ) / a;
        break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T_real> static T_real evalHexic(
    T_real const & x,
    T_real const & q3,
    T_real const & q2,
    T_real const & q1,
    T_real const & q0 )
  {
    T_real       t1   = x + q3;
    T_real       t2   = x + t1;
    T_real const t3   = x + t2;
    t1                = t1 * x + q2;
    t2                = t2 * x + t1;
    t1                = t1 * x + q1;
    T_real const Q    = t1 * x + q0;                     // Q(x)
    T_real const dQ   = t2 * x + t1;                     // Q'(x)
    T_real const ddQ  = t3 * x + t2;                     // Q''(x) / 2
    T_real const dddQ = x + t3;                          // Q'''(x) / 6
    return Q * dddQ * dddQ - dQ * ddQ * dddQ + dQ * dQ;  // H(x), usually < 0
  }

  template <typename T_real> static void evalHexic(
    T_real const & x,
    T_real const & q3,
    T_real const & q2,
    T_real const & q1,
    T_real const & q0,
    T_real &       p,
    T_real &       dp )
  {
    T_real       t1   = x + q3;
    T_real       t2   = x + t1;
    T_real const t3   = x + t2;
    t1                = t1 * x + q2;
    t2                = t2 * x + t1;
    t1                = t1 * x + q1;
    T_real const Q    = t1 * x + q0;                                   // Q(x)
    T_real const dQ   = t2 * x + t1;                                   // Q'(x)
    T_real const ddQ  = t3 * x + t2;                                   // Q''(x) / 2
    T_real const dddQ = x + t3;                                        // Q'''(x) / 6
    p                 = Q * dddQ * dddQ - dQ * ddQ * dddQ + dQ * dQ;   // H(x), usually < 0
    dp                = 2 * dddQ * ( 4 * Q - dQ * dddQ - ddQ * ddQ );  // H'(x)
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = (x-r)*(a4*x^3+b2*x^2+b1*x+b0)
  template <typename T_real> static void deflateQuarticPolynomial(
    T_real const & a4,
    T_real const & a3,
    T_real const & a2,
    T_real const & a1,
    T_real const & a0,
    T_real const & r,
    T_real &       b2,
    T_real &       b1,
    T_real &       b0 )
  {
    /*
    ..  The cross-over index is the k maximizing |a_k| |r|^k.  Forming those
    ..  products directly overflows as soon as |r| is large - for the quartic
    ..  with roots 1e154, 1e152, 10, 1 every one of them becomes infinite, they
    ..  all compare equal and the first, wrong, branch is taken.  Dividing by
    ..  |r|^4 when |r| > 1 (and multiplying otherwise) leaves the ordering
    ..  unchanged and cannot overflow either way.
    */
    T_real const ar{ abs( r ) };
    T_real       w[5];
    if ( ar > 1 )
    {
      w[4] = abs( a4 );
      w[3] = abs( a3 ) / ar;
      w[2] = ( abs( a2 ) / ar ) / ar;
      w[1] = ( ( abs( a1 ) / ar ) / ar ) / ar;
      w[0] = ( ( ( abs( a0 ) / ar ) / ar ) / ar ) / ar;
    }
    else
    {
      T_real const r2{ ar * ar };
      w[0] = abs( a0 );
      w[1] = abs( a1 ) * ar;
      w[2] = abs( a2 ) * r2;
      w[3] = abs( a3 ) * ar * r2;
      w[4] = abs( a4 ) * r2 * r2;
    }
    integer i_cross{ 0 };
    for ( integer k{ 1 }; k <= 4; ++k )
      if ( w[k] > w[i_cross] ) i_cross = k;
    switch ( i_cross )
    {
      case 0:
        b2 = a3 + a4 * r;
        b1 = a2 + r * b2;
        b0 = a1 + r * b1;
        break;
      case 1:
        b2 = a3 + a4 * r;
        b1 = a2 + r * b2;
        b0 = -a0 / r;
        break;
      case 2:
        b2 = a3 + a4 * r;
        b0 = -a0 / r;
        b1 = ( b0 - a1 ) / r;
        break;
      case 3:
      case 4:
        b0 = -a0 / r;
        b1 = ( b0 - a1 ) / r;
        b2 = ( b1 - a2 ) / r;
        break;
    }
  }

  /*
  ||   _   _               _              ____  _               _   _
  ||  | \ | | _____      _| |_ ___  _ __ | __ )(_)___  ___  ___| |_(_) ___  _ __
  ||  |  \| |/ _ \ \ /\ / / __/ _ \| '_ \|  _ \| / __|/ _ \/ __| __| |/ _ \| '_ \
  ||  | |\  |  __/\ V  V /| || (_) | | | | |_) | \__ \  __/ (__| |_| | (_) | | | |
  ||  |_| \_|\___| \_/\_/  \__\___/|_| |_|____/|_|___/\___|\___|\__|_|\___/|_| |_|
  */

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Translate to C from Polynomial234RootSolvers
  template <typename T_real> static integer zeroQuarticByNewtonBisection(
    T_real const & a,
    T_real const & b,
    T_real const & c,
    T_real const & d,
    T_real &       x,
    bool &         ok )
  {
    ok = false;
    T_real p, dp;
    eval_monic_quartic( x, a, b, c, d, p, dp );
    T_real t = p;  // save p(x) for sign comparison
    x -= p / dp;   // 1st improved root

    integer iter       = 1;
    integer oscillate  = 0;
    integer nconverged = 0;
    bool    bisection  = false;
    bool    converged  = false;
    bool    have_s     = false;
    bool    have_u     = false;
    T_real  s          = 0;
    T_real  u          = 0;  // to mute warning
    while ( !( nconverged > 1 || bisection ) && iter < MAX_ITER_SAFETY )
    {
      ++iter;
      T_real ddp;
      eval_monic_quartic( x, a, b, c, d, p, dp, ddp );
      if ( p * t < 0 )
      {  // does Newton start oscillating ?
        if ( p < 0 )
        {
          ++oscillate;  // increment oscillation counter
          s      = x;   // save lower bisection bound
          have_s = true;
        }
        else
        {
          u      = x;  // save upper bisection bound
          have_u = true;
        }
        t = p;  // save current p(x)
      }
      // dp = (p*dp)/(dp*dp-p*ddp); // double zero correction
      dp = ( p * dp ) / ( dp * dp - p * ddp / 2 );  // Halley correction
      // dp = p/dp; // Newton correction
      x -= dp;                    // new Newton root
      bisection = oscillate > 2;  // activate bisection
      /*
      ..  RELATIVE convergence test.  The customary (1+|x|)*eps carries an
      ..  absolute floor of one epsilon: harmless while the roots are O(1), but
      ..  after Flocke's scaling a quartic with a tiny leading coefficient has
      ..  roots down at 1e-18, where that floor is eighty times the root itself
      ..  and the iteration stops with eight correct digits instead of sixteen.
      ..  The iteration cap above keeps the loop bounded when x tends to zero.
      */
      converged = abs( dp ) <= abs( x ) * machine_eps<T_real>();  // Newton convergence indicator
      if ( converged )
        ++nconverged;
      else
        nconverged = 0;
    }
    // The oscillation brackets are usable only when BOTH of them have been
    // recorded: s and u start at 0 and, without this guard, a Newton run that
    // never oscillated would "bisect" over the meaningless interval [0,0] and
    // silently accept the last, unconverged iterate.
    if ( bisection && have_s && have_u )
    {
      t = u - s;  // initial bisection interval
      while ( abs( t ) > abs( x ) * machine_eps<T_real>() && iter < MAX_ITER_SAFETY )
      {  // bisection iterates
        ++iter;
        p = eval_monic_quartic( x, a, b, c, d );
        if ( p < 0 )
          s = x;
        else
          u = x;            // keep bracket on root
        t = ( u - s ) / 2;  // new bisection interval
        x = s + t;          // new bisection root
      }
      ok = true;
    }
    else
    {
      ok = converged;
    }
    ok = ok && isfinite( x );
    return iter;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // Safeguarded Newton (Newton step when it stays inside the bracket and makes
  // progress, bisection otherwise) for the monic quartic
  //
  //    Q(x) = x^4 + q3 x^3 + q2 x^2 + q1 x + q0
  //
  // on an interval [lo,hi] where Q changes sign.  Unlike the heuristic routine
  // above this one is unconditionally convergent, and is used as a fallback
  // whenever the fast path fails.
  //
  template <typename T_real> static integer zeroQuarticInBracket(
    T_real const & q3,
    T_real const & q2,
    T_real const & q1,
    T_real const & q0,
    T_real const & lo,
    T_real const & hi,
    T_real &       x )
  {
    T_real const flo{ eval_monic_quartic( lo, q3, q2, q1, q0 ) };
    T_real const fhi{ eval_monic_quartic( hi, q3, q2, q1, q0 ) };
    if ( flo == 0 )
    {
      x = lo;
      return 0;
    }
    if ( fhi == 0 )
    {
      x = hi;
      return 0;
    }
    // xl is the end where Q < 0, xh the end where Q > 0
    T_real xl{ flo < 0 ? lo : hi };
    T_real xh{ flo < 0 ? hi : lo };

    T_real dx_old{ abs( hi - lo ) };
    T_real dx{ dx_old };
    x = ( lo + hi ) / 2;

    T_real p, dp;
    eval_monic_quartic( x, q3, q2, q1, q0, p, dp );

    integer iter{ 0 };
    for ( ; iter < 200; ++iter )
    {
      bool const   newton_out_of_range{ ( ( x - xh ) * dp - p ) * ( ( x - xl ) * dp - p ) > 0 };
      bool const   newton_too_slow{ abs( 2 * p ) > abs( dx_old * dp ) };
      T_real const xl_prev{ xl };
      T_real const xh_prev{ xh };
      dx_old = dx;
      if ( newton_out_of_range || newton_too_slow )
      {
        dx = ( xh - xl ) / 2;
        x  = xl + dx;
      }
      else
      {
        dx = p / dp;
        x -= dx;
      }
      /*
      ..  The stopping rule is RELATIVE.  The usual (1+|x|)*eps form hides an
      ..  absolute floor of one eps, which is harmless while the roots are O(1)
      ..  but fatal here: after Flocke's scaling a quartic with a tiny leading
      ..  coefficient has roots of the order of eps itself, and an absolute
      ..  floor would stop the iteration an order of magnitude away from them.
      ..  The bracket collapsing onto two adjacent floating point numbers is the
      ..  companion criterion, needed when the root is exactly zero.
      */
      if ( abs( dx ) <= abs( x ) * machine_eps<T_real>() ) break;
      if ( x == xl_prev || x == xh_prev ) break;
      eval_monic_quartic( x, q3, q2, q1, q0, p, dp );
      if ( p < 0 )
        xl = x;
      else
        xh = x;
    }
    return iter;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // Relative residual of the monic quartic at x, used to validate the root
  // returned by the fast path.
  //
  template <typename T_real> static T_real quarticRelativeResidual(
    T_real const & q3,
    T_real const & q2,
    T_real const & q1,
    T_real const & q0,
    T_real const & x )
  {
    T_real const ax{ abs( x ) };
    T_real const den{ ( ( ( ax + abs( q3 ) ) * ax + abs( q2 ) ) * ax + abs( q1 ) ) * ax + abs( q0 ) };
    if ( !( den > 0 ) ) return T_real( 0 );
    return abs( eval_monic_quartic( x, q3, q2, q1, q0 ) ) / den;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // Unconditionally correct search for ONE real root of the monic quartic.
  //
  // Q is monotone on each interval delimited by its stationary points, and all
  // its real roots lie in [-R,R] with R = 1 + max|q_k| (Cauchy bound).  Hence
  // scanning the sign of Q over the breakpoints -R, u, t, s, +R either exhibits
  // a bracket - and then `zeroQuarticInBracket` converges - or proves that no
  // real root exists.
  //
  template <typename T_real> static bool findRealRootOfQuartic(
    T_real const & q3,
    T_real const & q2,
    T_real const & q1,
    T_real const & q0,
    T_real const & u,
    T_real const & t,
    T_real const & s,
    bool const     three_stationary_points,
    T_real &       x,
    integer &      iter )
  {
    T_real R{ abs( q3 ) };
    if ( abs( q2 ) > R ) R = abs( q2 );
    if ( abs( q1 ) > R ) R = abs( q1 );
    if ( abs( q0 ) > R ) R = abs( q0 );
    R += 1;

    T_real  bp[5];
    integer n{ 0 };
    bp[n++] = -R;
    if ( three_stationary_points )
    {
      // u <= t <= s, all of them inside [-R,R]
      if ( u > bp[n - 1] && u < R ) bp[n++] = u;
      if ( t > bp[n - 1] && t < R ) bp[n++] = t;
      if ( s > bp[n - 1] && s < R ) bp[n++] = s;
    }
    else if ( s > bp[n - 1] && s < R ) { bp[n++] = s; }
    bp[n++] = R;

    /*
    ..  Every sign change is solved, and the root of LARGEST modulus is kept.
    ..  That is Flocke's own rationale for starting Newton at +-2: the extreme
    ..  root is the best conditioned one and the one the composite deflation
    ..  handles best.  Taking the first bracket instead is wrong whenever the
    ..  leading coefficient is tiny, because the scaled quartic then has one
    ..  root near +-1 and the others down at the level of the machine epsilon,
    ..  and the leftmost bracket is precisely the hopeless one.
    */
    bool    found{ false };
    T_real  best_x{ 0 };
    integer total_iter{ 0 };

    T_real f_prev{ eval_monic_quartic( bp[0], q3, q2, q1, q0 ) };
    if ( f_prev == 0 )
    {
      found  = true;
      best_x = bp[0];
    }
    for ( integer i{ 1 }; i < n; ++i )
    {
      T_real const f_cur{ eval_monic_quartic( bp[i], q3, q2, q1, q0 ) };
      T_real       cand{ 0 };
      bool         have_cand{ false };
      if ( f_cur == 0 )
      {
        cand      = bp[i];
        have_cand = true;
      }
      else if ( ( f_prev < 0 ) != ( f_cur < 0 ) )
      {
        total_iter += zeroQuarticInBracket( q3, q2, q1, q0, bp[i - 1], bp[i], cand );
        have_cand = true;
      }
      if ( have_cand && ( !found || abs( cand ) > abs( best_x ) ) )
      {
        found  = true;
        best_x = cand;
      }
      f_prev = f_cur;
    }
    if ( found )
    {
      x    = best_x;
      iter = total_iter;
      return true;
    }

    /*
    ..  No sign change.  That still does not prove the absence of real roots: a
    ..  root of EVEN multiplicity touches the axis without crossing it, so for
    ..  (x-r)^2 (x-s)^2 and for (x-r)^4 the scan above finds nothing while all
    ..  four roots are real.  Such a root sits on a stationary point and makes Q
    ..  numerically vanish there, so the minimum of the RELATIVE residual over
    ..  the stationary points settles the question.  The relative form is what
    ..  makes the test usable: the absolute residual depends on the scale of the
    ..  coefficients and no fixed threshold works for it.
    */
    T_real touch_x{ 0 };
    T_real touch_res{ 0 };
    bool   have_touch{ false };
    for ( integer i{ 1 }; i + 1 < n; ++i )
    {
      T_real const res{ quarticRelativeResidual( q3, q2, q1, q0, bp[i] ) };
      if ( !have_touch || res < touch_res )
      {
        touch_x    = bp[i];
        touch_res  = res;
        have_touch = true;
      }
    }
    if ( have_touch && touch_res <= QRF_TOUCHING_ROOT_FACTOR * machine_eps<T_real>() )
    {
      x    = touch_x;
      iter = 0;
      return true;
    }
    return false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Translate to C from Polynomial234RootSolvers
  template <typename T_real> static integer zeroHexicByNewtonBisection(
    T_real const & q3,
    T_real const & q2,
    T_real const & q1,
    T_real const & q0,
    T_real &       x )
  {
    T_real p, dp;
    evalHexic( x, q3, q2, q1, q0, p, dp );
    T_real t = p;  // save p(x) for sign comparison
    x -= p / dp;   // 1st improved root

    integer iter      = 1;
    integer oscillate = 0;
    bool    bisection = false;
    bool    converged = false;
    bool    have_s    = false;
    bool    have_u    = false;
    T_real  s         = 0;
    T_real  u         = 0;  // to mute warning
    // Both loops below are bounded by MAX_ITER_SAFETY and use the absolute plus
    // relative criterion (1+|x|)*eps.  The original formulation, |dx| <= |x|*eps
    // with no iteration cap, cannot terminate when x converges to zero.
    while ( !( converged || bisection ) && iter < MAX_ITER_SAFETY )
    {
      ++iter;
      evalHexic( x, q3, q2, q1, q0, p, dp );
      if ( p * t < 0 )
      {  // does Newton start oscillating ?
        if ( p < 0 )
        {
          ++oscillate;  // increment oscillation counter
          s      = x;   // save lower bisection bound
          have_s = true;
        }
        else
        {
          u      = x;  // save upper bisection bound
          have_u = true;
        }
        t = p;  // save current p(x)
      }
      if ( dp == 0 ) break;
      dp = p / dp;  // Newton correction
      if ( !isfinite( dp ) ) break;
      x -= dp;                                                    // new Newton root
      bisection = oscillate > 2;                                  // activate bisection
      converged = abs( dp ) <= abs( x ) * machine_eps<T_real>();  // Newton convergence indicator
    }
    if ( bisection && have_s && have_u )
    {
      t = u - s;  // initial bisection interval
      while ( abs( t ) > abs( x ) * machine_eps<T_real>() && iter < MAX_ITER_SAFETY )
      {  // bisection iterates
        ++iter;
        p = evalHexic( x, q3, q2, q1, q0 );
        if ( p < 0 )
          s = x;
        else
          u = x;            // keep bracket on root
        t = ( u - s ) / 2;  // new bisection interval
        x = s + t;          // new bisection root
      }
    }
    return iter;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*
  ||       _       __ _       _
  ||    __| | ___ / _| | __ _| |_ ___
  ||   / _` |/ _ \ |_| |/ _` | __/ _ \
  ||  | (_| |  __/  _| | (_| | ||  __/
  ||   \__,_|\___|_| |_|\__,_|\__\___|
  */
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // a3*x^3 + a2*x^2 + a1*x + a0 = (x-r)*(a3*x^2+b1*x+b0)


  /*\
   *  Calculate the zeros of the quartic a*z^4 + b*z^3 + c*z^2 + d*x + e.
  \*/
  template <typename T_real, typename T_complex> void QuarticT<T_real, T_complex>::find_roots()
  {
    T_real const & A = m_ABCDE[0];
    T_real const & B = m_ABCDE[1];
    T_real const & C = m_ABCDE[2];
    T_real const & D = m_ABCDE[3];
    T_real const & E = m_ABCDE[4];

    m_iter = m_nreal = m_ncplx = 0;

    // special cases
    if ( A == 0 )
    {
      CubicT<T_real, T_complex> csolve( B, C, D, E );
      m_nreal = csolve.num_roots();
      m_r0 = m_r1 = m_r2 = m_r3 = 0;
      switch ( m_nreal )
      {
        case 3: m_r2 = csolve.real_root2(); [[fallthrough]];
        case 2: m_r1 = csolve.real_root1(); [[fallthrough]];
        case 1: m_r0 = csolve.real_root0();
      }
      if ( csolve.complex_root() )
      {
        m_ncplx = 2;
        m_nreal -= 2;
      }
      return;
    }
    if ( E == 0 )
    {
      CubicT<T_real, T_complex> csolve( A, B, C, D );
      m_nreal = csolve.num_roots();
      m_r0 = m_r1 = m_r2 = m_r3 = 0;
      switch ( m_nreal )
      {
        case 3: m_r2 = csolve.real_root2(); [[fallthrough]];
        case 2: m_r1 = csolve.real_root1(); [[fallthrough]];
        case 1: m_r0 = csolve.real_root0();
      }
      ++m_nreal;
      if ( csolve.complex_root() )
      {
        m_ncplx = 2;
        m_nreal -= 2;
      }
      if ( m_nreal == 4 )
      {  // caso regolare, 4 radici reali maintain order
        if ( m_r3 < m_r2 ) std::swap( m_r2, m_r3 );
        if ( m_r2 < m_r1 ) std::swap( m_r1, m_r2 );
        if ( m_r1 < m_r0 ) std::swap( m_r0, m_r1 );
      }
      else if ( m_ncplx == 2 && m_nreal == 2 && m_r2 > m_r3 )
      {  // (m_r0,m_r1) complex pair, m_r2 and m_r3 real: keep them ordered
        std::swap( m_r2, m_r3 );
      }
      return;
    }
    if ( B == 0 && D == 0 )
    {  // biquadratic case
      // A x^4 + C x^2 + E
      QuadraticT<T_real, T_complex> qsolve( A, C, E );
      T_real                        x = qsolve.real_root0();
      T_real                        y = qsolve.real_root1();
      if ( qsolve.complex_root() )
      {
        // complex conjugate pair biquadratic roots x +/- iy.
        m_ncplx = 4;
        x /= 2;  // re
        y /= 2;  // im
        T_real z = hypot( x, y );
        y        = sqrt( z - x );
        x        = sqrt( z + x );
        m_r0     = -x;
        m_r1     = y;
        m_r2     = x;
        m_r3     = y;
      }
      else
      {
        // real roots of quadratic are ordered x <= y
        if ( x >= 0 )
        {  // y >= 0
          m_nreal = 4;
          x       = sqrt( x );
          y       = sqrt( y );
          m_r0    = -y;
          m_r1    = -x;
          m_r2    = x;
          m_r3    = y;
        }
        else if ( y >= 0 )
        {  // x < 0 && y >= 0
          m_nreal = m_ncplx = 2;
          x                 = sqrt( -x );
          y                 = sqrt( y );
          m_r0              = 0;
          m_r1              = x;  // (real,imaginary)
          m_r2              = -y;
          m_r3              = y;
        }
        else
        {  // x < 0 && y < 0
          m_ncplx = 4;
          x       = sqrt( -x );
          y       = sqrt( -y );
          m_r0    = 0;
          m_r1    = x;
          m_r2    = 0;
          m_r3    = y;  // 2 x (real,imaginary)
        }
      }
      return;
    }

    /*
    .. The general case. Rescale quartic polynomial, such that largest absolute
    .. coefficient is (exactly!) equal to 1. Honor the presence of a special
    .. quartic case that might have been obtained during the rescaling process
    .. (due to underflow in the coefficients).
    */

    T_real A3 = B / A;
    T_real A2 = C / A;
    T_real A1 = D / A;
    T_real A0 = E / A;

    /*
    ..  Extremely unbalanced coefficients, such as |e| / |a| > DBL_MAX, make the
    ..  normalized quartic overflow and the solver used to return `inf` roots.
    ..  In that case the unknown is rescaled, x = root_scale * y with
    ..  root_scale a power of two chosen so that |A0| ~ 1; the whole solution is
    ..  then carried out on y and every root is multiplied back at the end.
    */
    T_real root_scale{ 1 };
    bool   rescaled{ false };
    if ( !( isfinite( A3 ) && isfinite( A2 ) && isfinite( A1 ) && isfinite( A0 ) ) )
    {
      using std::frexp;
      using std::ldexp;
      int eA{ 0 }, eE{ 0 };
      frexp( A, &eA );
      if ( E != 0 ) frexp( E, &eE );
      int const m{ ( eE - eA ) / 4 };
      root_scale = ldexp( T_real( 1 ), m );
      A3         = scaled_ratio( B, A, -m );
      A2         = scaled_ratio( C, A, -2 * m );
      A1         = scaled_ratio( D, A, -3 * m );
      A0         = scaled_ratio( E, A, -4 * m );
      rescaled   = true;
      root_assert(
        isfinite( A3 ) && isfinite( A2 ) && isfinite( A1 ) && isfinite( A0 ) && isfinite( root_scale ),
        "QuarticT: coefficients a={}, b={}, c={}, d={}, e={} are too unbalanced to be normalized",
        A,
        B,
        C,
        D,
        E );
    }

    T_real  q0, q1, q2, q3, scale;
    integer i_case;
    scaleQuarticMonicPolynomial( A3, A2, A1, A0, q3, q2, q1, q0, i_case, scale );

    // Q(x)    = q0 + q1 * x + q2 * x^2 + q3 * x^3 + x^4
    // Q'(x)/4 = q1/4 + (q2/2) * x + (3/4)*q3 * x^2 + x^3

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

    T_real c2 = 3 * ( q3 / 4 );
    T_real c1 = q2 / 2;
    T_real c0 = q1 / 4;

    CubicT<T_real, T_complex> qsolve( T_real( 1 ), c2, c1, c0 );
    // When the derivative has three real roots they come out sorted as
    // u <= t <= s.  When it has a conjugate pair, real_root0()/real_root1() are
    // the real and imaginary part of that pair and the single real stationary
    // point is real_root2(): u and t are then meaningless and must not be used.
    bool const three_stationary_points = !qsolve.complex_root();
    T_real     u                       = qsolve.real_root0();  // root according to paper
    T_real     t                       = qsolve.real_root1();
    T_real     s                       = qsolve.real_root2();
    if ( !three_stationary_points ) u = t = s;

    bool must_refine_r3 = true;
    if ( three_stationary_points )
    {
      T_real Qs    = eval_monic_quartic( s, q3, q2, q1, q0 );
      T_real Qu    = eval_monic_quartic( u, q3, q2, q1, q0 );
      bool   q0pos = q0 > 0;
      m_nreal      = 1;
      if ( Qs < 0 && Qu < 0 )
      {
        if ( Qs < Qu )
        {
          m_r3 = 2;
          if ( q0pos && s < 0 ) m_r3 = 0;
        }
        else
        {
          m_r3 = -2;
          if ( q0pos && u > 0 ) m_r3 = 0;
        }
      }
      else if ( Qs < 0 )
      {
        if ( 4 * s < -q3 )
        {
          m_r3 = -2;
          if ( q0pos && s > 0 ) m_r3 = 0;
        }
        else
        {
          m_r3 = 2;
          if ( q0pos && s < 0 ) m_r3 = 0;
        }
      }
      else if ( Qu < 0 )
      {
        if ( 4 * u < -q3 )
        {
          m_r3 = -2;
          if ( q0pos && u > 0 ) m_r3 = 0;
        }
        else
        {
          m_r3 = 2;
          if ( q0pos && u < 0 ) m_r3 = 0;
        }
      }
      else
      {
        // check for astrological combination when s or u are root of the quartic
        /*
        DO NOT WORK
        real_type epsi = 10 * ( 1 + std::abs(q3) +
                                    std::abs(q2) +
                                    std::abs(q1) +
                                    std::abs(q0) ) * machepsi;
        if      ( Qs <= epsi ) r3 = s;
        else if ( Qu <= epsi ) r3 = u;
        else                   nreal = 0;
        */
        if ( Qs == 0 )
        {
          must_refine_r3 = false;
          m_r3           = s;
        }
        else if ( Qu == 0 )
        {
          must_refine_r3 = false;
          m_r3           = u;
        }
        else
        {
          m_nreal = 0;
        }
      }
    }
    else
    {
      // Q has a single stationary point (a minimum) located at x = s.  The test
      // below must therefore look at s: using real_root0() here compared the
      // REAL PART of the complex stationary pair and picked the wrong Newton
      // starting point.
      if ( T_real const Qs = eval_monic_quartic( s, q3, q2, q1, q0 ); Qs <= 0 )
      {
        T_real const tmp = q0 >= 0 ? 0 : 2;
        m_r3             = s > 0 ? -tmp : -2;
        m_nreal          = 1;
      }
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
    bool real_root_ok{ false };
    if ( m_nreal > 0 )
    {
      if ( must_refine_r3 )
      {
        bool conv{ false };
        m_iter += zeroQuarticByNewtonBisection( q3, q2, q1, q0, m_r3, conv );
        real_root_ok = conv;
      }
      else
      {
        real_root_ok = true;
      }
      // Accept the fast path only if the returned value really is a root.  The
      // bound is eps^(3/4): a correctly rounded root leaves a relative residual
      // of a few eps, so this is generous, while a root known to half precision
      // only - residual around sqrt(eps) - is rejected and recomputed below.
      if ( real_root_ok ) real_root_ok = quarticRelativeResidual( q3, q2, q1, q0, m_r3 ) <= QRF_ACCEPT_RESIDUAL;
    }

    /*
    ..  Guaranteed fallback.  The heuristic above can both miss an existing real
    ..  root and return an unconverged iterate (the Newton loop gives up after
    ..  MAX_ITER_SAFETY steps and the oscillation brackets may never have been
    ..  recorded).  Deciding the question on the monotonicity intervals of Q is
    ..  cheap - at most four extra evaluations - and always conclusive.
    */
    if ( !real_root_ok )
    {
      T_real  x_real{ 0 };
      integer extra_iter{ 0 };
      // the fallback result is validated exactly like the fast path one: a
      // quartic whose real root cannot be isolated to the working precision is
      // better served by the all-complex branch than by a bogus deflation
      if (
        findRealRootOfQuartic( q3, q2, q1, q0, u, t, s, three_stationary_points, x_real, extra_iter ) &&
        quarticRelativeResidual( q3, q2, q1, q0, x_real ) <= QRF_ACCEPT_RESIDUAL )
      {
        m_r3    = x_real;
        m_nreal = 1;
        m_iter += extra_iter;
        real_root_ok = true;
      }
      else
      {
        m_nreal = 0;
      }
    }

    if ( real_root_ok )
    {
      m_r3 *= scale;

      /*
      ..  Find remaining roots -> reduce to cubic. The reduction to a cubic polynomial
      ..  is done using composite deflation to minimize rounding errors. Also, while
      ..  the composite deflation analysis is done on the reduced quartic, the actual
      ..  deflation is being performed on the original quartic again to avoid enhanced
      ..  propagation of root errors.
      */
      // Flocke deflates on the ORIGINAL quartic to limit error propagation; in
      // the rescaled case the original coefficients live in a different
      // variable, so the normalized ones must be used instead.
      T_real const d4{ rescaled ? T_real( 1 ) : A };
      if ( rescaled )
        deflateQuarticPolynomial( T_real( 1 ), A3, A2, A1, A0, m_r3, q2, q1, q0 );
      else
        deflateQuarticPolynomial( A, B, C, D, E, m_r3, q2, q1, q0 );

      CubicT<T_real, T_complex> csolve( d4, q2, q1, q0 );
      m_r0 = csolve.real_root0();
      m_r1 = csolve.real_root1();
      m_r2 = csolve.real_root2();
      if ( csolve.complex_root() )
      {
        m_nreal = m_ncplx = 2;
        if ( m_r2 > m_r3 ) std::swap( m_r2, m_r3 );
      }
      else
      {
        m_nreal = 4;
        if ( m_r2 > m_r3 ) std::swap( m_r2, m_r3 );
        if ( m_r1 > m_r2 ) std::swap( m_r1, m_r2 );
        if ( m_r0 > m_r1 ) std::swap( m_r0, m_r1 );
      }
    }
    else
    {
      /*
      .. If no real roots have been found by now, only complex roots are
      .. possible. Find real parts of roots first, followed by imaginary
      .. components.
      */
      s            = q3 / 2;
      t            = s * s - q2;
      u            = s * t + q1;                             // value of Q'(-q3/4) at stationary point -q3/4
      bool notZero = ( abs( u ) >= machine_eps<T_real>() );  // H(-a3/4) is considered > 0 at stationary point
      bool minimum;
      if ( q3 != 0 )
      {
        s       = q1 / q3;
        minimum = q0 > s * s;  // H''(-q3/4) > 0 -> minimum
      }
      else
      {
        minimum = 4 * q0 > q2 * q2;  // H''(-q3/4) > 0 -> minimum
      }

      bool iterate = notZero || minimum;

      T_real a, b, c, d;
      if ( iterate )
      {
        T_real x = q3 >= 0 ? 2 : -2;  // initial root -> target = smaller mag root
        m_iter += zeroHexicByNewtonBisection( q3, q2, q1, q0, x );

        a = x * scale;    // 1st real component -> a
        b = -A3 / 2 - a;  // 2nd real component -> b

        x        = 4 * a + A3;  // Q'''(a)
        T_real y = x + 2 * A3;
        y        = y * a + 2 * A2;  // Q'(a)
        y        = y * a + A1;
        y /= x;
        if ( y < 0 ) y = 0;  // ensure Q'(a) / Q'''(a) >= 0

        x        = 4 * b + A3;  // Q'''(b)
        T_real z = x + 2 * A3;
        z        = z * b + 2 * A2;  // Q'(b)
        z        = z * b + A1;
        z /= x;
        if ( z < 0 ) z = 0;  // ensure Q'(b) / Q'''(b) >= 0

        c         = a * a;  // store a^2 for later
        d         = b * b;  // store b^2 for later
        T_real ss = c + y;  // magnitude^2 of (a + iy) root
        T_real tt = d + z;  // magnitude^2 of (b + iz) root

        // the arguments below are non negative in exact arithmetic but can turn
        // slightly negative because of rounding: clamp them as already done in
        // the non-iterative branch
        if ( ss > tt )
        {                                   // minimize imaginary error
          T_real const cc2{ y };            // 1st imaginary component -> c
          T_real const dd2{ A0 / ss - d };  // 2nd imaginary component -> d
          c = cc2 > 0 ? sqrt( cc2 ) : T_real( 0 );
          d = dd2 > 0 ? sqrt( dd2 ) : T_real( 0 );
        }
        else
        {
          T_real const cc2{ A0 / tt - c };  // 1st imaginary component -> c
          T_real const dd2{ z };            // 2nd imaginary component -> d
          c = cc2 > 0 ? sqrt( cc2 ) : T_real( 0 );
          d = dd2 > 0 ? sqrt( dd2 ) : T_real( 0 );
        }
      }
      else
      {  // no bisection -> real components equal

        a = -A3 / 4;  // 1st real component -> a
        b = a;        // 2nd real component -> b = a

        T_real x = a + A3;
        x        = x * a + A2;  // Q(a)
        x        = x * a + A1;
        x        = x * a + A0;
        T_real y = A2 / 2 - 3 * ( a * a );  // Q''(a) / 2
        T_real z = y * y - x;
        z        = z > 0 ? sqrt( z ) : 0;  // force discriminant to be >= 0
                                           // square root of discriminant
        y        = y > 0 ? y + z : y - z;  // larger magnitude root
        if ( y == 0 )
        {
          /*
          ..  Both coefficients of the resolvent vanish: with t = x - a the
          ..  quartic has degenerated to  Q(a+t) = t^4 + x,  x >= 0.  Its roots
          ..  lie on the diagonals through a,
          ..        a +- m +- i m,   m = ( x^(1/4) ) / sqrt(2)
          ..  and x == 0 gives the quadruple root a.  Without this branch the
          ..  division below was 0/0 and produced infinite imaginary parts, as
          ..  happens for every (x-r)^4 whose Q(a) rounds to a value >= 0.
          */
          T_real const m = x > 0 ? sqrt( sqrt( x ) / 2 ) : T_real( 0 );
          b              = a - m;  // 2nd real component
          a              = a + m;  // 1st real component
          c = d = m;               // both imaginary components
        }
        else
        {
          x /= y;                     // smaller magnitude root
          c = y < 0 ? 0 : sqrt( y );  // ensure root of biquadratic > 0
          d = x < 0 ? 0 : sqrt( x );  // large magnitude imaginary component
        }
      }

      m_ncplx = 4;
      if ( a > b )
      {
        m_r0 = a;
        m_r1 = c;
        m_r2 = b;
        m_r3 = d;
      }
      else if ( a < b )
      {
        m_r0 = b;
        m_r1 = d;
        m_r2 = a;
        m_r3 = c;
      }
      else
      {
        m_r0 = a;
        m_r1 = c;
        m_r2 = a;
        m_r3 = d;
      }
    }

    if ( rescaled )
    {
      // undo x = root_scale * y (root_scale is a positive power of two, so the
      // ordering established above is preserved)
      m_r0 *= root_scale;
      m_r1 *= root_scale;
      m_r2 *= root_scale;
      m_r3 *= root_scale;
    }
  }

  template <typename T_real, typename T_complex> void QuarticT<T_real, T_complex>::info( ostream_type & s ) const
  {
    T_real const & A = m_ABCDE[0];
    T_real const & B = m_ABCDE[1];
    T_real const & C = m_ABCDE[2];
    T_real const & D = m_ABCDE[3];
    T_real const & E = m_ABCDE[4];

    s << std::format(
      "poly a={} b={} c={} d={} e={}\n"
      "n. complex = {}\n"
      "n. real    = {}\n",
      A,
      B,
      C,
      D,
      E,
      m_ncplx,
      m_nreal );
    if ( m_ncplx > 0 )
      s << std::format(
        "x0 = ({},{})\n"
        "x1 = ({},{})\n",
        m_r0,
        m_r1,
        m_r0,
        -m_r1 );
    else
    {
      if ( m_nreal > 0 ) s << std::format( "x0 = {}\n", m_r0 );
      if ( m_nreal > 1 ) s << std::format( "x1 = {}\n", m_r1 );
    }
    if ( m_ncplx > 2 )
      s << std::format(
        "x2 = ({},{})\n"
        "x3 = ({},{})\n",
        m_r2,
        m_r3,
        m_r2,
        -m_r3 );
    else
    {
      if ( m_nreal > 2 || ( m_ncplx > 0 && m_nreal > 0 ) ) s << std::format( "x2 = {}\n", m_r2 );
      if ( m_nreal > 3 || ( m_ncplx > 0 && m_nreal > 1 ) ) s << std::format( "x3 = {}\n", m_r3 );
    }
  }

  template <typename T_real, typename T_complex> T_real QuarticT<T_real, T_complex>::backward_error() const
  {
    T_real        err{ 0 };
    integer const n{ num_roots() };
    for ( integer i{ 0 }; i < n; ++i )
    {
      T_real const e{ root_backward_error<T_real, T_complex>( m_ABCDE, 4, root( i ) ) };
      if ( e > err ) err = e;
    }
    return err;
  }

  template <typename T_real, typename T_complex> bool QuarticT<T_real, T_complex>::check( ostream_type & s ) const
  {
    integer const n{ num_roots() };
    if ( n <= 0 )
    {
      s << "no root computed\n";
      return true;
    }
    T_real const tol{ residual_tolerance<T_real>() };
    bool         ok{ true };
    for ( integer i{ 0 }; i < n; ++i )
    {
      T_complex const r{ root( i ) };
      T_real const    be{ root_backward_error<T_real, T_complex>( m_ABCDE, 4, r ) };
      s << std::format( "|p(r{})| = {}, backward error = {}\n", i, abs( eval( r ) ), be );
      ok = ok && be <= tol;
    }
    return ok;
  }

  template class QuarticT<real_type, real_complex>;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template class QuarticT<quad_real, quad_complex>;
#endif

}  // namespace PolynomialRoots

// EOF: PolynomialRoots-3-Quartic.cc
