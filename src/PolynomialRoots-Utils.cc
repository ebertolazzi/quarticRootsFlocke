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

//! Largest degree handled by the multiplicity helpers below.
#define QRF_MAX_DEGREE 7

namespace PolynomialRoots
{

  template <> real_type machine_eps<real_type>()
  { return std::numeric_limits<real_type>::epsilon(); }
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template <> quad_real machine_eps<quad_real>()
  { return std::numeric_limits<quad_real>::epsilon(); }
#endif
  template <> real_type residual_tolerance<real_type>()
  { return 1e-12; }
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template <> quad_real residual_tolerance<quad_real>()
  { return quad_real( "1e-20" ); }
#endif

  using std::abs;
  using std::cbrt;
  using std::pow;
  using std::sqrt;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of polinomial
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // (op[0] + op[1]/x.... + op[n]/x^n)*x^n
  //
  template <typename T_real> T_real eval_poly( T_real const op[], integer const Degree, T_real const & x )
  {
    if ( Degree < 0 ) return T_real( 0 );
    if ( abs( x ) > 1 )
    {
      T_real res = op[Degree];
      T_real xn  = 1;
      for ( integer i = 1; i <= Degree; ++i )
      {
        res /= x;
        res += op[Degree - i];
        xn *= x;
      }
      res *= xn;
      return res;
    }
    T_real res( op[0] );
    for ( integer i = 1; i <= Degree; ++i )
    {
      res *= x;
      res += op[i];
    }
    return res;
  }

  template double eval_poly( double const op[], integer const Degree, double const & x );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template quad_real eval_poly( quad_real const op[], integer const Degree, quad_real const & x );
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of polinomial and its derivative
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // (op[0] + op[1]/x.... + op[n]/x^n)*x^n
  //
  template <typename T_real>
  void eval_poly_Dpoly( T_real const op[], integer const Degree, T_real const & x, T_real & p, T_real & dp )
  {
    // Degree 0 (and the ill-formed negative case) must be handled apart: the
    // Horner loops below would otherwise read op[1] out of bounds.
    if ( Degree < 0 )
    {
      p = dp = 0;
      return;
    }
    if ( Degree == 0 )
    {
      p  = op[0];
      dp = 0;
      return;
    }
    if ( abs( x ) > 1 )
    {
      T_real const   y  = T_real( 1 ) / x;
      T_real const * pc = op + Degree;

      p  = *pc--;
      dp = 0;
      for ( integer i{ 1 }; i <= Degree; ++i )
      {
        dp = dp * y + p;
        p  = p * y + *pc--;
      }

      T_real xn = 1;
      for ( integer i = 0; i < Degree; ++i ) xn *= x;

      dp = static_cast<T_real>( Degree ) * p * ( xn / x ) - dp * ( xn / ( x * x ) );
      p *= xn;
    }
    else
    {
      // Horner for p(x) and, in lockstep, for p'(x):
      //   b_0 = op[0], b_i = b_{i-1} * x + op[i]   -->  p (x) = b_n
      //   c_0 = 0    , c_i = c_{i-1} * x + b_{i-1} -->  p'(x) = c_n
      p  = op[0];
      dp = 0;
      for ( integer i = 1; i <= Degree; ++i )
      {
        dp = dp * x + p;
        p  = p * x + op[i];
      }
    }
  }

  template void eval_poly_Dpoly( double const op[], integer const Degree, double const & x, double & p, double & dp );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template void eval_poly_Dpoly(
    quad_real const   op[],
    integer const     Degree,
    quad_real const & x,
    quad_real &       p,
    quad_real &       dp );
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // stable computation of a Newton step
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  //
  // Both p(x) and p'(x) come from `eval_poly_Dpoly`, which already switches to
  // the reciprocal Horner scheme when |x| > 1.  The step is refused, and x left
  // untouched, when the derivative vanishes or the correction overflows.
  //
  template <typename T_real> bool Newton_step( T_real const op[], integer const Degree, T_real & x )
  {
    if ( Degree < 1 ) return false;
    T_real p, dp;
    eval_poly_Dpoly( op, Degree, x, p, dp );
    if ( dp == 0 ) return false;
    T_real const dx{ p / dp };
    if ( !isfinite( dx ) ) return false;
    x -= dx;
    return true;
  }

  template bool Newton_step( double const op[], integer const Degree, double & x );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template bool Newton_step( quad_real const op[], integer const Degree, quad_real & x );
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // op[0] * x^n + .... + op[n-1]*x + op[n]
  template <typename T_real, typename T_complex>
  T_complex eval_poly_complex( T_real const op[], integer const Degree, T_complex const & x )
  {
    if ( Degree < 0 ) return T_complex( 0, 0 );
    if ( abs( x ) > 1 )
    {
      T_complex res( op[Degree], 0 );
      T_complex xn( 1, 0 );
      for ( integer i = 1; i <= Degree; ++i )
      {
        res /= x;
        res += op[Degree - i];
        xn *= x;
      }
      res *= xn;
      return res;
    }
    T_complex res( op[0], 0 );
    for ( integer i = 1; i <= Degree; ++i )
    {
      res *= x;
      res += op[i];
    }
    return res;
  }

  template real_complex eval_poly_complex( real_type const op[], integer const Degree, real_complex const & x );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template quad_complex eval_poly_complex( quad_real const op[], integer const Degree, quad_complex const & x );
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // Relative backward error  |p(x)| / sum_k |op[k]| |x|^(n-k).
  //
  // When |x| > 1 both the numerator and the denominator are evaluated already
  // divided by |x|^n, so the ratio is obtained without any intermediate
  // overflow even for roots of huge modulus.
  //
  template <typename T_real, typename T_complex>
  T_real root_backward_error( T_real const op[], integer const Degree, T_complex const & x )
  {
    if ( Degree < 0 ) return T_real( 0 );

    T_real const ax{ abs( x ) };
    T_real       den{ 0 };
    T_real       pw{ 1 };
    T_complex    val( op[0], 0 );

    if ( ax > 1 )
    {
      // reversed polynomial evaluated at y = 1/x, that is p(x) / x^n
      val = T_complex( op[Degree], 0 );
      for ( integer i = 1; i <= Degree; ++i )
      {
        val /= x;
        val += op[Degree - i];
      }
      // sum_k |op[k]| |x|^(-k)
      for ( integer k = 0; k <= Degree; ++k )
      {
        den += abs( op[k] ) * pw;
        pw /= ax;
      }
    }
    else
    {
      for ( integer i = 1; i <= Degree; ++i )
      {
        val *= x;
        val += op[i];
      }
      // sum_k |op[k]| |x|^(n-k)
      for ( integer k = Degree; k >= 0; --k )
      {
        den += abs( op[k] ) * pw;
        pw *= ax;
      }
    }

    if ( !( den > 0 ) ) return T_real( 0 );
    return abs( val ) / den;
  }

  template real_type root_backward_error( real_type const op[], integer const Degree, real_complex const & x );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template quad_real root_backward_error( quad_real const op[], integer const Degree, quad_complex const & x );
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // Taylor coefficients of op at c, by repeated synthetic division:
  // after the j-th pass the remainder is p^(j)(c)/j!
  //
  template <typename T_value, typename T_point>
  static void taylor_shift( T_value const op[], integer const Degree, T_point const & c, T_value out[] )
  {
    T_value t[QRF_MAX_DEGREE + 1];
    for ( integer j{ 0 }; j <= Degree; ++j ) t[j] = op[j];
    for ( integer j{ 0 }; j <= Degree; ++j )
    {
      for ( integer i{ 1 }; i <= Degree - j; ++i ) t[i] += c * t[i - 1];
      out[j] = t[Degree - j];
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // A root of multiplicity m at c satisfies p(c) = p'(c) = ... = p^(m-1)(c) = 0.
  // Each of those Taylor coefficients is compared with the one of the polynomial
  // built from the absolute values of the coefficients, scaled by the unit
  // roundoff: that is exactly the size a relative perturbation eps of the data
  // can produce in it.
  //
  //
  // j = 0 condition alone: |p(c)| <= K eps A_0(|c|).  It is necessary for a
  // root of ANY multiplicity >= 1 at c, costs two Horner evaluations instead of
  // two full Taylor shifts, and lets the common case - no multiple root at all -
  // be rejected without ever building the shifts.
  //
  template <typename T_real, typename T_complex>
  static bool residual_vanishes( T_real const op[], integer const Degree, T_complex const & centre )
  {
    T_complex    v( op[0], 0 );
    T_real       a{ abs( op[0] ) };
    T_real const ac{ abs( centre ) };
    for ( integer j{ 1 }; j <= Degree; ++j )
    {
      v = v * centre + op[j];
      a = a * ac + abs( op[j] );
    }
    return abs( v ) <= T_real( multiple_root_safety_factor ) * machine_eps<T_real>() * a;
  }

  template <typename T_real, typename T_complex>
  bool is_multiple_root( T_real const op[], integer const Degree, integer const multiplicity, T_complex const & centre )
  {
    if ( multiplicity < 2 || Degree < multiplicity || Degree > QRF_MAX_DEGREE ) return false;
    if ( !residual_vanishes<T_real, T_complex>( op, Degree, centre ) ) return false;

    T_complex tay[QRF_MAX_DEGREE + 1];
    T_complex cop[QRF_MAX_DEGREE + 1];
    T_real    abs_op[QRF_MAX_DEGREE + 1];
    T_real    abs_tay[QRF_MAX_DEGREE + 1];
    for ( integer j{ 0 }; j <= Degree; ++j )
    {
      cop[j]    = T_complex( op[j], 0 );
      abs_op[j] = abs( op[j] );
    }
    taylor_shift( cop, Degree, centre, tay );
    taylor_shift( abs_op, Degree, abs( centre ), abs_tay );

    T_real const k_eps{ T_real( multiple_root_safety_factor ) * machine_eps<T_real>() };
    for ( integer j{ 0 }; j < multiplicity; ++j )
      if ( !( abs( tay[j] ) <= k_eps * abs_tay[j] ) ) return false;
    return true;
  }

  template bool is_multiple_root<real_type, real_complex>(
    real_type const      op[],
    integer const        Degree,
    integer const        multiplicity,
    real_complex const & centre );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template bool is_multiple_root<quad_real, quad_complex>(
    quad_real const      op[],
    integer const        Degree,
    integer const        multiplicity,
    quad_complex const & centre );
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  //   delta ~ ( eps * S(c) / |p^(m)(c)|/m! )^(1/m),   S(c) = sum_j |a_j| |c|^(n-j)
  //
  template <typename T_real, typename T_complex> T_real multiple_root_separation(
    T_real const      op[],
    integer const     Degree,
    integer const     multiplicity,
    T_complex const & centre )
  {
    if ( multiplicity < 2 || Degree < multiplicity || Degree > QRF_MAX_DEGREE ) return T_real( 0 );

    T_real const ac{ abs( centre ) };
    T_real       S{ 0 };
    T_real       pw{ 1 };
    for ( integer j{ Degree }; j >= 0; --j )
    {
      S += abs( op[j] ) * pw;
      pw *= ac;
    }

    T_complex tay[QRF_MAX_DEGREE + 1];
    T_complex cop[QRF_MAX_DEGREE + 1];
    for ( integer j{ 0 }; j <= Degree; ++j ) cop[j] = T_complex( op[j], 0 );
    taylor_shift( cop, Degree, centre, tay );

    T_real const Dk{ abs( tay[multiplicity] ) };
    if ( !( Dk > 0 ) ) return std::numeric_limits<T_real>::max();

    T_real const ratio{ machine_eps<T_real>() * S / Dk };
    if ( !( ratio > 0 ) ) return T_real( 0 );

    switch ( multiplicity )
    {
      case 2: return sqrt( ratio );
      case 3: return cbrt( ratio );
      case 4: return sqrt( sqrt( ratio ) );
      default: return pow( ratio, T_real( 1 ) / T_real( multiplicity ) );
    }
  }

  template real_type multiple_root_separation<real_type, real_complex>(
    real_type const      op[],
    integer const        Degree,
    integer const        multiplicity,
    real_complex const & centre );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template quad_real multiple_root_separation<quad_real, quad_complex>(
    quad_real const      op[],
    integer const        Degree,
    integer const        multiplicity,
    quad_complex const & centre );
#endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //
  // Multiplicity of each computed root.
  //
  // The search goes TOP DOWN: the largest group of mutually close roots is
  // looked for first, because the admissible separation grows with the
  // multiplicity.  A bottom-up merge would be wrong: the four roots of (x-r)^4
  // spread over eps^(1/4), far beyond the eps^(1/2) allowed for a double root,
  // so no pair of them would ever be merged.
  //
  // With at most four roots the subsets are enumerated exhaustively; among the
  // candidates of a given size the one with the smallest diameter is taken, so
  // the result does not depend on the order of the input.
  //
  template <typename T_real, typename T_complex> void root_multiplicities(
    T_real const    op[],
    integer const   Degree,
    T_complex const roots[],
    integer const   n,
    integer         mult[] )
  {
    if ( n <= 0 ) return;

    constexpr integer MAX_N{ 4 };
    integer const     nn{ n < MAX_N ? n : MAX_N };
    if ( Degree > QRF_MAX_DEGREE ) return;

    for ( integer i{ 0 }; i < n; ++i ) mult[i] = 1;

    /*
    ..  Cheap rejection of the common case.  Any cluster contains a pair, and the
    ..  midpoint of a pair inside a cluster is itself a point where p vanishes to
    ..  the working precision.  So if no pair midpoint passes the j = 0 test
    ..  there is no multiple root at all, and the subset enumeration - which
    ..  costs two Taylor shifts per candidate - is skipped entirely.
    */
    {
      bool possible{ false };
      for ( integer i{ 0 }; i < nn && !possible; ++i )
      {
        for ( integer j{ i + 1 }; j < nn && !possible; ++j )
        {
          T_complex centre{ roots[i] };
          centre += roots[j];
          centre /= T_real( 2 );
          possible = residual_vanishes<T_real, T_complex>( op, Degree, centre );
        }
      }
      if ( !possible ) return;
    }

    unsigned pending{ ( 1u << static_cast<unsigned>( nn ) ) - 1u };  // unassigned roots

    while ( pending != 0 )
    {
      integer left{ 0 };
      for ( integer i{ 0 }; i < nn; ++i )
        if ( ( pending >> static_cast<unsigned>( i ) ) & 1u ) ++left;

      unsigned best_set{ 0 };
      integer  best_k{ 0 };
      T_real   best_diam{ 0 };

      for ( integer k{ left }; k >= 2 && best_k == 0; --k )
      {
        for ( unsigned sub{ 1 }; sub < ( 1u << static_cast<unsigned>( nn ) ); ++sub )
        {
          if ( ( sub & ~pending ) != 0 ) continue;  // only unassigned roots
          integer card{ 0 };
          for ( integer i{ 0 }; i < nn; ++i )
            if ( ( sub >> static_cast<unsigned>( i ) ) & 1u ) ++card;
          if ( card != k ) continue;

          T_complex centre( 0, 0 );
          T_real    diam{ 0 };
          for ( integer i{ 0 }; i < nn; ++i )
          {
            if ( !( ( sub >> static_cast<unsigned>( i ) ) & 1u ) ) continue;
            centre += roots[i];
            for ( integer j{ i + 1 }; j < nn; ++j )
            {
              if ( !( ( sub >> static_cast<unsigned>( j ) ) & 1u ) ) continue;
              T_real const d{ abs( roots[i] - roots[j] ) };
              if ( d > diam ) diam = d;
            }
          }
          centre /= T_real( k );

          if ( !is_multiple_root<T_real, T_complex>( op, Degree, k, centre ) ) continue;
          if ( best_k == 0 || diam < best_diam )
          {
            best_set  = sub;
            best_k    = k;
            best_diam = diam;
          }
        }
      }

      if ( best_k == 0 )
      {
        // no cluster left: the lowest unassigned root is simple
        for ( integer i{ 0 }; i < nn; ++i )
        {
          if ( ( pending >> static_cast<unsigned>( i ) ) & 1u )
          {
            mult[i] = 1;
            pending &= ~( 1u << static_cast<unsigned>( i ) );
            break;
          }
        }
      }
      else
      {
        for ( integer i{ 0 }; i < nn; ++i )
          if ( ( best_set >> static_cast<unsigned>( i ) ) & 1u ) mult[i] = best_k;
        pending &= ~best_set;
      }
    }
  }

  template void root_multiplicities<real_type, real_complex>(
    real_type const    op[],
    integer const      Degree,
    real_complex const roots[],
    integer const      n,
    integer            mult[] );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template void root_multiplicities<quad_real, quad_complex>(
    quad_real const    op[],
    integer const      Degree,
    quad_complex const roots[],
    integer const      n,
    integer            mult[] );
#endif

}  // namespace PolynomialRoots

// EOF: PolynomialRoots-Utils.cc
