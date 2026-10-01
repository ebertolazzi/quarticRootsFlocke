//
// Regression suite for the defects fixed in the solvers.
//
// Every case below reproduces a bug that was present in the library: each one
// states the expected behaviour, so a future refactoring that reintroduces the
// defect makes this executable fail.
//

#include "PolynomialRoots.hh"
#include "TestHelpers.hh"
#include "TestReporter.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

  using PolynomialRoots::Cubic;
  using PolynomialRoots::integer;
  using PolynomialRoots::Quadratic;
  using PolynomialRoots::Quartic;
  using PolynomialRoots::real_type;
  using PolynomialRoots::RootsError;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  using PolynomialRoots::QuarticHQ;
#endif

  // --------------------------------------------------------------------------
  // Relative backward error of x as a root of op[0..Degree]
  // --------------------------------------------------------------------------
  double backward_error( std::vector<double> const & op, double re, double im )
  {
    std::size_t const n{ op.size() - 1 };
    double            vr{ op[0] }, vi{ 0 };
    for ( std::size_t k = 1; k <= n; ++k )
    {
      double const nr{ vr * re - vi * im + op[k] };
      vi = vr * im + vi * re;
      vr = nr;
    }
    double const ax{ std::hypot( re, im ) };
    double       den{ 0 }, pw{ 1 };
    for ( std::size_t k = n + 1; k-- > 0; )
    {
      den += std::abs( op[k] ) * pw;
      pw *= ax;
    }
    if ( !( den > 0 ) ) return 0;
    return std::hypot( vr, vi ) / den;
  }

  template <typename Solver> double solver_backward_error( Solver const & s, std::vector<double> const & op )
  {
    double err{ 0 };
    for ( PolynomialRoots::integer i{ 0 }; i < s.num_roots(); ++i )
    {
      auto const r{ s.root( i ) };
      err = std::max( err, backward_error( op, r.real(), r.imag() ) );
    }
    return err;
  }

  // --------------------------------------------------------------------------
  void report( TestReporter::Summary & summary, int & id, char const * label, char const * source, bool ok, double v )
  {
    summary.case_header( id++, label, source );
    summary.comparison_table_header();
    summary.comparison_table_row( "double", ok, v );
    summary.comparison_table_footer();
    if ( ok )
      summary.pass();
    else
      summary.fail();
  }

  // ==========================================================================
  // 1) root_assert must THROW: it used to build a std::runtime_error and drop
  //    it on the floor, so every isfinite() guard was a no-op.
  // ==========================================================================
  bool case_root_assert_throws()
  {
    double const nan{ std::numeric_limits<double>::quiet_NaN() };
    double const inf{ std::numeric_limits<double>::infinity() };
    int          thrown{ 0 };

    auto expect_throw = [&thrown]( auto && fun )
    {
      try
      {
        fun();
      }
      catch ( std::runtime_error const & )
      {
        ++thrown;
      }
      catch ( ... )
      {
      }
    };

    expect_throw(
      [&]
      {
        Quadratic const q( nan, 1, 2 );
        (void) q;
      } );
    expect_throw(
      [&]
      {
        Cubic const c( 1, inf, 2, 3 );
        (void) c;
      } );
    expect_throw(
      [&]
      {
        Quartic const p( 1, 2, nan, 3, 4 );
        (void) p;
      } );
    expect_throw(
      [&]
      {
        Quadratic q;
        q.setup( 1, nan, 2 );
      } );
    expect_throw(
      [&]
      {
        Cubic c;
        c.setup( 1, 2, inf, 3 );
      } );
    expect_throw(
      [&]
      {
        Quartic p;
        p.setup( nan, 1, 2, 3, 4 );
      } );

    return thrown == 6;
  }

  // ==========================================================================
  // 2) The quartic used to return four completely wrong roots on roughly one
  //    random polynomial out of 5000: the Newton loop gave up after
  //    MAX_ITER_SAFETY iterations and the unconverged iterate was accepted
  //    because the bisection brackets had never been recorded.
  // ==========================================================================
  double case_quartic_known_failures()
  {
    std::vector<std::vector<double>> const bad{
      { 2.36261, -7.99568, 0.388084, -2.45634, 5.74222 },
      { 0.520492, -7.86652, -4.99027, -6.56348, 8.50988 },
      { -0.000115087, -6.21172, 3.41514, 0.28747, 4.12569 },
      { 2.76848, 1.75009, -3.77253, -6.75626, 5.11305 },
    };
    double worst{ 0 };
    for ( auto const & c : bad )
    {
      Quartic const p( c[0], c[1], c[2], c[3], c[4] );
      worst = std::max( worst, solver_backward_error( p, c ) );
    }
    return worst;
  }

  //! Randomized sweep: no quartic may exceed the backward error bound.
  double case_quartic_random_sweep( int n_samples )
  {
    std::mt19937_64                        gen{ 20260101u };
    std::uniform_real_distribution<double> U( -10, 10 );
    double                                 worst{ 0 };
    for ( int t{ 0 }; t < n_samples; ++t )
    {
      std::vector<double> c( 5 );
      for ( auto & v : c ) v = U( gen );
      if ( std::abs( c[0] ) < 1e-6 ) continue;
      Quartic const p( c[0], c[1], c[2], c[3], c[4] );
      worst = std::max( worst, solver_backward_error( p, c ) );
    }
    return worst;
  }

  // ==========================================================================
  // 3) get_real_roots() and friends must honour num_real_roots(): a degenerate
  //    quartic (a == 0) used to report a fourth, never computed, root.
  // ==========================================================================
  bool case_quartic_degenerate_getters()
  {
    bool ok{ true };
    {  // (x-1)(x-2)(x-3), leading coefficient zero
      Quartic const p( 0, 1, -6, 11, -6 );
      double        r[4]{ -9, -9, -9, -9 };
      ok = ok && p.num_real_roots() == 3 && p.num_complex_roots() == 0;
      ok = ok && p.num_roots() == 3;
      ok = ok && p.get_real_roots( r ) == 3;
      ok = ok && std::abs( r[0] - 1 ) < 1e-12 && std::abs( r[1] - 2 ) < 1e-12 && std::abs( r[2] - 3 ) < 1e-12;
      ok = ok && r[3] == -9;  // the fourth slot must be left untouched
      ok = ok && p.get_positive_roots( r ) == 3;
      ok = ok && p.get_negative_roots( r ) == 0;
      ok = ok && p.get_roots_in_range( 1.5, 3.5, r ) == 2;
      ok = ok && p.get_roots_in_open_range( 1.5, 3.5, r ) == 2;
      ok = ok && p.get_roots_in_open_range( 2.5, 3.5, r ) == 1;
    }
    {  // (x-1)(x^2+1), one real root plus a conjugate pair
      Quartic const p( 0, 1, -1, 1, -1 );
      double        r[4]{ -9, -9, -9, -9 };
      ok = ok && p.num_real_roots() == 1 && p.num_complex_roots() == 2;
      ok = ok && p.get_real_roots( r ) == 1 && std::abs( r[0] - 1 ) < 1e-14;
    }
    {  // identically zero polynomial: nothing has been computed
      Quartic const p( 0, 0, 0, 0, 0 );
      double        r[4]{ -9, -9, -9, -9 };
      ok = ok && p.num_roots() == 0 && p.get_real_roots( r ) == 0;
    }
    {  // regular case must keep working
      Quartic const p( 1, 0, -5, 0, 4 );
      double        r[4]{};
      ok = ok && p.get_real_roots( r ) == 4;
      ok = ok && r[0] == -2 && r[1] == -1 && r[2] == 1 && r[3] == 2;
    }
    return ok;
  }

  // ==========================================================================
  // 4) Quadratic::get_real_roots() used to emit one root even when none had
  //    been computed (a == b == 0).
  // ==========================================================================
  bool case_quadratic_getters()
  {
    bool ok{ true };
    {
      Quadratic const q( 0, 0, 1 );  // no root at all
      double          r[2]{ -9, -9 };
      ok = ok && q.num_roots() == 0 && q.get_real_roots( r ) == 0 && r[0] == -9;
    }
    {
      Quadratic const q( 0, 2, -4 );  // degenerates to 2x-4
      double          r[2]{ -9, -9 };
      ok = ok && q.num_roots() == 1 && q.get_real_roots( r ) == 1 && r[0] == 2;
    }
    {
      Quadratic const q( 1, 0, 1 );  // complex pair
      double          r[2]{ -9, -9 };
      ok = ok && q.num_real_roots() == 0 && q.num_complex_roots() == 2 && q.get_real_roots( r ) == 0;
    }
    {
      Quadratic const q( 1, -3, 2 );
      double          r[2]{};
      ok = ok && q.get_real_roots( r ) == 2 && r[0] == 1 && r[1] == 2;
      ok = ok && q.num_real_roots() == 2 && q.num_complex_roots() == 0;
    }
    return ok;
  }

  // ==========================================================================
  // 5) Newton_step() computed a wrong derivative in the |x| > 1 branch: it
  //    started the Horner recurrence from op[Degree] and used the weights
  //    (Degree-i) instead of (i+1).
  // ==========================================================================
  double case_newton_step()
  {
    struct Case
    {
      int    n;
      double op[5];
    };
    Case const cases[]{
      { 2, { 1, -3, 2, 0, 0 } },
      { 3, { 1, -6, 11, -6, 0 } },
      { 4, { 1, -10, 35, -50, 24 } },
      { 4, { 2, 1, -3, 0.5, -7 } },
    };
    double worst{ 0 };
    for ( auto const & c : cases )
    {
      for ( double x0 : { -3.5, -1.25, -0.4, 0.3, 0.9, 1.5, 2.5, 7.25 } )
      {
        // reference: plain Horner for p and p'
        double p{ c.op[0] }, dp{ 0 };
        for ( int i{ 1 }; i <= c.n; ++i )
        {
          dp = dp * x0 + p;
          p  = p * x0 + c.op[i];
        }
        if ( dp == 0 ) continue;
        double const expected{ x0 - p / dp };

        double x{ x0 };
        if ( !PolynomialRoots::Newton_step<double>( c.op, c.n, x ) ) return 1;
        worst = std::max( worst, std::abs( x - expected ) / ( 1 + std::abs( expected ) ) );
      }
    }
    return worst;
  }

  // ==========================================================================
  // 6) eval_poly_Dpoly() read op[1] out of bounds for Degree == 0 and returned
  //    dp = op[0] instead of 0.
  // ==========================================================================
  bool case_eval_poly_dpoly_degree0()
  {
    bool ok{ true };
    {
      double const op[1]{ 7 };  // a single element: op[1] must never be touched
      double       p{ -1 }, dp{ -1 };
      PolynomialRoots::eval_poly_Dpoly<double>( op, 0, 0.5, p, dp );
      ok = ok && p == 7 && dp == 0;
      PolynomialRoots::eval_poly_Dpoly<double>( op, 0, 25.0, p, dp );
      ok = ok && p == 7 && dp == 0;
      ok = ok && PolynomialRoots::eval_poly<double>( op, 0, 25.0 ) == 7;
      double x{ 3 };
      ok = ok && !PolynomialRoots::Newton_step<double>( op, 0, x ) && x == 3;
    }
    {
      double const op[2]{ 3, -6 };  // 3x - 6
      double       p{}, dp{};
      PolynomialRoots::eval_poly_Dpoly<double>( op, 1, 0.5, p, dp );
      ok = ok && p == -4.5 && dp == 3;
      PolynomialRoots::eval_poly_Dpoly<double>( op, 1, 4.0, p, dp );
      ok = ok && std::abs( p - 6 ) < 1e-14 && std::abs( dp - 3 ) < 1e-14;
    }
    return ok;
  }

  // ==========================================================================
  // 7) The smaller root of the quadratic was obtained as (C/r0)/A, where C/r0
  //    underflows to zero when |r0| is huge and |C| tiny.
  // ==========================================================================
  bool case_quadratic_underflow()
  {
    bool ok{ true };
    {
      Quadratic const q( 1e-200, 1.0, 1e-200 );
      ok = ok && std::abs( q.real_root0() + 1e200 ) <= 1e185;
      ok = ok && q.real_root1() != 0 && std::abs( q.real_root1() + 1e-200 ) <= 1e-215;
    }
    {
      Quadratic const q( 1e-180, -1.0, 1e-180 );
      ok = ok && q.real_root0() != 0 && q.real_root1() != 0;
      ok = ok && q.backward_error() <= 1e-14;
    }
    return ok;
  }

  // ==========================================================================
  // 8) The hexic solver had no iteration cap and used |dx| <= |x|*eps, which
  //    cannot be satisfied when x converges to zero.  These quartics have only
  //    complex roots of very small modulus and used to be able to hang.
  // ==========================================================================
  double case_hexic_termination()
  {
    double worst{ 0 };
    for ( double eps : { 1e-2, 1e-4, 1e-6, 1e-8, 1e-10 } )
    {
      // (x^2 + eps^2)(x^2 + 2 eps^2) : four purely imaginary, tiny roots
      std::vector<double> const c{ 1, 0, 3 * eps * eps, 0, 2 * eps * eps * eps * eps };
      Quartic const             p( c[0], c[1], c[2], c[3], c[4] );
      worst = std::max( worst, solver_backward_error( p, c ) );
    }
    return worst;
  }

  // ==========================================================================
  // 9) find_roots() did not reset the stored roots, so a degenerate setup()
  //    leaked the previous result through real_root1()/real_root2().
  // ==========================================================================
  bool case_state_reset()
  {
    bool ok{ true };
    {
      Cubic c( 1, -6, 11, -6 );  // roots 1, 2, 3
      c.setup( 0, 1, -3, 2 );    // degenerates to a quadratic
      ok = ok && c.num_roots() == 2 && c.real_root2() == 0 && c.root2().real() == 0;
    }
    {
      Quartic p( 1, 0, -5, 0, 4 );  // roots -2,-1,1,2
      p.setup( 0, 1, -6, 11, -6 );  // degenerates to a cubic
      double r[4]{ -9, -9, -9, -9 };
      ok = ok && p.get_real_roots( r ) == 3;
    }
    return ok;
  }

  // ==========================================================================
  // 10) root(i) / get_root(i) used to return a silent zero for an out of range
  //     index.
  // ==========================================================================
  bool case_root_index_checked()
  {
    int  thrown{ 0 };
    auto expect_throw = [&thrown]( auto && fun )
    {
      try
      {
        fun();
      }
      catch ( std::runtime_error const & )
      {
        ++thrown;
      }
      catch ( ... )
      {
      }
    };
    Quadratic const q( 1, -3, 2 );
    Cubic const     c( 1, -6, 11, -6 );
    Quartic const   p( 1, 0, -5, 0, 4 );
    expect_throw( [&] { (void) q.root( 2 ); } );
    expect_throw( [&] { (void) c.root( 3 ); } );
    expect_throw( [&] { (void) p.root( 4 ); } );
    expect_throw( [&] { (void) p.root( -1 ); } );
    double re{}, im{};
    expect_throw( [&] { q.get_root( 7, re, im ); } );
    // valid indices must not throw
    bool no_throw{ true };
    try
    {
      for ( int i{ 0 }; i < 4; ++i ) (void) p.root( i );
    }
    catch ( ... )
    {
      no_throw = false;
    }
    return thrown == 5 && no_throw;
  }

  // ==========================================================================
  // 11) check()/backward_error must not flag a perfectly computed root just
  //     because the polynomial is large: the old criterion compared |p(r)|
  //     with a tolerance scaled only by the coefficients.
  // ==========================================================================
  bool case_check_scale_invariance()
  {
    bool ok{ true };
    {
      Quadratic const q( 1, -2e8, 1e16 );  // double root at 1e8
      ok = ok && q.check() && q.backward_error() <= 1e-15;
    }
    {
      Cubic const c( 1, -3e6, 3e12, -1e18 );  // triple root at 1e6
      ok = ok && c.backward_error() <= 1e-10;
    }
    {
      Quartic const p( 1, 0, -5e-20, 0, 4e-40 );  // roots of order 1e-10
      ok = ok && p.check() && p.backward_error() <= 1e-15;
    }
    return ok;
  }

  // ==========================================================================
  // 12) Extremely unbalanced coefficients used to overflow the normalization
  //     and produce infinite roots.
  // ==========================================================================
  double case_unbalanced_coefficients()
  {
    std::vector<std::vector<double>> const cases{
      { 1e-300, 1, 1, 1, 1e300 },
      { 1e-300, 1, 1, 1, -1e300 },
      { 1e300, 1e300, 1e300, 1e300, 1e300 },
      { 1e-200, 1, 0, 0, 1e200 },
    };
    double worst{ 0 };
    for ( auto const & c : cases )
    {
      Quartic const p( c[0], c[1], c[2], c[3], c[4] );
      for ( PolynomialRoots::integer i{ 0 }; i < p.num_roots(); ++i )
      {
        auto const r{ p.root( i ) };
        if ( !std::isfinite( r.real() ) || !std::isfinite( r.imag() ) ) return 1;
      }
      worst = std::max( worst, solver_backward_error( p, c ) );
    }
    return worst;
  }

  // ==========================================================================
  // 13) Jenkins-Traub entry point: distinct error codes.
  // ==========================================================================
  bool case_jenkins_traub_error_codes()
  {
    bool                   ok{ true };
    std::vector<real_type> zr( PolynomialRoots::max_degree + 2, 0 );
    std::vector<real_type> zi( PolynomialRoots::max_degree + 2, 0 );
    {
      std::vector<real_type> const c{ 1, 0, -1 };
      ok = ok &&
           PolynomialRoots::jenkins_traub_roots( c.data(), 0, zr.data(), zi.data() ) == RootsError::DEGREE_TOO_SMALL;
      ok = ok && PolynomialRoots::jenkins_traub_roots( c.data(), 2, zr.data(), zi.data() ) == RootsError::SUCCESS;
    }
    {
      std::vector<real_type> const c{ 0, 1, -1 };
      ok = ok && PolynomialRoots::jenkins_traub_roots( c.data(), 2, zr.data(), zi.data() ) == RootsError::LEADING_ZERO;
    }
    {
      std::vector<real_type> c( PolynomialRoots::max_degree + 3, 0 );
      c.front() = 1;
      ok = ok &&
           PolynomialRoots::jenkins_traub_roots( c.data(), PolynomialRoots::max_degree + 1, zr.data(), zi.data() ) ==
             RootsError::DEGREE_TOO_LARGE;
    }
    // the legacy integer wrapper must agree
    {
      std::vector<real_type> const c{ 0, 1, -1 };
      ok = ok && PolynomialRoots::roots( c.data(), 2, zr.data(), zi.data() ) == -2;
    }
    ok = ok && std::string( PolynomialRoots::to_string( RootsError::SUCCESS ) ) == "success";
    return ok;
  }

  // ==========================================================================
  // 14) The deprecated camelCase spellings must keep forwarding to the new API.
  // ==========================================================================
#if defined( __GNUC__ ) || defined( __clang__ )
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#elif defined( _MSC_VER )
#pragma warning( push )
#pragma warning( disable : 4996 )
#endif
  bool case_deprecated_aliases()
  {
    bool ok{ true };
    {
      Quadratic const q( 1, -3, 2 );
      double          a[2]{}, b[2]{};
      ok = ok && q.numRoots() == q.num_roots();
      ok = ok && q.complexRoot() == q.complex_root();
      ok = ok && q.doubleRoot() == q.double_root();
      ok = ok && q.getRealRoots( a ) == q.get_real_roots( b ) && a[0] == b[0] && a[1] == b[1];
      ok = ok && q.getPositiveRoots( a ) == q.get_positive_roots( b );
      ok = ok && q.getNegativeRoots( a ) == q.get_negative_roots( b );
      ok = ok && q.getRootsInRange( 0, 1.5, a ) == q.get_roots_in_range( 0, 1.5, b );
      ok = ok && q.getRootsInOpenRange( 0, 1.5, a ) == q.get_roots_in_open_range( 0, 1.5, b );
      double r0{}, i0{}, r1{}, i1{};
      q.getRoot0( r0, i0 );
      q.get_root0( r1, i1 );
      ok = ok && r0 == r1 && i0 == i1;
    }
    {
      Cubic const c( 1, -6, 11, -6 );
      ok = ok && c.tripleRoot() == c.triple_root() && c.numRoots() == c.num_roots();
      double re{}, im{};
      c.getRoot2( re, im );
      ok = ok && re == c.real_root2();
    }
    {
      Quartic const p( 1, 0, -5, 0, 4 );
      ok = ok && p.numRealRoots() == p.num_real_roots();
      ok = ok && p.numComplexRoots() == p.num_complex_roots();
      double re{}, im{};
      p.getRoot3( re, im );
      ok = ok && re == p.real_root3();
    }
    {
      // free functions
      double const x{ 2.0 }, a{ 1.0 }, b{ 2.0 }, c{ 3.0 }, d{ 4.0 };
      ok = ok && PolynomialRoots::evalMonicCubic( x, a, b, c ) == PolynomialRoots::eval_monic_cubic( x, a, b, c );
      ok = ok &&
           PolynomialRoots::evalMonicQuartic( x, a, b, c, d ) == PolynomialRoots::eval_monic_quartic( x, a, b, c, d );
      ok = ok && PolynomialRoots::machepsiT<double>() == PolynomialRoots::machine_eps<double>();
      ok = ok && PolynomialRoots::toleranceT<double>() == PolynomialRoots::residual_tolerance<double>();
      ok = ok && PolynomialRoots::MAXDEGREE == PolynomialRoots::max_degree;
    }
    return ok;
  }
#if defined( __GNUC__ ) || defined( __clang__ )
#pragma GCC diagnostic pop
#elif defined( _MSC_VER )
#pragma warning( pop )
#endif

#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  // ==========================================================================
  // 15) The same fixes must hold for the multiprecision instantiation.
  // ==========================================================================
  double case_quartic_hq()
  {
    using quad_real = PolynomialRoots::quad_real;
    std::vector<std::vector<double>> const bad{
      { 2.36261, -7.99568, 0.388084, -2.45634, 5.74222 },
      { 0.520492, -7.86652, -4.99027, -6.56348, 8.50988 },
      { 2.76848, 1.75009, -3.77253, -6.75626, 5.11305 },
      { 1, 0, -5, 0, 4 },
    };
    double worst{ 0 };
    for ( auto const & c : bad )
    {
      quad_real const a0{ c[0] }, a1{ c[1] }, a2{ c[2] }, a3{ c[3] }, a4{ c[4] };
      QuarticHQ const p( a0, a1, a2, a3, a4 );
      worst = std::max( worst, static_cast<double>( p.backward_error() ) );
      // the double and the HQ solver must agree on the nature of the roots
      Quartic const d( c[0], c[1], c[2], c[3], c[4] );
      if ( d.num_real_roots() != p.num_real_roots() ) return 1;
      if ( d.num_complex_roots() != p.num_complex_roots() ) return 1;
    }
    {  // degenerate quartic: the getters must honour num_real_roots()
      quad_real const b0{ 0 }, b1{ 1 }, b2{ -6 }, b3{ 11 }, b4{ -6 };
      QuarticHQ const p( b0, b1, b2, b3, b4 );
      quad_real const minus9{ -9 };
      quad_real       r[4]{ minus9, minus9, minus9, minus9 };
      if ( p.get_real_roots( r ) != 3 ) return 1;
      if ( r[3] != minus9 ) return 1;
    }
    {  // multiplicity, exactly as for the double instantiation
      quad_real const c0{ 1 }, c1{ -4 }, c2{ 6 }, c3{ -4 }, c4{ 1 };
      QuarticHQ const p( c0, c1, c2, c3, c4 );  // (x-1)^4
      if ( !p.quadruple_root() || p.max_root_multiplicity() != 4 ) return 1;
      quad_real const d0{ 1 }, d1{ -8 }, d2{ 22 }, d3{ -24 }, d4{ 9 };
      QuarticHQ const q( d0, d1, d2, d3, d4 );  // (x-1)^2 (x-3)^2
      if ( !q.double_root() ) return 1;
    }
    {
      // The resolution of the multiplicity test follows the scalar type with no
      // tuning: two roots 1e-10 apart are below the double precision limit and
      // are reported as a double root, while the multiprecision instantiation,
      // whose epsilon is around 1e-100, keeps them apart.
      double const    r1{ 1.0 }, r2{ 1.0 + 1e-10 }, r3{ 4.0 }, r4{ -2.5 };
      double const    e1{ -( r1 + r2 + r3 + r4 ) };
      double const    e2{ r1 * r2 + r1 * r3 + r1 * r4 + r2 * r3 + r2 * r4 + r3 * r4 };
      double const    e3{ -( r1 * r2 * r3 + r1 * r2 * r4 + r1 * r3 * r4 + r2 * r3 * r4 ) };
      double const    e4{ r1 * r2 * r3 * r4 };
      Quartic const   d( 1, e1, e2, e3, e4 );
      quad_real const q0{ 1 }, q1{ e1 }, q2{ e2 }, q3{ e3 }, q4{ e4 };
      QuarticHQ const q( q0, q1, q2, q3, q4 );
      if ( !d.double_root() ) return 1;
      if ( q.double_root() ) return 1;
    }
    return worst;
  }
#endif

  // ==========================================================================
  // 16) Quartics with a real root of EVEN multiplicity used to be classified as
  //     having no real root at all: Q touches the axis without crossing it, so
  //     a plain sign scan does not see it.  Backward error was ~0.1.
  // ==========================================================================
  double case_quartic_touching_roots()
  {
    double       worst{ 0 };
    double const rr[]{ -0.764603, 2.9, 13.093, -1.60251, 1e-3, 7.5e2 };
    double const ss[]{ -953.224, 1.4, -3.25, 5.0, 11.0, -2.0 };
    for ( double r : rr )
    {
      {  // (x-r)^4 : the resolvent of the all-complex branch degenerates to 0/0
        std::vector<double> const c{ 1, -4 * r, 6 * r * r, -4 * r * r * r, r * r * r * r };
        Quartic const             p( c[0], c[1], c[2], c[3], c[4] );
        for ( PolynomialRoots::integer i{ 0 }; i < p.num_roots(); ++i )
        {
          auto const z{ p.root( i ) };
          if ( !std::isfinite( z.real() ) || !std::isfinite( z.imag() ) ) return 1;
        }
        if ( !p.quadruple_root() ) return 1;
        worst = std::max( worst, solver_backward_error( p, c ) );
      }
      for ( double s : ss )
      {  // (x-r)^2 (x-s)^2 : four real roots, no sign change anywhere
        if ( std::abs( r - s ) < 1e-3 ) continue;
        double const              e1{ -( 2 * r + 2 * s ) };
        double const              e2{ r * r + 4 * r * s + s * s };
        double const              e3{ -( 2 * r * s * s + 2 * r * r * s ) };
        double const              e4{ r * r * s * s };
        std::vector<double> const c{ 1, e1, e2, e3, e4 };
        Quartic const             p( c[0], c[1], c[2], c[3], c[4] );
        if ( p.max_root_multiplicity() != 2 ) return 1;
        worst = std::max( worst, solver_backward_error( p, c ) );
      }
    }
    return worst;
  }

  // ==========================================================================
  // 17) Multiplicity detection.  The tolerance is derived from the perturbation
  //     analysis, so no hand tuned constant is involved and the same code works
  //     for `double` and for the multiprecision instantiation.
  // ==========================================================================
  bool case_root_multiplicity()
  {
    bool ok{ true };

    {                                   // quadratic
      Quadratic const dbl( 1, -2, 1 );  // (x-1)^2
      ok = ok && dbl.double_root() && dbl.max_root_multiplicity() == 2;
      ok = ok && dbl.root_multiplicity( 0 ) == 2 && dbl.root_multiplicity( 1 ) == 2;
      Quadratic const simple( 1, -3, 2 );  // (x-1)(x-2)
      ok = ok && !simple.double_root() && simple.max_root_multiplicity() == 1;
      Quadratic const cplx( 1, 0, 1 );  // x^2+1
      ok = ok && !cplx.double_root();
      // roots separated well above the resolution limit must stay distinct
      Quadratic const close( 1, -2.0000001, 1.00000005 );
      ok = ok && !close.double_root();
    }
    {  // cubic: this is the case the old implementation got wrong, because it
      // only looked at the quadratic left by the deflation
      Cubic const a( 1, -4, 5, -2 );  // (x-1)^2 (x-2)
      ok = ok && a.double_root() && !a.triple_root() && a.max_root_multiplicity() == 2;
      Cubic const b( 1, -5, 8, -4 );  // (x-1)(x-2)^2
      ok = ok && b.double_root() && !b.triple_root();
      Cubic const c( 1, -3, 3, -1 );  // (x-1)^3
      ok = ok && c.triple_root() && !c.double_root() && c.max_root_multiplicity() == 3;
      Cubic const d( 1, -6, 11, -6 );  // (x-1)(x-2)(x-3)
      ok = ok && !d.double_root() && !d.triple_root() && d.max_root_multiplicity() == 1;
      Cubic const e( 1, -1, 1, -1 );  // (x-1)(x^2+1)
      ok = ok && !e.double_root() && !e.triple_root();
    }
    {                                      // quartic
      Quartic const a( 1, -4, 6, -4, 1 );  // (x-1)^4
      ok = ok && a.quadruple_root() && !a.triple_root() && !a.double_root();
      Quartic const b( 1, -5, 9, -7, 2 );  // (x-1)^3 (x-2)
      ok = ok && b.triple_root() && !b.quadruple_root() && b.max_root_multiplicity() == 3;
      Quartic const c( 1, -8, 22, -24, 9 );  // (x-1)^2 (x-3)^2
      ok = ok && c.double_root() && c.max_root_multiplicity() == 2;
      Quartic const d( 1, 0, 2, 0, 1 );  // (x^2+1)^2, a repeated CONJUGATE pair
      ok = ok && d.double_root() && d.num_complex_roots() == 4;
      Quartic const e( 1, -2, 2, -2, 1 );  // (x-1)^2 (x^2+1)
      ok = ok && e.double_root() && e.max_root_multiplicity() == 2;
      Quartic const f( 1, 0, -5, 0, 4 );  // (x^2-1)(x^2-4), four simple roots
      ok = ok && !f.double_root() && f.max_root_multiplicity() == 1;
    }
    {  // out of range index must throw
      Quadratic const q( 1, -3, 2 );
      bool            thrown{ false };
      try
      {
        (void) q.root_multiplicity( 2 );
      }
      catch ( std::runtime_error const & )
      {
        thrown = true;
      }
      ok = ok && thrown;
    }
    return ok;
  }

  // ==========================================================================
  // 18) eval_poly_complex: both branches |x| <= 1 and |x| > 1.
  // ==========================================================================
  bool case_eval_poly_complex()
  {
    using PolynomialRoots::real_complex;
    bool ok{ true };

    double const op[5]{ 2, -3, 0, 1, -5 };  // 2x^4 - 3x^3 + x - 5
    auto         reference = [&]( double re, double im )
    {
      double vr{ op[0] }, vi{ 0 };
      for ( int k{ 1 }; k < 5; ++k )
      {
        double const nr{ vr * re - vi * im + op[k] };
        vi = vr * im + vi * re;
        vr = nr;
      }
      return std::pair<double, double>{ vr, vi };
    };
    for ( auto [re, im] : { std::pair<double, double>{ 0.25, 0.5 },
                            { -0.75, 0.125 },
                            { 3.0, -2.0 },
                            { -4.5, 7.25 },
                            { 0.0, 0.0 },
                            { 1.0, 0.0 } } )
    {
      real_complex const z( re, im );
      real_complex const v{ PolynomialRoots::eval_poly_complex<double, real_complex>( op, 4, z ) };
      auto const [er, ei]{ reference( re, im ) };
      double const scale{ 1 + std::hypot( er, ei ) };
      ok = ok && std::abs( v.real() - er ) <= 1e-12 * scale && std::abs( v.imag() - ei ) <= 1e-12 * scale;
    }
    // a real root must evaluate to (almost) zero through the complex path too
    {
      Quartic const      p( 1, 0, -5, 0, 4 );
      real_complex const v{ p.eval( p.root( 0 ) ) };
      ok = ok && std::abs( v.real() ) < 1e-13 && std::abs( v.imag() ) < 1e-13;
    }
    // degree 0 and the empty polynomial must not read out of bounds
    {
      double const       k[1]{ 7 };
      real_complex const v{ PolynomialRoots::eval_poly_complex<double, real_complex>( k, 0, real_complex( 9, 9 ) ) };
      ok = ok && v.real() == 7 && v.imag() == 0;
    }
    return ok;
  }

  // ==========================================================================
  // 19) Degenerate quartic with a vanishing constant term (e == 0): x is a
  //     factor and the solver falls back on the cubic.
  // ==========================================================================
  bool case_quartic_zero_constant()
  {
    bool ok{ true };
    {  // x (x-1)(x-2)(x-3)
      Quartic const p( 1, -6, 11, -6, 0 );
      double        r[4]{};
      ok = ok && p.num_real_roots() == 4 && p.num_complex_roots() == 0;
      ok = ok && p.get_real_roots( r ) == 4;
      for ( int i{ 0 }; i < 4; ++i ) ok = ok && std::abs( r[i] - i ) < 1e-12;
      ok = ok && p.get_positive_roots( r ) == 3 && p.get_negative_roots( r ) == 0;
    }
    {  // x (x+2)(x^2+1) : one real root plus the pair, plus the root at zero
      Quartic const p( 1, 2, 1, 2, 0 );
      double        r[4]{ -9, -9, -9, -9 };
      ok = ok && p.num_complex_roots() == 2 && p.num_real_roots() == 2;
      ok = ok && p.get_real_roots( r ) == 2;
      ok = ok && std::abs( r[0] + 2 ) < 1e-12 && std::abs( r[1] ) < 1e-12;
      ok = ok && r[2] == -9;  // nothing written beyond the real roots
    }
    {  // x^2 (x-1)(x+1) : double root at the origin
      Quartic const p( 1, 0, -1, 0, 0 );
      ok = ok && p.num_real_roots() == 4 && p.double_root();
    }
    {  // both a == 0 and e == 0 : x (b x^2 + c x + d)
      Quartic const p( 0, 1, -3, 2, 0 );
      double        r[4]{ -9, -9, -9, -9 };
      ok = ok && p.num_real_roots() == 3 && p.get_real_roots( r ) == 3 && r[3] == -9;
    }
    return ok;
  }

  // ==========================================================================
  // 20) Cubic getters: the family where the Quadratic and the Quartic ones were
  //     broken, checked here for the remaining class.
  // ==========================================================================
  bool case_cubic_getters()
  {
    bool ok{ true };
    {  // (x+2)(x-1)(x-3)
      Cubic const c( 1, -2, -5, 6 );
      double      r[3]{ -9, -9, -9 };
      ok = ok && c.num_roots() == 3 && c.num_real_roots() == 3 && c.num_complex_roots() == 0;
      ok = ok && c.get_real_roots( r ) == 3;
      ok = ok && std::abs( r[0] + 2 ) < 1e-12 && std::abs( r[1] - 1 ) < 1e-12 && std::abs( r[2] - 3 ) < 1e-12;
      ok = ok && c.get_positive_roots( r ) == 2;
      ok = ok && c.get_negative_roots( r ) == 1;
      ok = ok && c.get_roots_in_range( -2.0, 1.0, r ) == 2;
      ok = ok && c.get_roots_in_open_range( -2.0, 1.0, r ) == 0;
      ok = ok && c.get_roots_in_range( 0.0, 10.0, r ) == 2;
    }
    {  // (x-1)(x^2+1) : only one real root is reported
      Cubic const c( 1, -1, 1, -1 );
      double      r[3]{ -9, -9, -9 };
      ok = ok && c.num_real_roots() == 1 && c.num_complex_roots() == 2;
      ok = ok && c.get_real_roots( r ) == 1 && std::abs( r[0] - 1 ) < 1e-12 && r[1] == -9;
      ok = ok && c.get_positive_roots( r ) == 1;
      ok = ok && c.get_negative_roots( r ) == 0;
    }
    {  // degenerate: a == 0 reduces to a quadratic
      Cubic const c( 0, 1, -3, 2 );
      double      r[3]{ -9, -9, -9 };
      ok = ok && c.num_roots() == 2 && c.get_real_roots( r ) == 2 && r[2] == -9;
    }
    {  // nothing to report
      Cubic const c( 0, 0, 0, 1 );
      double      r[3]{ -9, -9, -9 };
      ok = ok && c.num_roots() == 0 && c.get_real_roots( r ) == 0 && r[0] == -9;
    }
    return ok;
  }

  // ==========================================================================
  // 21) A leading coefficient many orders of magnitude below the others pushes
  //     one root out to 1/a and squeezes the other three down to the level of
  //     the machine epsilon ONCE Flocke's scaling is applied.  Two separate
  //     defects showed up there:
  //       - the stopping rules were written as (1+|x|)*eps, whose absolute
  //         floor of one epsilon is larger than the roots themselves;
  //       - the bracketing fallback kept the FIRST sign change it met, which is
  //         the hopeless one, instead of the best conditioned outermost root.
  // ==========================================================================
  double case_tiny_leading_coefficient()
  {
    std::vector<std::vector<double>> const cases{
      { 1.07256e-16, -2.50049, 3.51302, -0.504858, -9.08732 },
      { 1.24908e-19, -4.20306, -6.54116, -6.74844, -0.94316 },
      { 7.15249e-18, -3.81214, 4.60778, -5.13646, 8.38667 },
      { 2.94582e-08, -0.000473582, 5.60742, -4.2555, -3.05728 },
      { 1e-16, 1, 1, 1, 1 },
      { -1e-18, 2, -3, 5, -7 },
    };
    double worst{ 0 };
    for ( auto const & c : cases )
    {
      Quartic const p( c[0], c[1], c[2], c[3], c[4] );
      for ( PolynomialRoots::integer i{ 0 }; i < p.num_roots(); ++i )
      {
        auto const z{ p.root( i ) };
        if ( !std::isfinite( z.real() ) || !std::isfinite( z.imag() ) ) return 1;
      }
      worst = std::max( worst, solver_backward_error( p, c ) );
    }
    // randomized sweep over the same regime
    std::mt19937_64                        gen{ 90210u };
    std::uniform_real_distribution<double> U( -10, 10 );
    for ( int t{ 0 }; t < 200000; ++t )
    {
      double const a{ U( gen ) * std::pow( 10.0, -std::abs( U( gen ) ) * 2 ) };
      if ( a == 0 ) continue;
      std::vector<double> const c{ a, U( gen ), U( gen ), U( gen ), U( gen ) };
      Quartic const             p( c[0], c[1], c[2], c[3], c[4] );
      worst = std::max( worst, solver_backward_error( p, c ) );
    }
    return worst;
  }

  // ==========================================================================
  // 22) Cubic deflation.  Choosing between the forward and the backward formula
  //     on |r2| < 1 alone makes b1 = -(cc/r2 + bb)/r2 cancel exactly for
  //     (x-R)^2 (x-1) with R large: b1 collapsed to zero and the real double
  //     root came back out as a purely imaginary pair.
  // ==========================================================================
  double case_cubic_deflation()
  {
    double worst{ 0 };
    for ( double R : { 1e4, 1e8, 1e12, 1e15, 1e20, 1e25, 1e30, -1e18, 1e6 } )
    {
      for ( double s : { 1.0, -1.0, 1e-6, 7.0 } )
      {  // (x-R)^2 (x-s)
        std::vector<double> const c{ 1, -( 2 * R + s ), R * R + 2 * R * s, -R * R * s };
        Cubic const               p( c[0], c[1], c[2], c[3] );
        // The double root may legitimately come out as a conjugate pair whose
        // imaginary part sits at the sqrt(eps) level: one ulp on the
        // coefficients is enough to push it off the real axis.  What must hold
        // is that it is recognised as a double root, with the right location,
        // and not turned into the purely imaginary pair the old deflation gave.
        if ( !p.double_root() ) return 1;
        integer n_double{ 0 };
        for ( integer i{ 0 }; i < p.num_roots(); ++i )
        {
          if ( p.root_multiplicity( i ) != 2 ) continue;
          ++n_double;
          auto const z{ p.root( i ) };
          if ( std::abs( z.real() - R ) > 1e-6 * std::abs( R ) ) return 1;
        }
        if ( n_double != 2 ) return 1;
        worst = std::max( worst, solver_backward_error( p, c ) );
      }
    }
    return worst;
  }

  // ==========================================================================
  // 23) Composite deflation of the quartic.  The cross-over weights
  //     |a_k| |r|^k were formed as plain products, so for a root of modulus
  //     1e154 all five overflowed to infinity, compared equal, and the first
  //     (wrong) branch was selected.  Orellana & De Michele 2020, case 22.
  // ==========================================================================
  double case_quartic_deflation_overflow()
  {
    // roots 1e154, 1e152, 10, 1 -> the two small roots used to come out as
    // -1.7e139 and -5.8e-139
    std::vector<double> const c{ 1.0, -1.0100000000000001e154, 1e306, -1.1e307, 9.9999999999999999e306 };
    Quartic const             p( c[0], c[1], c[2], c[3], c[4] );
    double                    r[4]{};
    integer const             nr{ p.get_real_roots( r ) };
    if ( nr != 4 ) return 1;
    if ( std::abs( r[0] - 1 ) > 1e-12 ) return 1;
    if ( std::abs( r[1] - 10 ) > 1e-11 ) return 1;
    return solver_backward_error( p, c );
  }

}  // namespace

int main()
{
  TestReporter::Summary summary( std::cout, "Regression suite for the fixed defects" );

  int id{ 1 };

  report( summary, id, "root_assert throws", "non finite coefficients must be rejected", case_root_assert_throws(), 0 );

  {
    double const be{ case_quartic_known_failures() };
    report( summary, id, "quartic: known catastrophic cases", "unconverged Newton accepted", be <= 1e-13, be );
  }
  {
    double const be{ case_quartic_random_sweep( 200000 ) };
    report( summary, id, "quartic: 200k random polynomials", "max relative backward error", be <= 1e-13, be );
  }

  report(
    summary,
    id,
    "quartic: degenerate getters",
    "get_*_roots must honour num_real_roots()",
    case_quartic_degenerate_getters(),
    0 );

  report( summary, id, "quadratic: getters", "get_real_roots must honour num_roots()", case_quadratic_getters(), 0 );

  {
    double const err{ case_newton_step() };
    report( summary, id, "Newton_step", "correct derivative for |x| > 1 too", err <= 1e-14, err );
  }

  report(
    summary,
    id,
    "eval_poly_Dpoly: Degree 0 and 1",
    "no out of bounds read, dp = 0",
    case_eval_poly_dpoly_degree0(),
    0 );

  report(
    summary,
    id,
    "quadratic: smaller root underflow",
    "(C/A)/r0 instead of (C/r0)/A",
    case_quadratic_underflow(),
    0 );

  {
    double const be{ case_hexic_termination() };
    report( summary, id, "hexic solver: termination", "roots clustered at the origin", be <= 1e-12, be );
  }

  report( summary, id, "solver state reset", "setup() must not leak the previous roots", case_state_reset(), 0 );

  report( summary, id, "root(i): index checked", "out of range index must throw", case_root_index_checked(), 0 );

  report(
    summary,
    id,
    "check(): scale invariance",
    "relative backward error, not bare residual",
    case_check_scale_invariance(),
    0 );

  {
    double const be{ case_unbalanced_coefficients() };
    report( summary, id, "quartic: unbalanced coefficients", "|e|/|a| beyond DBL_MAX", be <= 1e-13, be );
  }

  report(
    summary,
    id,
    "Jenkins-Traub: error codes",
    "distinct codes for each failure",
    case_jenkins_traub_error_codes(),
    0 );

  report( summary, id, "deprecated aliases", "camelCase spellings still forward", case_deprecated_aliases(), 0 );

  {
    double const be{ case_quartic_touching_roots() };
    report( summary, id, "quartic: even multiplicity real roots", "Q touches the axis", be <= 1e-13, be );
  }

  report( summary, id, "root multiplicity", "tolerance from the perturbation analysis", case_root_multiplicity(), 0 );

  report( summary, id, "eval_poly_complex", "both |x| <= 1 and |x| > 1", case_eval_poly_complex(), 0 );

  report(
    summary,
    id,
    "quartic: vanishing constant term",
    "e == 0 degenerate branch",
    case_quartic_zero_constant(),
    0 );

  report( summary, id, "cubic getters", "get_*_roots on the remaining class", case_cubic_getters(), 0 );

  {
    double const be{ case_tiny_leading_coefficient() };
    report( summary, id, "quartic: tiny leading coefficient", "roots at the epsilon level", be <= 1e-12, be );
  }

  {
    double const be{ case_cubic_deflation() };
    report( summary, id, "cubic: deflation by a large root", "(x-R)^2 (x-s), R up to 1e30", be <= 1e-13, be );
  }

  {
    double const be{ case_quartic_deflation_overflow() };
    report( summary, id, "quartic: cross-over overflow", "Orellana & De Michele 2020, case 22", be <= 1e-13, be );
  }

#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  {
    double const be{ case_quartic_hq() };
    report( summary, id, "QuarticHQ: same fixes", "multiprecision instantiation", be <= 1e-30, be );
  }
#else
  summary.case_header( id++, "QuarticHQ: same fixes", "multiprecision instantiation" );
  summary.note( "multiprecision backend disabled, case skipped" );
#endif

  return summary.finish();
}
