#include "PolynomialRoots.hh"
#include "TestReporter.hh"
#include "TestHelpers.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <utility>
#include <array>
#include <tuple>

namespace
{

  using PolynomialRoots::Cubic;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  using PolynomialRoots::CubicHQ;
  using PolynomialRoots::MAXDEGREE;
  using PolynomialRoots::quad_real;
#else
  using PolynomialRoots::MAXDEGREE;
#endif
  using PolynomialRoots::Quadratic;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  using PolynomialRoots::QuadraticHQ;
#endif
  using PolynomialRoots::Quartic;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  using PolynomialRoots::QuarticHQ;
#endif
  using PolynomialRoots::real_type;

  // -------------------------------------------------------------------------
  // Trait per ottenere il tipo scalare usato da un solutore
  // -------------------------------------------------------------------------
  template <typename Solver> struct solver_scalar_type;
  template <> struct solver_scalar_type<Quadratic>
  {
    using type = real_type;
  };
  template <> struct solver_scalar_type<Cubic>
  {
    using type = real_type;
  };
  template <> struct solver_scalar_type<Quartic>
  {
    using type = real_type;
  };
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  template <> struct solver_scalar_type<QuadraticHQ>
  {
    using type = quad_real;
  };
  template <> struct solver_scalar_type<CubicHQ>
  {
    using type = quad_real;
  };
  template <> struct solver_scalar_type<QuarticHQ>
  {
    using type = quad_real;
  };
#endif

  template <typename Solver> using solver_scalar_t = typename solver_scalar_type<Solver>::type;

  // -------------------------------------------------------------------------
  // Confronto generico per qualsiasi tipo numerico
  // -------------------------------------------------------------------------
  template <typename T> bool close_enough( T a, T b )
  {
    using std::abs;
    using std::max;
    T const scale = max( T( 1 ), max( abs( a ), abs( b ) ) );
    return abs( a - b ) <= 100 * std::numeric_limits<T>::epsilon() * scale;
  }

  // -------------------------------------------------------------------------
  // Struttura per un singolo test: x, valore atteso p(x), valore atteso p'(x)
  // -------------------------------------------------------------------------
  template <typename Real> struct EvalTest
  {
    Real x;
    Real expected_p;
    Real expected_dp;
  };

  // -------------------------------------------------------------------------
  // Esecuzione dei test per un dato solutore e una lista di test
  // -------------------------------------------------------------------------
  template <typename Solver, typename CoeffArray, typename TestList>
  std::pair<bool, double> run_eval_case( const CoeffArray & coeffs, const TestList & tests )
  {
    using Real = solver_scalar_t<Solver>;

    // Costruisce il solutore
    Solver solver( std::apply( []( auto... args ) { return Solver( args... ); }, coeffs ) );

    bool   ok      = true;
    double max_err = 0.0;

    for ( const auto & test : tests )
    {
      Real p{ 0 }, dp{ 0 };
      solver.eval( test.x, p, dp );

      double err_p  = std::abs( static_cast<double>( p - test.expected_p ) );
      double err_dp = std::abs( static_cast<double>( dp - test.expected_dp ) );
      max_err       = std::max( max_err, err_p );
      max_err       = std::max( max_err, err_dp );

      ok = ok && close_enough( p, test.expected_p );
      ok = ok && close_enough( dp, test.expected_dp );
    }

    return { ok, max_err };
  }

  // Helper per convertire una lista di tuple in una sequenza di test
  template <typename Real> auto make_tests( const std::vector<std::tuple<Real, Real, Real>> & data )
  {
    std::vector<EvalTest<Real>> tests;
    tests.reserve( data.size() );
    for ( const auto & [x, ep, edp] : data ) tests.push_back( { x, ep, edp } );
    return tests;
  }

  // -------------------------------------------------------------------------
  // Casi di test specifici
  // -------------------------------------------------------------------------
  template <typename Solver> std::pair<bool, double> run_quadratic_eval_case()
  {
    using Real = solver_scalar_t<Solver>;
    // coefficienti per x^2 - 3x + 2
    std::array<Real, 3> coeffs = { Real( 1 ), Real( -3 ), Real( 2 ) };
    // test: (x, p(x), p'(x))
    std::vector<EvalTest<Real>> tests = { { Real( 0.5 ), Real( 0.75 ), Real( -2.0 ) },
                                          { Real( 2.0 ), Real( 0.0 ), Real( 1.0 ) } };
    return run_eval_case<Solver>( coeffs, tests );
  }

  template <typename Solver> std::pair<bool, double> run_cubic_eval_case()
  {
    using Real = solver_scalar_t<Solver>;
    // x^3 - 6x^2 + 11x - 6
    std::array<Real, 4>         coeffs = { Real( 1 ), Real( -6 ), Real( 11 ), Real( -6 ) };
    std::vector<EvalTest<Real>> tests  = { { Real( 0.5 ), Real( -1.875 ), Real( 5.75 ) },
                                           { Real( 3.0 ), Real( 0.0 ), Real( 2.0 ) } };
    return run_eval_case<Solver>( coeffs, tests );
  }

  template <typename Solver> std::pair<bool, double> run_quartic_eval_case()
  {
    using Real = solver_scalar_t<Solver>;
    // x^4 - 5x^2 + 4
    std::array<Real, 5>         coeffs = { Real( 1 ), Real( 0 ), Real( -5 ), Real( 0 ), Real( 4 ) };
    std::vector<EvalTest<Real>> tests  = { { Real( 0.25 ), Real( 3.69140625 ), Real( -2.4375 ) },
                                           { Real( 2.5 ), Real( 11.8125 ), Real( 37.5 ) } };
    return run_eval_case<Solver>( coeffs, tests );
  }

  void run_mode_comparison(
    TestReporter::Summary & summary,
    int                     index,
    char const *            label,
    char const *            source,
    std::pair<bool, double> ( *double_runner )()
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
    , std::pair<bool, double> ( *hq_runner )()
#endif
  )
  {
    summary.case_header( index, label, source );
    auto [ok_double, err_double] = double_runner();
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
    auto [ok_quad, err_quad]     = hq_runner();
    summary.comparison_table_header();
    summary.comparison_table_row( "double", ok_double, err_double );
    summary.comparison_table_row( "quad", ok_quad, err_quad );
    summary.comparison_table_footer();
    if ( ok_double && ok_quad )
      summary.pass( "double and quad passed" );
    else
      summary.fail( "double/quad mismatch" );
#else
    summary.comparison_table_header();
    summary.comparison_table_row( "double", ok_double, err_double );
    summary.comparison_table_footer();
    if ( ok_double )
      summary.pass( "double passed" );
    else
      summary.fail( "double failed" );
#endif
  }

}  // namespace

int main()
{
  TestReporter::Summary summary( std::cout, "Polynomial/derivative regression suite" );

  int case_id{ 1 };
  run_mode_comparison(
    summary,
    case_id++,
    "quadratic evalPolyDPoly",
    "regression for |x| <= 1 and |x| > 1",
    run_quadratic_eval_case<Quadratic>
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
    ,
    run_quadratic_eval_case<QuadraticHQ>
#endif
  );
  run_mode_comparison(
    summary,
    case_id++,
    "cubic evalPolyDPoly",
    "regression for |x| <= 1 and |x| > 1",
    run_cubic_eval_case<Cubic>
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
    ,
    run_cubic_eval_case<CubicHQ>
#endif
  );
  run_mode_comparison(
    summary,
    case_id++,
    "quartic evalPolyDPoly",
    "regression for |x| <= 1 and |x| > 1",
    run_quartic_eval_case<Quartic>
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
    ,
    run_quartic_eval_case<QuarticHQ>
#endif
  );

  {
    summary.case_header( case_id, "MAXDEGREE guard", "Jenkins-Traub entry-point" );
    std::vector<real_type> coeffs( MAXDEGREE + 2, 0 );
    std::vector<real_type> zr( MAXDEGREE + 1, 0 );
    std::vector<real_type> zi( MAXDEGREE + 1, 0 );
    coeffs.front() = 1;

    if ( PolynomialRoots::roots( coeffs.data(), MAXDEGREE + 1, zr.data(), zi.data() ) == -3 )
    {
      summary.pass( "roots(MAXDEGREE overflow)" );
    }
    else
    {
      summary.fail( "roots(MAXDEGREE overflow)" );
    }
  }

  return summary.finish();
}
