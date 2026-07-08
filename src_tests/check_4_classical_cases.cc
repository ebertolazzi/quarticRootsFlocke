#include "PolynomialRoots.hh"
#include "TestReporter.hh"
#include "TestHelpers.hh"

#include <cstddef>
#include <iostream>
#include <string>
#include <array>
#include <utility>

namespace
{

  using PolynomialRoots::Cubic;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  using PolynomialRoots::CubicHQ;
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
  using std::array;

  struct QuadraticCase
  {
    char const * label;
    char const * source;
    real_type    a, b, c;
  };

  struct CubicCase
  {
    char const * label;
    char const * source;
    real_type    a, b, c, d;
  };

  struct QuarticCase
  {
    char const * label;
    char const * source;
    real_type    a, b, c, d, e;
  };

  QuadraticCase quadratic_cases[] = { { "(x-1)^2", "classical multiple-root family", 1, -2, 1 },
                                      { "x^2+1", "classical conjugate-pair family", 1, 0, 1 },
                                      { "W2(x)=(x-1)(x-2)", "Wilkinson 1963", 1, -3, 2 } };

  CubicCase cubic_cases[] = { { "(x-1)^3", "classical multiple-root family", 1, -3, 3, -1 },
                              { "(x-1)^2(x+2)", "classical multiple-root family", 1, 0, -3, 2 },
                              { "W3(x)=(x-1)(x-2)(x-3)", "Wilkinson 1963", 1, -6, 11, -6 },
                              { "x^3-3x+1", "classical casus irreducibilis", 1, 0, -3, 1 } };

  QuarticCase quartic_cases[] = { { "(x-1)^4", "Chavez-Pichardo et al. 2022, M4 family", 1, -4, 6, -4, 1 },
                                  { "(x-1)^3(x+2)", "Chavez-Pichardo et al. 2022, M3 family", 1, -1, -3, 5, -2 },
                                  { "(x-1)^2(x-3)^2", "Chavez-Pichardo et al. 2022, DM2 family", 1, -8, 22, -24, 9 },
                                  { "W4(x)=(x-1)(x-2)(x-3)(x-4)", "Wilkinson 1963", 1, -10, 35, -50, 24 },
                                  { "(x^2-1)(x^2-4)", "classical biquadratic family", 1, 0, -5, 0, 4 },
                                  { "(x^2+1)^2", "classical double complex pair", 1, 0, 2, 0, 1 } };

  template <typename Solver> std::pair<bool, double> run_quadratic( QuadraticCase const & tc )
  {
    Solver solver( tc.a, tc.b, tc.c );
    solver.info( std::cout );
    bool   ok  = solver.check( std::cout );
    double res = TestHelpers::max_residual( solver );
    return { ok, res };
  }

  template <typename Solver> std::pair<bool, double> run_cubic( CubicCase const & tc )
  {
    Solver solver( tc.a, tc.b, tc.c, tc.d );
    solver.info( std::cout );
    bool   ok  = solver.check( std::cout );
    double res = TestHelpers::max_residual( solver );
    return { ok, res };
  }

  template <typename Solver> std::pair<bool, double> run_quartic( QuarticCase const & tc )
  {
    Solver solver( tc.a, tc.b, tc.c, tc.d, tc.e );
    solver.info( std::cout );
    bool   ok  = solver.check( std::cout );
    double res = TestHelpers::max_residual( solver );
    return { ok, res };
  }

  template <typename CaseType, typename Runner> void run_case(
    TestReporter::Summary & summary,
    int &                   idx,
    CaseType const &        tc,
    Runner                  runner_double,
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
    Runner                  runner_quad,
#endif
    bool                    known_difficult = false )
  {
    summary.case_header( idx++, tc.label, tc.source );
    auto [ok_double, res_double] = runner_double( tc );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
    auto [ok_quad, res_quad]     = runner_quad( tc );
    summary.comparison_table_header();
    summary.comparison_table_row( "double", ok_double, res_double );
    summary.comparison_table_row( "quad", ok_quad, res_quad );
    summary.comparison_table_footer();
    if ( ok_double && ok_quad ) { summary.pass( "double and quad passed" ); }
    else if ( known_difficult ) { summary.warn( "known difficult DM2 family" ); }
    else
    {
      summary.fail( "double/quad mismatch" );
    }
#else
    summary.comparison_table_header();
    summary.comparison_table_row( "double", ok_double, res_double );
    summary.comparison_table_footer();
    if ( ok_double ) { summary.pass( "double passed" ); }
    else if ( known_difficult ) { summary.warn( "known difficult DM2 family" ); }
    else
    {
      summary.fail( "double failed" );
    }
#endif
  }

}  // namespace

int main()
{
  std::cout.precision( 18 );
  TestReporter::Summary summary( std::cout, "Classical and literature-driven polynomial cases" );

  int idx{ 1 };
  for ( auto const & tc : quadratic_cases )
  {
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
    run_case( summary, idx, tc, run_quadratic<Quadratic>, run_quadratic<QuadraticHQ> );
#else
    run_case( summary, idx, tc, run_quadratic<Quadratic> );
#endif
  }

  for ( auto const & tc : cubic_cases ) {
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
    run_case( summary, idx, tc, run_cubic<Cubic>, run_cubic<CubicHQ> );
#else
    run_case( summary, idx, tc, run_cubic<Cubic> );
#endif
  }

  for ( std::size_t i{ 0 }; i < sizeof( quartic_cases ) / sizeof( quartic_cases[0] ); ++i )
  {
    auto const & tc = quartic_cases[i];
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
    run_case( summary, idx, tc, run_quartic<Quartic>, run_quartic<QuarticHQ>, i == 2 );
#else
    run_case( summary, idx, tc, run_quartic<Quartic>, i == 2 );
#endif
  }

  return summary.finish();
}
