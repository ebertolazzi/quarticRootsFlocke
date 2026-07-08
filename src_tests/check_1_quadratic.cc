/*
.. This program solves a set of 8 cubic and 13 quartic polynomials using
.. the cubic and quartic solvers as described in the manuscript.
*/

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
#include "TestReporter.hh"
#include "TestHelpers.hh"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <utility>

using PolynomialRoots::integer;
using PolynomialRoots::Quadratic;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
using PolynomialRoots::QuadraticHQ;
#endif
using std::array;

// Set the polynomial coefficients corresponding and the exact roots.
static double qq[][3] = { { 50, 50, 0 },      { 50, 50, 1 },     { 50, 50, 50 },    { 50, 50, 99 },
                          { 50, 50, 100 },    { 50, 0, 50 },     { 50, 1, 50 },     { 50, 99, 50 },
                          { 50, 100, 50 },    { 0, 50, 50 },     { 1, 50, 50 },     { 99, 50, 50 },
                          { 1e10, 50, 50 },   { 100, 1e10, 50 }, { 100, 50, 1e10 }, { 1e-10, 50, 50 },
                          { 100, 1e-10, 50 }, { 100, 50, 1e-10 } };

template <typename Solver> static std::pair<bool, double> do_test( double const p[3] )
{
  Solver const qsolve( p[0], p[1], p[2] );
  std::cout.precision( 8 );
  qsolve.info( std::cout );
  bool   ok  = qsolve.check( std::cout );
  double res = TestHelpers::max_residual( qsolve );  // nessun template parameter
  return { ok, res };
}

static void run_case( TestReporter::Summary & summary, int index, double const p[3] )
{
  summary.case_header( index, "scaled quadratic coefficients", "Flocke-style regression set" );

  auto [ok_double, res_double] = do_test<Quadratic>( p );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  auto [ok_quad, res_quad]     = do_test<QuadraticHQ>( p );

  summary.comparison_table_header();
  summary.comparison_table_row( "double", ok_double, res_double );
  summary.comparison_table_row( "quad", ok_quad, res_quad );
  summary.comparison_table_footer();

  if ( ok_double && ok_quad )
    summary.pass( "double and quad passed" );
  else
    summary.fail( "double/quad mismatch" );
#else
  summary.comparison_table_header();
  summary.comparison_table_row( "double", ok_double, res_double );
  summary.comparison_table_footer();

  if ( ok_double )
    summary.pass( "double passed" );
  else
    summary.fail( "double failed" );
#endif
}

int main()
{
  std::cout.precision( 20 );
  TestReporter::Summary summary( std::cout, "Quadratic solver regression suite" );
  constexpr integer     N = static_cast<integer>( sizeof( qq ) / sizeof( qq[0] ) );
  for ( integer k = 0; k < N; ++k ) run_case( summary, k + 1, qq[k] );
  // run_case( summary, 14, qq[13] );
  return summary.finish();
}
