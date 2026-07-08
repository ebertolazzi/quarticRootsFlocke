/*
.. This program solves a set of 8 cubic and 13 quartic polynomials using
.. the cubic and quartic solvers as described in the manuscript.
*/

#include "PolynomialRoots.hh"
#include "TestReporter.hh"
#include "TestHelpers.hh"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <utility>
#include <array>

using PolynomialRoots::Quartic;
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
using PolynomialRoots::QuarticHQ;
#endif
using std::array;

// Set the polynomial coefficients corresponding and the exact roots.
// Set the exact quartic roots.

static double rootQuarticReal0[] = { -1.017229096421223459824821857711825869530e15,
                                     0.5313292845913352486945402695550732313951e-1,
                                     -0.2656646422956676243472701347775328370522e-1,
                                     -0.2656646422956676243472701347775328370522e-1 };
static double rootQuarticImag0[] = { 0,
                                     0,
                                     0.4601446582307080400727162429104527153041e-1,
                                     -0.4601446582307080400727162429104527153041e-1 };
static double quartic0[]         = { -5.898376305898030e-09,
                                     -6.000000000001007e+06,
                                     3.664603343001005e-27,
                                     0,
                                     9.000000000003021e+02 };

static double rootQuarticReal1[] = { 1, 1e3, 1e6, 1e+9 };
static double rootQuarticImag1[] = { 0, 0, 0, 0 };
static double quartic1[] = { 1, -1001001001.0, 1001002001001000.0, -1001001001000000000.0, 1000000000000000000.0 };

static double rootQuarticReal2[] = { 1.003, 1.002, 1.001, 1.000 };
static double rootQuarticImag2[] = { 0, 0, 0, 0 };
static double quartic2[]         = { 1, -4.006, +6.018011, -4.018022006, +1.006011006 };

static double rootQuarticReal3[] = { -1.0e+77, -1.0e+76, -1.0e+74, 1.0e+80 };
static double rootQuarticImag3[] = { 0, 0, 0, 0 };
static double quartic3[]         = { 1, -9.988990e79, -1.100898900e157, -1.010999e233, -1e307 };

static double rootQuarticReal4[] = { 1.0e+14, 2, 1, -1 };
static double rootQuarticImag4[] = { 0, 0, 0, 0 };
static double quartic4[]         = { 1, -1.00000000000002e+14, +1.99999999999999e+14, +1.00000000000002e+14, -2.e+14 };

static double rootQuarticReal5[] = { 1.0e+7, 1, -1, -2.0e+7 };
static double rootQuarticImag5[] = { 0, 0, 0, 0 };
static double quartic5[]         = { 1, +1.e+7, -2.00000000000001e+14, -1.e+7, +2.e+14 };

static double rootQuarticReal6[] = { 1.0e+7, -1.0e+6, 1, 1 };
static double rootQuarticImag6[] = { 0, 0, 1, -1 };
static double quartic6[]         = { 1, -9.000002e+6, -0.9999981999998e+13, +1.9999982e+13, -2.e+13 };

static double rootQuarticReal7[] = { -1.0e+6, -1.0e+6, -7, -4 };
static double rootQuarticImag7[] = { 1.0e+5, -1.0e+5, 0, 0 };
static double quartic7[]         = { 1, +2.000011e+6, +1.010022000028e+12, +1.1110056e+13, +2.828e+13 };

static double rootQuarticReal8[] = { 1000, 1000, 11.0, 1.0e+8 };
static double rootQuarticImag8[] = { 1, -1, 0, 0 };
static double quartic8[]         = { 1, -1.00002011e8, +2.01101022001e11, -1.02200111000011e14, 1.1000011e15 };
// 1,-1.00002011e+8,+2.01101022001e+11,-1.02200111000011e+14,+1.1000011e+15};

static double rootQuarticReal9[] = { 1.0e+7, 1.0e+7, 1, 1 };
static double rootQuarticImag9[] = { 1.0e+6, -1.0e+6, 2, -2 };
static double quartic9[]         = { 1, -2.0000002e+7, +1.01000040000005e+14, -2.020001e+14, +5.05e+14 };

static double rootQuarticReal10[] = { 1.0e+4, 1.0e+4, -7.0, -7.0 };
static double rootQuarticImag10[] = { 3.0, -3.0, 1.0e+3, -1.0e+3 };
static double quartic10[]         = { 1, -1.9986e+4, +1.00720058e+8, -1.8600979874e+10, +1.00004909000441e+14 };

static double rootQuarticReal11[] = { 1.002, 1.002, 1.001, 1.001 };
static double rootQuarticImag11[] = { 4.998, -4.998, 5.001, -5.001 };
static double quartic11[]         = { 1, -4.006, +5.6008018e+1, -1.04148036024e+2, +6.75896068064016e+2 };

static double rootQuarticReal12[] = { 1.0e+3, 1.0e+3, 1.0e+3, 1.0e+3 };
static double rootQuarticImag12[] = { 3.0, -3.0, 1.0, -1.0 };
static double quartic12[]         = { 1, -4.0e+3, +6.00001e+6, -4.00002e+9, +1.000010000009e+12 };

static double rootQuarticReal13[] = { 2.0, 2.0, 1.0, 1.0 };
static double rootQuarticImag13[] = { 1.0e+4, -1.0e+4, 1.0e+3, -1.0e+3 };
static double quartic13[]         = { 1, -6.0, +1.01000013e+8, -2.04000012e+8, +1.00000104000004e+14 };

static double rootQuarticReal14[] = { 0.064423196657504, 0.064423196657504, 0.046605762356559, -0.035452155671566 };
static double rootQuarticImag14[] = { 0.122215870576901, -0.122215870576901, 0, 0 };
static double quartic14[]         = { 714285.71428571432,
                                      -100000,
                                      13479.924671428573,
                                      -2.2737367544323206e-13,
                                      -22.526485131747265 };

static double r                   = std::sqrt( 2.0 );
static double rootQuarticReal15[] = { 0, r, r, r };
static double rootQuarticImag15[] = { 0, 0, 0, 0 };
static double quartic15[]         = { 1, -3 * r, 6, -2 * r, 0 };
// double quartic15[]         = { 0, 1, -3*r, 6, -2*r };

static double rootQuarticReal16[] = { -r, -r, -r, 0 };
static double rootQuarticImag16[] = { 0, 0, 0, 0 };
static double quartic16[]         = { 1, 3 * r, 6, 2 * r, 0 };
// double quartic16[]         = { 0, 1, 3*r, 6, 2*r };

static double rootQuarticReal17[] = { 19.02897754136850147151903,
                                      19.02897754136850147151903,
                                      0.5103896497990552414378477e-1,
                                      3.264106471395880584257508 };
static double rootQuarticImag17[] = { 3.948839665289545498594455, -3.948839665289545498594455, 0, 0 };
static double quartic17[]         = { 0.000158925, -0.00657522, 0.0801029, -0.2, 0.01 };

static double rootQuarticReal18[] = { -8.767990511, -8.767990511, 0.7679905093, 0.7679905118 };
static double rootQuarticImag18[] = { 0, 0, 0, 0 };
static double quartic18[]         = { 2.25,
                                      36,
                                      113.698199211166496525038382969796657562255859375,
                                      -242.41440631066808464311179704964160919189453125,
                                      102.0221256717945692571447580121457576751708984375 };

// Emanuele Siego example
static double rootQuarticReal19[] = { 0.14271266144792, 0.14271266144792, 0.36842201250793, 0.36842202205209 };
static double rootQuarticImag19[] = { 0.060284865351686, -0.060284865351686, 0, 0 };
static double quartic19[]         = { 5.0 * ( -1462.2326183798939 ),
                                      4.0 * ( 1868.4944990527565 ),
                                      3.0 * ( -901.83171888631921 ),
                                      2.0 * ( 206.27445843560130 ),
                                      1.0 * ( -23.818259047289210 ) };

// Emanuele Siego example #2
static double rootQuarticReal20[] = { 0.14271266144792, 0.14271266144792, 0.36842201250793, 0.36842202205209 };
static double rootQuarticImag20[] = { 0.060284865351686, -0.060284865351686, 0, 0 };
static double quartic20[]         = { 5.0 * ( -4.1254841444253390e+02 ),
                                      4.0 * ( 5.9401042551598744e+02 ),
                                      3.0 * ( -3.1743361677460911e+02 ),
                                      2.0 * ( 7.7477065082592929e+01 ),
                                      1.0 * ( -8.9469482641160880e+00 ) };

static double rootQuarticReal21[] = { 50.29261772409999764477107,
                                      50.29261772409999764477107,
                                      -109.5695612171199968543078,
                                      -32.79238302114399772335673 };
static double rootQuarticImag21[] = { 3.050397502900299873118684, -3.050397502900299873118684, 0, 0 };
static double quartic21[]         = { -1.27321767367025728318135329174e-09,
                                      -5.31908439792814016019367533931e-08,
                                      1.04248648319060507416809052694e-05,
                                      0,
                                      -0.0116136513934813456055961111701 };

// Emanuele Siego example #3
static double rootQuarticReal22[] = { -1.5226962261856962e-08,
                                      1.5226962261856962e-08,
                                      0.0029426967763176593,
                                      0.0029426967763176593 };
static double rootQuarticImag22[] = { 0, 0, 0, 0 };
static double quartic22[]         = { 86594.643180459301,
                                      -509.64355468067714,
                                      0.74986322275510664,
                                      0.0000000000000000,
                                      -1.7386537088520670e-16 };

static char const * caseSource[] = { "Flocke 2015",          "Flocke 2015",          "Flocke 2015",
                                     "Flocke 2015",          "Flocke 2015",          "Flocke 2015",
                                     "Flocke 2015",          "Flocke 2015",          "Flocke 2015",
                                     "Flocke 2015",          "Flocke 2015",          "Flocke 2015",
                                     "Flocke 2015",          "extended stress case", "extended stress case",
                                     "extended stress case", "extended stress case", "extended stress case",
                                     "extended stress case", "Emanuele Siego",       "Emanuele Siego",
                                     "extended stress case", "Emanuele Siego" };

template <typename Solver>
static std::pair<bool, double> do_test( double const p[5], double const re[4], double const im[4] )
{
  Solver solver( p[0], p[1], p[2], p[3], p[4] );
  solver.info( std::cout );
  bool ok = solver.check( std::cout );
  if ( !ok )
  {
    std::cout << "\nExpected\n"
              << "x0 = (" << re[0] << ", " << im[0] << ")\n"
              << "x1 = (" << re[1] << ", " << im[1] << ")\n"
              << "x2 = (" << re[2] << ", " << im[2] << ")\n"
              << "x3 = (" << re[3] << ", " << im[3] << ")\n\n";
  }
  double res = TestHelpers::max_residual( solver );  // nessun template parameter
  return { ok, res };
}

static void run_case(
  TestReporter::Summary & summary,
  int                     index,
  char const *            source,
  double const            p[5],
  double const            re[4],
  double const            im[4],
  bool                    known_difficult = false )
{
  summary.case_header( index, "quartic stress case", source );
  auto [ok_double, res_double] = do_test<Quartic>( p, re, im );
#if POLYNOMIAL_ROOTS_HAS_MULTIPRECISION
  auto [ok_quad, res_quad]     = do_test<QuarticHQ>( p, re, im );
  summary.comparison_table_header();
  summary.comparison_table_row( "double", ok_double, res_double );
  summary.comparison_table_row( "quad", ok_quad, res_quad );
  summary.comparison_table_footer();
  if ( ok_double && ok_quad ) { summary.pass( "double and quad passed" ); }
  else if ( known_difficult ) { summary.warn( "known difficult legacy case" ); }
  else
  {
    summary.fail( "double/quad mismatch" );
  }
#else
  summary.comparison_table_header();
  summary.comparison_table_row( "double", ok_double, res_double );
  summary.comparison_table_footer();
  if ( ok_double ) { summary.pass( "double passed" ); }
  else if ( known_difficult ) { summary.warn( "known difficult legacy case" ); }
  else
  {
    summary.fail( "double failed" );
  }
#endif
}

int main()
{
  std::cout.precision( 25 );
  TestReporter::Summary summary( std::cout, "Quartic solver regression suite" );
  run_case( summary, 1, caseSource[0], quartic0, rootQuarticReal0, rootQuarticImag0 );
  run_case( summary, 2, caseSource[1], quartic1, rootQuarticReal1, rootQuarticImag1 );
  run_case( summary, 3, caseSource[2], quartic2, rootQuarticReal2, rootQuarticImag2 );
  run_case( summary, 4, caseSource[3], quartic3, rootQuarticReal3, rootQuarticImag3, true );
  run_case( summary, 5, caseSource[4], quartic4, rootQuarticReal4, rootQuarticImag4 );
  run_case( summary, 6, caseSource[5], quartic5, rootQuarticReal5, rootQuarticImag5 );
  run_case( summary, 7, caseSource[6], quartic6, rootQuarticReal6, rootQuarticImag6 );
  run_case( summary, 8, caseSource[7], quartic7, rootQuarticReal7, rootQuarticImag7, true );
  run_case( summary, 9, caseSource[8], quartic8, rootQuarticReal8, rootQuarticImag8 );
  run_case( summary, 10, caseSource[9], quartic9, rootQuarticReal9, rootQuarticImag9, true );
  run_case( summary, 11, caseSource[10], quartic10, rootQuarticReal10, rootQuarticImag10 );
  run_case( summary, 12, caseSource[11], quartic11, rootQuarticReal11, rootQuarticImag11 );
  run_case( summary, 13, caseSource[12], quartic12, rootQuarticReal12, rootQuarticImag12 );
  run_case( summary, 14, caseSource[13], quartic13, rootQuarticReal13, rootQuarticImag13 );
  run_case( summary, 15, caseSource[14], quartic14, rootQuarticReal14, rootQuarticImag14 );
  run_case( summary, 16, caseSource[15], quartic15, rootQuarticReal15, rootQuarticImag15 );
  run_case( summary, 17, caseSource[16], quartic16, rootQuarticReal16, rootQuarticImag16 );
  run_case( summary, 18, caseSource[17], quartic17, rootQuarticReal17, rootQuarticImag17 );
  run_case( summary, 19, caseSource[18], quartic18, rootQuarticReal18, rootQuarticImag18 );
  run_case( summary, 20, caseSource[19], quartic19, rootQuarticReal19, rootQuarticImag19 );
  run_case( summary, 21, caseSource[20], quartic20, rootQuarticReal20, rootQuarticImag20 );
  run_case( summary, 22, caseSource[21], quartic21, rootQuarticReal21, rootQuarticImag21 );
  run_case( summary, 23, caseSource[22], quartic22, rootQuarticReal22, rootQuarticImag22 );
  return summary.finish();
}

// 3, 8, 17
