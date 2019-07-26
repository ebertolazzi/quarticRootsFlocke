/*
.. This program solves a set of 8 cubic polynomials using
.. the cubic solvers as described in
..
.. N. FLOCKE
.. Algorithm 954: An Accurate and Efficient Cubic and Quartic Equation Solver
.. for Physical Applications
.. ACM TOMS, Vol. 41, No. 4, 2015.
.. DOI: http://dx.doi.org/10.1145/2699468
*/

#include "PolynomialRoots.hh"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace PolynomialRoots;

// Set the polynomial coefficients corresponding and the exact roots.

double rootCubicReal1[] = { 1, 1.0e+7, 1.0e+14 };
double rootCubicImag1[] = { 0, 0, 0 };
double cubic1[]         = { 1, -1.00000010000001e+14, +1.00000010000001e+21, -1.e+21 };

double rootCubicReal2[] = { 1, 1.000001, 1.000002 };
double rootCubicImag2[] = { 0, 0, 0 };
double cubic2[]         = { 1, -3.000003, +3.000006000002, -1.000003000002 };

double rootCubicReal3[] = { -1.0e+81, 1.0e+77, 1.0e+80 };
double rootCubicImag3[] = { 0, 0, 0 };
double cubic3[]         = { 1, 8.999e+80, -1.0009e+161, +1.e+238 };

double rootCubicReal4[] = { 1, -1, -1.0e+24 };
double rootCubicImag4[] = { 0, 0, 0 };
double cubic4[]         = { 1, +1.e+24, -1.0, -1.e+24 };

double rootCubicReal5[] = { -1, 1.0e+14, 1.0e+14 };
double rootCubicImag5[] = { 0, 0, 0 };
double cubic5[]         = { 1, -1.99999999999999e+14, +9.9999999999998e+27, +1.e+28 };

double rootCubicReal6[] = { 1.0e+5, 1.0e+5, 1.0e+5 };
double rootCubicImag6[] = { 0, 1, -1 };
double cubic6[]         = { 1, -3.e+5, +3.0000000001e+10, -1.0000000001e+15 };

double rootCubicReal7[] = { 1, 1, 1 };
double rootCubicImag7[] = { 0, 1.0e+7, -1.0e+7 };
double cubic7[]         = { 1, -3.0, +1.00000000000003e+14, -1.00000000000001e+14 };

double rootCubicReal8[] = { 1, 1.0e+7, 1.0e+7 };
double rootCubicImag8[] = { 0, 1, -1 };
double cubic8[]         = { 1, -2.0000001e+7, +1.00000020000001e+14, -1.00000000000001e+14 };

double rootCubicReal9[] = { -1.0e+14, 1, 1 };
double rootCubicImag9[] = { 0, 1, -1 };
double cubic9[]         = { 1, 9.9999999999998e13, -1.99999999999998e14, 2e+14 };

double r = sqrt(2.0);
double rootCubicReal10[] = { r, r, r };
double rootCubicImag10[] = { 0, 0, 0 };
double cubic10[]         = { 1, -3*r, 6, -2*r };

double rootCubicReal11[] = { -r, -r, -r };
double rootCubicImag11[] = { 0, 0, 0 };
double cubic11[]         = { 1, 3*r, 6, 2*r };

double rootCubicReal12[] = { -r, -r, -r };
double rootCubicImag12[] = { 0, 0, 0 };
double cubic12[]         = { 2.25, 37.727978648722122, 142.67292879536021, -132.84295088443082 };

static
void
do_test( double const p[4] ) {
  //Cubic csolve( p[3], p[2], p[1], p[0] );
  Cubic csolve( p[0], p[1], p[2], p[3] );
  csolve.info( cout );
  if ( !csolve.check( cout ) ) {
    cout << "\nFailed!\n\n\n";
    //std::exit(0);
  }
}

#define DO_TEST( N ) \
  cout << "\n\nText N." << N << '\n'; do_test( cubic##N )

int
main() {
  cout.precision(18);
  DO_TEST(1);
  DO_TEST(2);
  DO_TEST(3); // failed check
  DO_TEST(4);
  DO_TEST(5); // failed check
  DO_TEST(6);
  DO_TEST(7);
  DO_TEST(8);
  DO_TEST(9);
  DO_TEST(10);
  DO_TEST(11);
  DO_TEST(12);
  cout << "\n\nALL DONE!\n";
  return 0;
}
