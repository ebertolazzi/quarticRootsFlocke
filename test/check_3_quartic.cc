/*
.. This program solves a set of 8 cubic and 13 quartic polynomials using
.. the cubic and quartic solvers as described in the manuscript.
*/

#include "PolynomialRoots.hh"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace PolynomialRoots;

// Set the polynomial coefficients corresponding and the exact roots.
// Set the exact quartic roots.

double rootQuarticReal1[] = { 1.0e+9, 1.0e+6, 1.0e+3, 1.0 };
double rootQuarticImag1[] = { 0, 0, 0, 0 };
double quartic1[]         = { 1,-1001001001.0, 1001002001001000.0,-1001001001000000000.0,1000000000000000000.0 };

double rootQuarticReal2[] = { 1.003, 1.002, 1.001, 1.000 };
double rootQuarticImag2[] = { 0, 0, 0, 0 };
double quartic2[]         = { 1,-4.006,+6.018011,-4.018022006,+1.006011006 };

double rootQuarticReal3[] = { 1.0e+80, -1.0e+74, -1.0e+76, -1.0e+77 };
double rootQuarticImag3[] = { 0, 0, 0, 0 };
double quartic3[]         = { 1,-9.98899e+79,-1.1008989e+157,-1.010999e+233,-1.e+307 };

double rootQuarticReal4[] = { 1.0e+14, 2, 1, -1 };
double rootQuarticImag4[] = { 0, 0, 0, 0 };
double quartic4[]         = {1,-1.00000000000002e+14,+1.99999999999999e+14,+1.00000000000002e+14,-2.e+14};

double rootQuarticReal5[] = { 1.0e+7, 1, -1, -2.0e+7 };
double rootQuarticImag5[] = { 0, 0, 0, 0 };
double quartic5[]         = { 1,+1.e+7,-2.00000000000001e+14,-1.e+7,+2.e+14 };

double rootQuarticReal6[] = { 1.0e+7, -1.0e+6, 1, 1 };
double rootQuarticImag6[] = { 0, 0, 1, -1 };
double quartic6[]         = { 1,-9.000002e+6,-0.9999981999998e+13,+1.9999982e+13,-2.e+13 };

double rootQuarticReal7[] = { -4.0, -7.0, -1.0e+6, -1.0e+6 };
double rootQuarticImag7[] = { 0, 0, 1.0e+5, -1.0e+5 };
double quartic7[]         = { 1,+2.000011e+6,+1.010022000028e+12,+1.1110056e+13,+2.828e+13 };

double rootQuarticReal8[] = { 1.0e+8, 11.0, 1.0e+3, 1.0e+3 };
double rootQuarticImag8[] = { 0, 0, 1, -1 };
double quartic8[]         = {1,-1.00002011e+8,+2.01101022001e+11,-1.02200111000011e+14,+1.1000011e+15};

double rootQuarticReal9[] = { 1.0e+7, 1.0e+7, 1, 1 };
double rootQuarticImag9[] = { 1.0e+6, -1.0e+6, 2, -2 };
double quartic9[]         = { 1,-2.0000002e+7,+1.01000040000005e+14,-2.020001e+14,+5.05e+14 };

double rootQuarticReal10[] = { 1.0e+4, 1.0e+4, -7.0, -7.0 };
double rootQuarticImag10[] = { 3.0, -3.0, 1.0e+3, -1.0e+3 };
double quartic10[]         = { 1,-1.9986e+4,+1.00720058e+8,-1.8600979874e+10,+1.00004909000441e+14 };

double rootQuarticReal11[] = { 1.002, 1.002, 1.001, 1.001 };
double rootQuarticImag11[] = { 4.998, -4.998, 5.001, -5.001 };
double quartic11[]         = { 1,-4.006,+5.6008018e+1,-1.04148036024e+2,+6.75896068064016e+2 };

double rootQuarticReal12[] = { 1.0e+3, 1.0e+3, 1.0e+3, 1.0e+3 };
double rootQuarticImag12[] = { 3.0, -3.0, 1.0, -1.0 };
double quartic12[]         = { 1,-4.0e+3,+6.00001e+6,-4.00002e+9,+1.000010000009e+12 };

double rootQuarticReal13[] = { 2.0, 2.0, 1.0, 1.0 };
double rootQuarticImag13[] = { 1.0e+4, -1.0e+4, 1.0e+3, -1.0e+3 };
double quartic13[]         = { 1,-6.0,+1.01000013e+8,-2.04000012e+8,+1.00000104000004e+14 };

double rootQuarticReal14[] = { 0.064423196657504, 0.064423196657504, 0.046605762356559, -0.035452155671566 };
double rootQuarticImag14[] = { 0.122215870576901, -0.122215870576901, 0, 0 };
double quartic14[]         = { 714285.71428571432,
                               -100000,
                               13479.924671428573,
                               -2.2737367544323206e-13,
                               -22.526485131747265 };

#define DO_TEST( N ) do_test( quartic##N, rootQuarticReal##N, rootQuarticImag##N );

static
void
do_test( double const p[5],
         double const re[4],
         double const im[4] ) {
  Quartic qsolve( p[0], p[1], p[2], p[3], p[4] );
  qsolve.info( cout );
  if ( !qsolve.check( cout ) ) {
    cout << "\n\nFailed!\n\nExpected\n"
         << "x0 = (" << re[0] << ", " << im[0] << ")\n"
         << "x1 = (" << re[1] << ", " << im[1] << ")\n"
         << "x2 = (" << re[2] << ", " << im[2] << ")\n"
         << "x3 = (" << re[3] << ", " << im[3] << ")\n\n\n\n\n";
    //std::exit(0);
  } else {
    cout << "\nOK!\n";
  }
}

#define TESTCUBIC(A) do_test3( cubic##A, rootCubicReal##A, rootCubicImag##A )

int
main() {
  cout.precision(14);
  DO_TEST(1);
  DO_TEST(2);
  DO_TEST(3);
  DO_TEST(4);
  DO_TEST(5);
  DO_TEST(6);
  DO_TEST(7);
  DO_TEST(8);
  DO_TEST(9);
  DO_TEST(10);
  DO_TEST(11);
  DO_TEST(12);
  DO_TEST(13);
  DO_TEST(14);
  cout << "\n\nALL DONE!\n";
  return 0;
}

