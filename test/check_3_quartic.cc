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

static double rootQuarticReal0[] = { -1.017229096421223459824821857711825869530e15, 0.5313292845913352486945402695550732313951e-1, -0.2656646422956676243472701347775328370522e-1, -0.2656646422956676243472701347775328370522e-1 };
static double rootQuarticImag0[] = { 0, 0, 0.4601446582307080400727162429104527153041e-1, -0.4601446582307080400727162429104527153041e-1 };
static double quartic0[]         = { -5.898376305898030e-09, -6.000000000001007e+06, 3.664603343001005e-27, 0, 9.000000000003021e+02 };

static double rootQuarticReal1[] = { 1, 1e3, 1e6, 1e+9 };
static double rootQuarticImag1[] = { 0, 0, 0, 0 };
static double quartic1[]         = { 1, -1001001001.0, 1001002001001000.0, -1001001001000000000.0, 1000000000000000000.0 };

static double rootQuarticReal2[] = { 1.003, 1.002, 1.001, 1.000 };
static double rootQuarticImag2[] = { 0, 0, 0, 0 };
static double quartic2[]         = { 1,-4.006,+6.018011,-4.018022006,+1.006011006 };

static double rootQuarticReal3[] = { -1.0e+77, -1.0e+76, -1.0e+74, 1.0e+80 };
static double rootQuarticImag3[] = { 0, 0, 0, 0 };
static double quartic3[]         = { 1, -9.988990e79, -1.100898900e157, -1.010999e233, -1e307 };

static double rootQuarticReal4[] = { 1.0e+14, 2, 1, -1 };
static double rootQuarticImag4[] = { 0, 0, 0, 0 };
static double quartic4[]         = {1,-1.00000000000002e+14,+1.99999999999999e+14,+1.00000000000002e+14,-2.e+14};

static double rootQuarticReal5[] = { 1.0e+7, 1, -1, -2.0e+7 };
static double rootQuarticImag5[] = { 0, 0, 0, 0 };
static double quartic5[]         = { 1,+1.e+7,-2.00000000000001e+14,-1.e+7,+2.e+14 };

static double rootQuarticReal6[] = { 1.0e+7, -1.0e+6, 1, 1 };
static double rootQuarticImag6[] = { 0, 0, 1, -1 };
static double quartic6[]         = { 1,-9.000002e+6,-0.9999981999998e+13,+1.9999982e+13,-2.e+13 };

static double rootQuarticReal7[] = { -4.0, -7.0, -1.0e+6, -1.0e+6 };
static double rootQuarticImag7[] = { 0, 0, 1.0e+5, -1.0e+5 };
static double quartic7[]         = { 1,+2.000011e+6,+1.010022000028e+12,+1.1110056e+13,+2.828e+13 };

static double rootQuarticReal8[] = { 1.0e+3, 1.0e+3, 11.0, 1.0e+8 };
static double rootQuarticImag8[] = {      1,     -1,    0,    0   };
static double quartic8[]         = { 1, -1.00002011e8, +2.01101022001e11, -1.02200111000011e14, 1.1000011e15 };
// 1,-1.00002011e+8,+2.01101022001e+11,-1.02200111000011e+14,+1.1000011e+15};

static double rootQuarticReal9[] = { 1.0e+7, 1.0e+7, 1, 1 };
static double rootQuarticImag9[] = { 1.0e+6, -1.0e+6, 2, -2 };
static double quartic9[]         = { 1,-2.0000002e+7,+1.01000040000005e+14,-2.020001e+14,+5.05e+14 };

static double rootQuarticReal10[] = { 1.0e+4, 1.0e+4, -7.0, -7.0 };
static double rootQuarticImag10[] = { 3.0, -3.0, 1.0e+3, -1.0e+3 };
static double quartic10[]         = { 1,-1.9986e+4,+1.00720058e+8,-1.8600979874e+10,+1.00004909000441e+14 };

static double rootQuarticReal11[] = { 1.002, 1.002, 1.001, 1.001 };
static double rootQuarticImag11[] = { 4.998, -4.998, 5.001, -5.001 };
static double quartic11[]         = { 1,-4.006,+5.6008018e+1,-1.04148036024e+2,+6.75896068064016e+2 };

static double rootQuarticReal12[] = { 1.0e+3, 1.0e+3, 1.0e+3, 1.0e+3 };
static double rootQuarticImag12[] = { 3.0, -3.0, 1.0, -1.0 };
static double quartic12[]         = { 1,-4.0e+3,+6.00001e+6,-4.00002e+9,+1.000010000009e+12 };

static double rootQuarticReal13[] = { 2.0, 2.0, 1.0, 1.0 };
static double rootQuarticImag13[] = { 1.0e+4, -1.0e+4, 1.0e+3, -1.0e+3 };
static double quartic13[]         = { 1,-6.0,+1.01000013e+8,-2.04000012e+8,+1.00000104000004e+14 };

static double rootQuarticReal14[] = { 0.064423196657504, 0.064423196657504, 0.046605762356559, -0.035452155671566 };
static double rootQuarticImag14[] = { 0.122215870576901, -0.122215870576901, 0, 0 };
static double quartic14[]         = {
  714285.71428571432,
  -100000,
  13479.924671428573,
  -2.2737367544323206e-13,
  -22.526485131747265
};

static double r = sqrt(2.0);
static double rootQuarticReal15[] = { 0, r, r, r };
static double rootQuarticImag15[] = { 0, 0, 0, 0 };
static double quartic15[]         = { 1, -3*r, 6, -2*r, 0 };
//double quartic15[]         = { 0, 1, -3*r, 6, -2*r };

static double rootQuarticReal16[] = { -r, -r, -r, 0 };
static double rootQuarticImag16[] = { 0, 0, 0, 0 };
static double quartic16[]         = { 1, 3*r, 6, 2*r, 0 };
//double quartic16[]         = { 0, 1, 3*r, 6, 2*r };

static double rootQuarticReal17[] = { 19.01804207, 3.926187087, 0.5001001646e-3, 3.336516286 };
static double rootQuarticImag17[] = { 0, 0, 0, 0 };
static double quartic17[]         = { 0.000158925, -0.00657522, 0.0801029, -0.2, 0.01 };

static double rootQuarticReal18[] = { -8.767990511, -8.767990511, 0.7679905093, 0.7679905118 };
static double rootQuarticImag18[] = { 0, 0, 0, 0 };
static double quartic18[]         = {
  2.25, 36,
  113.698199211166496525038382969796657562255859375,
  -242.41440631066808464311179704964160919189453125,
  102.0221256717945692571447580121457576751708984375
};

#define DO_TEST( N ) \
  cout << "\n\nText N." << N << '\n'; \
  do_test( quartic##N, rootQuarticReal##N, rootQuarticImag##N )

static
void
do_test(
  double const p[5],
  double const re[4],
  double const im[4]
) {
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

int
main() {
  cout.precision(14);
  //DO_TEST(9);
  //cout << "\n\nALL DONE!\n";
  DO_TEST(0);
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
  DO_TEST(15);
  DO_TEST(16);
  DO_TEST(17);
  DO_TEST(18);
  cout << "\n\nALL DONE!\n";
  return 0;
}

// 3, 8, 17
