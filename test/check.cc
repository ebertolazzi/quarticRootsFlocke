/*
.. This program solves a set of 8 cubic and 13 quartic polynomials using
.. the cubic and quartic solvers as described in the manuscript.
*/

#include "PolynomialRoots.hh"
#include <iostream>
#include <iomanip>

using namespace std ;
using namespace PolynomialRoots ;

// Set the polynomial coefficients corresponding and the exact roots.

double rootCubicReal1[] = { 1, 1.0e+7, 1.0e+14 } ;
double rootCubicImag1[] = { 0, 0, 0 } ;
double cubic1[]         = { 1,-1.00000010000001e+14,+1.00000010000001e+21,-1.e+21 } ;

double rootCubicReal2[] = { 1, 1.000001, 1.000002 } ;
double rootCubicImag2[] = { 0, 0, 0 } ;
double cubic2[]         = { 1,-3.000003,+3.000006000002,-1.000003000002 } ;

double rootCubicReal3[] = { 1.0e+80, 1.0e+77, -1.0e+81 } ;
double rootCubicImag3[] = { 0, 0, 0 } ;
double cubic3[]         = { 1,+8.999e+80,-1.0009e+161,+1.e+238 } ;

double rootCubicReal4[] = { 1, -1, -1.0e+24 } ;
double rootCubicImag4[] = { 0, 0, 0 } ;
double cubic4[]         = {1,+1.e+24,-1.0,-1.e+24} ;

double rootCubicReal5[] = { 1.0e+14, 1.0e+14, -1 } ;
double rootCubicImag5[] = { 0, 0, 0 } ;
double cubic5[]         = {1,-1.99999999999999e+14,+0.99999999999998e+28,+1.e+28} ;

double rootCubicReal6[] = { 1.0e+5, 1.0e+5, 1.0e+5 } ;
double rootCubicImag6[] = { 0, 1, -1 } ;
double cubic6[]         = {1,-3.e+5,+3.0000000001e+10,-1.0000000001e+15} ;

double rootCubicReal7[] = { 1, 1, 1 } ;
double rootCubicImag7[] = { 0, 1.0e+7, -1.0e+7 } ;
double cubic7[]         = {1,-3.0,+1.00000000000003e+14,-1.00000000000001e+14} ;

double rootCubicReal8[] = { 1, 1.0e+7, 1.0e+7 } ;
double rootCubicImag8[] = { 0, 1, -1 } ;
double cubic8[]         = {1,-2.0000001e+7,+1.00000020000001e+14,-1.00000000000001e+14} ;

double rootCubicReal9[] = { -1.0e+14, 1, 1 } ;
double rootCubicImag9[] = { 0, 1, -1 } ;
double cubic9[]         = {1,+0.99999999999998e+14,-1.99999999999998e+14,2.e+14} ;

// Set the exact quartic roots.

double rootQuarticReal1[] = { 1.0e+9, 1.0e+6, 1.0e+3, 1.0 } ;
double rootQuarticImag1[] = { 0, 0, 0, 0 } ;
double quartic1[]         = {1,-1.001001001e+9,+1.001002001001e+15,-1.001001001e+18,+1.e+18} ;

double rootQuarticReal2[] = { 1.003, 1.002, 1.001, 1.000 } ;
double rootQuarticImag2[] = { 0, 0, 0, 0 } ;
double quartic2[]         = {1,-4.006,+6.018011,-4.018022006,+1.006011006} ;

double rootQuarticReal3[] = { 1.0e+80, -1.0e+74, -1.0e+76, -1.0e+77 } ;
double rootQuarticImag3[] = { 0, 0, 0, 0 } ;
double quartic3[]         = {1,-9.98899e+79,-1.1008989e+157,-1.010999e+233,-1.e+307} ;

double rootQuarticReal4[] = { 1.0e+14, 2, 1, -1 } ;
double rootQuarticImag4[] = { 0, 0, 0, 0 } ;
double quartic4[]         = {1,-1.00000000000002e+14,+1.99999999999999e+14,+1.00000000000002e+14,-2.e+14} ;

double rootQuarticReal5[] = { 1.0e+7, 1, -1, -2.0e+7 } ;
double rootQuarticImag5[] = { 0, 0, 0, 0 } ;
double quartic5[]         = {1,+1.e+7,-2.00000000000001e+14,-1.e+7,+2.e+14} ;

double rootQuarticReal6[] = { 1.0e+7, -1.0e+6, 1, 1 } ;
double rootQuarticImag6[] = { 0, 0, 1, -1 } ;
double quartic6[]         = {1,-9.000002e+6,-0.9999981999998e+13,+1.9999982e+13,-2.e+13} ;

double rootQuarticReal7[] = { -4.0, -7.0, -1.0e+6, -1.0e+6 } ;
double rootQuarticImag7[] = { 0, 0, 1.0e+5, -1.0e+5 } ;
double quartic7[]         = {1,+2.000011e+6,+1.010022000028e+12,+1.1110056e+13,+2.828e+13} ;

double rootQuarticReal8[] = { 1.0e+8, 11.0, 1.0e+3, 1.0e+3 } ;
double rootQuarticImag8[] = { 0, 0, 1, -1 } ;
double quartic8[]         = {1,-1.00002011e+8,+2.01101022001e+11,-1.02200111000011e+14,+1.1000011e+15} ;

double rootQuarticReal9[] = { 1.0e+7, 1.0e+7, 1, 1 } ;
double rootQuarticImag9[] = { 1.0e+6, -1.0e+6, 2, -2 } ;
double quartic9[]         = {1,-2.0000002e+7,+1.01000040000005e+14,-2.020001e+14,+5.05e+14} ;

double rootQuarticReal10[] = { 1.0e+4, 1.0e+4, -7.0, -7.0 } ;
double rootQuarticImag10[] = { 3.0, -3.0, 1.0e+3, -1.0e+3 } ;
double quartic10[]         = {1,-1.9986e+4,+1.00720058e+8,-1.8600979874e+10,+1.00004909000441e+14} ;

double rootQuarticReal11[] = { 1.002, 1.002, 1.001, 1.001 } ;
double rootQuarticImag11[] = { 4.998, -4.998, 5.001, -5.001 } ;
double quartic11[]         = {1,-4.006,+5.6008018e+1,-1.04148036024e+2,+6.75896068064016e+2} ;

double rootQuarticReal12[] = { 1.0e+3, 1.0e+3, 1.0e+3, 1.0e+3 } ;
double rootQuarticImag12[] = { 3.0, -3.0, 1.0, -1.0 } ;
double quartic12[]         = {1,-4.0e+3,+6.00001e+6,-4.00002e+9,+1.000010000009e+12} ;

double rootQuarticReal13[] = { 2.0, 2.0, 1.0, 1.0 } ;
double rootQuarticImag13[] = { 1.0e+4, -1.0e+4, 1.0e+3, -1.0e+3 } ;
double quartic13[]         = {1,-6.0,+1.01000013e+8,-2.04000012e+8,+1.00000104000004e+14} ;


static
void
do_test( double const p[5] ) {
  int    nr, nc ;
  double r1, r2, r3, r4 ;
  int iter = solveQuartic( p[0], p[1], p[2], p[3], p[4], r1, r2, r3, r4, nr, nc ) ;
  cout.precision(8) ;
  cout << "\n\n"
       << p[0] << " * x^4 + "
       << p[1] << " * x^3 + "
       << p[2] << " * x^2 + "
       << p[3] << " * x + "
       << p[4] << " --> zeros" ;
  if ( nr == 4 ) {
    valueType res1 = evalPoly( p, 4, r1, true ) ;
    valueType res2 = evalPoly( p, 4, r2, true ) ;
    valueType res3 = evalPoly( p, 4, r3, true ) ;
    valueType res4 = evalPoly( p, 4, r4, true ) ;
    cout << "\nr1 = " << setw(8) << r1 << "\tp(r1) = " << res1
         << "\nr2 = " << setw(8) << r2 << "\tp(r2) = " << res2
         << "\nr3 = " << setw(8) << r3 << "\tp(r3) = " << res3
         << "\nr4 = " << setw(8) << r4 << "\tp(r4) = " << res4
         << "\niter = " << iter << "\n\n" ;
  } else if ( nr == 2 ) {
    valueType res1 = evalPoly( p, 4, r1, true ) ;
    valueType res2 = evalPoly( p, 4, r2, true ) ;
    complex<valueType> z1(r3,r4) ;
    complex<valueType> z2(r3,-r4) ;
    complex<valueType> res3 = evalPolyC( p, 4, z1, true ) ;
    complex<valueType> res4 = evalPolyC( p, 4, z2, true ) ;
    cout << "\nr1 = " << setw(8) << r1 << "\tp(r1) = " << res1
         << "\nr2 = " << setw(8) << r2 << "\tp(r2) = " << res2
         << "\nr3 = " << setw(8) << z1 << "\tp(r3) = " << abs(res3)
         << "\nr4 = " << setw(8) << z2 << "\tp(r4) = " << abs(res4)
         << "\niter = " << iter << "\n\n" ;
  } else {
    complex<valueType> z1(r1,r2) ;
    complex<valueType> z2(r1,-r2) ;
    complex<valueType> z3(r3,r4) ;
    complex<valueType> z4(r3,-r4) ;
    complex<valueType> res1 = evalPolyC( p, 4, z1, true ) ;
    complex<valueType> res2 = evalPolyC( p, 4, z2, true ) ;
    complex<valueType> res3 = evalPolyC( p, 4, z3, true ) ;
    complex<valueType> res4 = evalPolyC( p, 4, z4, true ) ;
    cout << "\nr1 = " << setw(8) << z1 << "\tp(r1) = " << abs(res1)
         << "\nr2 = " << setw(8) << z2 << "\tp(r2) = " << abs(res2)
         << "\nr3 = " << setw(8) << z3 << "\tp(r3) = " << abs(res3)
         << "\nr4 = " << setw(8) << z4 << "\tp(r4) = " << abs(res4)
         << "\niter = " << iter << "\n\n" ;
  }
}

static
void
do_test3( double const p[4], double const re[3], double const im[3] ) {
  int    nr, nc ;
  double r1, r2, r3 ;
  double re1, re2, re3 ;
  double im1, im2, im3 ;
  int iter = solveCubic( p[0], p[1], p[2], p[3], r1, r2, r3, nr, nc ) ;
  cout.precision(8) ;
  cout << "\n\n"
       << p[0] << " * x^3 + "
       << p[1] << " * x^2 + "
       << p[2] << " * x + "
       << p[3] << " --> zeros" ;
  if ( nr == 3 ) {
    valueType res1 = evalPoly( p, 3, r1, true ) ;
    valueType res2 = evalPoly( p, 3, r2, true ) ;
    valueType res3 = evalPoly( p, 3, r3, true ) ;
    cout << "\nr1 = " << setw(8) << r1 << "\tp(r1) = " << res1
         << "\nr2 = " << setw(8) << r2 << "\tp(r2) = " << res2
         << "\nr3 = " << setw(8) << r3 << "\tp(r3) = " << res3
         << "\niter = " << iter << "\n\n" ;
    re1 = r1 ;
    re2 = r2 ;
    re3 = r3 ;
    im1 = 0 ;
    im2 = 0 ;
    im3 = 0 ;
  } else {
    complex<valueType> z2(r2,r3) ;
    complex<valueType> z3(r2,-r3) ;
    complex<valueType> res1 = evalPoly( p, 3, r1, true ) ;
    complex<valueType> res2 = evalPolyC( p, 3, z2, true ) ;
    complex<valueType> res3 = evalPolyC( p, 3, z3, true ) ;
    cout << "\nr1 = " << setw(8) << r1 << "\tp(r1) = " << abs(res1)
         << "\nr2 = " << setw(8) << z2 << "\tp(r2) = " << abs(res2)
         << "\nr3 = " << setw(8) << z3 << "\tp(r3) = " << abs(res3)
         << "\niter = " << iter << "\n\n" ;
    re1 = r1 ;
    re2 = r2 ;
    re3 = r2 ;
    im1 = 0 ;
    im2 = r3 ;
    im3 = -r3 ;
  }
  cout << "(1) = real " << setw(14) << re[0] << " err = " << re1-re[0] << '\n' ;
  cout << "(1) = imag " << setw(14) << im[0] << " err = " << im1-im[0] << '\n' ;
  cout << "(2) = real " << setw(14) << re[1] << " err = " << re2-re[1] << '\n' ;
  cout << "(2) = imag " << setw(14) << im[1] << " err = " << im2-im[1] << '\n' ;
  cout << "(3) = real " << setw(14) << re[2] << " err = " << re3-re[2] << '\n' ;
  cout << "(4) = imag " << setw(14) << im[2] << " err = " << im3-im[2] << '\n' ;
}

#define TESTCUBIC(A) do_test3( cubic##A, rootCubicReal##A, rootCubicImag##A )

int
main() {
  //TESTCUBIC(1) ;
  TESTCUBIC(2) ;
  //TESTCUBIC(3) ;
  //TESTCUBIC(4) ;
  //TESTCUBIC(5) ;
  //TESTCUBIC(6) ;
  //TESTCUBIC(7) ;
  //TESTCUBIC(8) ;
  //TESTCUBIC(9) ;
  return 0 ;
}
