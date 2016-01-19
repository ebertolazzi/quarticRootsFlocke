#include "PolynomialRoots.hh"
#include <iostream>
#include <iomanip>

using namespace std ;
using namespace PolynomialRoots ;

double p1[5] = { 1, -10, 35,-50, 24 } ; // 1, 2, 3, 4
double p2[5] = { 1,-3, 3,-3, 2 } ; // 1, 2, I, -I
double p3[5] = { 1,0,5,0,4 } ; // I, -I, 2*I, -2*I
double p4[5] = { 1,-2,6,-2,5 } ; // I, -I, 1+2*I, 1-2*I
double p5[5] = { 1,-802,+161606,-324002,800005 } ; // 400+I, 400-I, 1+2*I, 1-2*I
double p6[5] = { 1, -2.000001000, 1.00000001e8,-100.000001 } ; // 1/1000000, 1+-10000*I
double p7[5] = { 1, -20000000001.0/1000000.0, 100000000019999/1000000.0, -99999999999999/1000000000000.0 } ; // 10000+-1/1000,-1/1000000

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


int
main() {
  do_test( p1 ) ;
  do_test( p2 ) ;
  do_test( p3 ) ;
  do_test( p4 ) ;
  do_test( p5 ) ;
  return 0 ;
}

