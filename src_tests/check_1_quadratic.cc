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

static double qq[][3] = {
  {  50,  50,   0 },
  {  50,  50,   1 },
  {  50,  50,  50 },
  {  50,  50,  99 },
  {  50,  50, 100 },
  {  50,   0,  50 },
  {  50,   1,  50 },
  {  50,  99,  50 },
  {  50, 100,  50 },
  {   0,  50,  50 },
  {   1,  50,  50 },
  {  99,  50,  50 },
  { 1e10, 50,  50 },
  { 100,  1e10,  50 },
  { 100,  50,  1e10 },
  { 1e-10, 50,  50 },
  { 100,  1e-10,  50 },
  { 100,  50,  1e-10 }
};

static
void
do_test( double const p[3] ) {
  Quadratic qsolve(p[0],p[1],p[2]);
  cout.precision(8);
  qsolve.info(cout);
  if ( !qsolve.check(cout) ) {
    cout << "\nFailed!\n\n\n";
    std::exit(0);
  }
}

int
main() {
  cout.precision(20);
  for ( int k = 0; k < sizeof(qq)/sizeof(qq[0]); ++k ) {
    cout << "\n\nExample N." << k << '\n';
    do_test(qq[k]);
  }
  cout << "\n\nALL DONE!\n";
  return 0;
}
