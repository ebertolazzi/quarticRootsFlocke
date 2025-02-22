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
  Quadratic const qsolve(p[0],p[1],p[2]);
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
  constexpr integer N{static_cast<integer>(sizeof(qq) / sizeof(qq[0]))};
  for ( integer k{0}; k < N; ++k ) {
    cout << "\n\nExample N." << k << '\n';
    do_test(qq[k]);
  }
  cout << "\n\nALL DONE!\n";
  return 0;
}
