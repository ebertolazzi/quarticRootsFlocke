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

#include <iomanip>
#include <iostream>

using PolynomialRoots::Quadratic;
using PolynomialRoots::integer;

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
bool
do_test( double const p[3] ) {
  Quadratic const qsolve(p[0],p[1],p[2]);
  std::cout.precision(8);
  qsolve.info(std::cout);
  return qsolve.check(std::cout);
}

int
main() {
  std::cout.precision(20);
  TestReporter::Summary summary(
    std::cout,
    "Quadratic solver regression suite"
  );
  constexpr integer N{static_cast<integer>(sizeof(qq) / sizeof(qq[0]))};
  for ( integer k{0}; k < N; ++k ) {
    summary.case_header(
      static_cast<int>(k+1),
      "scaled quadratic coefficients",
      "Flocke-style regression set"
    );
    if ( do_test(qq[k]) ) summary.pass();
    else                  summary.fail();
  }
  return summary.finish();
}
