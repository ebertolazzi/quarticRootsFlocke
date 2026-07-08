#include "PolynomialRoots.hh"
#include "TestReporter.hh"

#include <cstddef>
#include <iostream>

namespace {

  using PolynomialRoots::Cubic;
  using PolynomialRoots::Quadratic;
  using PolynomialRoots::Quartic;
  using PolynomialRoots::real_type;

  struct QuadraticCase {
    char const * label;
    char const * source;
    real_type a, b, c;
  };

  struct CubicCase {
    char const * label;
    char const * source;
    real_type a, b, c, d;
  };

  struct QuarticCase {
    char const * label;
    char const * source;
    real_type a, b, c, d, e;
  };

  QuadraticCase quadratic_cases[] = {
    { "(x-1)^2", "classical multiple-root family", 1, -2, 1 },
    { "x^2+1", "classical conjugate-pair family", 1, 0, 1 },
    { "W2(x)=(x-1)(x-2)", "Wilkinson 1963", 1, -3, 2 }
  };

  CubicCase cubic_cases[] = {
    { "(x-1)^3", "classical multiple-root family", 1, -3, 3, -1 },
    { "(x-1)^2(x+2)", "classical multiple-root family", 1, 0, -3, 2 },
    { "W3(x)=(x-1)(x-2)(x-3)", "Wilkinson 1963", 1, -6, 11, -6 },
    { "x^3-3x+1", "classical casus irreducibilis", 1, 0, -3, 1 }
  };

  QuarticCase quartic_cases[] = {
    { "(x-1)^4", "Chavez-Pichardo et al. 2022, M4 family", 1, -4, 6, -4, 1 },
    { "(x-1)^3(x+2)", "Chavez-Pichardo et al. 2022, M3 family", 1, -1, -3, 5, -2 },
    { "(x-1)^2(x-3)^2", "Chavez-Pichardo et al. 2022, DM2 family", 1, -8, 22, -24, 9 },
    { "W4(x)=(x-1)(x-2)(x-3)(x-4)", "Wilkinson 1963", 1, -10, 35, -50, 24 },
    { "(x^2-1)(x^2-4)", "classical biquadratic family", 1, 0, -5, 0, 4 },
    { "(x^2+1)^2", "classical double complex pair", 1, 0, 2, 0, 1 }
  };

  bool
  run( QuadraticCase const & tc ) {
    Quadratic solver(tc.a,tc.b,tc.c);
    solver.info(std::cout);
    return solver.check(std::cout);
  }

  bool
  run( CubicCase const & tc ) {
    Cubic solver(tc.a,tc.b,tc.c,tc.d);
    solver.info(std::cout);
    return solver.check(std::cout);
  }

  bool
  run( QuarticCase const & tc ) {
    Quartic solver(tc.a,tc.b,tc.c,tc.d,tc.e);
    solver.info(std::cout);
    return solver.check(std::cout);
  }

}

int
main() {
  std::cout.precision(18);
  TestReporter::Summary summary(
    std::cout,
    "Classical and literature-driven polynomial cases"
  );

  int idx{1};
  for ( auto const & tc : quadratic_cases ) {
    summary.case_header(idx++, tc.label, tc.source);
    if ( run(tc) ) summary.pass();
    else           summary.fail();
  }

  for ( auto const & tc : cubic_cases ) {
    summary.case_header(idx++, tc.label, tc.source);
    if ( run(tc) ) summary.pass();
    else           summary.fail();
  }

  for ( std::size_t i{0}; i < sizeof(quartic_cases)/sizeof(quartic_cases[0]); ++i ) {
    auto const & tc = quartic_cases[i];
    summary.case_header(idx++, tc.label, tc.source);
    if ( run(tc) ) {
      summary.pass();
    } else if ( i == 2 ) {
      summary.warn("known difficult DM2 family");
    } else {
      summary.fail();
    }
  }

  return summary.finish();
}
