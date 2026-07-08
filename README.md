# Quartic Roots Solver

## Overview

`quarticRootsFlocke` is a C++ port of Norbert Flocke's cubic/quartic solver
for real-coefficient polynomials up to degree 4. The repository also contains
an experimental real Jenkins-Traub implementation for higher-degree real
polynomials.

## Usage

```cpp
#include "PolynomialRoots.hh"

double coeffs[] = { 8, -8, 16, -16, 8, -8 };
double zeror[5], zeroi[5];
int degree = 5;

int ok = PolynomialRoots::roots(coeffs, degree, zeror, zeroi);
std::cout << "ok = " << ok << '\n';
for ( int i = 0; i < degree; ++i ) {
  std::cout << zeror[i] << " + I*" << zeroi[i] << '\n';
}
```

For specialized low-degree solvers:

```cpp
Quadratic q2(a,b,c);
q2.info(std::cout);

Cubic q3(a,b,c,d);
q3.info(std::cout);

Quartic q4(a,b,c,d,e);
q4.info(std::cout);
```

See the class interfaces in `src/PolynomialRoots.hh` and
`src/PolynomialRootsHQ.hh` for the available accessors.

The library now exposes separate solver families:

- `Quadratic`, `Cubic`, `Quartic` for the standard `double` implementation
- `QuadraticHQ`, `CubicHQ`, `QuarticHQ` for the high-quality multiprecision path

```cpp
Quartic q4(a,b,c,d,e);      // double solver
QuarticHQ q4hq(a,b,c,d,e);  // high-quality multiprecision solver
```

The public root interface remains in `double`, while the `HQ` classes use an
internal multiprecision backend to improve the solve/refinement phase.

Internally the current `HQ` classes use `stillwater-sc/universal` and its
`sw::universal::quad` type. The dependency is resolved in this order:

- local sibling checkout: `../universal`
- `FetchContent` fallback from the upstream repository

## Build And Test

```bash
rake build
rake test
```

The test executables now use Unicode markers and `termcolor.hh` for a more
readable terminal output. When stdout is attached to a terminal, each suite
prints colored case headers plus a final `✓` / `✗` summary. Each regression
case also prints a small comparison table with the `double` and `HQ` result
status side by side, including the maximum residual for each class. Legacy
numerically difficult cases are shown as warnings (`⚠`) instead of being
silently ignored. Every quadratic, cubic and quartic regression case is run
with separate solver classes:

- standard `double` solver classes
- `HQ` multiprecision solver classes

Current test programs:

- `check_1_quadratic`
- `check_2_cubic`
- `check_3_quartic`
- `check_4_classical_cases`
- `check_eval_derivative`

## Test Families

The current regression suite covers:

- Flocke's published cubic and quartic stress cases.
- Extended scaled and clustered-root examples already present in this repository.
- Regression checks for `evalPolyDPoly` in both branches `|x| <= 1` and `|x| > 1`.
- The `MAXDEGREE` guard in the Jenkins-Traub entry point.
- Classical literature-driven families:
  - Wilkinson-style factorized polynomials `W2`, `W3`, `W4`.
  - Multiple-root families such as `(x-1)^3`, `(x-1)^4`,
    `(x-1)^3(x+2)`, `(x-1)^2(x-3)^2`.
  - Biquadratic and repeated conjugate-pair quartics such as
    `(x^2-1)(x^2-4)` and `(x^2+1)^2`.

## Third-Party Components

- `src/termcolor.hh`
  - upstream: [ikalnytskyi/termcolor](https://github.com/ikalnytskyi/termcolor)
  - license: BSD 3-Clause
  - local copy: `licenses3rd/LICENSE-termcolor.txt`

- internal high-precision backend
  - upstream: [stillwater-sc/universal](https://github.com/stillwater-sc/universal)
  - license: MIT
  - local copy: `licenses3rd/LICENSE-universal.txt`

## Build Dependencies

- C++20 compiler
- `stillwater-sc/universal`
  - used only for the internal high-precision backend
  - not exposed by the public headers
  - resolved from `../universal` when available
  - otherwise downloaded with `FetchContent` during CMake configuration

## Repository And Online Documentation

- [github](https://github.com/ebertolazzi/quarticRootsFlocke)
- [online doc](http://ebertolazzi.github.io/quarticRootsFlocke)

## References

- **N. Flocke**  
  Algorithm 954: An Accurate and Efficient Cubic and Quartic Equation Solver
  for Physical Applications  
  ACM TOMS, 41(4), 2015.

- **M. A. Jenkins and J. F. Traub**  
  A Three-Stage Algorithm for Real Polynomials Using Quadratic Iteration  
  SIAM Journal on Numerical Analysis, 7(4), 1970.

- **J. H. Wilkinson**  
  Rounding Errors in Algebraic Processes  
  Prentice Hall, 1963.

- **M. Chavez-Pichardo, M. A. Martinez-Cruz, A. Trejo-Martinez,
  D. Martinez-Carbajal, T. Arenas-Resendiz**  
  A Complete Review of the General Quartic Equation with Real Coefficients
  and Multiple Roots  
  Mathematics, 2022.

## Author

Enrico Bertolazzi  
Dipartimento di Ingegneria Industriale  
Universita degli Studi di Trento  
email: enrico.bertolazzi@unitn.it
