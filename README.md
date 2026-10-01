# Quartic Roots Solver

## Overview

`quarticRootsFlocke` is a C++ port of Norbert Flocke's cubic/quartic solver
for real-coefficient polynomials up to degree 4. The repository also contains
a real Jenkins-Traub implementation for higher-degree real polynomials.

The public header is `src/PolynomialRoots.hh`; `src/PolynomialRoots-complex.hxx`
is pulled in by it and provides the lightweight complex types used by the
solvers.

## Requirements

- a C++20 compiler **with `<format>`**: GCC >= 13, clang >= 17, MSVC >= 19.29
  (the build stops with an explicit message on older toolchains)
- CMake >= 3.16
- Boost headers, **only** when the optional multiprecision backend is enabled

## Usage

```cpp
#include "PolynomialRoots.hh"
#include <iostream>

using namespace PolynomialRoots;

double coeffs[] = { 8, -8, 16, -16, 8, -8 };
double zeror[5], zeroi[5];
int    degree = 5;

RootsError err = jenkins_traub_roots( coeffs, degree, zeror, zeroi );
if ( err != RootsError::SUCCESS ) {
  std::cerr << "failed: " << to_string(err) << '\n';
} else {
  for ( int i = 0; i < degree; ++i )
    std::cout << zeror[i] << " + I*" << zeroi[i] << '\n';
}
```

`roots()` is kept as a thin wrapper returning the same codes as plain integers.

For the specialized low-degree solvers:

```cpp
Quadratic q2(a,b,c);
Cubic     q3(a,b,c,d);
Quartic   q4(a,b,c,d,e);

q4.info( std::cout );          // human readable report
bool   ok = q4.check();        // silent quality check
double be = q4.backward_error();
```

### Reading the roots

Every solver stores its roots in a flat set of real components; the complex
pairs are reconstructed on demand:

```cpp
int n = q4.num_roots();              // real + complex
int r = q4.num_real_roots();
int c = q4.num_complex_roots();

for ( int i = 0; i < n; ++i ) {
  auto z = q4.root(i);               // complex_type, throws if i is out of range
  std::cout << z.real() << " + I*" << z.imag() << '\n';
}

double rr[4];
int nr = q4.get_real_roots( rr );    // only the real ones, ascending
nr = q4.get_positive_roots( rr );
nr = q4.get_negative_roots( rr );
nr = q4.get_roots_in_range( lo, hi, rr );
nr = q4.get_roots_in_open_range( lo, hi, rr );
```

Note that `real_root0()` / `real_root1()` return the *stored components*: when
the roots are a conjugate pair they are the real and the imaginary part of that
pair, not two real roots. Use `root(i)` or `get_real_roots()` when in doubt.

Non-finite coefficients and out-of-range root indices throw `std::runtime_error`.

### Multiple roots

```cpp
int m  = q4.root_multiplicity(i);     // multiplicity of the i-th root
int mx = q4.max_root_multiplicity();
bool d = q4.double_root();            // a root of multiplicity exactly 2
bool t = q4.triple_root();
bool f = q4.quadruple_root();         // Quartic only
```

In floating point two roots are "the same root" when no perturbation of the
last bits of the coefficients can tell them apart, so the test cannot be an
exact comparison nor a fixed tolerance. The library applies the classical
characterization of a root of multiplicity `m` at `c`, namely
`p(c) = p'(c) = ... = p^(m-1)(c) = 0`, and compares each of those Taylor
coefficients with the one of the polynomial built from the *absolute values* of
the coefficients, scaled by the unit roundoff:

```
|p^(j)(c)/j!|  <=  K * eps * A_j(|c|),    j = 0 ... m-1
```

Both sides come from the same Horner shift, so the criterion needs no
assumption on the magnitude of the roots and adapts by itself to the scalar
type: `machine_eps<T_real>()` is the only precision-dependent quantity, which is
why the very same code resolves roots `1e-10` apart with `QuarticHQ` while
`Quartic` reports them as a double root.

`is_multiple_root()` exposes the test and `multiple_root_separation()` the
companion estimate of how far such roots spread,
`(eps * S(c) / (|p^(m)(c)|/m!))^(1/m)` — the familiar rule that a double root
keeps half the significant digits and a triple root a third of them. The safety
factor is `multiple_root_safety_factor`; with its default value of 64, measured
over 100000 random polynomials spanning six decades of scale, every exact
multiple root of the tested families is recognized (worst case 99.5%) and no
multiplicity is reported for roots separated by more than `1e-5` in relative
terms. The band around `sqrt(eps)` is intrinsically ambiguous and no algorithm
can do better there.

### Accuracy

`check()` and `backward_error()` measure the **relative backward error**

```
|p(r)| / sum_k |a_k| |r|^(n-k)
```

rather than the bare residual `|p(r)|`, which grows with the magnitude of the
root and of the coefficients and is therefore meaningless as an absolute
quality threshold.

## Precision backends

The library is built around `double` by default:

- `Quadratic`, `Cubic`, `Quartic` - the standard `double` solvers
- `QuadraticHQ`, `CubicHQ`, `QuarticHQ` - the optional high-quality solvers

The `HQ` classes are the same class templates instantiated on
`boost::multiprecision::cpp_bin_float_100` (exposed as
`PolynomialRoots::quad_real`). They are compiled **only** when the project is
configured with `-DQRF_ENABLE_ARBITRARY_PRECISION=ON` *and* the Boost headers
are found; in that case the macro `POLYNOMIAL_ROOTS_HAS_MULTIPRECISION` is
defined to 1 for the library and for every consumer of the CMake target.

Boost is resolved, in order, from a local `../boost` checkout, from
`find_package(Boost)`, from `BOOST_ROOT` / `BOOST_INCLUDEDIR`, and finally from
the usual Homebrew prefixes (see `cmake/ResolveBoostMultiprecision.cmake`).

The Jenkins-Traub entry point is `double` only.

## Build And Test

```bash
rake build          # configure + build + install into ./lib
rake test           # build the tests and run CTest
```

or directly:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

### CMake options

| option | default | meaning |
|---|---|---|
| `QRF_ENABLE_TESTS` | `ON` | build the executables in `src_tests/` |
| `QRF_BUILD_SHARED` | `OFF` | also build a shared library |
| `QRF_INSTALL_LOCAL` | `ON` | install into `./lib` instead of the system prefix |
| `QRF_ENABLE_ARBITRARY_PRECISION` | `OFF` | enable the Boost-based `HQ` solvers |
| `QRF_BUILD_FAT` | `ON` on Apple | build universal binaries |
| `QRF_GIT_VERSION` | `ON` | take the project version from the Git tags |

The historical `UTILS_*` spellings are still accepted and seed the default of
the matching `QRF_*` option, with a deprecation message.

## Tests

The test executables use Unicode markers and `src_tests/termcolor.hh` for a
readable terminal output. When stdout is attached to a terminal, each suite
prints coloured case headers plus a final `✓` / `✗` summary; each regression
case also prints a comparison table with the `double` and, when enabled, the
`HQ` result side by side.

Current test programs:

- `check_1_quadratic`
- `check_2_cubic`
- `check_3_quartic`
- `check_4_classical_cases`
- `check_5_regressions`
- `check_6_literature`
- `check_Jenkins-Traub`
- `check_eval_derivative`

### Test families

- Flocke's published cubic and quartic stress cases.
- Extended scaled and clustered-root examples.
- Regression checks for `eval_poly_Dpoly` in both branches `|x| <= 1` and
  `|x| > 1`, and for `Newton_step`.
- The `max_degree` guard and the error codes of the Jenkins-Traub entry point.
- Classical literature-driven families:
  - Wilkinson-style factorized polynomials `W2` ... `W10`.
  - Multiple-root families such as `(x-1)^3`, `(x-1)^4`,
    `(x-1)^3(x+2)`, `(x-1)^2(x-3)^2`.
  - Biquadratic and repeated conjugate-pair quartics such as
    `(x^2-1)(x^2-4)` and `(x^2+1)^2`.
- `check_6_literature` collects the numerically difficult polynomials published
  in the literature, each with its reference next to the coefficients:

  | source | what it contributes |
  |---|---|
  | Forsythe, CS-TR-66-40, 1966 | the classical quadratics: cancellation of the small root, overflow of `b^2`, underflow of `c/(a r)`, nearly double root |
  | Kahan, *On the Cost of Floating-Point Computation…*, 2004 | the quadratics whose discriminant `B^2-AC` is lost to cancellation, including the integer one with discriminant exactly `-1` out of `8.1e31` |
  | Wilkinson, *Rounding Errors in Algebraic Processes*, 1963 | `W2`, `W3` |
  | Strobach, *The Low-Rank LDLT Quartic Solver*, 2015, §6.5 | the extreme benchmarks, reduced to cubics where applicable |
  | Orellana & De Michele, Algorithm 1010, ACM TOMS 46(2), 2020, Table 1 | the 24 reference quartics, from four roots spread over 150 decades to the quadruple root and the clustered complex pairs |

  The pass criterion is the relative backward error, the quantity a root finder
  is responsible for; the forward error against the nominal roots is reported
  for information, since several of those polynomials are so ill conditioned
  that their small roots are not recoverable from double precision coefficients
  by any algorithm. Those cases are marked with a warning rather than a failure.

- `check_5_regressions` pins every defect fixed in the solvers: the polynomials
  that used to make the quartic return four wrong roots, a 200k-sample
  randomized backward-error sweep, the quartics with a real root of even
  multiplicity (which used to be reported as having no real root at all), the
  `(x-r)^4` family (whose resolvent degenerated to `0/0` and produced infinite
  imaginary parts), the degenerate getters, multiplicity detection,
  `eval_poly_complex` in both branches, and the behaviour of the deprecated
  aliases.

## Deprecated API

The camelCase spellings (`getRealRoots`, `numRoots`, `getRoot0`, `evalMonicCubic`,
`MAXDEGREE`, `machepsiT`, `toleranceT`, ...) are marked `[[deprecated]]` and
forward to their snake_case counterparts (`get_real_roots`, `num_roots`,
`get_root0`, `eval_monic_cubic`, `max_degree`, `machine_eps`,
`residual_tolerance`). They still work; they will emit a compiler warning.

Three further changes may affect existing code:

- `double_root()` and `triple_root()` no longer come from the exact tests used
  by the solvers internally, which missed most multiple roots; they now follow
  the criterion described above. Code relying on the previous, essentially
  always false, behaviour will see different answers — for the better.

- `PolynomialRoots.hh` no longer includes `<iostream>`. Translation units that
  relied on that transitive include must include it themselves.
- `roots()` used to return `-2` both for a zero leading coefficient and for a
  convergence failure; the latter is now `-4` (`RootsError::NO_CONVERGENCE`).

## Third-Party Components

- `src_tests/termcolor.hh`
  - upstream: [ikalnytskyi/termcolor](https://github.com/ikalnytskyi/termcolor)
  - license: BSD 3-Clause
  - local copy: `licenses3rd/LICENSE-termcolor.txt`

- optional high-precision backend
  - upstream: [Boost.Multiprecision](https://github.com/boostorg/multiprecision)
  - license: Boost Software License 1.0
  - local copy: `licenses3rd/LICENSE-boost.txt`
  - header only, used only when `QRF_ENABLE_ARBITRARY_PRECISION=ON`

## Repository And Online Documentation

- [github](https://github.com/ebertolazzi/quarticRootsFlocke)
- [online doc](http://ebertolazzi.github.io/quarticRootsFlocke)

The `docs/` directory holds the generated Doxygen output published through
GitHub Pages; it is regenerated from `doxygen/`.

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
