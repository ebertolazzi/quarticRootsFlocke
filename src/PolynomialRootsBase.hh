/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2014                                                      |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef POLYNOMIAL_ROOTS_BASE_HH
#define POLYNOMIAL_ROOTS_BASE_HH

#include <complex>
#include <cstdint>
#include <iosfwd>
#include <limits>

namespace PolynomialRoots {

  using real_type    = double;
  using integer      = int;
  using complex_type = std::complex<real_type>;
  using ostream_type = std::basic_ostream<char>;
  using istream_type = std::basic_istream<char>;

  inline constexpr integer MAXDEGREE{100};
  inline constexpr int bitsValueType = std::numeric_limits<real_type>::digits;
  inline constexpr real_type splitFactor =
    static_cast<real_type>((std::uint64_t(1) << (bitsValueType - 2)) + 1);

} // namespace PolynomialRoots

#endif
