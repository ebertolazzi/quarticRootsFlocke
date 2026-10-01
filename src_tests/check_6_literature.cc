//
// Numerically difficult quadratics, cubics and quartics taken from the
// literature.  Each case carries the reference it comes from next to its data.
//
// References used here:
//
//   [F66]  G. E. Forsythe, "How do you solve a quadratic equation?",
//          Technical Report CS-TR-66-40, Computer Science Department,
//          Stanford University, 1966.
//
//   [K04]  W. Kahan, "On the Cost of Floating-Point Computation Without
//          Extra-Precise Arithmetic", 2004.
//          https://people.eecs.berkeley.edu/~wkahan/Qdrtcs.pdf
//
//   [W63]  J. H. Wilkinson, "Rounding Errors in Algebraic Processes",
//          Prentice Hall, 1963.
//
//   [F15]  N. Flocke, "Algorithm 954: An Accurate and Efficient Cubic and
//          Quartic Equation Solver for Physical Applications",
//          ACM TOMS 41(4), 2015.  DOI 10.1145/2699468
//
//   [S15]  P. Strobach, "The Low-Rank LDLT Quartic Solver", AST Consulting
//          technical report, 2015, section 6.5 (extreme benchmark tests).
//
//   [ODM]  A. G. Orellana, C. De Michele, "Algorithm 1010: Boosting Efficiency
//          in Solving Quartic Equations with No Compromise in Accuracy",
//          ACM TOMS 46(2), 2020, Table 1.  DOI 10.1145/3386241
//
// The pass criterion is the relative BACKWARD error
//
//     |p(r)| / sum_k |a_k| |r|^(n-k)
//
// which is the quantity a root finder is responsible for.  The forward error
// against the nominal roots is reported as information only: several of these
// polynomials are so ill conditioned that their small roots are simply not
// recoverable from double precision coefficients, no matter the algorithm.
//

#include "PolynomialRoots.hh"
#include "TestHelpers.hh"
#include "TestReporter.hh"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{

  using PolynomialRoots::Cubic;
  using PolynomialRoots::Quadratic;
  using PolynomialRoots::Quartic;
  using Cplx = std::complex<double>;

  // -------------------------------------------------------------------------
  // relative backward error of z as a root of op[0..n]
  // -------------------------------------------------------------------------
  //! Both the numerator and the denominator are divided by |z|^n when |z| > 1,
  //! so that roots of modulus 1e150 do not overflow the comparison.
  double backward_error( std::vector<double> const & op, Cplx const & z )
  {
    std::size_t const n{ op.size() - 1 };
    double const      az{ std::abs( z ) };
    double            den{ 0 }, pw{ 1 };
    Cplx              v;
    if ( az > 1 )
    {
      v = Cplx{ op[n], 0.0 };  // reversed polynomial at 1/z, that is p(z)/z^n
      for ( std::size_t i{ 1 }; i <= n; ++i ) v = v / z + op[n - i];
      for ( std::size_t k{ 0 }; k <= n; ++k )
      {
        den += std::abs( op[k] ) * pw;
        pw /= az;
      }
    }
    else
    {
      v = Cplx{ op[0], 0.0 };
      for ( std::size_t k{ 1 }; k <= n; ++k ) v = v * z + op[k];
      for ( std::size_t k{ n + 1 }; k-- > 0; )
      {
        den += std::abs( op[k] ) * pw;
        pw *= az;
      }
    }
    if ( !( den > 0 ) ) return 0;
    return std::abs( v ) / den;
  }

  //! Greedy matching of the computed roots against the nominal ones; returns the
  //! largest relative distance.
  double forward_error( std::vector<Cplx> const & exact, std::vector<Cplx> const & got )
  {
    if ( exact.size() != got.size() ) return std::numeric_limits<double>::infinity();
    std::vector<bool> used( got.size(), false );
    double            worst{ 0 };
    for ( auto const & e : exact )
    {
      double      best{ std::numeric_limits<double>::infinity() };
      std::size_t bi{ 0 };
      for ( std::size_t j{ 0 }; j < got.size(); ++j )
      {
        if ( used[j] ) continue;
        double const d{ std::abs( e - got[j] ) };
        if ( d < best )
        {
          best = d;
          bi   = j;
        }
      }
      used[bi] = true;
      double const scale{ std::max( 1.0, std::abs( e ) ) };
      worst = std::max( worst, best / scale );
    }
    return worst;
  }

  template <typename Solver> std::vector<Cplx> roots_of( Solver const & s )
  {
    std::vector<Cplx> v;
    for ( PolynomialRoots::integer i{ 0 }; i < s.num_roots(); ++i )
    {
      auto const z{ s.root( i ) };
      v.emplace_back( z.real(), z.imag() );
    }
    return v;
  }

  //! Expand prod (x - r_k) into real coefficients, exactly as [ODM] eq. (109).
  std::vector<double> coefficients_from_roots( std::vector<Cplx> const & r )
  {
    std::vector<Cplx> c{ Cplx{ 1.0, 0.0 } };
    for ( auto const & z : r )
    {
      std::vector<Cplx> n( c.size() + 1, Cplx{ 0.0, 0.0 } );
      for ( std::size_t i{ 0 }; i < c.size(); ++i )
      {
        n[i] += c[i];
        n[i + 1] -= c[i] * z;
      }
      c.swap( n );
    }
    std::vector<double> op;
    op.reserve( c.size() );
    for ( auto const & z : c ) op.push_back( z.real() );
    return op;
  }

  // -------------------------------------------------------------------------
  // reporting
  // -------------------------------------------------------------------------
  constexpr double BACKWARD_TOL{ 1e-13 };

  struct Totals
  {
    int total{ 0 };
    int failed{ 0 };
  };

  void report(
    TestReporter::Summary &     summary,
    int &                       id,
    Totals &                    tot,
    char const *                label,
    char const *                reference,
    std::vector<double> const & op,
    std::vector<Cplx> const &   computed,
    std::vector<Cplx> const &   nominal )
  {
    summary.case_header( id++, label, reference );

    double be{ 0 };
    for ( auto const & z : computed ) be = std::max( be, backward_error( op, z ) );
    double const fe{ nominal.empty() ? -1.0 : forward_error( nominal, computed ) };

    bool const ok{ be <= BACKWARD_TOL };
    ++tot.total;
    if ( !ok ) ++tot.failed;

    summary.comparison_table_header();
    summary.comparison_table_row( "backward", ok, be );
    if ( fe >= 0 ) summary.comparison_table_row( "forward", fe < 1e-8, fe );
    summary.comparison_table_footer();

    if ( ok )
    {
      if ( fe >= 1e-8 )
        summary.warn( "backward stable, forward error limited by the conditioning" );
      else
        summary.pass( "backward and forward error within tolerance" );
    }
    else
    {
      summary.fail( "backward error above tolerance" );
    }
  }

  // =========================================================================
  // QUADRATICS
  //
  // Several references state the quadratic as A x^2 - 2 B x + C, so that the
  // discriminant is simply B^2 - A C.  The coefficient of x is then -2B.
  // =========================================================================
  struct QuadCase
  {
    char const * label;
    char const * reference;
    double       a, b, c;  // a x^2 + b x + c
    Cplx         r0, r1;   // nominal roots, (0,0) twice when unknown
    bool         have_roots;
  };

  QuadCase const quadratic_cases[]{
    // [F66] p.3, the example that "causes no difficulty for a computer"
    { "6 x^2 + 5 x - 4", "Forsythe 1966, easy baseline", 6, 5, -4, Cplx{ -4.0 / 3.0, 0 }, Cplx{ 0.5, 0 }, true },

    // [F66] the classical cancellation example: computing the small root as
    // (-b + sqrt(b^2-4ac)) / (2a) loses it completely
    { "x^2 - 1e5 x + 1",
      "Forsythe 1966, catastrophic cancellation",
      1,
      -1e5,
      1,
      Cplx{ 1.0000000000100000000002e-5, 0 },
      Cplx{ 99999.99999, 0 },
      true },

    // [F66] same roots as the baseline, but b^2 overflows in single precision;
    // scaled here to 1e150 so that b^2 overflows in double as well
    { "6e150 x^2 + 5e150 x - 4e150",
      "Forsythe 1966, overflow of b^2 (rescaled)",
      6e150,
      5e150,
      -4e150,
      Cplx{ -4.0 / 3.0, 0 },
      Cplx{ 0.5, 0 },
      true },

    // [F66] the two roots differ by 1e200 in magnitude: computing the smaller
    // one as c/(a*larger) underflows unless the divisions are reordered
    { "1e-200 x^2 + x - 1e-200",
      "Forsythe 1966, underflow of the small root (rescaled)",
      1e-200,
      1.0,
      -1e-200,
      Cplx{ -1e200, 0 },
      Cplx{ 1e-200, 0 },
      true },

    // [F66] the roots agree to eight digits: the discriminant 4e-7 is the
    // difference of two numbers of size 16
    { "x^2 - 4 x + 3.9999999",
      "Forsythe 1966, nearly double root",
      1,
      -4.0,
      3.9999999,
      Cplx{ 2.0 - 3.1622776601683794e-4, 0 },
      Cplx{ 2.0 + 3.1622776601683794e-4, 0 },
      true },

    // [K04] Table 1: D = B^2 - A C = 0.0022, lost to cancellation in double
    { "10.27 x^2 - 2*29.61 x + 85.37",
      "Kahan 2004, Table 1 (first)",
      10.27,
      -2 * 29.61,
      85.37,
      Cplx{ 2.8785871686049835, 0 },
      Cplx{ 2.8877206336472626, 0 },
      true },

    // [K04] Table 1: D = 0.0492
    { "10.28 x^2 - 2*29.62 x + 85.34",
      "Kahan 2004, Table 1 (second)",
      10.28,
      -2 * 29.62,
      85.34,
      Cplx{ 2.8607185104906377, 0 },
      Cplx{ 2.9028988317431196, 0 },
      true },

    // [K04] Table 2: A, B, C are exactly representable, D = B^2 - A C = 121/64
    // and one root is exactly 1.  Computing B*B - A*C in double gives 2 instead
    // of 1.890625 unless the product difference is evaluated carefully.
    { "94906265.625 x^2 - 189812534 x + 94906268.375",
      "Kahan 2004, Table 2",
      94906265.625,
      -189812534.0,
      94906268.375,
      Cplx{ 1.0, 0 },
      Cplx{ 1.0000000289759583, 0 },
      true },

    // [K04] p.6: integer coefficients, exactly representable, with the tiniest
    // possible negative discriminant B^2 - A C = -1 out of about 8.1e31
    { "16675184804162890 x^2 - 2*9007199254740987 x + 4865291711450953",
      "Kahan 2004, discriminant = -1 exactly",
      16675184804162890.0,
      -2.0 * 9007199254740987.0,
      4865291711450953.0,
      Cplx{ 0.5401558879570785, 5.996935036967951e-17 },
      Cplx{ 0.5401558879570785, -5.996935036967951e-17 },
      true },

    // [W63] the degree two Wilkinson polynomial
    { "W2 = (x-1)(x-2)", "Wilkinson 1963", 1, -3, 2, Cplx{ 1, 0 }, Cplx{ 2, 0 }, true },
  };

  // =========================================================================
  // CUBICS, given by their roots.  The eight [F15] cases already live in
  // check_2_cubic.cc, so the families below complete the picture.
  // =========================================================================
  struct CubicCase
  {
    char const * label;
    char const * reference;
    Cplx         r[3];
  };

  constexpr double S15{ 1e15 };

  CubicCase const cubic_cases[]{
    // cubic reduction of the Kahan benchmark of [S15] sec. 6.5.3, whose quartic
    // version has roots [-S, 1/S, 1, S]: thirty orders of magnitude of spread
    { "-S, 1/S, S with S = 1e15",
      "Strobach 2015 sec. 6.5.3, cubic reduction",
      { Cplx{ -S15, 0 }, Cplx{ 1.0 / S15, 0 }, Cplx{ S15, 0 } } },

    // cubic reduction of [ODM] Table 1 case 1
    { "1e9, 1e6, 1e3",
      "Orellana & De Michele 2020, case 1 reduced",
      { Cplx{ 1e9, 0 }, Cplx{ 1e6, 0 }, Cplx{ 1e3, 0 } } },

    // cubic reduction of [ODM] Table 1 case 2
    { "2.002, 2.001, 2",
      "Orellana & De Michele 2020, case 2 reduced",
      { Cplx{ 2.002, 0 }, Cplx{ 2.001, 0 }, Cplx{ 2.0, 0 } } },

    // triple root, cubic version of the [S15] sec. 6.5.1 benchmark
    { "1000, 1000, 1000",
      "Strobach 2015 sec. 6.5.1, cubic reduction",
      { Cplx{ 1000, 0 }, Cplx{ 1000, 0 }, Cplx{ 1000, 0 } } },

    // double root plus a tiny one, cubic version of [S15] sec. 6.5.2
    { "1000, 1000, 1e-15",
      "Strobach 2015 sec. 6.5.2, cubic reduction",
      { Cplx{ 1000, 0 }, Cplx{ 1000, 0 }, Cplx{ 1e-15, 0 } } },

    // cubic reduction of [ODM] Table 1 case 12: a conjugate pair whose
    // imaginary part is eleven orders of magnitude below the real one
    { "1e3+i, 1e3-i, 1e3",
      "Orellana & De Michele 2020, case 12 reduced",
      { Cplx{ 1e3, 1 }, Cplx{ 1e3, -1 }, Cplx{ 1e3, 0 } } },

    // cubic reduction of [ODM] Table 1 case 13
    { "2+1e4 i, 2-1e4 i, 1",
      "Orellana & De Michele 2020, case 13 reduced",
      { Cplx{ 2, 1e4 }, Cplx{ 2, -1e4 }, Cplx{ 1, 0 } } },

    // cubic reduction of [ODM] Table 1 case 6
    { "1e7, 1+i, 1-i", "Orellana & De Michele 2020, case 6 reduced", { Cplx{ 1e7, 0 }, Cplx{ 1, 1 }, Cplx{ 1, -1 } } },

    // the casus irreducibilis: three real roots that Cardano's formula can only
    // reach through complex arithmetic
    { "x^3 - 3x + 1",
      "classical casus irreducibilis",
      { Cplx{ 2 * std::cos( 2 * 3.14159265358979323846 / 9 ), 0 },
        Cplx{ 2 * std::cos( 8 * 3.14159265358979323846 / 9 ), 0 },
        Cplx{ 2 * std::cos( 14 * 3.14159265358979323846 / 9 ), 0 } } },

    // [W63] degree three
    { "W3 = (x-1)(x-2)(x-3)", "Wilkinson 1963", { Cplx{ 1, 0 }, Cplx{ 2, 0 }, Cplx{ 3, 0 } } },
  };

  // =========================================================================
  // QUARTICS: Table 1 of [ODM], given by their roots as in the paper.
  // =========================================================================
  struct QuarticCase
  {
    char const * label;
    char const * reference;
    Cplx         r[4];
  };

  constexpr double S30{ 1e30 };

  QuarticCase const quartic_cases[]{
    // Orellana & De Michele 2020, Table 1, case 1: 1e9, 1e6, 1e3, 1
    { "ODM case 01", "four widely spaced real", { { 1000000000.0, 0 }, { 1000000.0, 0 }, { 1000.0, 0 }, { 1, 0 } } },
    // Orellana & De Michele 2020, Table 1, case 2: 2.003, 2.002, 2.001, 2
    { "ODM case 02", "four closely spaced small real", { { 2.003, 0 }, { 2.002, 0 }, { 2.001, 0 }, { 2.0, 0 } } },
    // Orellana & De Michele 2020, Table 1, case 3: 1e53, 1e50, 1e49, 1e47
    { "ODM case 03", "four large real", { { 1e+53, 0 }, { 1e+50, 0 }, { 1e+49, 0 }, { 1e+47, 0 } } },
    // Orellana & De Michele 2020, Table 1, case 4: 1e14, 2, 1, -1
    { "ODM case 04", "one large, three small real", { { 100000000000000.0, 0 }, { 2, 0 }, { 1, 0 }, { -1, 0 } } },
    // Orellana & De Michele 2020, Table 1, case 5: -2e7, 1e7, 1, -1
    { "ODM case 05", "two large, two small real", { { -20000000.0, 0 }, { 10000000.0, 0 }, { 1, 0 }, { -1, 0 } } },
    // Orellana & De Michele 2020, Table 1, case 6: 1e7, -1e6, 1+i, 1-i
    { "ODM case 06",
      "two large real, two small complex",
      { { 10000000.0, 0 }, { -1000000.0, 0 }, { 1, 1 }, { 1, -1 } } },
    // Orellana & De Michele 2020, Table 1, case 7: -7, -4, -1e6+1e5 i, -1e6-1e5 i
    { "ODM case 07",
      "two small real, two large complex",
      { { -7, 0 }, { -4, 0 }, { -1000000.0, 100000.0 }, { -1000000.0, -100000.0 } } },
    // Orellana & De Michele 2020, Table 1, case 8: 1e8, 11, 1e3+i, 1e3-i
    { "ODM case 08",
      "one large and one small real, two medium complex",
      { { 100000000.0, 0 }, { 11, 0 }, { 1000.0, 1 }, { 1000.0, -1 } } },
    // Orellana & De Michele 2020, Table 1, case 9: 1e7+1e6 i, 1e7-1e6 i, 1+2i, 1-2i
    { "ODM case 09",
      "two large, two small complex",
      { { 10000000.0, 1000000.0 }, { 10000000.0, -1000000.0 }, { 1, 2 }, { 1, -2 } } },
    // Orellana & De Michele 2020, Table 1, case 10: 1e4+3i, 1e4-3i, -7+1e3 i, -7-1e3 i
    { "ODM case 10",
      "four complex, mixed magnitudes",
      { { 10000.0, 3 }, { 10000.0, -3 }, { -7, 1000.0 }, { -7, -1000.0 } } },
    // Orellana & De Michele 2020, Table 1, case 11: 1.001+4.998i, 1.001-4.998i, 1+5.001i, 1-5.001i
    { "ODM case 11",
      "four closely spaced complex",
      { { 1.001, 4.998 }, { 1.001, -4.998 }, { 1, 5.001 }, { 1, -5.001 } } },
    // Orellana & De Michele 2020, Table 1, case 12: 1e3+3i, 1e3-3i, 1e3+i, 1e3-i
    { "ODM case 12",
      "four complex, equal real parts, small imaginary parts",
      { { 1000.0, 3 }, { 1000.0, -3 }, { 1000.0, 1 }, { 1000.0, -1 } } },
    // Orellana & De Michele 2020, Table 1, case 13: 2+1e4 i, 2-1e4 i, 1+1e3 i, 1-1e3 i
    { "ODM case 13",
      "four complex, small real, large imaginary parts",
      { { 2, 10000.0 }, { 2, -10000.0 }, { 1, 1000.0 }, { 1, -1000.0 } } },
    // Orellana & De Michele 2020, Table 1, case 14: 1000 (quadruple)
    { "ODM case 14", "quadruple root", { { 1000, 0 }, { 1000, 0 }, { 1000, 0 }, { 1000, 0 } } },
    // Orellana & De Michele 2020, Table 1, case 15: 1000, 1000, 1000, 1e-15
    { "ODM case 15", "triple root plus one tiny root", { { 1000, 0 }, { 1000, 0 }, { 1000, 0 }, { 1e-15, 0 } } },
    // Orellana & De Michele 2020, Table 1, case 16: 1+0.1i, 1-0.1i, 1e16+1e7 i, 1e16-1e7 i
    { "ODM case 16",
      "four complex, widely spaced",
      { { 1, 0.1 }, { 1, -0.1 }, { 1e+16, 10000000.0 }, { 1e+16, -10000000.0 } } },
    // Orellana & De Michele 2020, Table 1, case 17: 10000, 10001, 10010, 10100
    { "ODM case 17", "four closely spaced large real", { { 10000, 0 }, { 10001, 0 }, { 10010, 0 }, { 10100, 0 } } },
    // Orellana & De Michele 2020, Table 1, case 18: 40000+300i, 40000-300i, 30000+7000i, 30000-7000i
    { "ODM case 18",
      "four quite large complex",
      { { 40000, 300 }, { 40000, -300 }, { 30000, 7000 }, { 30000, -7000 } } },
    // Orellana & De Michele 2020, Table 1, case 19: 1e44, 1e30, 1e30, 1
    { "ODM case 19",
      "four very widely spaced real, one double",
      { { 1e+44, 0 }, { 1e+30, 0 }, { 1e+30, 0 }, { 1, 0 } } },
    // Orellana & De Michele 2020, Table 1, case 20: 1e14, 1e7, 1e7, 1
    { "ODM case 20",
      "four widely spaced real, one double",
      { { 100000000000000.0, 0 }, { 10000000.0, 0 }, { 10000000.0, 0 }, { 1, 0 } } },
    // Orellana & De Michele 2020, Table 1, case 21: 1e15, 1e7, 1e7, 1
    { "ODM case 21",
      "four widely spaced real, one double",
      { { 1000000000000000.0, 0 }, { 10000000.0, 0 }, { 10000000.0, 0 }, { 1, 0 } } },
    // Orellana & De Michele 2020, Table 1, case 22: 1e154, 1e152, 10, 1
    { "ODM case 22", "two very large real roots", { { 1e+154, 0 }, { 1e+152, 0 }, { 10, 0 }, { 1, 0 } } },
    // [ODM] Table 1, case 24 = the Kahan benchmark of [S15] sec. 6.5.3 with
    // S = 1e30: roots -S, 1/S, 1, S, sixty orders of magnitude of spread
    { "ODM case 24",
      "widely spaced coefficients, S = 1e30",
      { { -S30, 0 }, { 1.0 / S30, 0 }, { 1.0, 0 }, { S30, 0 } } },
  };

  // [ODM] Table 1, case 23 is given by its coefficients rather than its roots:
  // x^4 + x^3 + x^2 + (3/8) x + 1e-3, "four distinct roots with d2 about 0"
  double const odm_case_23[]{ 1.0, 1.0, 1.0, 3.0 / 8.0, 1e-3 };

  // -------------------------------------------------------------------------
  template <typename Solver> Solver make_solver( std::vector<double> const & op );

  template <> Quadratic make_solver<Quadratic>( std::vector<double> const & op )
  { return Quadratic( op[0], op[1], op[2] ); }

  template <> Cubic make_solver<Cubic>( std::vector<double> const & op )
  { return Cubic( op[0], op[1], op[2], op[3] ); }

  template <> Quartic make_solver<Quartic>( std::vector<double> const & op )
  { return Quartic( op[0], op[1], op[2], op[3], op[4] ); }

  template <typename Solver, typename CaseType>
  void run_from_roots( TestReporter::Summary & summary, int & id, Totals & tot, CaseType const & tc, int nr )
  {
    std::vector<Cplx> const nominal( tc.r, tc.r + nr );
    auto const              op{ coefficients_from_roots( nominal ) };
    Solver const            solver{ make_solver<Solver>( op ) };
    report( summary, id, tot, tc.label, tc.reference, op, roots_of( solver ), nominal );
  }

}  // namespace

int main()
{
  TestReporter::Summary summary( std::cout, "Numerically difficult cases from the literature" );
  Totals                tot;
  int                   id{ 1 };

  summary.note( "quadratics" );
  for ( auto const & tc : quadratic_cases )
  {
    std::vector<double> const op{ tc.a, tc.b, tc.c };
    Quadratic const           q{ make_solver<Quadratic>( op ) };
    std::vector<Cplx>         nominal;
    if ( tc.have_roots ) nominal = { tc.r0, tc.r1 };
    report( summary, id, tot, tc.label, tc.reference, op, roots_of( q ), nominal );
  }

  summary.note( "cubics" );
  for ( auto const & tc : cubic_cases ) run_from_roots<Cubic>( summary, id, tot, tc, 3 );

  summary.note( "quartics" );
  for ( auto const & tc : quartic_cases ) run_from_roots<Quartic>( summary, id, tot, tc, 4 );

  {
    std::vector<double> const op( odm_case_23, odm_case_23 + 5 );
    Quartic const             p( op[0], op[1], op[2], op[3], op[4] );
    report( summary, id, tot, "ODM case 23", "four distinct roots with d2 about 0", op, roots_of( p ), {} );
  }

  std::cout << "\n  literature cases: " << tot.total - tot.failed << " / " << tot.total
            << " within the backward error tolerance\n";
  return summary.finish();
}
