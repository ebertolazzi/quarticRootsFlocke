// check_Jenkins-Traub.cc - Test Jenkins-Traub solver (roots()) for degrees 4..10
// Compares against known exact roots where possible and checks residuals.

#include "PolynomialRoots.hh"
#include "TestReporter.hh"
#include "TestHelpers.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using PolynomialRoots::real_type;
using std::abs;
using std::array;
using std::complex;
using std::vector;

// ---------------------------------------------------------------------------
// Helper: build monic polynomial coefficients from roots.
// Returns vector [a_n, a_{n-1}, ..., a_0] with a_n = 1 (leading coefficient first).
// roots: vector of complex roots.
// If roots is empty, returns constant polynomial 1.
// ---------------------------------------------------------------------------
vector<real_type> build_polynomial( const vector<complex<real_type>> & roots )
{
  std::size_t n = roots.size();
  // coeffs[0] = constant term, coeffs[n] = leading coefficient during construction
  vector<complex<real_type>> coeffs( n + 1, 0.0 );
  coeffs[0] = 1.0;  // start with constant polynomial 1

  for ( std::size_t i = 0; i < n; ++i )
  {
    // multiply current polynomial (coeffs[0..i]) by (x - roots[i])
    vector<complex<real_type>> new_coeffs( n + 1, 0.0 );
    for ( std::size_t j = 0; j <= i; ++j )
    {
      new_coeffs[j] += coeffs[j] * ( -roots[i] );  // constant term
      new_coeffs[j + 1] += coeffs[j];              // x term
    }
    coeffs = std::move( new_coeffs );
  }

  // Convert to real_type and reverse order so that real_coeffs[0] = leading coefficient
  // This matches the expected input for PolynomialRoots::roots()
  vector<real_type> real_coeffs( n + 1 );
  for ( std::size_t i = 0; i <= n; ++i )
  {
    real_coeffs[i] = coeffs[n - i].real();  // reverse order: leading first
  }
  return real_coeffs;
}

// ---------------------------------------------------------------------------
// Compute the maximum RELATIVE backward error over the computed roots:
//
//            |p(r)| / sum_k |a_k| |r|^(n-k)
//
// The bare residual |p(r)| grows with |r| and with the size of the
// coefficients, so an absolute threshold is unreachable for Wilkinson-like
// polynomials even when every root is correct to the last bit.  The quantity
// below is scale free and is the standard accuracy indicator for root finders.
//
// coeffs: [a_n, a_{n-1}, ..., a_0]
// roots:  vectors of real and imaginary parts (same length)
// Returns -1.0 if input vectors have inconsistent sizes.
// ---------------------------------------------------------------------------
double max_backward_error(
  const vector<real_type> & coeffs,
  const vector<real_type> & zeror,
  const vector<real_type> & zeroi )
{
  if ( zeror.size() != zeroi.size() ) return -1.0;

  std::size_t n = coeffs.size() - 1;
  if ( zeror.size() != n ) return -1.0;

  double maxRes = 0.0;
  for ( std::size_t i = 0; i < n; ++i )
  {
    complex<real_type> const r( zeror[i], zeroi[i] );
    complex<real_type>       val( coeffs[0], 0.0 );
    double const             ar  = std::abs( r );
    double                   den = 0.0;
    double                   pw  = 1.0;
    for ( std::size_t k = 1; k <= n; ++k ) val = val * r + coeffs[k];  // Horner
    for ( std::size_t k = n + 1; k-- > 0; )
    {
      den += std::abs( coeffs[k] ) * pw;
      pw *= ar;
      if ( !std::isfinite( pw ) ) break;
    }
    if ( den > 0 ) maxRes = std::max( maxRes, std::abs( val ) / den );
  }
  return maxRes;
}

// ---------------------------------------------------------------------------
// Match computed roots to exact roots (greedy assignment by nearest distance)
// Returns the maximum Euclidean distance after matching.
// If vectors are empty or sizes differ, returns -1.0.
// ---------------------------------------------------------------------------
double matching_error( const vector<complex<real_type>> & exact, const vector<complex<real_type>> & computed )
{
  if ( exact.empty() || computed.empty() || exact.size() != computed.size() ) return -1.0;

  std::size_t  n = exact.size();
  vector<bool> used( n, false );
  double       maxErr = 0.0;

  for ( std::size_t i = 0; i < n; ++i )
  {
    double      bestDist = std::numeric_limits<double>::max();
    std::size_t bestIdx  = 0;
    for ( std::size_t j = 0; j < n; ++j )
    {
      if ( !used[j] )
      {
        double d = std::abs( exact[i] - computed[j] );
        if ( d < bestDist )
        {
          bestDist = d;
          bestIdx  = j;
        }
      }
    }
    used[bestIdx] = true;
    maxErr        = std::max( maxErr, bestDist );
  }
  return maxErr;
}

// ---------------------------------------------------------------------------
// Test case structures
// ---------------------------------------------------------------------------
struct QuadraticCase
{
  const char * label;
  const char * source;
  real_type    a, b, c;
};

struct CubicCase
{
  const char * label;
  const char * source;
  real_type    a, b, c, d;
};

struct QuarticCase
{
  const char * label;
  const char * source;
  real_type    a, b, c, d, e;
};

struct GeneralCase
{
  std::string                label;        // FIXED: now owns its string, no dangling pointer
  const char *               source;       // static literal
  vector<real_type>          coeffs;       // [a_n, ..., a_0] (leading first)
  vector<complex<real_type>> exact_roots;  // known roots
};

// ---------------------------------------------------------------------------
// Predefined test cases (copied from check_4_classical_cases.cc)
// ---------------------------------------------------------------------------
QuadraticCase quadratic_cases[] = { { "(x-1)^2", "classical multiple-root family", 1, -2, 1 },
                                    { "x^2+1", "classical conjugate-pair family", 1, 0, 1 },
                                    { "W2(x)=(x-1)(x-2)", "Wilkinson 1963", 1, -3, 2 } };

CubicCase cubic_cases[] = { { "(x-1)^3", "classical multiple-root family", 1, -3, 3, -1 },
                            { "(x-1)^2(x+2)", "classical multiple-root family", 1, 0, -3, 2 },
                            { "W3(x)=(x-1)(x-2)(x-3)", "Wilkinson 1963", 1, -6, 11, -6 },
                            { "x^3-3x+1", "classical casus irreducibilis", 1, 0, -3, 1 } };

QuarticCase quartic_cases[] = { { "(x-1)^4", "Chavez-Pichardo et al. 2022, M4 family", 1, -4, 6, -4, 1 },
                                { "(x-1)^3(x+2)", "Chavez-Pichardo et al. 2022, M3 family", 1, -1, -3, 5, -2 },
                                { "(x-1)^2(x-3)^2", "Chavez-Pichardo et al. 2022, DM2 family", 1, -8, 22, -24, 9 },
                                { "W4(x)=(x-1)(x-2)(x-3)(x-4)", "Wilkinson 1963", 1, -10, 35, -50, 24 },
                                { "(x^2-1)(x^2-4)", "classical biquadratic family", 1, 0, -5, 0, 4 },
                                { "(x^2+1)^2", "classical double complex pair", 1, 0, 2, 0, 1 } };

// ---------------------------------------------------------------------------
// Build exact roots for quartic cases (for comparison)
// ---------------------------------------------------------------------------
vector<complex<real_type>> exact_roots_quartic( const QuarticCase & tc )
{
  // These are known from the labels; we hard-code them for simplicity.
  // (x-1)^4 -> 1 (multiplicity 4)
  if ( std::string( tc.label ) == "(x-1)^4" ) { return { { 1, 0 }, { 1, 0 }, { 1, 0 }, { 1, 0 } }; }
  // (x-1)^3(x+2) -> 1 (x3), -2
  if ( std::string( tc.label ) == "(x-1)^3(x+2)" ) { return { { 1, 0 }, { 1, 0 }, { 1, 0 }, { -2, 0 } }; }
  // (x-1)^2(x-3)^2 -> 1 (x2), 3 (x2)
  if ( std::string( tc.label ) == "(x-1)^2(x-3)^2" ) { return { { 1, 0 }, { 1, 0 }, { 3, 0 }, { 3, 0 } }; }
  // W4 -> 1,2,3,4
  if ( std::string( tc.label ) == "W4(x)=(x-1)(x-2)(x-3)(x-4)" ) { return { { 1, 0 }, { 2, 0 }, { 3, 0 }, { 4, 0 } }; }
  // (x^2-1)(x^2-4) -> -1,1,-2,2
  if ( std::string( tc.label ) == "(x^2-1)(x^2-4)" ) { return { { -1, 0 }, { 1, 0 }, { -2, 0 }, { 2, 0 } }; }
  // (x^2+1)^2 -> i, -i, i, -i  (order not important)
  if ( std::string( tc.label ) == "(x^2+1)^2" ) { return { { 0, 1 }, { 0, -1 }, { 0, 1 }, { 0, -1 } }; }
  return {};
}

// ---------------------------------------------------------------------------
// Generate general cases for degrees 5..10
// ---------------------------------------------------------------------------
vector<GeneralCase> generate_general_cases()
{
  vector<GeneralCase> cases;

  // Wilkinson-like: roots 1..n
  for ( int deg = 5; deg <= 10; ++deg )
  {
    vector<complex<real_type>> roots;
    for ( int i = 1; i <= deg; ++i ) roots.emplace_back( static_cast<real_type>( i ), 0.0 );
    vector<real_type> coeffs = build_polynomial( roots );
    // FIXED: use std::string for label, so we can safely use .c_str() later
    cases.push_back(
      { "W" + std::to_string( deg ),  // now a std::string, not a temporary
        "Wilkinson-like roots 1..n",
        coeffs,
        roots } );
  }

  // Alternating roots: 1, -2, 3, -4, ...
  for ( int deg = 5; deg <= 10; ++deg )
  {
    vector<complex<real_type>> roots;
    for ( int i = 1; i <= deg; ++i )
    {
      real_type r = ( i % 2 == 1 ) ? static_cast<real_type>( i ) : -static_cast<real_type>( i );
      roots.emplace_back( r, 0.0 );
    }
    vector<real_type> coeffs = build_polynomial( roots );
    cases.push_back( { "Alt" + std::to_string( deg ), "Alternating sign roots", coeffs, roots } );
  }

  // Complex conjugate pairs, limited to degree 8 to avoid numerical issues,
  // and with smaller moduli (real part = i/2, imag = 0.5)
  for ( int deg = 5; deg <= 8; ++deg )
  {
    vector<complex<real_type>> roots;
    for ( int i = 1; i <= deg / 2; ++i )
    {
      roots.emplace_back( static_cast<real_type>( i ) * 0.5, 0.5 );
      roots.emplace_back( static_cast<real_type>( i ) * 0.5, -0.5 );
    }
    if ( deg % 2 == 1 ) { roots.emplace_back( 0.0, 0.0 ); }
    vector<real_type> coeffs = build_polynomial( roots );
    cases.push_back( { "Complex" + std::to_string( deg ), "Complex conjugate pairs (scaled)", coeffs, roots } );
  }

  return cases;
}

// ---------------------------------------------------------------------------
// Run a single test case using Jenkins-Traub roots()
// ---------------------------------------------------------------------------
std::pair<bool, double> run_JT( const vector<real_type> & coeffs, vector<complex<real_type>> & computed_roots )
{
  int n = static_cast<int>( coeffs.size() ) - 1;
  if ( n <= 0 )
  {
    computed_roots.clear();
    return { false, std::numeric_limits<double>::infinity() };
  }

  vector<real_type> zeror( static_cast<std::size_t>( n ) ), zeroi( static_cast<std::size_t>( n ) );
  int               status = PolynomialRoots::roots( coeffs.data(), n, zeror.data(), zeroi.data() );
  if ( status != 0 )
  {
    computed_roots.clear();
    return { false, std::numeric_limits<double>::infinity() };
  }

  computed_roots.resize( static_cast<std::size_t>( n ) );
  for ( int i = 0; i < n; ++i )
  {
    computed_roots[static_cast<std::size_t>( i )] = complex<real_type>(
      zeror[static_cast<std::size_t>( i )],
      zeroi[static_cast<std::size_t>( i )] );
  }

  double res = max_backward_error( coeffs, zeror, zeroi );
  // If max_backward_error returns a negative value, something went wrong.
  if ( res < 0.0 ) { return { false, std::numeric_limits<double>::infinity() }; }

  // Scale free threshold: a backward stable solver stays within a small
  // multiple of the unit roundoff whatever the degree and the conditioning.
  bool ok = ( res < 1e-13 );
  return { ok, res };
}

// ---------------------------------------------------------------------------
// Report a test case
// ---------------------------------------------------------------------------
template <typename CaseType> void report_case(
  TestReporter::Summary &            summary,
  int &                              idx,
  const CaseType &                   tc,
  const vector<complex<real_type>> & computed,
  double                             res,
  bool                               ok,
  const vector<complex<real_type>> & exact = {} )
{
  // FIXED: use c_str() only when tc.label is a std::string (or a const char*)
  if constexpr ( std::is_same_v<CaseType, GeneralCase> ) { summary.case_header( idx++, tc.label.c_str(), tc.source ); }
  else
  {
    summary.case_header( idx++, tc.label, tc.source );
  }

  summary.comparison_table_header();
  summary.comparison_table_row( "Jenkins-Traub", ok, res );

  // Compute matching error and use it to possibly override the 'ok' status.
  bool   final_ok = ok;
  double err      = -1.0;
  if ( !exact.empty() && !computed.empty() && ok )
  {
    err = matching_error( exact, computed );
    if ( err >= 0.0 )
    {
      summary.comparison_table_row( "matching error", ok, err );
      // The forward error cannot be smaller than (unit roundoff) x (condition
      // number of the root), which for Wilkinson-like families grows very
      // quickly with the degree: scale the tolerance with the magnitude of the
      // exact roots instead of using a fixed absolute bound.
      double scale = 1.0;
      for ( auto const & z : exact ) scale = std::max( scale, std::abs( z ) );
      double const matching_tolerance = 1e-9 * scale * static_cast<double>( exact.size() * exact.size() );
      if ( err > matching_tolerance ) { final_ok = false; }
    }
    else
    {
      summary.comparison_table_row( "matching error", false, err );
      final_ok = false;
    }
  }

  summary.comparison_table_footer();
  if ( final_ok ) { summary.pass( "Backward error and matching within tolerance" ); }
  else
  {
    summary.fail( "Backward error too high or matching error too large" );
  }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main()
{
  std::cout.precision( 18 );
  TestReporter::Summary summary( std::cout, "Jenkins-Traub solver test (degrees 4..10)" );

  int idx = 1;

  // ----- Test quartic cases (degree 4) -----
  for ( const auto & tc : quartic_cases )
  {
    // Build coefficient vector
    vector<real_type>          coeffs = { tc.a, tc.b, tc.c, tc.d, tc.e };
    vector<complex<real_type>> computed;
    auto [ok, res] = run_JT( coeffs, computed );
    auto exact     = exact_roots_quartic( tc );
    report_case( summary, idx, tc, computed, res, ok, exact );
  }

  // ----- Test cubic cases (degree 3) -----
  for ( const auto & tc : cubic_cases )
  {
    vector<real_type>          coeffs = { tc.a, tc.b, tc.c, tc.d };
    vector<complex<real_type>> computed;
    auto [ok, res] = run_JT( coeffs, computed );
    report_case( summary, idx, tc, computed, res, ok );
  }

  // ----- Test quadratic cases (degree 2) -----
  for ( const auto & tc : quadratic_cases )
  {
    vector<real_type>          coeffs = { tc.a, tc.b, tc.c };
    vector<complex<real_type>> computed;
    auto [ok, res] = run_JT( coeffs, computed );
    report_case( summary, idx, tc, computed, res, ok );
  }

  // ----- Test general cases (degrees 5..10) -----
  auto general_cases = generate_general_cases();
  for ( const auto & tc : general_cases )
  {
    vector<complex<real_type>> computed;
    auto [ok, res] = run_JT( tc.coeffs, computed );
    report_case( summary, idx, tc, computed, res, ok, tc.exact_roots );
  }

  std::cout << summary.finish() << '\n';

  return 0;
}
