#include "PolynomialRoots.hh"
#include "TestReporter.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace {

  using PolynomialRoots::MAXDEGREE;
  using PolynomialRoots::Quadratic;
  using PolynomialRoots::Cubic;
  using PolynomialRoots::Quartic;
  using PolynomialRoots::real_type;

  TestReporter::Summary *
  summary_ptr{nullptr};

  bool
  close_enough( real_type a, real_type b ) {
    real_type const scale = std::max(
      real_type(1),
      std::max( std::abs(a), std::abs(b) )
    );
    return std::abs(a-b) <= 100*std::numeric_limits<real_type>::epsilon()*scale;
  }

  void
  check_value(
    char const * label,
    real_type      value,
    real_type      expected
  ) {
    if ( close_enough(value,expected) ) {
      summary_ptr->pass(label);
      return;
    }
    std::cerr
      << "    value    = " << value << '\n'
      << "    expected = " << expected << '\n';
    summary_ptr->fail(label);
  }

}

int
main() {
  TestReporter::Summary summary(
    std::cout,
    "Polynomial/derivative regression suite"
  );
  summary_ptr = &summary;

  {
    summary.case_header(1,"quadratic evalPolyDPoly","regression for |x| <= 1 and |x| > 1");
    Quadratic q(1,-3,2);
    real_type p{0};
    real_type dp{0};

    q.eval(0.5,p,dp);
    check_value("quadratic p(|x|<=1)",p,0.75);
    check_value("quadratic dp(|x|<=1)",dp,-2.0);

    q.eval(2.0,p,dp);
    check_value("quadratic p(|x|>1)",p,0.0);
    check_value("quadratic dp(|x|>1)",dp,1.0);
  }

  {
    summary.case_header(2,"cubic evalPolyDPoly","regression for |x| <= 1 and |x| > 1");
    Cubic c(1,-6,11,-6);
    real_type p{0};
    real_type dp{0};

    c.eval(0.5,p,dp);
    check_value("cubic p(|x|<=1)",p,-1.875);
    check_value("cubic dp(|x|<=1)",dp,5.75);

    c.eval(3.0,p,dp);
    check_value("cubic p(|x|>1)",p,0.0);
    check_value("cubic dp(|x|>1)",dp,2.0);
  }

  {
    summary.case_header(3,"quartic evalPolyDPoly","regression for |x| <= 1 and |x| > 1");
    Quartic q(1,0,-5,0,4);
    real_type p{0};
    real_type dp{0};

    q.eval(0.25,p,dp);
    check_value("quartic p(|x|<=1)",p,3.69140625);
    check_value("quartic dp(|x|<=1)",dp,-2.4375);

    q.eval(2.5,p,dp);
    check_value("quartic p(|x|>1)",p,11.8125);
    check_value("quartic dp(|x|>1)",dp,37.5);
  }

  {
    summary.case_header(4,"MAXDEGREE guard","Jenkins-Traub entry-point");
    std::vector<real_type> coeffs(MAXDEGREE+2,0);
    std::vector<real_type> zr(MAXDEGREE+1,0);
    std::vector<real_type> zi(MAXDEGREE+1,0);
    coeffs.front() = 1;

    if ( PolynomialRoots::roots(coeffs.data(),MAXDEGREE+1,zr.data(),zi.data()) == -3 ) {
      summary.pass("roots(MAXDEGREE overflow)");
    } else {
      summary.fail("roots(MAXDEGREE overflow)");
    }
  }

  return summary.finish();
}
