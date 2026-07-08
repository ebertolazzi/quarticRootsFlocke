#ifndef TEST_REPORTER_HH
#define TEST_REPORTER_HH

#include "termcolor.hh"

#include <iomanip>
#include <iostream>
#include <string_view>

namespace TestReporter {

  class Summary {
    std::ostream & m_out;
    std::string_view m_suite;
    int m_total{0};
    int m_failed{0};
    int m_warned{0};

  public:

    explicit
    Summary(
      std::ostream     & out,
      std::string_view   suite
    )
    : m_out(out)
    , m_suite(suite)
    {
      m_out << termcolor::colorize;
      m_out
        << '\n'
        << termcolor::bold << termcolor::cyan
        << "┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        << "┃ " << m_suite << '\n'
        << "┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        << termcolor::reset;
    }

    void
    case_header(
      int               index,
      std::string_view  label,
      std::string_view  source = {}
    ) {
      m_out
        << '\n'
        << termcolor::bold << termcolor::blue
        << "▶ Case " << std::setw(2) << index << "  " << label
        << termcolor::reset << '\n';
      if ( !source.empty() ) {
        m_out
          << termcolor::dark << "  source: " << source << termcolor::reset
          << '\n';
      }
    }

    void
    note( std::string_view msg ) {
      m_out
        << termcolor::yellow
        << "  • " << msg
        << termcolor::reset
        << '\n';
    }

    void
    pass( std::string_view msg = "check() passed" ) {
      ++m_total;
      m_out
        << termcolor::green
        << "  ✓ " << msg
        << termcolor::reset
        << '\n';
    }

    void
    fail( std::string_view msg = "check() failed" ) {
      ++m_total;
      ++m_failed;
      m_out
        << termcolor::red << termcolor::bold
        << "  ✗ " << msg
        << termcolor::reset
        << '\n';
    }

    void
    warn( std::string_view msg = "known difficult case" ) {
      ++m_total;
      ++m_warned;
      m_out
        << termcolor::yellow << termcolor::bold
        << "  ⚠ " << msg
        << termcolor::reset
        << '\n';
    }

    int
    finish() {
      m_out << '\n';
      if ( m_failed == 0 ) {
        if ( m_warned == 0 ) {
          m_out
            << termcolor::bold << termcolor::green
            << "✓ Summary: all " << m_total << " cases passed"
            << termcolor::reset
            << "\n\n";
        } else {
          m_out
            << termcolor::bold << termcolor::green
            << "✓ Summary: no regressions in " << m_total << " cases"
            << termcolor::reset
            << '\n'
            << termcolor::bold << termcolor::yellow
            << "⚠ Warnings: " << m_warned << " known difficult cases"
            << termcolor::reset
            << "\n\n";
        }
        return 0;
      }

      m_out
        << termcolor::bold << termcolor::red
        << "✗ Summary: " << m_failed << " / " << m_total << " cases failed"
        << termcolor::reset
        << "\n\n";
      return 1;
    }
  };

}

#endif
