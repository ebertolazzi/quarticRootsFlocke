#ifndef TEST_REPORTER_HH
#define TEST_REPORTER_HH

#include "termcolor.hh"
#include <format>  // C++20
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>

namespace TestReporter
{

  class Summary
  {
    std::ostream &   m_out;
    std::string_view m_suite;
    int              m_total{ 0 };
    int              m_failed{ 0 };
    int              m_warned{ 0 };

  public:
    explicit Summary( std::ostream & out, std::string_view suite ) : m_out( out ), m_suite( suite )
    {
      m_out << termcolor::colorize;
      m_out << '\n'
            << termcolor::bold << termcolor::cyan << "┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
            << "┃ " << m_suite << '\n'
            << "┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
            << termcolor::reset;
    }

    void case_header( int index, std::string_view label, std::string_view source = {} )
    {
      m_out << '\n'
            << termcolor::bold << termcolor::blue << "▶ Case " << std::setw( 2 ) << index << "  " << label
            << termcolor::reset << '\n';
      if ( !source.empty() ) { m_out << termcolor::dark << "  source: " << source << termcolor::reset << '\n'; }
    }

    void note( std::string_view msg ) { m_out << termcolor::yellow << "  • " << msg << termcolor::reset << '\n'; }

    void comparison_table_header( std::string_view title = "mode comparison" )
    {
      m_out << termcolor::bold << termcolor::magenta << "  " << title << '\n'
            << termcolor::reset << "  ┌────────┬────────┬──────────────────────────────────────────┐\n"
            << "  │ mode   │ result │ residual & note                          │\n"
            << "  ├────────┼────────┼──────────────────────────────────────────┤\n";
    }

    // Nuova versione con residuo numerico e colore tramite termcolor
    void comparison_table_row( std::string_view mode, bool ok, double residual, std::string_view note = {} )
    {
      std::string status{ ok ? "pass" : "fail" };
      if ( note.empty() ) note = ok ? "check() passed" : "check() failed";

      m_out << "  │ " << std::left << std::setw( 6 ) << mode << " │ ";
      if ( ok )
        m_out << termcolor::green;
      else
        m_out << termcolor::red << termcolor::bold;
      m_out << std::setw( 6 ) << status << termcolor::reset;
      m_out << " │ ";

      // Usa std::format per formattare il residuo in notazione scientifica
      std::string residual_str = std::format( "{:.4e}", residual );
      if ( !note.empty() ) residual_str += std::format( " ({})", note );

      // Scegli il colore in base al residuo usando termcolor
      if ( residual < 1e-12 )
        m_out << termcolor::green;
      else if ( residual < 1e-6 )
        m_out << termcolor::yellow;
      else if ( residual < 1.0 )
        m_out << termcolor::magenta;
      else
        m_out << termcolor::red;

      m_out << std::format( "{:<40}", residual_str ) << termcolor::reset;

      m_out << " │\n" << std::right;
    }

    void comparison_table_footer() { m_out << "  └────────┴────────┴──────────────────────────────────────────┘\n"; }

    void pass( std::string_view msg = "check() passed" )
    {
      ++m_total;
      m_out << termcolor::green << "  ✓ " << msg << termcolor::reset << '\n';
    }

    void fail( std::string_view msg = "check() failed" )
    {
      ++m_total;
      ++m_failed;
      m_out << termcolor::red << termcolor::bold << "  ✗ " << msg << termcolor::reset << '\n';
    }

    void warn( std::string_view msg = "known difficult case" )
    {
      ++m_total;
      ++m_warned;
      m_out << termcolor::yellow << termcolor::bold << "  ⚠ " << msg << termcolor::reset << '\n';
    }

    int finish()
    {
      m_out << '\n';
      if ( m_failed == 0 )
      {
        if ( m_warned == 0 )
        {
          m_out << termcolor::bold << termcolor::green << "✓ Summary: all " << m_total << " cases passed"
                << termcolor::reset << "\n\n";
        }
        else
        {
          m_out << termcolor::bold << termcolor::green << "✓ Summary: no regressions in " << m_total << " cases"
                << termcolor::reset << '\n'
                << termcolor::bold << termcolor::yellow << "⚠ Warnings: " << m_warned << " known difficult cases"
                << termcolor::reset << "\n\n";
        }
        return 0;
      }

      m_out << termcolor::bold << termcolor::red << "✗ Summary: " << m_failed << " / " << m_total << " cases failed"
            << termcolor::reset << "\n\n";
      return 1;
    }
  };

}  // namespace TestReporter
#endif
