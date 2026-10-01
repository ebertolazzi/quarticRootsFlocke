#ifndef TEST_REPORTER_HH
#define TEST_REPORTER_HH

#include "termcolor.hh"
#include <algorithm>
#include <format>  // C++20
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>

namespace TestReporter
{

  //
  // Column widths of the comparison table, in characters of CONTENT; each cell
  // is printed with one blank on either side, so the rule segment above it is
  // two characters wider.  Keeping the three widths in one place is what stops
  // the header, the rows and the footer from drifting apart.
  //
  inline constexpr std::size_t MODE_WIDTH{ 14 };   // longest label in use is "matching error"
  inline constexpr std::size_t RESULT_WIDTH{ 6 };  // "pass" / "fail"
  inline constexpr std::size_t NOTE_WIDTH{ 40 };   // residual plus its note

  namespace detail
  {
    //! Repeat a (possibly multi byte) string \p n times.
    inline std::string repeat( std::string_view piece, std::size_t n )
    {
      std::string s;
      s.reserve( piece.size() * n );
      for ( std::size_t i{ 0 }; i < n; ++i ) s.append( piece );
      return s;
    }

    //!
    //! Fit \p text into exactly \p width characters: pad on the right when it is
    //! shorter, cut it and mark the cut with a trailing `~` when it is longer.
    //! Without the cut a long label - "Jenkins-Traub" in a six character column -
    //! pushes the right border of the table out of line.
    //!
    inline std::string fit( std::string_view text, std::size_t width )
    {
      if ( text.size() == width ) return std::string{ text };
      if ( text.size() < width ) return std::string{ text } + std::string( width - text.size(), ' ' );
      if ( width == 0 ) return {};
      return std::string{ text.substr( 0, width - 1 ) } + '~';
    }

    //! One horizontal rule, built from the column widths.
    inline std::string rule( std::string_view left, std::string_view mid, std::string_view right )
    {
      return std::string{ left } + repeat( "─", MODE_WIDTH + 2 ) + std::string{ mid } +
             repeat( "─", RESULT_WIDTH + 2 ) + std::string{ mid } + repeat( "─", NOTE_WIDTH + 2 ) +
             std::string{ right };
    }
  }  // namespace detail

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
      // the banner follows the title instead of being a fixed length bar that a
      // longer suite name would stick out of
      std::size_t const bar{ std::max<std::size_t>( 62, m_suite.size() + 2 ) };
      m_out << termcolor::colorize;
      m_out << '\n'
            << termcolor::bold << termcolor::cyan << "┏" << detail::repeat( "━", bar ) << '\n'
            << "┃ " << m_suite << '\n'
            << "┗" << detail::repeat( "━", bar ) << '\n'
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
      // the reset goes before the newline, otherwise the escape sequence lands
      // at the start of the next line and shows up in captured output
      m_out << termcolor::bold << termcolor::magenta << "  " << title << termcolor::reset << '\n'
            << "  " << detail::rule( "┌", "┬", "┐" ) << '\n'
            << "  │ " << detail::fit( "mode", MODE_WIDTH ) << " │ " << detail::fit( "result", RESULT_WIDTH ) << " │ "
            << detail::fit( "residual & note", NOTE_WIDTH ) << " │\n"
            << "  " << detail::rule( "├", "┼", "┤" ) << '\n';
    }

    //! One row: label, outcome and the residual with its note.
    void comparison_table_row( std::string_view mode, bool ok, double residual, std::string_view note = {} )
    {
      std::string_view const status{ ok ? "pass" : "fail" };
      if ( note.empty() ) note = ok ? "check() passed" : "check() failed";

      std::string text{ std::format( "{:.4e}", residual ) };
      text += std::format( " ({})", note );

      m_out << "  │ " << detail::fit( mode, MODE_WIDTH ) << " │ ";

      if ( ok )
        m_out << termcolor::green;
      else
        m_out << termcolor::red << termcolor::bold;
      m_out << detail::fit( status, RESULT_WIDTH ) << termcolor::reset << " │ ";

      // colour by the size of the residual
      if ( residual < 1e-12 )
        m_out << termcolor::green;
      else if ( residual < 1e-6 )
        m_out << termcolor::yellow;
      else if ( residual < 1.0 )
        m_out << termcolor::magenta;
      else
        m_out << termcolor::red;

      m_out << detail::fit( text, NOTE_WIDTH ) << termcolor::reset << " │\n";
    }

    void comparison_table_footer() { m_out << "  " << detail::rule( "└", "┴", "┘" ) << '\n'; }

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
                << termcolor::bold << termcolor::yellow << "⚠ Warnings: " << m_warned << " known difficult case"
                << ( m_warned == 1 ? "" : "s" ) << termcolor::reset << "\n\n";
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
