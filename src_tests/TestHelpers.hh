#ifndef TEST_HELPERS_HH
#define TEST_HELPERS_HH

#include <cstddef>
#include <type_traits>
#include <utility>

#include "PolynomialRoots.hh"

namespace TestHelpers
{

  using PolynomialRoots::abs;
  using std::abs;

  // ---------------------------------------------------------------------------
  // Detection idiom: non tutti i solver hanno 3 o 4 radici
  // (Quadratic si ferma a root1(), Cubic a root2(), Quartic arriva a root3()).
  // ---------------------------------------------------------------------------

  template <typename...> using void_t = void;

  template <typename T, typename = void> struct has_root2 : std::false_type
  {
  };
  template <typename T> struct has_root2<T, void_t<decltype( std::declval<T const &>().root2() )>> : std::true_type
  {
  };

  template <typename T, typename = void> struct has_root3 : std::false_type
  {
  };
  template <typename T> struct has_root3<T, void_t<decltype( std::declval<T const &>().root3() )>> : std::true_type
  {
  };

  // ---------------------------------------------------------------------------
  // Ogni solver espone `value_type` (reale) e `complex_type` (root0()...root3()
  // restituiscono già complex_type, ed eval(complex_type) è sempre presente):
  // non serve altro macchinario di conversione.
  // ---------------------------------------------------------------------------

  template <typename Solver> using value_type_t   = typename Solver::value_type;
  template <typename Solver> using complex_type_t = typename Solver::complex_type;

  inline double to_double( double value )
  { return value; }

  template <typename T> double to_double( T const & value )
  { return static_cast<double>( value ); }

  // ---------------------------------------------------------------------------
  // Valutazione del polinomio e residui
  // ---------------------------------------------------------------------------

  //! Horner generico, utile per valutare un array di coefficienti esterno
  //! (o per solver privi di un proprio eval()).
  template <typename CoeffContainer, typename Complex>
  Complex eval_poly( CoeffContainer const & coeff, Complex const & z )
  {
    using Real = std::remove_cvref_t<decltype( coeff[0] )>;

    if ( coeff.size() == 0 ) return Complex{};

    Complex p{ coeff[0], Real{ 0 } };
    for ( std::size_t k = 1; k < coeff.size(); ++k ) p = p * z + Complex{ coeff[k], Real{ 0 } };
    return p;
  }

  template <typename Solver> double max_residual( Solver const & solver )
  {
    double max_res = 0.0;
    for ( int i = 0; i < solver.num_roots(); ++i )
    {
      auto const   r   = solver.root( i );
      auto const   pv  = solver.eval( r );
      double const res = to_double( abs( pv ) );
      if ( res > max_res ) max_res = res;
    }
    return max_res;
  }

}  // namespace TestHelpers

#endif
