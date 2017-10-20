/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2014                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                | 
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef RPOLY_HH
#define RPOLY_HH

#include <utility>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <iostream>

/*
..
.. N. FLOCKE
.. Algorithm 954: An Accurate and Efficient Cubic and Quartic Equation Solver
.. for Physical Applications
.. ACM TOMS, Vol. 41, No. 4, 2015.
.. DOI: http://dx.doi.org/10.1145/2699468
..
*/

namespace PolynomialRoots {

  typedef double valueType ;
  typedef int    indexType ;
  typedef std::complex<valueType> complexType ;

  // a * x^2 + b * x + c ;
  class Quadratic {
    valueType a, b, c ;
    valueType r0, r1 ;
    indexType nrts ;
    bool      cplx ;
    bool      dblx ;

    void findRoots();

  public:

    Quadratic() : a(0), b(0), c(0), nrts(0), cplx(false), dblx(false) {}
    Quadratic( valueType _a, valueType _b, valueType _c )
    : a(_a), b(_b), c(_c), nrts(0), cplx(false), dblx(false)  {
      findRoots();
    }

    void
    setup( valueType _a, valueType _b, valueType _c ) {
      a = _a ; b = _b ; c = _c ;
      findRoots();
    }

    bool numRoots()     const { return nrts ; }
    bool complexRoots() const { return cplx ; }
    bool dobleRoot()    const { return dblx ; }

    valueType real_root0() const { return r0 ; }
    valueType real_root1() const { return r1 ; }

    complexType
    root0() const
    { return cplx ? complexType(r0,r1) : complexType(r0,0)  ; }

    complexType
    root1() const
    { return cplx ? complexType(r0,-r1) : complexType(r1,0)  ; }

    void
    getRoot0( valueType & re, valueType & im ) const {
      if ( cplx ) {
        re = r0 ; im = r1 ;
      } else {
        re = r0 ; im = 0 ;
      }
    }

    void
    getRoot1( valueType & re, valueType & im ) const {
      if ( cplx ) {
        re = r0 ; im = -r1 ;
      } else {
        re = r1 ; im = 0 ;
      }
    }

    void
    getRoot0( complexType & r ) const {
      if ( cplx ) {
        r = complexType(r0,r1) ;
      } else {
        r = complexType(r0,0) ;
      }
    }

    void
    getRoot1( complexType & r ) const {
      if ( cplx ) {
        r = complexType(r0,-r1) ;
      } else {
        r = complexType(r1,0) ;
      }
    }

    valueType
    eval( valueType x ) const
    { return (a*x+b)*x+c ; }

    complexType
    eval( complexType x ) const
    { return (a*x+b)*x+c ; }

    void
    info( std::ostream & s ) const ;

    bool
    check( std::ostream & s ) const ;

  } ;
/*
  void
  solveQuadratic( valueType   a,
                  valueType   b,
                  valueType   c,
                  valueType & r1,
                  valueType & r2,
                  indexType & nr,
                  indexType & nc ) ;
*/
  indexType
  solveCubic( valueType   A,
              valueType   B,
              valueType   C,
              valueType   D,
              valueType & r1,
              valueType & r2,
              valueType & r3,
              indexType & nr,
              indexType & nc ) ;

  indexType
  solveQuartic( valueType   A,
                valueType   B,
                valueType   C,
                valueType   D,
                valueType   E,
                valueType & r1,
                valueType & r2,
                valueType & r3,
                valueType & r4,
                indexType & nr,
                indexType & nc ) ;

  valueType
  evalPoly( valueType const op[],
            indexType       Degree,
            valueType       x,
            bool            reverse ) ;

  std::complex<valueType>
  evalPolyC( valueType const         op[],
             indexType               Degree,
             std::complex<valueType> x,
             bool                    reverse ) ;

  valueType
  CompHorner( valueType const p[],
              indexType       Degree,
              valueType       x,
              bool            reverse ) ;

  int
  roots( valueType const op[],
         indexType       Degree,
         valueType       zeror[],
         valueType       zeroi[] ) ;
}

extern "C" void quartic_solver_( double op[], int & degree, double zeror[], double zeroi[] ) ;

#endif
