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

  void
  solveQuadratic( valueType   a,
                  valueType   b,
                  valueType   c,
                  valueType & r1,
                  valueType & r2,
                  indexType & nr,
                  indexType & nc ) ;

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

  int
  roots( valueType const op[],
         indexType       Degree,
         valueType       zeror[],
         valueType       zeroi[] ) ;
}

extern "C" void quartic_solver_( double op[], int & degree, double zeror[], double zeroi[] ) ;

#endif
