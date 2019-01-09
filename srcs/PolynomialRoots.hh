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

#ifndef POLYNOMIAL_ROOTS_HH
#define POLYNOMIAL_ROOTS_HH

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

  typedef double valueType;
  typedef int    indexType;
  typedef std::complex<valueType> complexType;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // a * x^2 + b * x + c;
  class Quadratic {
    valueType a, b, c;
    valueType r0, r1;
    indexType nrts;
    bool      cplx;
    bool      dblx;

    void findRoots();

  public:

    Quadratic() : a(0), b(0), c(0), nrts(0), cplx(false), dblx(false) {}
    Quadratic( valueType _a, valueType _b, valueType _c )
    : a(_a), b(_b), c(_c), nrts(0), cplx(false), dblx(false)  {
      findRoots();
    }

    void
    setup( valueType _a, valueType _b, valueType _c ) {
      a = _a; b = _b; c = _c;
      findRoots();
    }

    indexType numRoots()     const { return nrts; }
    bool      complexRoots() const { return cplx; }
    bool      doubleRoot()   const { return dblx; }

    indexType getRealRoots( valueType r[] ) const;
    indexType getPositiveRoots( valueType r[] ) const;
    indexType getNegativeRoots( valueType r[] ) const;

    valueType real_root0() const { return r0; }
    valueType real_root1() const { return r1; }

    complexType
    root0() const
    { return cplx ? complexType(r0,r1) : complexType(r0,0); }

    complexType
    root1() const
    { return cplx ? complexType(r0,-r1) : complexType(r1,0); }

    void
    getRoot0( valueType & re, valueType & im ) const {
      if ( cplx ) { re = r0; im = r1; }
      else        { re = r0; im = 0;  }
    }

    void
    getRoot0( complexType & r ) const {
      if ( cplx ) r = complexType(r0,r1);
      else        r = complexType(r0,0);
    }

    void
    getRoot1( valueType & re, valueType & im ) const {
      if ( cplx ) { re = r0; im = -r1; }
      else        { re = r1; im = 0;   }
    }

    void
    getRoot1( complexType & r ) const {
      if ( cplx ) r = complexType(r0,-r1);
      else        r = complexType(r1,0);
    }

    valueType   eval( valueType x ) const;
    complexType eval( complexType const & x ) const;

    void eval( valueType x, valueType & p, valueType & dp ) const;

    void
    info( std::ostream & s ) const;

    bool
    check( std::ostream & s ) const;

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // A * x^3 + B * x^2 + C * x + D;
  class Cubic {
    valueType A, B, C, D;
    valueType r0, r1, r2;
    indexType nrts, iter;
    bool      cplx; // complex root
    bool      dblx; // double root
    bool      trpx; // triple root

    void findRoots();

  public:

    Cubic() : A(0), B(0), C(0), D(0), nrts(0), iter(0), cplx(false), trpx(false) {}
    Cubic( valueType _a, valueType _b, valueType _c, valueType _d )
    : A(_a), B(_b), C(_c), D(_d), nrts(0), iter(0), cplx(false), trpx(false) {
      findRoots();
    }

    void
    setup( valueType _a, valueType _b, valueType _c, valueType _d ) {
      A = _a; B = _b; C = _c; D = _d;
      findRoots();
    }

    indexType numRoots()     const { return nrts; }
    bool      complexRoots() const { return cplx; }
    bool      doubleRoot()   const { return dblx; }
    bool      tripleRoot()   const { return trpx; }

    indexType getRealRoots( valueType r[] ) const;
    indexType getPositiveRoots( valueType r[] ) const;
    indexType getNegativeRoots( valueType r[] ) const;

    valueType real_root0() const { return r0; }
    valueType real_root1() const { return r1; }
    valueType real_root2() const { return r2; }

    complexType
    root0() const
    { return cplx ? complexType(r0,r1) : complexType(r0,0); }

    complexType
    root1() const
    { return cplx ? complexType(r0,-r1) : complexType(r1,0); }

    complexType
    root2() const
    { return complexType(r2,0); }

    void
    getRoot0( valueType & re, valueType & im ) const {
      if ( cplx ) { re = r0; im = r1; }
      else        { re = r0; im = 0;  }
    }

    void
    getRoot0( complexType & r ) const {
      if ( cplx ) r = complexType(r0,r1);
      else        r = complexType(r0,0);
    }

    void
    getRoot1( valueType & re, valueType & im ) const {
      if ( cplx ) { re = r0; im = -r1; }
      else        { re = r1; im = 0;   }
    }

    void
    getRoot1( complexType & r ) const {
      if ( cplx ) r = complexType(r0,-r1);
      else        r = complexType(r1,0);
    }

    void
    getRoot2( valueType & re, valueType & im ) const
    { re = r2; im = 0; }

    void
    getRoot2( complexType & r ) const
    { r = complexType(r2,0); }

    valueType   eval( valueType x ) const;
    complexType eval( complexType const & x ) const;

    void eval( valueType x, valueType & p, valueType & dp ) const;

    void
    info( std::ostream & s ) const;

    bool
    check( std::ostream & s ) const;

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // A * x^3 + B * x^2 + C * x + D;
  class Quartic {
    valueType A, B, C, D, E;
    valueType r0, r1, r2, r3;
    indexType iter, nreal, ncplx;

    void findRoots();

    inline bool cplx0() const { return ncplx > 0; }
    inline bool cplx1() const { return ncplx > 0; }
    inline bool cplx2() const { return ncplx > 2; }
    inline bool cplx3() const { return ncplx > 2; }

  public:

    Quartic() : A(0), B(0), C(0), D(0), E(0), iter(0), nreal(0), ncplx(0) {}
    Quartic( valueType _a, valueType _b, valueType _c, valueType _d, valueType _e )
    : A(_a), B(_b), C(_c), D(_d), E(_e), iter(0), nreal(0), ncplx(0) {
      findRoots();
    }

    void
    setup( valueType _a, valueType _b, valueType _c, valueType _d, valueType _e ) {
      A = _a; B = _b; C = _c; D = _d; E = _e;
      findRoots();
    }

    indexType numRealRoots()    const { return nreal; }
    indexType numComplexRoots() const { return ncplx; }

    indexType getRealRoots( valueType r[] ) const;
    indexType getPositiveRoots( valueType r[] ) const;
    indexType getNegativeRoots( valueType r[] ) const;

    valueType real_root0() const { return r0; }
    valueType real_root1() const { return r1; }
    valueType real_root2() const { return r2; }
    valueType real_root3() const { return r3; }

    complexType
    root0() const
    { return cplx0() ? complexType(r0,r1) : complexType(r0,0); }

    complexType
    root1() const
    { return cplx1() ? complexType(r0,-r1) : complexType(r1,0); }

    complexType
    root2() const
    { return cplx2() ? complexType(r2,r3) : complexType(r2,0); }

    complexType
    root3() const
    { return cplx3() ? complexType(r2,-r3) : complexType(r3,0); }

    void
    getRoot0( valueType & re, valueType & im ) const {
      if ( cplx0() ) { re = r0; im = r1; }
      else           { re = r0; im = 0;  }
    }

    void
    getRoot0( complexType & r ) const {
      if ( cplx0() ) r = complexType(r0,r1);
      else           r = complexType(r0,0);
    }

    void
    getRoot1( valueType & re, valueType & im ) const {
      if ( cplx1() ) { re = r0; im = -r1; }
      else           { re = r1; im = 0;   }
    }

    void
    getRoot1( complexType & r ) const {
      if ( cplx1() ) r = complexType(r0,-r1);
      else           r = complexType(r1,0);
    }

    void
    getRoot2( valueType & re, valueType & im ) const {
      if ( cplx2() ) { re = r2; im = r3; }
      else           { re = r2; im = 0;  }
    }

    void
    getRoot2( complexType & r ) const {
      if ( cplx2() ) r = complexType(r2,r3);
      else           r = complexType(r2,0);
    }

    void
    getRoot3( valueType & re, valueType & im ) const {
      if ( cplx3() ) { re = r2; im = -r3; }
      else           { re = r3; im = 0;   }
    }

    void
    getRoot3( complexType & r ) const {
      if ( cplx3() ) r = complexType(r2,-r3);
      else           r = complexType(r3,0);
    }

    valueType   eval( valueType x ) const;
    complexType eval( complexType const & x ) const;

    void
    info( std::ostream & s ) const;

    bool
    check( std::ostream & s ) const;

  };

  valueType
  evalPoly( valueType const op[],
            indexType       Degree,
            valueType       x,
            bool            reverse );

  std::complex<valueType>
  evalPolyC( valueType const         op[],
             indexType               Degree,
             std::complex<valueType> x,
             bool                    reverse );

  valueType
  CompHorner( valueType const p[],
              indexType       Degree,
              valueType       x,
              bool            reverse );

  // find roots of a generic polinomial using Jenkins-Traub method
  int
  roots( valueType const op[],
         indexType       Degree,
         valueType       zeror[],
         valueType       zeroi[] );
}

extern "C" void quartic_solver_( double op[], int & degree, double zeror[], double zeroi[] );

#endif
