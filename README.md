Quartic Roots Solver
====================

[![Build Status](https://travis-ci.org/ebertolazzi/quarticRootsFlocke.svg?branch=master)](https://travis-ci.org/ebertolazzi/quarticRootsFlocke)

Polynomial Roots Solver
-----------------------

Port to C++ of Jenkins-Traub real
polynomial root finder and the algorithm
of Norbert Flocke for polynomial roots
up to degree 4.

[Doxygen documentation: http://ebertolazzi.github.io/quarticRootsFlocke/](http://ebertolazzi.github.io/quarticRootsFlocke/)

[Repository: https://github.com/ebertolazzi/quarticRootsFlocke](https://github.com/ebertolazzi/quarticRootsFlocke)

Usage
-----

```cpp
  #include "PolynomialRoots.hh"
  ....

  double coeffs[] = { 8, -8, 16,-16, 8,-8 }; // polynomial coeffs

  double zeror[5], zeroi[5];
  int    info[5];
  int    degree = 5;

  int ok = PolynomialRoots::roots( coeffs, degree, zeror, zeroi ); // ok < 0 failed
  cout << " ok = " << ok << '\n' ;
  for ( int i = 0 ; i < degree ; ++i )
    cout << zeror[i] << " + I* " << zeroi[i] << '\n';
```

To solve quadratic, cubic or quartic use specialized classes

```cpp
  Quadratic qsolve(a,b,c);
  qsolve.info(cout);

  Cubic csolve(a,b,c,d);
  csolve.info(cout);

  Quartic q4(a,b,c,d,e);
  q4.info(cout);
```

look at the class definition to see how to access the computed roots.

Online Documentation
--------------------

Available at: [http://ebertolazzi.github.io/quarticRootsFlocke](http://ebertolazzi.github.io/quarticRootsFlocke)

References
----------

- **N.Flocke**  
  Algorithm 954: An Accurate and Efficient Cubic and Quartic  
  Equation Solver for Physical Applications  
  ACM TOMS, vol 41, n.4, 2015

- **M.A. Jenkins and J.F.Traub**  
  A Three-Stage Algorithm for Real Polynomials
  Using Quadratic Iteration    
  SIAM Journal on Numerical Analysis  
  Vol. 7, No. 4 (Dec., 1970), pp. 545-566

Author
------

    Enrico Bertolazzi  
    Dipartimento di Ingegneria Industriale  
    Universita` degli Studi di Trento  
    email: enrico.bertolazzi@unitn.it
