Rpoly
=====

Port to C++ of Jenkins-Traub real polynomial root finder

Usage
-----

~~~~~~~~~~~~~~~~~~~~~
  #include "Rpoly.hh"
  ....

  Rpoly::Rpoly rpoly ;
  double coeffs[] = { 8, -8, 16,-16, 8,-8 } ; // polynomial coeffs
  
  double zeror[5], zeroi[5] ;
  int    info[5] ;
  int    degree = 5 ;

  int ok = rpoly.eval( coeffs, degree, zeror, zeroi, info ) ; // ok = number or zeros or -1
  cout << " ok = " << ok << '\n' ;
  for ( int i = 0 ; i < degree ; ++i )
    cout << zeror[i] << " + I* " << zeroi[i] << '\n' ;
~~~~~~~~~~~~~~~~~~~~~

* * *

Enrico Bertolazzi<br>
Dipartimento di Ingegneria Industriale<br>
Universita` degli Studi di Trento<br>
email: enrico.bertolazzi@unitn.it