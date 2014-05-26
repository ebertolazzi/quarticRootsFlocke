#include "rpoly.hh"
#include <iostream>

using namespace std ;

int
main() {

  Rpoly::Rpoly rpoly ;

  double op[] = {1,-55,1320,-18150,157773,-902055,3416930,-8409500,
                 12753576,-10628640,3628800} ;
  
  double zeror[10], zeroi[10] ;
  int    info[10] ;

  int ok = rpoly.eval( op, 10, zeror, zeroi, info ) ;
  cout << " ok = " << ok << '\n' ;
  for ( int i = 0 ; i < 10 ; ++i )
    cout << zeror[i] << " + I* " << zeroi[i] << '\n' ;

  double op1[] = { 8, -8, 16,-16, 8,-8};

  ok = rpoly.eval( op1, 5, zeror, zeroi, info ) ;
  cout << " ok = " << ok << '\n' ;
  for ( int i = 0 ; i < 5 ; ++i )
    cout << zeror[i] << " + I* " << zeroi[i] << '\n' ;

  double op2[] = { 8e9, -8e9, 16e9,-16e9, 8e9,-8e9};

  ok = rpoly.eval( op2, 5, zeror, zeroi, info ) ;
  cout << " ok = " << ok << '\n' ;
  for ( int i = 0 ; i < 5 ; ++i )
    cout << zeror[i] << " + I* " << zeroi[i] << '\n' ;

  return 0 ;
}

