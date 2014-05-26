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

namespace Rpoly {

  typedef double valueType ;

  class Rpoly {

    void
    quad( valueType a,
          valueType b1,
          valueType c,
          valueType & sr,
          valueType & si,
          valueType & lr,
          valueType & li) ;

    void
    quadsd( int         n,
            valueType & u,
            valueType & v,
            valueType   p[],
            valueType   q[],
            valueType & a,
            valueType & b) ;

    void
    fxshfr( int l2, int & nz ) ;

    void
    quadit( valueType & uu,
            valueType & vv,
            int       & nz) ;

    void
    realit( valueType sss, int  & nz, int & iflag) ;

    void calcsc(int & type) ;
    void nextk(int type) ;

    void
    newest( int type, valueType & uu, valueType & vv ) ;

    valueType const base ;
    valueType const eta ;
    valueType const infin ;
    valueType const smalno ;

    valueType *p,*qp,*k,*qk,*svk;
    valueType sr,si,u,v ;
    valueType a1,a2,a3,a6,a7 ;
    valueType a,b,c,d,e,f,g,h ;
    valueType szr,szi,lzr,lzi ;
    valueType are,mre;

    int n, itercnt;

  public:
  
    Rpoly() ;

    int
    eval( valueType const op[],
          int       const degree,
          valueType       zeror[],
          valueType       zeroi[],
          int             info[] ) ;
  };

}

#endif
