#include <cstdio>
#include "catch.hpp"
#include <ibex.h>
#include <flpl.h>

using namespace Catch;
using namespace Detail;
using namespace std;
using namespace ibex;
using namespace intflpl;

TEST_CASE("LU")
{
  SECTION("Test LU with ID left")
  {
     IntervalMatrix M(3,5);
     M[0][0]=1.0; M[0][1]=0.0; M[0][2]=0.0; M[0][3]=-2.0; M[0][4]=-2.5;
     M[1][0]=0.0; M[1][1]=1.0; M[1][2]=0.0; M[1][3]=2.0; M[1][4]=-2.5;
     M[2][0]=0.0; M[2][1]=0.0; M[2][2]=1.0; M[2][3]=-2.0; M[2][4]=1.5;

     IntLU lu(M,true,true); 
     Interval d = lu.getDeterminant();
     CHECK (d == Interval::one());

     IntervalVector A(3,1.0);
     IntervalVector X = lu.MXeqA(A);    /* X = (1,1,1,0,0) */
     CHECK (X[0] == Interval::one());
     CHECK (X[1] == Interval::one());
     CHECK (X[3] == Interval::zero());
     CHECK (X[4] == Interval::zero());
     IntervalVector Xb = lu.MbXeqA(A);  /* X = A */
     CHECK (Xb[0] == Interval::one());
     CHECK (Xb[1] == Interval::one());
     CHECK (Xb[2] == Interval::one());
     IntervalVector Xc = lu.MbXeqCol(3); /* X = M[.][3] */
     CHECK (Xc[0] == Interval(-2.0));
     CHECK (Xc[1] == Interval(2.0));
     CHECK (Xc[2] == Interval(-2.0));
     IntervalVector Xd = lu.XMbeqA(A); /* X = A */
     CHECK (Xd[0] == Interval::one());
     CHECK (Xd[1] == Interval::one());
     CHECK (Xd[2] == Interval::one());
     IntervalVector Xe = lu.XMeqA(X); /* Xe = A */
     CHECK (Xe[0] == Interval::one());
     CHECK (Xe[1] == Interval::one());
     CHECK (Xe[2] == Interval::one());
     IntervalVector Xf = lu.extendXMeqA(X); /* Xf = (1,1,1,-2.0,-3.5) */
     CHECK (Xf[0] == Interval::one());
     CHECK (Xf[1] == Interval::one());
     CHECK (Xf[2] == Interval::one());
     CHECK (Xf[3] == Interval(-2.0));
     CHECK (Xf[4] == Interval(-3.5));
  }

  SECTION("Test LU with ID left and exchange col")
  {
     IntervalMatrix M(3,5);
     M[0][0]=1.0; M[0][1]=0.0; M[0][2]=0.0; M[0][3]=-2.0; M[0][4]=-2.5;
     M[1][0]=0.0; M[1][1]=1.0; M[1][2]=0.0; M[1][3]=2.0; M[1][4]=-2.5;
     M[2][0]=0.0; M[2][1]=0.0; M[2][2]=1.0; M[2][3]=-2.0; M[2][4]=1.5;

     IntLU lu(M,true,true); 
     lu.exchangeColsBasis(3,0);
      /* L-1 = {{1,0,0},{1,1,0},{-1,0,1}}
         M = {{1,0,0,-2,-2.5},{1,1,0,0,-5.0},{-1,0,1,0,4}} */
     IntervalVector A(3,1.0);
     Interval d2 = lu.getDeterminant();
     CHECK (d2 == Interval(-2.0));
     IntervalVector X2 = lu.MXeqA(A);    /* X = (0,2,0,-0.5,0) */
     CHECK (X2[0] == Interval::zero());
     CHECK (X2[1] == Interval(2.0));
     CHECK (X2[2] == Interval::zero());
     CHECK (X2[3] == Interval(-0.5));
     CHECK (X2[4] == Interval::zero());
     IntervalVector X2b = lu.MbXeqA(A);  /* X = (-0.5,2,0) */
     CHECK (X2b[0] == Interval(-0.5));
     CHECK (X2b[1] == Interval(2.0));
     CHECK (X2b[2] == Interval::zero());
     IntervalVector X2c = lu.MbXeqCol(0); /* X = (-0.5,1,-1) */
     CHECK (X2c[0] == Interval(-0.5));
     CHECK (X2c[1] == Interval(1.0));
     CHECK (X2c[2] == Interval(-1.0));
     IntervalVector X2d = lu.XMbeqA(A); /* X = (-0.5,1,1) */
     CHECK (X2d[0] == Interval(-0.5));
     CHECK (X2d[1] == Interval(1.0));
     CHECK (X2d[2] == Interval(1.0));
     IntervalVector X2e = lu.XMeqA(X2); /* X = (2.25,2,0) */
     CHECK (X2e[0] == Interval(2.25));
     CHECK (X2e[1] == Interval(2.0));
     CHECK (X2e[2] == Interval::zero());
     IntervalVector X2f = lu.extendXMeqA(X2); /* Xf = (2.25,2,0,-0.5,-10.625) */
     CHECK (X2f[0] == Interval(2.25));
     CHECK (X2f[1] == Interval(2.0));
     CHECK (X2f[2] == Interval::zero());
     CHECK (X2f[3] == Interval(-0.5));
     CHECK (X2f[4] == Interval(-10.625));
  }
}

