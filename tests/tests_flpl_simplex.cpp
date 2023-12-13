#include <cstdio>
#include "catch.hpp"
#include <ibex.h>
#include <flpl.h>

using namespace Catch;
using namespace Detail;
using namespace std;
using namespace ibex;
using namespace intflpl;

TEST_CASE("Simplex")
{
  SECTION("Test simplex")
  {
     /* start with a kind of octahedron */
     /* initial box */
     IntervalVector B(3);
     B[0]=Interval(-2.5,2); B[1]=Interval(-2,2); B[2]=Interval(-2,2); 
     std::vector<cstrrhs_status> clstats(3,0);
     Intsimplex simp(B,4,clstats);

     Vector cst(3,1.0); Interval val(-3.0,3.0);
     CHECK(simp.load_constraint(cst,val,0)==3);
     cst[0]=-1.0; 
     CHECK(simp.load_constraint(cst,val,0)==4);
     cst[1]=-1.0; 
     CHECK(simp.load_constraint(cst,val,0)==5);
     cst[0]=1.0; 
     CHECK(simp.load_constraint(cst,val,0)==6);

     /* max X */
     simp.generate_init_basis(0,true);
     simplex_ret ret = simp.simplex_mat();
     CHECK(ret==0);
     double opt = simp.get_objective_value();
     CHECK(opt==2.0);
     /* min X */
     simp.generate_init_basis(0,false);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==2.5);
     /* min X+Y+Z */
     simp.generate_init_basis(3,false);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==3.0);
     /* max 2X-Y+Z */
     IntervalVector Ob(3,0);
     Ob[0]=2.0; Ob[1]=-1.0; Ob[2]=1.0;
     simp.generate_init_basis(Ob,true);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==5.0);
     /* min 2X-Y+Z */
     simp.generate_init_basis(Ob,false);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==5.5);

     /* change objective row */
     simp.changeObjRowCol(3,Interval(1,5),0);
     simp.changeObjRowCol(4,Interval(-5,-1),0);
     /* max X */
     simp.generate_init_basis(0,true);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==2.0);
     /* min X */
     simp.generate_init_basis(0,false);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==-1.0);
     /* min X+Y+Z */
     simp.generate_init_basis(3,false);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==-1.0);
     /* max 2X-Y+Z */
     Ob[0]=2.0; Ob[1]=-1.0; Ob[2]=1.0;
     simp.generate_init_basis(Ob,true);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==5.0);
     /* min 2X-Y+Z */
     simp.generate_init_basis(Ob,false);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==0.0);

/* change objective row (solutions (2,1,0) and (2,0,1) */
     simp.changeObjRowCol(3,Interval(3,5),0);
/* max X */
     simp.generate_init_basis(0,true);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==2.0);
     /* min X */
     simp.generate_init_basis(0,false);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==-2.0);
     /* min X+Y+Z */
     simp.generate_init_basis(3,false);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==-3.0);
     /* max 2X-Y+Z */
     Ob[0]=2.0; Ob[1]=-1.0; Ob[2]=1.0;
     simp.generate_init_basis(Ob,true);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==5.0);
     /* min 2X-Y+Z */
     simp.generate_init_basis(Ob,false);
     ret = simp.simplex_mat();
     CHECK(ret==0);
     opt = simp.get_objective_value();
     CHECK(opt==-3.0);

/* change objective row (empty) */
     simp.changeObjRowCol(3,Interval(3.01,5),0);
/* max X */
     simp.generate_init_basis(0,true);
     ret = simp.simplex_mat();
     CHECK(ret==(1<<INFEASIBLE));
     /* min X */
     simp.generate_init_basis(0,false);
     ret = simp.simplex_mat();
     CHECK(ret==(1<<INFEASIBLE));
     /* min X+Y+Z */
     simp.generate_init_basis(3,false);
     ret = simp.simplex_mat();
     CHECK(ret==(1<<INFEASIBLE));
     /* max 2X-Y+Z */
     Ob[0]=2.0; Ob[1]=-1.0; Ob[2]=1.0;
     simp.generate_init_basis(Ob,true);
     ret = simp.simplex_mat();
     CHECK(ret==(1<<INFEASIBLE));
     /* min 2X-Y+Z */
     simp.generate_init_basis(Ob,false);
     ret = simp.simplex_mat();
     CHECK(ret==(1<<INFEASIBLE));

  }
}

