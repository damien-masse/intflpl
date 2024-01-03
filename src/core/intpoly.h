/**
 * Polyhedra using intervals
 * ----------------------------------------------------------------------------
 *  \date       2023
 *  \author     Damien Massé
 *  \copyright  Copyright 2023
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */


#ifndef __INTPOLY_H__
#define __INTPOLY_H__

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <vector>
#include <map>
#include <list>
#include <cmath>

#include "flpl_def.h"
#include "intcstrvect.h"

namespace intflpl {

class IntPoly;

IntPoly sum_tau(const IntPoly& iv, const IntervalVector& V,
                                                bool keep=false);

/** representation of polyhedron with ivbox and interval constraints */
class IntPoly
{
   public :
      /** empty constructor, cannot be used except with a resize()
        * or an intersection with an ivbox */
      IntPoly ();

      /** constructor from dimension (full domain)
       */
      IntPoly (int dim);

      /** build an empty or full polyhedron **/
      IntPoly (int dim,bool empty);

      /* with a box */
      IntPoly (const IntervalVector &Box);

      /* copy */
      IntPoly (const IntPoly &P);
      
      /* with a box and a set of (interval vector) constraints */
      IntPoly (const IntervalVector &Box, const std::vector<std::pair<IntervalVector,Interval>> &Csts);

      /* from a slanted box  ( M.V, with reverse known )
      */
      IntPoly(const IntervalMatrix& M, const IntervalMatrix &rM,
                        const IntervalVector &V);

      /* from a slanted box  ( M.V, reverse computed )
      */
      IntPoly(const IntervalMatrix& M, const IntervalVector &V);


      /***************** ACCESS **************************/
      /**
       * number of variables
       */
      int get_dim() const;
      int size() const; /* equivalent to get_dim */

      /**
       * number of constraints (not including bounding box)
       */
      int get_nbcsts() const; 

      /**
       * non-flat dim (in the bounding box)
       */
      int get_not_flat_dim() const; 
      /** is_empty
       */
      bool is_empty() const;
      /** bounding box
       */
      const IntervalVector& getBox() const;
      const IntervalVector& box() const;

      /** complete map of constraints */
      const CstrVectMap &getCsts() const;

      /** is unbounded
       */
      bool is_unbounded() const;
      /** is bounded
       */
      bool is_bounded() const;

      /** has a flat dimension */
      bool is_flat() const;

      /** is a box */
      bool is_box() const;

      /** ``component'' : interval of the bounding box 
        * (only const !!!) */
      const Interval& operator[](unsigned int i) const;
      /**
       * ``middle'' (e.g. center of the bounding box ?)
       */
      Vector mid() const;


      /* satisfaction of constraint :
       * true : yes      false : maybe */
      bool satisfy_cst (const std::pair<const CstrVect,CstrRhs> &cst) const;

      /** contains point (approximate...) more or less yes but not completely
       */
      bool contains(const Vector& iv) const;

      /** intersects (approximate) : 
       *  true = maybe           false = no */
      bool intersects(const IntervalVector& x) const;
      /** intersects (approximate) : 
       *  true = maybe 		 false = no */
      bool intersects(const IntPoly& x) const;

      /**
       * relative distance, fast?
       * (e.g. use the bounding box to make the difference between
       *    close constraints ?) 
       */
      double rel_distanceFast(const IntPoly& iv) const;

      /************* Modification **********************/
      /** empty */
      void set_empty();
      /** to 0 */
      void clear();
      /** setComponent : (very) approximate equivalent to 
       *  x[i] = interval  */
      void setComponent(int i, const Interval& a);
      /** inflation by a cube
       *  this <- this + [-r,r]^d
       *  return this
       */
      IntPoly& inflate(double rad);
      /** inflation by a ball ?
       *  this <- this + ball(rad)
       *  return this
       */
      /* expansion of a flat dimension */
      void unflat(int dm, Interval offset);
      /* IntPoly& inflate(double rad); */
      /** centered homothety
       *  x <- [c] + delta*([x]-[c]) ( (1-delta)[c] + delta*[x] )
       *  return this
       */
      IntPoly& homothety(IntervalVector c, double delta);
      /** set this to x
       */
      IntPoly& operator=(const IntervalVector& x);
      /** "project" x on this, get an overapproximation
       */
      /* IntPoly& assign(const IPoly& iv); */

      /** fast expansion from a bigger set (i.e. iv > this)
           this |= this+(iv-this)*fact */
      void inflate_from_baseFast(const IntPoly& iv, double fact);

      /** apply a (non-singular) endomorphism ;
          this <- M*this ;    IM : inverse of M */
      IntPoly& linMult(const IntervalMatrix& M, const IntervalMatrix& IM);

      /***** Intersections *****/

      /** intersection with a box */
      IntPoly& operator&=(const IntervalVector& x);

      /** intersection with a box, new value */
      friend IntPoly operator&(const IntPoly& iv, const IntervalVector& x);
      /** fast intersection, with the same constraint base */
      IntPoly& meetFast(const IntPoly& iv);
      /** meet with keeping the constraint base */
      IntPoly& meetKeep(const IntPoly& iv);

      /** optimal intersection */
      IntPoly &operator&=(const IntPoly &Q); 
      friend IntPoly operator&(const IntPoly &C1, const IntPoly &C2);
      /** different intersections ? */
      IntPoly &meet(const IntPoly& iv, bool ctcG=false);

      /** intersections with interval linear constraints,
       *  possibly with contractors */
      /** intersection with a interval linear constraint,
       *      V x \in b 
       *      ctc = true => contractance is guaranteed 
       *      V and b may also be contracted */
      bool meetLN(IntervalVector& V, Interval& b, bool ctc);
      /** intersection with a centered interval linear constraint,
       *      V (x-c) \in b 
       *      keep = true => contractance is guaranteed */
      bool meetLN(IntervalVector& V, IntervalVector& C,
                                Interval& b, bool ctc);
      /** with a set of linear constraints */
      bool meetLM(IntervalMatrix& S, IntervalVector& b, bool ctc);
      /** with a centered set of linear constraints */
      bool meetLM(IntervalMatrix& S, IntervalVector& C,
                                IntervalVector& b, bool ctc);
      /** with a centered set of linear constraints */
      bool meetLM(IntervalMatrix& S, const Vector& C,
                                IntervalVector& b, bool ctc);
      /** other version... */
      IntPoly &operator&=(const std::vector<std::pair<IntervalVector, Interval>> &Res);

      /** unions */

      /** union with a box 
       */
      IntPoly& operator|=(const IntervalVector& iv);
      friend IntPoly operator|(const IntPoly& iv, const IntervalVector& x);
      /** approximate (for now) union for polyhedra */
      IntPoly& operator|=(const IntPoly& Q);
      friend IntPoly operator|(const IntPoly &C1, const IntPoly &C2);

      /** union + widening */
      IntPoly &widen(const IntPoly &Q);      /* union + widening */

      /** functions returning boxes */

      friend IntervalVector operator&(const IntervalVector& x,
						 const IntPoly& iv);
      friend IntervalVector& operator&=(IntervalVector &V, const IntPoly &C);
      friend IntervalVector& operator|=(IntervalVector &V, const IntPoly &C);
      friend IntervalVector operator+(const IntervalVector &V, const IntPoly &C);
      /** product : compared with linMult :
           M is (generally) small and contains singular matrices 
           we just want an IntervalVector */
      friend IntervalVector operator*(const IntervalMatrix& M, 
                                const IntPoly& iv);

      /***** operations  *****/

      /** sum, difference */
      IntPoly& operator+=(const IntervalVector& V);
      IntPoly& operator-=(const IntervalVector& V);
      friend IntPoly operator+(const IntPoly& iv, const IntervalVector& V);
      friend IntPoly operator-(const IntPoly& iv, const IntervalVector& V);
      /* fast operations */
      IntPoly& sumFast(const IntPoly& iv);
      IntPoly& diffFast(const IntPoly& iv);

      /*** functions for differential inclusions */

      /** centered multiply and add 
          this <- center + M *(this-center) + V
       */
      void cmult_and_add(const Vector&center,
                const IntervalMatrix& M, const IntervalMatrix &invM,
                const IntervalVector& V);
      /** tau-centered multiply and add (not modifying the matrices)
       * IV = IV + [0,1]*(M (IV-center) + V)
       * quick algorithm, maybe not precise (generally M is thick)
       */
      void ctau_mult_and_add(const Vector&center,
                const IntervalMatrix& M,
                const IntervalVector& V);

      friend IntPoly sum_tau(const IntPoly& iv, const IntervalVector& V,
                                                bool keep);
#if 0
      /**
       * union with     box[d->val] 
                                intersected with 
                        iv + ]0,1]*(M (iv-center)+V) */
      bool join_intersect_with_tau
                (const IntPoly& iv, const Vector &center,
                const IntervalMatrix& M,
                const IntervalVector& V,
                 const IntervalVector& box, int d, double val);
#endif

      void intersect_paral(const IntervalMatrix &M, const IntervalVector &V);

      /* crée une description des contraintes centrée sur C,
           avec si possible la première colinéaire à Z1 et les autres
           orthogonales */
      std::vector<std::pair<IntervalVector,Interval>> 
	build_constraints_for_propag(const Vector &C, const Vector &Z1) const;
      /* compute the intervals of the domain */
      IntervalVector 
	build_constraints_for_propag(const IntervalMatrix &Z) const;

      /** comparisons */
      bool is_subset (const IntPoly &Q) const;  /* true is guaranteed */
      bool is_subsetFast (const IntPoly &Q) const;  /* true is guaranteed */
      bool is_subset (const IntervalVector &IV) const; /* true is guaranteed */
      bool is_superset (const IntervalVector &IV) const; /* true is guaranteed */
      friend bool operator==(const IntPoly &C1, const IntPoly &C2); 

      /* modification */
      friend void diff_hull_box(IntPoly &B1, const IntervalVector &IV);

      /** display
       */
      friend std::ostream& operator<<(std::ostream& str, const IntPoly& C);
      
      /** representations */

      /** generate a ConvexPolygon which is an overapproximation of the
       * projection of the polyhedron (basic form)
       */
//      ConvexPolygon over_polygon(const Matrix& M) const;

      void vertices2D(std::vector<double>&x, std::vector<double>&y) const;
      std::list<Vector> facet3D() const;
      std::list<std::list<Vector>> getFacets3D(const IntervalVector& iv, bool with_doors=true) const;

      /* internal modification */
      void minimize(cstrrhs_status ct);     

  private :
      int dim; /* dimension of the space */
      int dim_not_flat;  /* number of non-flat dimension, empty=0 */

      IntervalVector Box; /* bounding box */
      std::vector<cstrrhs_status> BoxStat; /* are box csts redundants? */
      IntervalVector Inc; /* one box inside, if possible */
      cstrrhs_status minimized;    
	/* 0 : things are what they are...
           rhs_nrfilter : some NR may be in fact R (intersection)
           rhs_rfilter : must recheck R (union, as it is non-exact)
           rhs_rtfilter : must recompute RT (union and intersection)
	*/
      CstrVectMap csts;

      enum INTPOLY_STATUS {
           INTPOLY_INVALID, /* not usable: no dimension */
           INTPOLY_EMPTY,   /* empty */
           INTPOLY_UNBOUNDED,   /* unbounded */
           SIZE_INTPOLY_STAT
      };

      std::bitset<SIZE_INTPOLY_STAT> pstatus; /* private status */
      

      void compute_dim_not_flat();
      bool add_cst(const IntervalVector& V, const Interval& I);
      void intersect_box(const IntervalVector &iv); 
      /** fast (safe) bound on a constraint (in form CstrVect) */
      Interval bound_cstFast (const CstrVect &cvect) const;
};

#include "intpoly_inline.hpp"
}

#endif
      
     
