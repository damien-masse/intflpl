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

/** representation of polyhedron with ivbox and interval constraints */
class ExpPoly
{
   public :
      /** build an empty or full polyhedron **/
      ExpPoly (int dim,bool empty);
      /* with a box */
      ExpPoly (const IntervalVector &Box);
      /* with a box and a set of (interval vector) constraints */
      ExpPoly (const IntervalVector &Box, const std::vector<std::pair<IntervalVector,Interval>> &Csts);
      /* clone */
      ExpPoly (const ExpPoly &P);
      
      int get_dim() const;
      int get_not_flat_dim() const; 
      const IntervalVector &getBox() const;
      const CstrVectMap &getCsts() const;

      /* comparaison */
      bool is_empty() const;
      bool is_flat() const;
      bool is_bounded() const;
      bool is_box() const;
      bool satisfy_cst (const std::pair<const CstrVect,CstrRhs> &cst) const;
			  /* true is guaranteed */
      bool is_subset (const ExpPoly &Q) const;  /* true is guaranteed */
      bool is_subset (const IntervalVector &IV) const; /* true is guaranteed */
      bool is_superset (const IntervalVector &IV) const; /* true is guaranteed */

      /* internal modification */
      void minimize(cstrrhs_status ct);     

      /* modification */
      void set_empty();
      ExpPoly &operator&=(const IntervalVector &iv);
      ExpPoly &operator|=(const IntervalVector &iv); /* non optimal union */
      ExpPoly &operator&=(const ExpPoly &Q);
      ExpPoly &operator|=(const ExpPoly &Q); /* non optimal union */
      ExpPoly &operator&=(const std::vector<std::pair<IntervalVector, Interval>> &Res);
      void intersect_paral(const IntervalMatrix &M, const IntervalVector &V);
      ExpPoly &widen(const ExpPoly &Q);      /* union + widening */
      friend ExpPoly operator&(const ExpPoly &C1, const ExpPoly &C2);
      friend ExpPoly operator|(const ExpPoly &C1, const ExpPoly &C2);
      friend void diff_hull_box(ExpPoly &B1, const IntervalVector &IV);
      friend bool operator==(const ExpPoly &C1, const ExpPoly &C2); 


      void unflat(int dm, Interval offset);

      /* crée une description des contraintes centrée sur C,
           avec si possible la première colinéaire à Z1 et les autres
           orthogonales */
      std::vector<std::pair<IntervalVector,Interval>> 
	build_constraints_for_propag(const Vector &C, const Vector &Z1) const;
      /* compute the intervals of the domain */
      IntervalVector 
	build_constraints_for_propag(const IntervalMatrix &Z) const;


						/* true is guaranteed */
      friend IntervalVector& operator&=(IntervalVector &V, const ExpPoly &C);
      friend IntervalVector& operator|=(IntervalVector &V, const ExpPoly &C);
      friend IntervalVector operator+(const ExpPoly &C, const IntervalVector &V);
        
      friend std::ostream& operator<<(std::ostream& str, const ExpPoly& C);
      
      void vertices2D(std::vector<double>&x, std::vector<double>&y);
      std::list<Vector> facet3D() const;
      std::list<std::list<Vector>> getFacets3D(const IntervalVector& iv, bool with_doors=true) const;


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
      void compute_dim_not_flat();
      
      bool add_cst(const IntervalVector& V, const Interval& I);

      void intersect_box(const IntervalVector &iv); 
};

inline int ExpPoly::get_dim() const { return this->dim; }
inline int ExpPoly::get_not_flat_dim() const { return this->dim_not_flat; }
inline const IntervalVector &ExpPoly::getBox() const { return this->Box; }
inline const CstrVectMap &ExpPoly::getCsts() const { return this->csts; }
inline bool ExpPoly::is_empty() const { return this->Box.is_empty(); }
inline bool ExpPoly::is_flat() const { return this->dim_not_flat<this->dim; }
inline bool ExpPoly::is_bounded() const {return !this->Box.is_unbounded(); }
inline bool ExpPoly::is_box() const { return this->csts.size()==0; }
inline void ExpPoly::set_empty() { 
	this->dim_not_flat=-1; 
	this->Inc.set_empty(); 
	this->Box.set_empty(); this->csts.clear(); this->minimized=0; }
}

#endif
      
     
