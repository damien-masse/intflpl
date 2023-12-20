/**
 * Polyhedra using intervals
 * ----------------------------------------------------------------------------
 *  \date       2023
 *  \author     Damien Massé
 *  \copyright  Copyright 2023
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */


#ifndef __INTCSTRVECT_H__
#define __INTCSTRVECT_H__

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <vector>
#include <map>
#include <list>
#include <cmath>

#include "flpl_def.h"

namespace intflpl {

CstrVect traduit_vect(const IntervalVector &box,
		const Vector &v, Interval &bounds);
CstrVect traduit_vect(const IntervalVector &box,
		const IntervalVector &v, Interval &bounds);

inline bool operator==(const CstrVect &lhs, const CstrVect& rhs) {
   return (lhs.bdim==rhs.bdim && lhs.vdim==rhs.vdim &&
          lhs.vect==rhs.vect);
}
inline bool operator!=(const CstrVect &lhs, const CstrVect& rhs) {
   return (lhs.bdim!=rhs.bdim || lhs.vdim!=rhs.vdim ||
          lhs.vect!=rhs.vect);
}
      

struct CstrVectComp {
inline bool operator()(const CstrVect& lhs, const CstrVect& rhs) const
    {
        if (lhs.bdim<rhs.bdim) return true;
        else if (lhs.bdim>rhs.bdim) return false;
        if (lhs.vdim<rhs.vdim) return true;
        else if (lhs.vdim>rhs.vdim) return false;
        for (int i=0;i<lhs.vect.size();i++) 
           if (lhs.vect[i]<rhs.vect[i]) return true;
           else if (lhs.vect[i]>rhs.vect[i]) return false;
        return false;
    }
};

struct CstrVectMap : std::map<CstrVect,CstrRhs,CstrVectComp> {
//      CstrVectMap();
//      CstrVectMap(const CstrVectMap &cvm);
      
#if 0
      bool and_constraint(const IntervalVector &box,
			  const Vector &v, Interval &&i);
#endif
      POLYOP_RET and_constraint(const IntervalVector &box,
			  const CstrVect &v, const Interval &bounds,
			  const cstrrhs_status stat=0);
      POLYOP_RET and_constraint(const IntervalVector &box,
			  const CstrVect &v, const CstrRhs &rhs);
#if 0
      bool and_constraint(const IntervalVector &box,
			  const CstrVect &&v, Interval &&i);
      bool and_constraint(const IntervalVector &box,
			  CstrVect &&v, const Interval &i);
      bool and_constraint(const IntervalVector &box,
			  CstrVect &&v, Interval &&i);
#endif
//      friend class ExpPoly;
};

int simplify_polyhedron(int dim, IntervalVector &ivbox,
     std::vector<cstrrhs_status> &clstats,
     const std::vector<std::pair<IntervalVector,Interval>> &csts,
     CstrVectMap &csts_rs, IntervalVector &inc);

int simplify_polyhedron(int dim, IntervalVector &ivbox,
     std::vector<cstrrhs_status> &clstats, CstrVectMap &csts,
    IntervalVector &inc, cstrrhs_status recheck);

Interval simplex_form(int dim, const IntervalVector &ivbox,
     std::vector<cstrrhs_status> clstats,
     const CstrVectMap &csts,
     const Vector &obj);

IntervalVector generate_including_box(const IntervalVector &ivbox,
                const  Intsimplex &simp,const CstrVectMap &csts,
                const std::vector<CstrVectMap::iterator> &intmap);

std::list<Vector>
                generate_facet(const IntervalVector &ivbox,
                 std::vector<cstrrhs_status> clstats,
                 const CstrVectMap &csts, int dm, bool upper);

std::list<Vector>
                generate_facet(const IntervalVector &ivbox,
                 std::vector<cstrrhs_status> clstats,
                 const CstrVectMap &csts, const CstrVect &cstF, bool upper);

/** constraint satisfy vector */
inline bool cst_satisfy_vec(const std::pair<CstrVect,CstrRhs> &cst, const Vector &v) {
   return cst.second.val.contains(cst.first.vect*v);
}

/** constraint intersect box */
inline bool cst_intersect_box(const std::pair<CstrVect,CstrRhs> &cst, const IntervalVector &bx) {
   return cst.second.val.intersects(cst.first.vect*bx);
}

/** box is in constraint */
inline bool cst_satisfy_box(const std::pair<CstrVect,CstrRhs> &cst, const IntervalVector &bx) {
   return (cst.first.vect*bx).is_subset(cst.second.val);
}
/** compute the "distance" a constraint+box has wrt a constraint */
double dist_csts_box(const IntervalVector &ivbox,
       const Vector &v1, const Interval& b, const Vector &obj);

}

#endif
