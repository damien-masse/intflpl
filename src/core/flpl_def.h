/**
 * Common definitions for simplex and polyhedra with intervals
 * ----------------------------------------------------------------------------
 *  \date       2023
 *  \author     Damien Massé
 *  \copyright  Copyright 2023
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#ifndef __FLPL_DEF_H__
#define __FLPL_DEF_H__

#include <bitset>

/* IBEX is used to handle matrix of intervals,
   may be changed if another library is used */
#include <ibex.h>
using Vector = ibex::Vector;
using Matrix = ibex::Matrix;
using Interval =  ibex::Interval;
using IntervalVector = ibex::IntervalVector;
using IntervalMatrix = ibex::IntervalMatrix;

namespace intflpl {

class ExpPoly;

class Intsimplex;

/* constraint in the polyhedron */

struct CstrVect {
   int bdim; /* dim of greatest value (=1) */
   double vdim; /* value (between -1 and 1) of second greatest absolute val
                   note : 0 should not be possible */
   Vector vect;  /* the constraint */

   CstrVect(int bdim, double vdim, const Vector &vect);
   friend bool operator==(const CstrVect &lhs, const CstrVect& rhs);
   friend bool operator!=(const CstrVect &lhs, const CstrVect& rhs);
};

enum CSTRRHS_STATUS {
      LB_NR , /* lower bound is non-redundant */
      LB_R , /* lower bound is redundant */
      LB_RT , /* lower bound is redundant and tight (implies LB_R) */
      UB_NR , /* upper bound is non-redundant */
      UB_R , /* upper bound is redundant and tight (implies LB_R) */
      UB_RT , /* upper bound is redundant and tight */
      SIZECSTRRHS_STATUS
};
using cstrrhs_status = std::bitset<SIZECSTRRHS_STATUS>;
constexpr cstrrhs_status rhs_lbfilter = (1<<LB_NR) | (1<<LB_RT) | (1<<LB_R);
constexpr cstrrhs_status rhs_ubfilter = (1<<UB_NR) | (1<<UB_RT) | (1<<UB_R);
constexpr cstrrhs_status rhs_rfilter = (1<<LB_R) | (1<<UB_R);
constexpr cstrrhs_status rhs_rtfilter = (1<<LB_RT) | (1<<UB_RT);
constexpr cstrrhs_status rhs_nrfilter = (1<<LB_NR) | (1<<UB_NR);
constexpr cstrrhs_status rhs_lbnr = (1<<LB_NR);
constexpr cstrrhs_status rhs_ubnr = (1<<UB_NR);
constexpr cstrrhs_status rhs_lbr = (1<<LB_R);
constexpr cstrrhs_status rhs_ubr = (1<<UB_R);
constexpr cstrrhs_status rhs_lbrt = (1<<LB_R)|(1<<LB_RT);
constexpr cstrrhs_status rhs_ubrt = (1<<UB_R)|(1<<UB_RT);

/** status of polyhedral operation return */
enum POLYOP_RET {
	POL_UNCHANGED , /* no change to the polyhedron */
        POL_CHANGED ,   /* changed */
        POL_EMPTY       /* the polyhedron is now empty */
};

struct CstrRhs {
    Interval val;
    cstrrhs_status status;
};

/** status of simplex return */ 
enum SIMPLEX_RET { 
                  INFEASIBLE , /* true : the polyhedron is empty */ 
                  UNBOUNDED , /* true : unbounded (not possible) */ 
                  APPROXBASIS , /* true : the basis is approximated 
                                 false : the basis is always correct */ 
                  REACHEDLIMIT , /* true : the algorithme reached the 
                                    limit of iterations */ 
                  SIZESIMPLEX_RET  
              }; 
using simplex_ret = std::bitset<SIZESIMPLEX_RET>; 

struct CstrVectMap;

/* a few "scripts" for bitsets */

template <std::size_t N>
inline void assign_filter(std::bitset<N> &a1, const std::bitset<N> &a2, 
				const std::bitset<N> &filter) {
   a1 = (a1 & ~filter) | (a2 & filter);
}

template <std::size_t N>
inline void reset_filter(std::bitset<N> &a1, const std::bitset<N> &filter) {
   a1 &= ~filter;
}

template <std::size_t N>
inline bool all_filter(std::bitset<N> &a1, const std::bitset<N> &filter) {
   return ((~a1) & filter).none();
}

}

#endif
