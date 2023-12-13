/**
 * Simplex manipulations for interval matrices
 * ----------------------------------------------------------------------------
 *  \date       2023
 *  \author     Damien Massé
 *  \copyright  Copyright 2023
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#ifndef _INTSIMPLEX_H
#define _INTSIMPLEX_H

#include <vector>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <bitset>

#include "flpl_def.h"
#include "permut.h"
#include "intLU.h"


namespace intflpl {

class Intsimplex {
    public:

    /** status of simplex matrix */
    enum SIMPLEX_STAT {
                     EMPTY , /* true : the matrix has no usable constraint */
		     IDENTLEFT , /* the left part of the matrix is id mat,
				    but we have not built the basis still  */
                     INITBASIS , /* the matrix has a initial basis */
                     CONSBASIS , /* the matrix has a consistent basis */
                     OPTBASIS ,  /* the matrix has a "optimal" basis */
		     UNSAT ,     /* basis building has failed */
		     UNFEAS ,    /* empty polyhedron */
	             INVALCOL ,  /* basis built, but the last column is
					 inconsistent */
	             INVALROW ,  /* basis built, but the last row is
					 inconsistent */
		     SIZESIMPLEX_STAT
                 };
    using simplex_stat = std::bitset<SIZESIMPLEX_STAT>;

    /* build an "empty" simplex matrix, on point 0,0 */
    Intsimplex(int dim, int szalloc, bool useInterval=false);
    /* build an "empty" simplex matrix from an ivBox, allocating the
       place for the other constraints */
    Intsimplex(const IntervalVector &box, int szalloc, std::vector<cstrrhs_status> clstats, bool useInterval=false);
    /* duplicate the simplex */
//    Intsimplex(const Intsimplex &is);

    /* destructeur */
    ~Intsimplex();

    /* restore the simplex */
    void resetComputation();

    /* add constraints ; return the number of the column added */
    void load_box(const IntervalVector &box, std::vector<cstrrhs_status> clstats);
    int load_constraint(const Vector &cst, const Interval &val,
		       cstrrhs_status clstat);
    int load_constraint(const IntervalVector &cst, const Interval &val,
		       cstrrhs_status clstat);
    int load_constraint(const std::pair<CstrVect,CstrRhs> newcst);
    
    /* get the extremal (vertice) point, defined as the difference
       between the first dim objective columns and the ivbox */
    Vector getExtremalPoint(const IntervalVector &ivbox) const;
    

    /* modify the values of the last column, after computing the 
       simplex. returns true if the basis is still possible (no change
       of sign) */
    bool translateLastColumn(const IntervalVector &column);
    /* change one coef of the last row (e.g. to bound it more) */
    void changeObjRowCol(int col, const Interval& val, cstrrhs_status clstat);
    void extendObjRowCol(int col, double val);
    /* change one coef of the last row (e.g. to bound it more) */
    void activeCol(int col, bool activate);
    /* change the last row. must be done before computing the simplex */
    void loadObjRow(const IntervalVector &row, 
		     std::vector<cstrrhs_status> clstats);
    /* modify the values of the last column, after computing the 
       simplex. returns true if the basis is still consistent (no value
       outside 0) */
    bool translateObjRow(const IntervalVector &row);

    void lock_row(int row, bool locked);

    /** status of simplex return (simplex_ret) defined in explibdef.h */
    const simplex_ret simplexret_Infeasible = (1 << INFEASIBLE);
    const simplex_ret simplexret_MaxIter = 
			(1 << APPROXBASIS) | (1 << REACHEDLIMIT);

    /** constraint the matrix using ibex::bwd_mul. can be useful
     *  e.g. if the leftinitial basis is unbounded, as it may give a bound.
     *  if returns is false, it means that even ibex have seen that the
     *  system is empty */
    // bool constraint_matrix();

    /* load the objective column and generate the first basis.
       also, save the matrix in matinit */
    /* getub = true : get the max (use column)
       getub = false : get the min (use -column)
    */
    bool generate_init_basis(const IntervalVector &column, bool getub);
    /* when the objective is an existing column */
    bool generate_init_basis(int col, bool getub);

    /* simplex solving from a feasible matrix */
    simplex_ret simplex_mat();

    const Permut& get_basis() const;
    bool is_neg_basis(int col) const;
    const IntervalMatrix& get_mat() const;
    IntervalVector get_objrow() const;
    double get_objective_value() const;

    /** status of column in the simplex */
    enum COLSTATUS { HASLB ,    /* inf-bounded */
   		     HASUB ,    /* sup-bounded */
		     INBASIS ,  /* in basis & can go out */
		     OUTBASIS , /* out basis & can go in */
		     NEGBASIS , /* negative if in basis */
		     LOCKEDIN , /* locked in basis */
		     UNUSED ,   /* out basis, to ignore */
	             SIZECOLSTATUS };
   using colstatus = std::bitset<SIZECOLSTATUS>;
   const colstatus colstatus_Changebasis = (1 << INBASIS) | (1 << OUTBASIS);
   const colstatus colstatus_Unused = (1 << UNUSED);
   const colstatus colstatus_BasisFilter = (1 << UNUSED) | (1 << INBASIS) |  (1 << OUTBASIS) | (1<<LOCKEDIN);
   const colstatus colstatus_InBasis = (1 << INBASIS);
   const colstatus colstatus_OutBasis = (1 << OUTBASIS);

		    
/** generating the vertices of a facet, from an already optimal matrix.
 *  ivbox is the global box around the matrix.
 *  all basic columns except 2 must be "locked" */
    std::list<Vector> generate_vertices_2D(const IntervalVector &ivbox);

   private:
      int dim, sz, szalloc, numobjcol; /* dim : dimension (number of rows)
		       sz : columns used (+dim)
                       szalloc : allocated place in the matrix, +dim 
		       objcol : column used for the objective */
      int objcolInd;    /* =0 : lastcol is "independant" (=objcol)
                           =1 : lastcol is equal to one constraint 
                           =-1: lastcol is the negation of one constraint */
      simplex_stat status;
      IntervalMatrix matInit; /* initial matrix */
      IntervalMatrix matInitT; /* transpose of initial matrix */
      IntervalVector objcol;   /* objective column, when not last */

      IntervalVector objrowInit; /* initial last row */
      IntervalVector objrow_offset; /* offset for the last row */
      std::vector<colstatus> colstat;
      bool useInterval; /* interval in the initial matrix
			   are used for [-1,+1] columns */
      IntLU *LUform; /* LU form */
      static const int maxit = 1000;

/** change_basis : change the basis of the matrix
 *  oldcol and ncol are the respective leaving and entering column
 *  the pivot may be an interval, everything is applied directly
 */
      void change_basis(int oldcol, int ncol);

      void recompute_objoffset();

/** find entry col : check the successive rows to find 
 *  a leaving coefficient. Gap is the "distance" from the current
 *  solution and the constraint.
 *  return : -1 : all coef are the wrong sign, or the gap can't be bridged:
 *                the polyhedron is empty
 *           -2 : no adequate row has been found, but the polyhedra
 *                may not be empty. Suggest use a different column.
 *           >=0 : row of the possible entry. 
 *                        The pivot must be of constant sign
 */
     int check_entry_col_C(int colR, const IntervalVector &R,
            const IntervalVector &obj, const Interval &ratio,
            const Interval &bd, double asup,
            const Permut &basis, std::vector<int> &sign,
            Interval &result, int rowC) const;
     int check_entry_col(int col, Interval &res, const Permut &basis,
                std::vector<int> &bsign, bool &asup) const;

/** check entry : check if a coefficient can be pivot inside the polyhedra
 *  (vertice generation). We know that basis[row] can exit, we do not
 *  know if col enter
 *  asup : col enter the basis as positive coef
 *  gain : the maximum "gain", i.e. the move on the polyhedron
 *  return : true if definitely possible, false if maybe not possible */
      int check_entry_obj(int colout, const IntervalVector &RowB,
                bool &toup, double &bcoef) const;

    bool generate_init_basis_internal(double mult, const IntervalVector &obj);

};


/** treat interval : change an interval value to a punctual one using
 *  the actual basis
 *  row and col are the place of the interval
 */
//inline void Intsimplex::treat_interval(int row, int col, bool toZ) {
//    double gap = mat[row][col].diam();
////    if (gap==0.0) return;
//    double target;
//    if (toZ || mat[row][col].contains(0.0)) target=0.0;
//    else if (mat[row][col].ub()<0.0) target=mat[row][col].lb();
//    else target=mat[row][col].ub();
//
//    matInitT[col]+=(target-mat[row][col])*matInitT[basis[row]/2];
//    objrowInit[col]+=objrowInit[basis[row]/2]*(target-mat[row][col]);
//    for (int i=0;i<dim;i++) {
//      matInit[i][col]+=(target-mat[row][col])*matInit[i][basis[row]/2];
//    }
//
//    objrow[col]+=objrow[basis[row]/2]*(target-mat[row][col]);
//    mat[row][col] = target;
//}
//
inline void Intsimplex::lock_row(int row, bool locked) {
    assert(status[CONSBASIS]);
    colstat[LUform->getColBasis()[row]][LOCKEDIN]=locked;
}

inline const Permut& Intsimplex::get_basis() const {
    assert(status[CONSBASIS]);
    return this->LUform->getColBasis();
}
inline const IntervalMatrix& Intsimplex::get_mat() const {
    return this->matInit;
}
inline IntervalVector Intsimplex::get_objrow() const {
    return this->objrow_offset+this->objrowInit;
}

inline double Intsimplex::get_objective_value() const {
    /* FIXME : faire un cas particulier quand les n premières colonnes
       sont Id */
    Interval value(0.0);
    IntervalVector obj_local(dim);
    const Permut& basis = get_basis();
    if (objcolInd==0) {
       obj_local = LUform->MbXeqA(this->objcol);
    } else {
       if (basis.is_image(numobjcol)) {
          return (objcolInd==1 ? objrowInit[numobjcol].ub() :
		-objrowInit[numobjcol].lb());
       }
       obj_local = LUform->MbXeqCol(numobjcol);
       if (objcolInd==-1) obj_local = -obj_local;
    }
    for (int row=0;row<dim;row++) {
      value += objrowInit[basis[row]]*obj_local[row];
    }     
    return value.ub();
}

inline void Intsimplex::activeCol(int col, bool activate) {
    colstat[col][UNUSED]=!activate;
}

inline bool Intsimplex::is_neg_basis(int col) const {
    return colstat[col][NEGBASIS];
}
}

#endif
