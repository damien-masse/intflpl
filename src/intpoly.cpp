/**
 * Polyhedra using intervals
 * ----------------------------------------------------------------------------
 *  \date       2023
 *  \author     Damien Massé
 *  \copyright  Copyright 2023
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include <vector>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <ibex.h>
#include "intsimplex.h"
#include "intpoly.h"

namespace intsimplex {

#if 1
static int maxsize=0;
#endif

ExpPoly::ExpPoly(int dim, bool empty) :
  dim(dim), dim_not_flat(empty ? -1 : dim), 
  Box(dim, (empty ? Interval::empty_set() : Interval::all_reals())),
  Inc(Box), BoxStat(dim, 0),
  minimized(0), csts()
{  }

ExpPoly::ExpPoly(const IntervalVector &Box) : 
  dim(Box.size()), dim_not_flat(Box.is_empty() ? -1 : 0),
  BoxStat(dim, (1<<LB_NR) | (1<<UB_NR)),
  minimized(0), Box(Box), Inc(Box), csts() {
  if (dim_not_flat==0) 
    for (int i=0;i<dim;i++) 
         if (!this->Box[i].is_degenerated()) dim_not_flat++;
}

ExpPoly::ExpPoly(const IntervalVector &Box, const std::vector<std::pair<IntervalVector,Interval>> &Csts) :
     ExpPoly(Box) 
{
   this->Inc.set_empty();
   dim_not_flat = simplify_polyhedron(dim,this->Box,this->BoxStat,Csts,this->csts, this->Inc);
}

ExpPoly::ExpPoly (const ExpPoly &P) :
     dim(P.dim), dim_not_flat(P.dim_not_flat), BoxStat(P.BoxStat),
     Inc(P.Inc),
     Box(P.Box), csts(P.csts), minimized(P.minimized) {
}

bool ExpPoly::satisfy_cst (const std::pair<const CstrVect,CstrRhs> &cst) const {
    /* basic checks first */
    if ((cst.first.vect*this->Box).is_subset(cst.second.val)) return true;
    const CstrVectMap::const_iterator el = this->csts.find(cst.first);
    if (el!=this->csts.end()) {
       const Interval& res = (*el).second.val;
       return (res.is_subset(cst.second.val));
    }
    Interval result = simplex_form(dim,this->Box,this->BoxStat,
			this->csts,cst.first.vect);
    return result.is_subset(cst.second.val);
}

bool ExpPoly::is_subset (const ExpPoly &Q) const {
    /* A inclus dans Q si toutes les contraintes de Q sont inutiles */
    /* peut-être pas correct si (*this) n'est pas minimisé */
//    this->minimize(0);
    if (!this->Box.is_subset(Q.Box)) return false;
    for (const std::pair<CstrVect,CstrRhs> &ct : Q.csts) {
        if (!this->satisfy_cst(ct)) return false;
    }
    return true;
}

bool ExpPoly::is_subset (const IntervalVector &IV) const {
    return (this->Box.is_subset(IV));
}


bool operator==(const ExpPoly &C1, const ExpPoly &C2) {
    if (C1.Box != C2.Box) return false;
//    C1.minimize(0);
//    C2.minimize(0);
    return (C1.is_subset(C2) && C2.is_subset(C1));
}

bool ExpPoly::is_superset (const IntervalVector &IV) const {
    /* A inclus dans Q si toutes les contraintes de Q sont inutiles */
    if (!IV.is_subset(this->Box)) return false;
    for (const std::pair<CstrVect,CstrRhs> &ct : this->csts) {
        const CstrVect &C = ct.first;
        const Interval &Ret= ct.second.val; 
        if (!(C.vect*IV).is_subset(Ret)) return false;
    }
    return true;
}


void ExpPoly::compute_dim_not_flat() {
    if (!this->Box.is_empty()) {
      this->dim_not_flat=0;
      for (int i=0;i<dim;i++) {
         if (!this->Box[i].is_degenerated()) this->dim_not_flat++;
      }
    } else this->dim_not_flat=-1;
}

void ExpPoly::minimize(cstrrhs_status ct) {
    minimized |= ct;
    if (this->Box.is_empty()) { this->set_empty(); return; }
/*
    if (this->csts.size()==0) { 
        this->minimized=0; 
	return;
     }
*/
//    std::cerr << "minim before " << (*this) << "\n";
    this->dim_not_flat = simplify_polyhedron(dim,this->Box,this->BoxStat,
			this->csts, this->Inc, minimized);
    if (this->dim_not_flat==-1) { this->set_empty(); return; }
    for (int i=0;i<dim;i++) if (this->BoxStat[i]==0) std::cerr << "erreur minimize1 " << this->csts.size() << "\n";
    for (auto &a:this->csts) if (a.second.status==0) std::cerr << "erreur minimize2\n";
//    std::cerr << "minim after " << (*this) << "\n";
    this->minimized=0;
}

ExpPoly &ExpPoly::operator&=(const IntervalVector &iv) {
    if (this->Box.is_subset(iv)) return (*this);
    this->Box &= iv; /* note all constraints become not guaranteed */
		     /* (redundant one stays redundant but not tight) */
    this->Inc &= iv;
    if (this->Box.is_empty()) { this->dim_not_flat=-1; this->csts.clear();
				return (*this); }
    /* recompute dim_not_flat */
    this->compute_dim_not_flat();
    if (this->csts.size()>0) { 
        this->minimize(rhs_nrfilter | rhs_rtfilter);
    }
    return (*this);
}


ExpPoly &ExpPoly::operator&=(const ExpPoly &Q) {
    if (this->Box.is_unbounded()) {
       (*this)=Q;
       return (*this);
    } 
    bool modified = false;
    if (!(this->Box.is_subset(Q.Box))) { 
      modified=true;
      this->Box &= Q.Box;
    }
    if (this->Box.is_empty()) { this->set_empty();
				return (*this); }
    /* recompute dim_not_flat */
    this->compute_dim_not_flat();
    if (this->csts.size()+Q.csts.size()==0) 
		{ this->minimize(0); return (*this); }
    for (auto &q : Q.csts) { 
        if ((q.first.vect*this->Box).is_subset(q.second.val)) continue;
        POLYOP_RET p = this->csts.and_constraint(this->Box,q.first,
				 q.second.val, q.second.status);
        if (p==POL_EMPTY) { this->set_empty(); return (*this); }
        if (p==POL_CHANGED) modified=true;
    }
    if (!modified) return (*this);
    this->Inc &= Q.Inc;
    if (this->csts.size()>0) { 
       this->minimize(rhs_nrfilter | rhs_rtfilter);
    }
    return (*this);
}

bool ExpPoly::add_cst(const IntervalVector& V, const Interval& I) {
    Interval fcheck = V*this->Box;
    if (fcheck.is_subset(I)) return false;
    if (fcheck.is_disjoint(I)) {
        this->set_empty();
        return true;
    }
    if (!(V*this->Inc).is_subset(I)) this->Inc.set_empty(); /* we'll have to
            compute it */ /* FIXME : try to keep it */
//    std::cerr << "add_cst before " << this->Box << " " << V << " " << I << "\n";
    Interval bounds(I);
    CstrVect cv = traduit_vect(this->Box, V, bounds);
//    std::cerr << "add_cst " << cv.vdim << " " << cv.bdim << " " << cv.vect << " " << bounds << "\n";
    if (fabs(cv.vdim)<threshold) { /* almost a dimension */
       cv.vect[cv.bdim]=0.0; /* to get the "rest" */
       bounds+=cv.vect*this->Box;
       if (Box[cv.bdim].is_subset(bounds)) return false;
       Box[cv.bdim]&=bounds;
       if (Box[cv.bdim].is_empty()) { this->set_empty(); }
       return true;
    } 
    CstrRhs rhs = { .val = bounds , .status = 0 };
    Interval base = cv.vect*this->Box;
    bool issub=true;
    if (base.lb()>=bounds.lb()) {
       assign_filter(rhs.status,rhs_lbr,rhs_lbfilter);
    } else issub=false;
    if (base.ub()<=bounds.ub()) {
       assign_filter(rhs.status,rhs_ubr,rhs_ubfilter);
    } else issub=false;
    if (issub) return false;
    rhs.val &= base;
    POLYOP_RET p = this->csts.and_constraint(this->Box,cv,rhs);
    if (p==POL_EMPTY) { this->set_empty(); return true; }
    return (p==POL_CHANGED);
}

ExpPoly & ExpPoly::operator&=
	(const std::vector<std::pair<IntervalVector, Interval>> &Res) {
    bool modified=false;
    for (auto &q : Res) { 
             if (this->add_cst(q.first, q.second)) {
                if (this->Box.is_empty()) return (*this);
		modified=true;
	     }
    }
    if (!modified) return (*this);
//    std::cerr << "minimize3 " << this->csts.size() << "\n";
    this->minimize(rhs_nrfilter | rhs_rtfilter);
    return (*this);
}

void ExpPoly::unflat(int dm, Interval offset) {
    if (!Box[dm].is_degenerated()) return;
    if (!offset.is_degenerated()) this->dim_not_flat--;
    Box[dm] += offset;
//    std::cerr << "minimize4\n";
    this->minimize(rhs_nrfilter | rhs_rtfilter); /* really ? */
}


void ExpPoly::intersect_paral(const IntervalMatrix &M, const IntervalVector &V) {
    bool modified=false;
    for (int i=0;i<dim;i++) {
           Vector Z(M[i].mid());
           if (this->add_cst(M[i],V[i])) {
                if (this->Box.is_empty()) return;
		modified=true;
	   }
    }
    if (!modified) return;
//    std::cerr << "minimize5\n";
    this->minimize(rhs_nrfilter | rhs_rtfilter);
}

/* union de deux ivbox : construction des contraintes */
static CstrVectMap union_ivbox(const IntervalVector &iv1, const IntervalVector& iv2) {
   CstrVectMap ret;
   IntervalVector ivunion = iv1 | iv2;
   const int dim = iv1.size();
   for (int d1=0;d1<dim;d1++) {
      if (iv1[d1].lb()>iv2[d1].lb()) {
          for (int d2=0;d2<dim;d2++) {
             if (d1==d2) continue;
             if (iv1[d2].lb()<iv2[d2].lb()) {
               Vector vres(dim,0.0);
               vres[d1]=iv1[d2].lb()-iv2[d2].lb(); /* <0 */
               vres[d2]=iv2[d1].lb()-iv1[d1].lb(); /* <0 */
               Interval r1 = vres*iv1;
               Interval r2 = vres*iv2;
               Interval runion = r1 | r2;
               CstrVect cv = traduit_vect(ivunion,vres,runion);
               CstrRhs rhs = { .val = runion , .status = 0 };
 	       ret.and_constraint(ivunion, cv, rhs); /* en théorie, ne peut donner empty */
             } 
             if (iv1[d2].ub()>iv2[d2].ub()) {
               Vector vres(dim,0.0);
               vres[d1]=iv2[d2].ub()-iv1[d2].ub(); /* <0 */
               vres[d2]=iv1[d1].lb()-iv2[d1].lb(); /* >0 */
               Interval r1 = vres*iv1;
               Interval r2 = vres*iv2;
               Interval runion = r1 | r2;
               CstrVect cv = traduit_vect(ivunion,vres,runion);
               CstrRhs rhs = { .val = runion , .status = 0 };
 	       ret.and_constraint(ivunion, cv, rhs); /* en théorie, ne peut donner empty */
             } 
          }
      }
      if (iv1[d1].ub()<iv2[d1].ub()) {
          for (int d2=0;d2<dim;d2++) {
             if (d1==d2) continue;
             if (iv1[d2].lb()<iv2[d2].lb()) {
               Vector vres(dim,0.0);
               vres[d1]=iv2[d2].lb()-iv1[d2].lb(); /* >0 */
               vres[d2]=iv1[d1].ub()-iv2[d1].ub(); /* <0 */
               Interval r1 = vres*iv1;
               Interval r2 = vres*iv2;
               Interval runion = r1 | r2;
               CstrVect cv = traduit_vect(ivunion,vres,runion);
               CstrRhs rhs = { .val = runion , .status = 0 };
 	       ret.and_constraint(ivunion, cv, rhs); /* en théorie, ne peut donn
er empty */
             } 
             if (iv1[d2].ub()>iv2[d2].ub()) {
               Vector vres(dim,0.0);
               vres[d1]=iv1[d2].ub()-iv2[d2].ub(); /* >0 */
               vres[d2]=iv2[d1].ub()-iv1[d1].ub(); /* >0 */
               Interval r1 = vres*iv1;
               Interval r2 = vres*iv2;
               Interval runion = r1 | r2;
               CstrVect cv = traduit_vect(ivunion,vres,runion);
               CstrRhs rhs = { .val = runion , .status = 0 };
 	       ret.and_constraint(ivunion, cv, rhs); /* en théorie, ne peut donn
er empty */
             } 
          }
      }
   }
   return ret;
}
 
ExpPoly &ExpPoly::operator|=(const IntervalVector &iv) {
    if (iv.is_empty()) return (*this);
    if (this->Box.is_subset(iv)) {
        this->csts.clear(); this->Box=iv;
	this->Inc=iv;
        this->BoxStat.clear(); this->BoxStat.resize(dim,rhs_nrfilter);
        this->minimized=0;
        this->compute_dim_not_flat();
        return (*this);
    }
    if (this->Inc.is_empty()) this->Inc=iv;
    else this->Inc = 0.5*(this->Inc+iv);
    bool modifie=false;
    CstrVectMap c1;
    if (!iv.is_subset(this->Box)) {
      c1 = union_ivbox(iv, this->Box);
      this->Box |= iv;
      modifie=true;
    }
    this->compute_dim_not_flat();
//    if (this->csts.size()==0) return (*this);
    CstrVectMap::iterator ct_it = this->csts.begin();
    while (ct_it != this->csts.end()) {
      const CstrVect& C = ct_it->first;
      Interval& Ret= ct_it->second.val;
      Interval resultB = C.vect*iv;
      if (!resultB.is_subset(Ret)) { 
         Ret |= resultB;
         ct_it->second.status.reset();
         modifie=true;
      }
      ct_it++;
    }
    ct_it = c1.begin();
    while (ct_it != c1.end()) {
      const CstrVect& C = ct_it->first;
      const CstrRhs& Ret= ct_it->second; /* note : Ret.status = 0 */
      POLYOP_RET p = this->csts.and_constraint(this->Box,C,Ret);
      if (p==POL_EMPTY) /* should NOT happen */
         { std::cerr << "union gives empty pol ???\n"; 
	   this->set_empty(); return (*this); }
      if (p==POL_CHANGED) modifie=true;
      ct_it++;
    }
    if (modifie) {
//    std::cerr << "minimize6\n";
       this->minimize(rhs_rfilter | rhs_rtfilter); 
    }   
    return (*this);
}

ExpPoly &ExpPoly::operator|=(const ExpPoly &Q) {
    if (Q.Box.is_empty()) return (*this);
    this->minimize(0);
    if (this->Box.is_empty()) {
       this->Box = Q.Box;
       this->csts = Q.csts;
       this->Inc=Q.Inc;
       this->minimized=Q.minimized;
       this->dim_not_flat=Q.dim_not_flat;
       return (*this);
    }
    if (this->Inc.is_empty()) this->Inc=Q.Inc;
    else if (!Q.Inc.is_empty()) this->Inc = 0.5*(this->Inc+Q.Inc);
    CstrVectMap c1 = union_ivbox(Q.Box, this->Box);
#if 0
    if (this->csts.size()+Q.csts.size()==0) {
      this->Box |= Q.Box;
      this->compute_dim_not_flat();
      return (*this);
    }
#endif
#if 1
    int actsize=this->csts.size()+Q.csts.size()+c1.size();
    if (actsize>maxsize) {
       maxsize = actsize;
       std::cout << "new size : " << maxsize << " : " << (*this) << "\n" << Q << "\n";
//    debug_simplify=true;
    }
#endif
//    bool debug=false;
//    if (this->csts.size()+Q.csts.size()>10) {
//    std::cout << "operator| " << this->dim << " " << this->dim_not_flat << " " << this->Box << " " << this->csts.size() << " " << Q.csts.size() << "\n";
//    debug=true;
//    }
//    bool filtre=false;
//    if (this->csts.size()>4) filtre=true;
//    std::cerr << "minimize7\n";
//    Q.minimize(false);
    CstrVectMap nw;
    for (const std::pair<CstrVect,CstrRhs> &ct : this->csts) {
      const CstrVect& C = ct.first;
      const Interval& Ret= ct.second.val;
      const CstrVectMap::const_iterator el = Q.csts.find(C);
      if (el!=Q.csts.end()) {
         nw.and_constraint(this->Box, C,Ret | el->second.val); continue;
      }
      Interval result = simplex_form(dim,Q.Box,Q.BoxStat,Q.csts,C.vect);
      POLYOP_RET p = nw.and_constraint(this->Box,C,Ret | result);
      if (p==POL_EMPTY) /* should NOT happen */
         { std::cerr << "union gives empty pol ???\n"; 
	   this->set_empty(); return (*this); }
    }
//    if (filtre) {
    for (const std::pair<CstrVect,CstrRhs> &ct : Q.csts) {
      const CstrVect& C = ct.first;
      const Interval& Ret= ct.second.val;
      const CstrVectMap::const_iterator el = this->csts.find(C);
      if (el!=this->csts.end()) continue;
      Interval result = simplex_form(dim,this->Box,this->BoxStat,this->csts,C.vect);
      POLYOP_RET p = nw.and_constraint(this->Box,C,Ret | result);
      if (p==POL_EMPTY) /* should NOT happen */
         { std::cerr << "union gives empty pol ???\n"; 
	   this->set_empty(); return (*this); }
    }
    this->Box |= Q.Box;
    for (const std::pair<CstrVect,CstrRhs> &ct : c1) {
      const CstrVect& C = ct.first;
      const CstrRhs& Ret= ct.second;
      POLYOP_RET p = nw.and_constraint(this->Box,C,Ret);
      if (p==POL_EMPTY) /* should NOT happen */
         { std::cerr << "union gives empty pol ???\n"; 
	   this->set_empty(); return (*this); }
    }
    this->compute_dim_not_flat();
    this->csts=nw;
#if 0
    if (debug_simplify) std::cout << "before min " << (*this) << "\n";
#endif
//    std::cerr << "minimize8\n";
    this->minimize(rhs_rfilter | rhs_rtfilter | rhs_nrfilter); 
      /* FIXME : nrfilter is due to inadequate construction */
#if 0
    debug_simplify=false;
#endif
//    if (this->csts.size()>200) {
//      std::cout << "trop gros, on arrête\n";
//      std::cout << (*this) << "\n";
//      assert(false);
//    }
//    if (debug) std::cout << "après minimize : " << this->csts.size() << "\n";
    return (*this);
}

/* constraint-based => hard to produce a good widening... */
ExpPoly &ExpPoly::widen(const ExpPoly &Q) {
    std::cerr << "widen ?\n";
    this->minimize(0);
    if (Q.Box.is_empty()) return (*this);
    if (this->Box.is_empty()) {
       this->Box = Q.Box;
       this->csts = Q.csts;
       this->minimized=Q.minimized;
       this->dim_not_flat=Q.dim_not_flat;
       return (*this);
    }
//    Q.minimize(false);
    IntervalVector join = this->Box | Q.Box;
    for (int i=0;i<dim;i++) {
       double a = join[i].rad();
       if (join[i].lb()!=this->Box[i].lb()) join[i] |= (join[i].lb()-a);
       if (join[i].ub()!=this->Box[i].ub()) join[i] |= (join[i].ub()+a);
    }
    this->Box = join;
    this->compute_dim_not_flat();
    if (this->dim_not_flat<=1) {
       this->csts.clear();
       this->minimized=true;
       return (*this);
    }
    CstrVectMap::iterator ct_it = this->csts.begin();
    while (ct_it != this->csts.end()) {
      const CstrVect& C = ct_it->first;
      Interval& Ret= ct_it->second.val;
      Interval join = simplex_form(dim,Q.Box,Q.BoxStat,Q.csts,C.vect);
      join |= Ret;
      double a = join.rad();
      if (join.lb()!=Ret.lb()) join|=(join.lb()-a);
      if (join.ub()!=Ret.ub()) join|=(join.ub()+a);
      ct_it->second.val = join;
      ct_it++;
    }
    this->minimize(rhs_rfilter | rhs_rtfilter | rhs_nrfilter);
    return (*this);
}

ExpPoly operator&(const ExpPoly &C1, const ExpPoly &C2) {
   ExpPoly A(C1);
   A &= C2;
   return A;
}
ExpPoly operator|(const ExpPoly &C1, const ExpPoly &C2) {
   ExpPoly A(C1);
   A |= C2;
   return A;
}

std::vector<std::pair<IntervalVector,Interval>> 
		ExpPoly::build_constraints_for_propag
			(const Vector &C, const Vector &Z1) const {
   std::vector<std::pair<IntervalVector,Interval>> result;
   int flat_dim=-1;
   for (int i=0;i<dim;i++) {
       if (Box[i].is_degenerated() && Z1[i]!=0.0) { flat_dim=i; break; }
   }
   if (flat_dim==-1)  {
        for (const std::pair<CstrVect,CstrRhs>&q: this->csts) {
           result.push_back(std::pair<IntervalVector,Interval>
		(IntervalVector(q.first.vect),q.second.val));
        }
        for (int i=0;i<dim;i++) {
           IntervalVector V(dim,0.0);
           V[i]=1.0;
           result.push_back(std::pair<IntervalVector,Interval>(V,Box[i]));
        }
	return result;
   }
   result.reserve(dim+this->csts.size()+1);
   Interval resform = simplex_form(dim,this->Box,this->BoxStat,this->csts,Z1);
   IntervalVector Dir(Z1);
   Interval speed = Dir.norm2();
   resform /= speed;
   Dir *= (1.0/speed);
   result.push_back(std::pair<IntervalVector,Interval>(Dir,resform));
   for (int i=0;i<dim;i++) {
           if (i==flat_dim) continue;
           IntervalVector V(dim,0.0);
           V[i]=Dir[flat_dim];
           V[flat_dim]=-Dir[i];
           resform = Box[i]*V[i]+V[flat_dim]*Box[flat_dim];
           result.push_back(std::pair<IntervalVector,Interval>(V,resform));
   }
   for (const std::pair<CstrVect,CstrRhs>&q: this->csts) {
       IntervalVector V(q.first.vect);
       Interval Coef = V*Dir;
       Coef /= Dir[flat_dim];
       V[flat_dim] -= Coef;
       resform = q.second.val-Coef*Box[flat_dim];
       result.push_back(std::pair<IntervalVector,Interval>(V,resform));
   }
   return result;
}

IntervalVector ExpPoly::build_constraints_for_propag
			(const IntervalMatrix &Z) const {
   IntervalVector result(dim,0.0);
   for (int i=0;i<dim;i++) {
         Interval resform = simplex_form(dim,this->Box,this->BoxStat,this->csts,Z[i].mid());
         resform += (Z[i] - Z[i].mid())*Box;
         result[i] = resform;
   }
   return result;
}


/* returns [C1 - (C1 \cap IV)] (plus ou moins). */
void diff_hull_box(ExpPoly &C1, const IntervalVector &IV) {
   IntervalVector IV2 = IV & C1.Box;
   if (IV2.is_empty()) return;
   ExpPoly res(C1.get_dim(), true);
   for (int i=0;i<C1.dim;i++) {
       if (IV[i].lb()>C1.Box[i].lb()) {
          ExpPoly A(C1);
          A.Box[i]=Interval(C1.Box[i].lb(),IV[i].lb());
//	  std::cerr << "minimize11\n";
          A.minimize(0);
          if (!A.is_empty()) res|=A;
       }
       if (IV[i].ub()<C1.Box[i].ub()) {
          ExpPoly A(C1);
          A.Box[i]=Interval(IV[i].ub(),C1.Box[i].ub());
//	  std::cerr << "minimize12\n";
          A.minimize(0);
          if (!A.is_empty()) res|=A;
       }
   }
   C1 &= res;
}

IntervalVector& operator&=(IntervalVector &V, const ExpPoly &C) {
   ExpPoly A(C);
   A &= V;
   V &= A.Box;
   return V;
}
IntervalVector& operator|=(IntervalVector &V, const ExpPoly &C) {
   V |= C.Box;
   return V;
}
IntervalVector operator+(const ExpPoly &C, const IntervalVector &V) {
   return C.Box + V;
}

std::ostream& operator<< (std::ostream &str, const ExpPoly& C) {
   if (C.is_empty()) { str << "EmptyPoly "; return str; }
   str << "Poly(" << C.Box ;
   str << "{" ;
   for (int i=0;i<C.dim;i++) str << C.BoxStat[i] << ",";
   str << "}"; 
   for (const std::pair<CstrVect,CstrRhs>&c : C.csts) {
       str << "; " << c.first.vect << ":" << c.second.val << "{" << c.second.status << "}" ;
   }
   str << "). ";
   return str;
}

void ExpPoly::vertices2D(std::vector<double>&X, std::vector<double>&Y) {
   assert(dim==2);
//   std::cerr << "vertices2D " << *this << "\n";
   X.clear(); Y.clear();
   /* on retrie les contraintes selon un angle */
   std::map<double,double> rwcsts;
   /* d'abord les bornes */
   rwcsts.insert(std::pair<double,double>(atan2(0.0,1.0),Box[0].ub()));
   rwcsts.insert(std::pair<double,double>(atan2(1.0,0.0),Box[1].ub()));
   rwcsts.insert(std::pair<double,double>(atan2(0.0,-1.0),-Box[0].lb()));
   rwcsts.insert(std::pair<double,double>(atan2(-1.0,0.0),-Box[1].lb()));
   for (const std::pair<CstrVect,CstrRhs>&c : this->csts) {
      double ang = atan2(c.first.vect[1],c.first.vect[0]);
      double nrm = sqrt(c.first.vect[1]*c.first.vect[1]+
			c.first.vect[0]*c.first.vect[0]);
      rwcsts.insert(std::pair<double,double>(ang,c.second.val.ub()/nrm));
      ang=atan2(-c.first.vect[1],-c.first.vect[0]);
      rwcsts.insert(std::pair<double,double>(ang,-c.second.val.lb()/nrm));
   }
   std::map<double,double>::iterator it1=rwcsts.begin();
   std::map<double,double>::iterator it2=it1; it2++;
   double px, py;
   bool ok=true;
   while (ok) {
     double px = (it1->second*sin(it2->first)-it2->second*sin(it1->first))/
                    sin(it2->first-it1->first);
     double py = (it1->second*cos(it2->first)-it2->second*cos(it1->first))/
                    sin(it1->first-it2->first);
     X.push_back(px);
     Y.push_back(py);
//    std::cerr << "point : " << px << "," << py << "\n";
     ok=false;
     while (it2!=rwcsts.begin() && !ok) {
        it1++; it2++; if (it2==rwcsts.end()) it2=rwcsts.begin();
        if (cos(it2->first)*px+sin(it2->first)*py<it2->second) ok=true;
     }
   }
}

/** generate the facet of a 3D flat polyhedron */
std::list<Vector> ExpPoly::facet3D() const {
     if (!this->is_empty())
       for (int i=0;i<dim;i++) {
           if (this->Box[i].is_degenerated())
              return generate_facet(this->getBox(),this->BoxStat,this->getCsts(),i,true);
       }
     /* did not work ? */
     return std::list<Vector>();
}

std::list<std::list<Vector>> ExpPoly::getFacets3D(const IntervalVector &pos, bool with_doors) const {
     std::cerr << "getFacets3D " << (*this) << "\n";
     std::list<std::list<Vector>> res;
     std::list<Vector> to_add;
//     this->minimize(false);
     if (!this->is_empty()) {
       if (with_doors) {
         for (int i=0;i<dim;i++) {
             if (!this->BoxStat[i][UB_R]) {
               to_add = generate_facet(this->Box,this->BoxStat,this->csts,i,true);
#if 0
       std::cerr << "Pol(" << to_add.size() << ") : ";
       for (auto &point:to_add) {
           std::cerr << point << "," ;
       } 
       std::cerr << "\n";
#endif
                  if (to_add.size()>2)  res.push_back(to_add);
               }
               if (!this->BoxStat[i][LB_R]) {
                 to_add = generate_facet(this->Box,this->BoxStat,this->csts,i,false);
#if 0
       std::cerr << "Pol(" << to_add.size() << ") : ";
       for (auto &point:to_add) {
          std::cerr << point << "," ;
       }
       std::cerr << "\n";
#endif
                  if (to_add.size()>2) res.push_back(to_add);
              }
         }
       }
       for (auto &cst : this->csts) {
           if (!cst.second.status[UB_R]) {
             to_add = generate_facet(this->Box,this->BoxStat,this->csts,cst.first,true);
#if 0
     std::cerr << "Pol(" << to_add.size() << ") : ";
     for (auto &point:to_add) {
        std::cerr << point << "," ;
     }
     std::cerr << "\n";
#endif
             if (to_add.size()>2)  res.push_back(to_add);
          }
          if (!cst.second.status[LB_R]) {
            to_add = generate_facet(this->Box,this->BoxStat,this->csts,cst.first,false);
#if 0
     std::cerr << "Pol(" << to_add.size() << ") : ";
     for (auto &point:to_add) {
        std::cerr << point << "," ;
     }
     std::cerr << "\n";
#endif
            if (to_add.size()>2) res.push_back(to_add);
          }
     }
   }
   std::cerr << "fin getFacets3D, faces générées : " << res.size() << "\n";
   return res;
}

}

#define TEST_POLY	0
#if (TEST_POLY)
using namespace invariant;

/* test de simplex_form */
int main() {
   IntervalVector ivbox(3);
   ivbox[0]=Interval(0.625, 0.7660014045188112);
   ivbox[1]=Interval(-0.15625, 0);
   ivbox[2]=Interval(4.84375, 5);

 
   std::vector<std::pair<IntervalVector,Interval>> csts;
   IntervalVector cst(3); cst[0]=1; cst[1]=-3.49959e-17;  cst[2]=-0.094554;
   csts.push_back(std::pair<IntervalVector,Interval>(cst,Interval(0.1522296908946226, 0.2932310954134339)));
   cst[0]=1; cst[1]=0;  cst[2]=-0.0436611;
   csts.push_back(std::pair<IntervalVector,Interval>(cst,Interval(0.406694071125451,0.5476954756442622)));
   cst[0]=1; cst[1]=0.295813;  cst[2]=-0;
   csts.push_back(std::pair<IntervalVector,Interval>(cst,Interval(0.5787792319480085, 0.7660014045188111)));

   ExpPoly ep(ivbox,csts,false);
 
   ep.unflat(1,Interval(0,0.4));
   std::cout << ep << "\n";
   std::list<std::list<Vector>> facets = ep.getFacets3D(ivbox);
#if 0
   for (auto &facet:facets) {
     std::cout << "Pol(" << facet.size() << ") : ";
     for (auto &point:facet) {
        std::cout << point << "," ;
     }
     std::cout << "\n";
   }
#endif


#if 0
   std::list<Vector> facet;
   for (int j=0;j<3;j++) {
     facet = generate_facet(ep.getBox(),ep.getCsts(),j,true);
     std::cout << "facet : " << "\n";
     for (int i=0;i<facet.size();i++) {
        std::cout << "[" << facet[i][0] << "," << facet[i][1] << "," << facet[i][2] << "]," ;
     }
     std::cout << "\n";
     facet = generate_facet(ep.getBox(),ep.getCsts(),j,false);
     std::cout << "facet : " << "\n";
     for (int i=0;i<facet.size();i++) {
        std::cout << "[" << facet[i][0] << "," << facet[i][1] << "," << facet[i][2] << "]," ;
     }
     std::cout << "\n";
   }
   for (auto &cst : ep.getCsts()) {
     facet = generate_facet(ep.getBox(),ep.getCsts(),cst.first,true);
     std::cout << "facet : " << "\n";
     for (int i=0;i<facet.size();i++) {
        std::cout << "[" << facet[i][0] << "," << facet[i][1] << "," << facet[i][2] << "]," ;
     }
     std::cout << "\n";
     facet = generate_facet(ep.getBox(),ep.getCsts(),cst.first,false);
     std::cout << "facet : " << "\n";
     for (int i=0;i<facet.size();i++) {
        std::cout << "[" << facet[i][0] << "," << facet[i][1] << "," << facet[i][2] << "]," ;
     }
     std::cout << "\n";
   }
#endif
   return 0;

}

#endif
