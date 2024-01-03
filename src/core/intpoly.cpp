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
#include "intsimplex.h"
#include "intcstrvect.h"
#include "intpoly.h"

namespace intflpl {

static unsigned int maxsize=0;
static constexpr double threshold=1e-5;

IntPoly::IntPoly() :
   dim(-1), dim_not_flat(-1),
   Box(), BoxStat(), Inc(), minimized(-1),
   csts(), pstatus(1 << INTPOLY_INVALID)
{  }

IntPoly::IntPoly(int dim) :
     IntPoly(dim,false)
{  }

IntPoly::IntPoly(int dim, bool empty) :
  dim(dim), dim_not_flat(empty ? -1 : dim), 
  Box(dim, (empty ? Interval::empty_set() : Interval::all_reals())),
  BoxStat(dim, 0), Inc(Box),
  minimized(0), csts(), 
  pstatus(empty ? (1 << INTPOLY_EMPTY) : (1 << INTPOLY_UNBOUNDED))
{  }

IntPoly::IntPoly(const IntervalVector &Box) : 
  dim(Box.size()), dim_not_flat(Box.is_empty() ? -1 : 0),
  Box(Box), BoxStat(dim, (1<<LB_NR) | (1<<UB_NR)),
  Inc(Box), minimized(0), csts(),
  pstatus(Box.is_empty() ? (1 << INTPOLY_EMPTY)
	   : (Box.is_unbounded() ? (1 << INTPOLY_UNBOUNDED) : 0))
{
  if (dim_not_flat==0) 
    for (int i=0;i<dim;i++) 
         if (!this->Box[i].is_degenerated()) dim_not_flat++;
}

IntPoly::IntPoly(const IntervalVector &Box, const std::vector<std::pair<IntervalVector,Interval>> &Csts) :
     IntPoly(Box) 
{
   this->Inc.set_empty();
   dim_not_flat = simplify_polyhedron(dim,this->Box,this->BoxStat,
			Csts,this->csts, this->Inc);
   if (dim_not_flat==-1) {
         pstatus = (1 << INTPOLY_EMPTY); /* FIXME : unbounded case */
   }
}

IntPoly::IntPoly (const IntPoly &P) :
     dim(P.dim), dim_not_flat(P.dim_not_flat), Box(P.Box),
     BoxStat(P.BoxStat),
     Inc(P.Inc),
     minimized(P.minimized), csts(P.csts), pstatus(P.pstatus)  {
}

IntPoly::IntPoly(const IntervalMatrix& M, const IntervalMatrix &rM,
                        const IntervalVector &V) :
      IntPoly(M*V) 
{
    if (pstatus[INTPOLY_UNBOUNDED] || pstatus[INTPOLY_EMPTY]) return;
    bool modif = false;
    double radmin=0.0;
    for (int i=0;i<dim;i++) {
        modif |= this->add_cst(rM[i],V[i]);
        if (i==0 || V[i].rad()<radmin) radmin=V[i].rad();
    } 
    this->Inc=Box.mid();
    this->Inc.inflate(radmin/infinite_norm(rM));
    if (!modif) return;
    this->minimize(rhs_nrfilter | rhs_rtfilter);
}

IntPoly::IntPoly(const IntervalMatrix& M, 
                        const IntervalVector &V) :
      IntPoly(M*V) 
{
    IntLU ilu = IntLU(M,false,true);
    if (!ilu.inLUform() || ilu.getDeterminant().contains(0.0)) return;
 			/* FIXME ? M n'est pas inversible, mais 
			   on peut peut-être trouver des contraintes ? */
    IntervalMatrix rM = ilu.getInvB();
    if (pstatus[INTPOLY_UNBOUNDED] || pstatus[INTPOLY_EMPTY]) return;
    bool modif = false;
    double radmin=0.0;
    for (int i=0;i<dim;i++) {
        modif |= this->add_cst(rM[i],V[i]);
        if (i==0 || V[i].rad()<radmin) radmin=V[i].rad();
    } 
    this->Inc=Box.mid();
    this->Inc.inflate(radmin/infinite_norm(rM));
    if (!modif) return;
    this->minimize(rhs_nrfilter | rhs_rtfilter);
}


Interval IntPoly::bound_cstFast(const CstrVect &cvect) const {
    Interval ret = cvect.vect*this->Box;
    CstrVectMap::const_iterator el = this->csts.lower_bound(cvect);
    if (el==this->csts.end() || el->first.vdim!=cvect.vdim) {
         if (el==this->csts.begin()) return ret;
         el--;
         if (el->first.vdim!=cvect.vdim) return ret;
    }
//    std::cout << "bound cst: " << el->first.vect << " " << cvect.vect << " " 
// << (el->first==cvect) << "\n";
    if (el->first==cvect) { ret &= el->second.val; return ret; }
    ret &= el->second.val + (cvect.vect-el->first.vect)*this->Box;
    return ret;
}

bool IntPoly::satisfy_cst (const std::pair<const CstrVect,CstrRhs> &cst)
const {
    if (pstatus[INTPOLY_EMPTY]) return true;
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED]) return false;
    /* basic checks first */
    Interval fastbound = this->bound_cstFast(cst.first);
    if ((cst.second.status[LB_R] || fastbound.lb()>=cst.second.val.lb())
	&& (cst.second.status[UB_R] || fastbound.ub()<=cst.second.val.ub()))
       return true;
    Interval result = simplex_form(dim,this->Box,this->BoxStat,
			this->csts,cst.first.vect); 
    /* FIXME : on peut peut-être éviter de faire les deux bornes ? */
    return result.is_subset(cst.second.val);
}

bool IntPoly::contains(const Vector& iv) const {
    if (pstatus[INTPOLY_EMPTY] || pstatus[INTPOLY_INVALID]) return false;
    if (pstatus[INTPOLY_UNBOUNDED]) return true;
    if (!Box.contains(iv)) return false;
    if (Inc.contains(iv)) return true;
    for (const std::pair<CstrVect,CstrRhs>& cst : csts) { 
        if (!cst_satisfy_vec(cst,iv)) return false;
    }
    return true;
}

bool IntPoly::intersects(const IntervalVector& iv) const {
    if (pstatus[INTPOLY_EMPTY] || pstatus[INTPOLY_INVALID]) return false;
    if (pstatus[INTPOLY_UNBOUNDED]) return true;
    if (this->Inc.intersects(iv)) return true;
    if (!Box.intersects(iv)) return false;
    if (this->is_box()) return true;
    IntPoly A = (*this) & iv; /* FIXME : not efficient */
    return !A.is_empty();
}

bool IntPoly::intersects(const IntPoly& Q) const {
    if (pstatus[INTPOLY_EMPTY] || pstatus[INTPOLY_INVALID]) return false;
    if (Q.pstatus[INTPOLY_EMPTY] || Q.pstatus[INTPOLY_INVALID]) return false;
    if (pstatus[INTPOLY_UNBOUNDED] || Q.pstatus[INTPOLY_UNBOUNDED]) return true;
    if (this->Inc.intersects(Q.Inc)) return true;
    if (!Box.intersects(Q.Box)) return false;
    if (this->is_box() && Q.is_box()) return true;
    IntPoly A = (*this) & Q; /* FIXME : not efficient */
    return !A.is_empty();
}

double IntPoly::rel_distanceFast(const IntPoly& Q) const {
    if (Q.pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_INVALID]) return false;
    if (pstatus[INTPOLY_EMPTY] || Q.pstatus[INTPOLY_UNBOUNDED]) 
			return POS_INFINITY;
    if (pstatus[INTPOLY_UNBOUNDED]) return 0.0;
    double dist= this->Box.rel_distance(Q.Box);
//    std::cout << "rel_distance box : " << dist << "\n";
    for (const std::pair<CstrVect,CstrRhs>& cst : this->csts) { 
        Interval b = Q.bound_cstFast(cst.first);
        double d2 = cst.second.val.rel_distance(b);
        if (d2>dist) dist=d2;
    }
    return dist;
}

void IntPoly::clear() {
    assert(!pstatus[INTPOLY_INVALID]);
    pstatus=0;
    this->Box.clear();
    this->csts.clear();
    this->dim_not_flat=0;
    std::vector<cstrrhs_status> bStat(dim,(1 << LB_NR)|(1 << UB_NR));
    this->BoxStat.swap(bStat);
    Inc=Box;
    minimized=0;
}

void IntPoly::setComponent(int i, const Interval &a) {
    assert(!pstatus[INTPOLY_INVALID]);
    if (pstatus[INTPOLY_EMPTY]) return;
    if (a.is_empty()) { this->set_empty(); return; }
    Interval old = Box[i];
    this->Box[i]=a;
    this->Inc[i]=a;
    if (a.is_unbounded()) { pstatus[INTPOLY_UNBOUNDED]=true; }
    else if (!this->Box.is_unbounded()) {
            pstatus=0;
         }
    this->compute_dim_not_flat();
    if (old.is_unbounded()) return;
    CstrVectMap ocsts;
    ocsts.swap(csts);
    bool modif = false;
    for (const std::pair<CstrVect,CstrRhs>& cst : ocsts) { 
        Vector V(cst.first.vect);
        Interval I(cst.second.val);
        I -= V[i]*a;
        V[i]=0.0;
        if (infinite_norm(V)==0.0) continue;
        modif |= this->add_cst(V,I);
    }
    if (!modif) return;
    this->minimize(rhs_nrfilter | rhs_rtfilter);
}

IntPoly& IntPoly::inflate(double rad) {
    if (pstatus[INTPOLY_EMPTY] || pstatus[INTPOLY_INVALID]) return (*this);
    if (rad<=0.0) { return (*this); }
    this->Box.inflate(rad);
    std::vector<cstrrhs_status> bStat(dim,(1 << LB_NR)|(1 << UB_NR));
    this->BoxStat.swap(bStat);
    this->Inc.inflate(rad);
    for (std::pair<const CstrVect,CstrRhs>& cst : csts) { 
        const Vector &V(cst.first.vect);
        Interval &I(cst.second.val);
        double rV = 0.0;   /* FIXME : arrondi du double */
        for (int i=0;i<dim;i++) rV+=fabs(V[i])*rad;
        I.inflate(rV);
    }
    return (*this);
}

IntPoly& IntPoly::homothety(IntervalVector c, double delta) {
    if (pstatus[INTPOLY_EMPTY] || pstatus[INTPOLY_INVALID]) return (*this);
    if (pstatus[INTPOLY_UNBOUNDED]) return (*this);
    this->Box = (1-delta)*c+delta*this->Box;
    this->Inc = (1-delta)*c+delta*this->Inc;
    CstrVectMap ocsts;
    ocsts.swap(csts);
    bool modif = false;
    for (std::pair<const CstrVect,CstrRhs>& cst : ocsts) { 
        Interval I = delta*cst.second.val+
		(1-delta)*(cst.first.vect*c);
        POLYOP_RET p = this->csts.and_constraint(this->Box,cst.first,
                                 I, cst.second.status);
       if (p==POL_EMPTY) { this->set_empty(); return (*this); }
       if (p==POL_CHANGED) modif=true;
    }
    if (modif) this->minimize(rhs_rtfilter); /* FIXME  :
				check if c is punctual => 0 */
    return (*this);
}

IntPoly& IntPoly::operator=(const IntervalVector& x) {
    this->Box = x;
    this->dim=x.size();
    this->csts.clear();
    this->Inc=x;
    this->pstatus=(x.is_unbounded() ? (1<<INTPOLY_UNBOUNDED) : 0);
    if (x.is_empty()) { this->set_empty(); return (*this); }
    std::vector<cstrrhs_status> bStat(dim,(1 << LB_NR)|(1 << UB_NR));
    this->BoxStat.swap(bStat);
    this->compute_dim_not_flat();
    this->minimized=0;
    return (*this);
}

static void inflate_interval_from_base
		(Interval &a, const Interval &b, double fact) {
     double mn = fact*(b.lb() - a.lb());
     if (mn>0) mn=0.0;
     double mx = fact*(b.ub() - a.ub());
     if (mx<0) mx=0.0;
     Interval ev(mn,mx);
     a += ev;
}

void IntPoly::inflate_from_baseFast(const IntPoly &iv, double fact) {
    for (int i=0;i<dim;i++) {
	inflate_interval_from_base(this->Box[i],iv.Box[i],fact);
    }
    for (std::pair<const CstrVect,CstrRhs>& cst : csts) { 
        Interval b = iv.bound_cstFast(cst.first);
        Interval &I(cst.second.val);
	inflate_interval_from_base(I,b,fact);
    }
//    std::cout << "inflate_before_minimize " << (*this) << "\n";
    this->minimize(rhs_nrfilter | rhs_rtfilter); /* FIXME : mieux faire ? */
//    std::cout << "inflate_after_minimize " << (*this) << "\n";
}

IntPoly& IntPoly::linMult(const IntervalMatrix& M, const IntervalMatrix& IM) {
    if (pstatus[INTPOLY_EMPTY] || pstatus[INTPOLY_INVALID]) return (*this);
    if (pstatus[INTPOLY_UNBOUNDED]) return (*this);
    IntervalVector oldBox=this->Box;
    this->Box = M*this->Box;
    CstrVectMap ocsts;
    csts.swap(ocsts);
    for (int i=0;i<dim;i++) {
       const IntervalVector &V = IM[i];
       this->add_cst(V,oldBox[i]);
    }
    for (const std::pair<CstrVect,CstrRhs>& cst : ocsts) { 
        IntervalVector V(cst.first.vect*IM);
        const Interval &I(cst.second.val);
        this->add_cst(V,I);
    }
//    std::cout << "après linMult avant minimize : " << (*this) << "\n";
    this->minimize(rhs_nrfilter | rhs_rtfilter); /* FIXME : mieux faire ? */
    return (*this);
}


bool IntPoly::is_subset (const IntPoly &Q) const {
    /* A inclus dans Q si toutes les contraintes de Q sont inutiles */
    /* peut-être pas correct si (*this) n'est pas minimisé */
//    this->minimize(0);
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED]) return false;
    if (Q.pstatus[INTPOLY_INVALID]) return false;
    if (pstatus[INTPOLY_EMPTY]) return true;
    if (!this->Box.is_subset(Q.Box)) return false; /* revoir ! */
    for (const std::pair<CstrVect,CstrRhs> &ct : Q.csts) {
        if (!this->satisfy_cst(ct)) {
           return false;
        }
    }
    return true;
}

bool IntPoly::is_subsetFast (const  IntPoly &Q) const {
    return this->is_subset(Q);
}

bool IntPoly::is_subset (const IntervalVector &IV) const {
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED]) return false;
    return (this->Box.is_subset(IV));
}


bool operator==(const IntPoly &C1, const IntPoly &C2) {
    if (C1.pstatus[IntPoly::INTPOLY_INVALID] 
     || C2.pstatus[IntPoly::INTPOLY_INVALID])
	 return false;
    if (C1.pstatus[IntPoly::INTPOLY_UNBOUNDED]) {
        return C2.pstatus[IntPoly::INTPOLY_UNBOUNDED];
    }
    if (C2.pstatus[IntPoly::INTPOLY_UNBOUNDED]) return false;
    if (C1.Box != C2.Box) return false;
//    C1.minimize(0);
//    C2.minimize(0);
    return (C1.is_subset(C2) && C2.is_subset(C1));
}

bool IntPoly::is_superset (const IntervalVector &IV) const {
    /* A inclus dans Q si toutes les contraintes de Q sont inutiles */
    if (pstatus[INTPOLY_INVALID]) return false;
    if (pstatus[INTPOLY_UNBOUNDED]) return true;
    if (!IV.is_subset(this->Box)) return false;
    for (const std::pair<CstrVect,CstrRhs> &ct : this->csts) {
        const CstrVect &C = ct.first;
        const Interval &Ret= ct.second.val; 
        if (!(C.vect*IV).is_subset(Ret)) return false;
    }
    return true;
}


void IntPoly::compute_dim_not_flat() {
    if (pstatus[INTPOLY_INVALID]) return;
    if (!this->Box.is_empty()) {
      this->dim_not_flat=0;
      for (int i=0;i<dim;i++) {
         if (!this->Box[i].is_degenerated()) this->dim_not_flat++;
      }
    } else this->dim_not_flat=-1;
}

void IntPoly::minimize(cstrrhs_status ct) {
    minimized |= ct;
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED]) return;
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

IntPoly &IntPoly::operator&=(const IntervalVector &iv) {
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED]) {
       this->Box = iv;
       this->dim = iv.size();
       if (iv.is_empty()) pstatus=(1 << INTPOLY_EMPTY);
       else if (iv.is_unbounded()) pstatus=(1 << INTPOLY_UNBOUNDED);
       else pstatus=0;
       this->compute_dim_not_flat();
       return (*this);
    }
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

IntPoly &IntPoly::meetFast(const IntPoly &Q) {
    if (Q.pstatus[INTPOLY_INVALID]) return (*this);
    if (pstatus[INTPOLY_INVALID]) { (*this)=Q;  return *(this); }
    if (pstatus[INTPOLY_EMPTY] || Q.pstatus[INTPOLY_UNBOUNDED]) {
      return (*this);
    }
    bool modified = false;
    this->Box &= Q.Box;
    if (this->Box.is_empty()) { this->set_empty();
				return (*this); }
    for (std::pair<const CstrVect,CstrRhs>& cst : this->csts) { 
        Interval fastbound = Q.bound_cstFast(cst.first);
        if (!cst.second.val.is_subset(fastbound)) {
           cst.second.val&=fastbound;
           modified=true;
        }
    }
    if (!modified) return (*this);
    this->Inc &= Q.Inc;
    if (this->csts.size()>0) { 
       this->minimize(rhs_nrfilter | rhs_rtfilter);
    }
    return (*this);
}

IntPoly &IntPoly::meetKeep(const IntPoly &Q) {
    return this->meetFast(Q);
}

IntPoly &IntPoly::meet(const IntPoly& iv, bool ctcG) {
    if (ctcG) return (this->meetKeep(iv));
    return ((*this) &= iv);
}


IntPoly &IntPoly::operator&=(const IntPoly &Q) {
    if (Q.pstatus[INTPOLY_INVALID]) return (*this);
    if (pstatus[INTPOLY_INVALID]) { (*this)=Q;  return *(this); }
    if (pstatus[INTPOLY_EMPTY] || Q.pstatus[INTPOLY_UNBOUNDED]) {
      return (*this);
    }
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

bool IntPoly::meetLN(IntervalVector &V, Interval &b, bool ctc) {
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED]) return true;
    if (pstatus[INTPOLY_EMPTY]) return false;
    bool modif = this->add_cst(V,b);
    if (modif) 
        this->minimize(rhs_nrfilter | rhs_rtfilter);
    if (this->is_empty()) {
        V.set_empty(); b.set_empty();
        return false;
    }
    b &= V*this->Box;
    return true;
}

bool IntPoly::meetLN(IntervalVector &V, IntervalVector &C,
			Interval &b, bool ctc) {
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED]) return true;
    if (pstatus[INTPOLY_EMPTY]) return false;
    Interval bBis = b+V*C;
    bool modif = this->add_cst(V,bBis);
    if (modif) 
 
        this->minimize(rhs_nrfilter | rhs_rtfilter);
    if (this->is_empty()) {
        V.set_empty(); b.set_empty();
        return false;
    }
    b &= V*(this->Box-C);
    return true;
}


bool IntPoly::meetLM(IntervalMatrix &S, IntervalVector &b, bool ctc) {
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED]) return true;
    if (pstatus[INTPOLY_EMPTY]) return false;
    bool modif=false;
    for (int i=0;i<S.nb_rows();i++) {
        modif |= this->add_cst(S[i],b[i]);
    }
    if (modif) 
        this->minimize(rhs_nrfilter | rhs_rtfilter);
    if (this->is_empty()) {
        S.set_empty(); b.set_empty();
        return false;
    }
    b &= S*this->Box;
    return true;
}

bool IntPoly::meetLM(IntervalMatrix &S, const Vector& C,
			IntervalVector &b, bool ctc) {
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED]) return true;
    if (pstatus[INTPOLY_EMPTY]) return false;
    bool modif=false;
    IntervalVector bBis = b+S*C;
    for (int i=0;i<S.nb_rows();i++) {
        modif |= this->add_cst(S[i],bBis[i]);
    }
    if (modif) 
        this->minimize(rhs_nrfilter | rhs_rtfilter);
    if (this->is_empty()) {
        S.set_empty(); b.set_empty();
        return false;
    }
    b &= S*(this->Box-C);
    return true;
}


bool IntPoly::add_cst(const IntervalVector& V, const Interval& I) {
    Interval fcheck = V*this->Box;
    if (fcheck.is_subset(I)) return false;
    if (fcheck.is_disjoint(I)) {
//        std::cout << "is_disjoint : " << I << " " << fcheck << "\n";
        this->set_empty();
        return true;
    }
    if (!(V*this->Inc).is_subset(I)) this->Inc.set_empty(); /* we'll have to
            compute it */ /* FIXME : try to keep it */
//    std::cout << "add_cst before " << this->Box << " " << V << " " << I << "\n";
    Interval bounds(I);
    CstrVect cv = traduit_vect(this->Box, V, bounds);
//    std::cout << "add_cst " << cv.vdim << " " << cv.bdim << " " << cv.vect << " " << bounds << "\n";
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

IntPoly & IntPoly::operator&=
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

void IntPoly::unflat(int dm, Interval offset) {
    if (!Box[dm].is_degenerated()) return;
    if (!offset.is_degenerated()) this->dim_not_flat--;
    Box[dm] += offset;
//    std::cerr << "minimize4\n";
    this->minimize(rhs_nrfilter | rhs_rtfilter); /* really ? */
}


void IntPoly::intersect_paral(const IntervalMatrix &M, const IntervalVector &V) {
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
 
IntPoly &IntPoly::operator|=(const IntervalVector &iv) {
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

IntPoly &IntPoly::operator|=(const IntPoly &Q) {
    if (Q.Box.is_empty()) return (*this);
    this->minimize(0);
    if (this->Box.is_empty()) {
       this->Box = Q.Box;
       this->csts = Q.csts;
       this->Inc=Q.Inc;
       this->minimized=Q.minimized;
       this->dim_not_flat=Q.dim_not_flat;
       this->pstatus=Q.pstatus;
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
    unsigned int actsize=this->csts.size()+Q.csts.size()+c1.size();
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
IntPoly &IntPoly::widen(const IntPoly &Q) {
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

IntPoly operator&(const IntPoly &C1, const IntPoly &C2) {
   IntPoly A(C1);
   A &= C2;
   return A;
}
IntPoly operator&(const IntPoly &C1, const IntervalVector &IV) {
   IntPoly A(C1);
   A &= IV;
   return A;
}
IntervalVector operator&(const IntervalVector &C1, const IntPoly &C2) {
   if (C2.is_box()) return (C2.Box & C1);
   IntPoly A(C2);
   A &= C1;
   return A.Box;
}
IntPoly operator|(const IntPoly &C1, const IntPoly &C2) {
   IntPoly A(C1);
   A |= C2;
   return A;
}

IntervalVector operator*(const IntervalMatrix& M, const IntPoly& iv) {
   return M*iv.Box; /* TODO : do better ? */
}

IntPoly& IntPoly::operator+=(const IntervalVector& V) {
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED] 
		|| pstatus[INTPOLY_EMPTY]) return (*this);
    this->Box+=V;
    this->Inc+=V;
    if (this->is_box()) return (*this);
    for (std::pair<const CstrVect,CstrRhs>& cst : csts) { 
        Interval &I(cst.second.val);
        I += cst.first.vect*V;
    }
    this->minimize(rhs_rfilter); /* FIXME : mieux faire ? */
    return (*this);
}

IntPoly& IntPoly::operator-=(const IntervalVector& V) {
    return ((*this)+=(-V));
}

IntPoly operator+(const IntPoly& iv, const IntervalVector& V) {
    IntPoly A(iv);
    A+=V; return A;
}

IntPoly operator-(const IntPoly& iv, const IntervalVector& V) {
    IntPoly A(iv);
    A-=V; return A;
}

IntPoly sum_tau(const IntPoly& iv, const IntervalVector& V,
                                                bool keep) {
    IntPoly A(iv);
    if (!keep) { A += V; return (iv | A); }
    Interval tau(0.0,1.0);
    A.Box+=tau*V;
    A.Inc+=V;
    if (A.is_box()) return A;
    for (std::pair<const CstrVect,CstrRhs>& cst : A.csts) { 
        Interval &I(cst.second.val);
        I += (tau)*(cst.first.vect*V);
    }
    A.minimize(rhs_rtfilter); /* FIXME : mieux faire ? */
    return A;
}

IntPoly& IntPoly::sumFast(const IntPoly& iv) {
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED] 
	|| pstatus[INTPOLY_EMPTY] || iv.pstatus[INTPOLY_INVALID])
		 return (*this);
    if (iv.pstatus[INTPOLY_EMPTY]) { this->set_empty(); return (*this); }
    if (iv.pstatus[INTPOLY_UNBOUNDED]) { (*this)=iv; return (*this); }
    this->Box+=iv.Box;
    this->Inc+=iv.Inc;
    if (this->is_box()) return (*this);
    for (std::pair<const CstrVect,CstrRhs>& cst : csts) { 
        Interval &I(cst.second.val);
        Interval fastbound = iv.bound_cstFast(cst.first);
        I += fastbound;
    }
    this->minimize(rhs_rtfilter | rhs_nrfilter); /* FIXME : mieux faire ? */
    return (*this);
}

IntPoly& IntPoly::diffFast(const IntPoly& iv) {
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED] 
	|| pstatus[INTPOLY_EMPTY] || iv.pstatus[INTPOLY_INVALID])
		 return (*this);
    if (iv.pstatus[INTPOLY_EMPTY]) { this->set_empty(); return (*this); }
    if (iv.pstatus[INTPOLY_UNBOUNDED]) { (*this)=iv; return (*this); }
    this->Box-=iv.Box;
    this->Inc-=iv.Inc;
    if (this->is_box()) return (*this);
    for (std::pair<const CstrVect,CstrRhs>& cst : csts) { 
        Interval &I(cst.second.val);
        Interval fastbound = iv.bound_cstFast(cst.first);
        I -= fastbound;
    }
    this->minimize(rhs_rtfilter | rhs_nrfilter); /* FIXME : mieux faire ? */
    return (*this);
}

void IntPoly::cmult_and_add(const Vector& center,
                const IntervalMatrix& M, const IntervalMatrix &invM,
                const IntervalVector& V)
{
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED] 
	|| pstatus[INTPOLY_EMPTY]) return;
//    std::cout << "avant linMult : " << (*this) << "\n";
    this->linMult(M,invM);
//    std::cout << "après linMult : " << (*this) << "\n";
    (*this) += (-M*center + center +V);
}

void IntPoly::ctau_mult_and_add(const Vector& center,
                const IntervalMatrix& M, 
                const IntervalVector& V)
{
    if (pstatus[INTPOLY_INVALID] || pstatus[INTPOLY_UNBOUNDED] 
	|| pstatus[INTPOLY_EMPTY]) return;
    Interval tau(0.0,1.0);
    IntervalVector cBox = Box-center;
//    std::cout << "M*cBox" << M*cBox << "\nV : " << V << "\n";
    this->Box += tau*(M*cBox+V);
    if (this->is_box()) return;
    for (std::pair<const CstrVect,CstrRhs>& cst : csts) { 
        Interval &I(cst.second.val);
        I += tau*((cst.first.vect*M)*(Box-center)+(cst.first.vect)*V);
    }
    this->minimize(rhs_rtfilter | rhs_nrfilter);
}

/* this algorithm is not the most effective one for what we want... 
   => not used in invariant-lib nor codac, is it really useful */
#if 0
bool join_intersect_with_tau
           (const IntPoly& iv, const Vector &center,
            const IntervalMatrix& M,
            const IntervalVector& V,
            const IntervalVector& box, int d, double val) {
       

}
#endif

    
std::vector<std::pair<IntervalVector,Interval>> 
		IntPoly::build_constraints_for_propag
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

IntervalVector IntPoly::build_constraints_for_propag
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
void diff_hull_box(IntPoly &C1, const IntervalVector &IV) {
   IntervalVector IV2 = IV & C1.Box;
   if (IV2.is_empty()) return;
   IntPoly res(C1.get_dim(), true);
   for (int i=0;i<C1.dim;i++) {
       if (IV[i].lb()>C1.Box[i].lb()) {
          IntPoly A(C1);
          A.Box[i]=Interval(C1.Box[i].lb(),IV[i].lb());
//	  std::cerr << "minimize11\n";
          A.minimize(0);
          if (!A.is_empty()) res|=A;
       }
       if (IV[i].ub()<C1.Box[i].ub()) {
          IntPoly A(C1);
          A.Box[i]=Interval(IV[i].ub(),C1.Box[i].ub());
//	  std::cerr << "minimize12\n";
          A.minimize(0);
          if (!A.is_empty()) res|=A;
       }
   }
   C1 &= res;
}

IntervalVector& operator&=(IntervalVector &V, const IntPoly &C) {
   IntPoly A(C);
   A &= V;
   V &= A.Box;
   return V;
}
IntervalVector& operator|=(IntervalVector &V, const IntPoly &C) {
   V |= C.Box;
   return V;
}
IntervalVector operator+(const IntervalVector &V, const IntPoly &C) {
   return C.Box + V;
}

std::ostream& operator<< (std::ostream &str, const IntPoly& C) {
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

void IntPoly::vertices2D(std::vector<double>&X, std::vector<double>&Y) const {
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
std::list<Vector> IntPoly::facet3D() const {
     if (!this->is_empty())
       for (int i=0;i<dim;i++) {
           if (this->Box[i].is_degenerated())
              return generate_facet(this->getBox(),this->BoxStat,this->getCsts(),i,true);
       }
     /* did not work ? */
     return std::list<Vector>();
}

std::list<std::list<Vector>> IntPoly::getFacets3D(const IntervalVector &pos, bool with_doors) const {
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
using namespace intflpl;

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

   IntPoly ep(ivbox,csts);
 
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
