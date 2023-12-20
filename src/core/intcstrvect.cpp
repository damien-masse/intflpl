/**
 * Simplex manipulations for interval matrices
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
#include <bitset>
#include "flpl_def.h"
#include "intsimplex.h"
#include "intcstrvect.h"


namespace intflpl {


/* generate the facet of a polyhedron. with an extremal dimension */
/* TODO : partly merge with next function ? */ 
std::list<Vector> generate_facet(const IntervalVector &ivbox,
		 std::vector<cstrrhs_status> clstats,
                 const CstrVectMap &csts, int dm, bool upper) {
    assert(ivbox.size()==3);
    int sz = csts.size();
    Intsimplex simp(ivbox,sz,clstats);
    for (const std::pair<const CstrVect,CstrRhs> &cst : csts) {
        simp.load_constraint(cst);
    }
    if (!simp.generate_init_basis(dm,upper)) { return std::list<Vector>(); }
    simplex_ret sret = simp.simplex_mat();
    if (sret[INFEASIBLE]) return std::list<Vector>();
    if (simp.get_basis()[dm] != dm) { 
			/* FIXME : peut-être pas une bonne approche? */
               return std::list<Vector>(); }
    simp.lock_row(dm,true);
    return simp.generate_vertices_2D(ivbox);
}

/* generate the facet of a polyhedron. with a existing constraint */
std::list<Vector> 
		generate_facet(const IntervalVector &ivbox,
		 std::vector<cstrrhs_status> clstats,
                 const CstrVectMap &csts, const CstrVect &cstF, bool upper) {
    assert(ivbox.size()==3);
    int dim = ivbox.size();
    int sz = csts.size();
    Intsimplex simp(ivbox,sz,clstats);
    int colCst = -1;
    for (const std::pair<const CstrVect,CstrRhs> &cst : csts) {
        int a = simp.load_constraint(cst);
        if (cst.first==cstF) colCst = a;
    }
    if (!simp.generate_init_basis(colCst,upper)) { return std::list<Vector>(); }
    simplex_ret sret = simp.simplex_mat();
    if (sret[INFEASIBLE]) return std::list<Vector>();
    /* look for place */
    int dm = -1;
    for (int i=0;i<dim;i++) {
       if (simp.get_basis()[i]==colCst) { dm=i; simp.lock_row(i,true); break; } 
    }
    if (dm==-1) { /* pas une facette ? */
               return std::list<Vector>(); }
    return simp.generate_vertices_2D(ivbox);
}

/** get upper bound of a linear form on a polyhedron defined with intervals
   (simplex algorithm). The upper bound is guaranteed to be safe (>= the exact
   bound), but not optimal 
   @param dim dimension of the space
   @param ivbox IntervalVector bounding box (used as a "first" space)
   @param csts std::vector<std::pair<IntervalVector,Interval>> the constraints
   @param obj Vector vector of the linear form (objective) 
   @param maxit maximum number of iterations 
   @return result (empty si infaisable)
**/
Interval simplex_form(int dim, const IntervalVector &ivbox,
     std::vector<cstrrhs_status> clstats,
     const CstrVectMap &csts,
     const Vector &obj) {
   int sz = csts.size(); /* number of constraints, beside the ivbox */
   if (sz==0) return ivbox*obj;
   Intsimplex simp(ivbox,sz,clstats);
   for (const std::pair<const CstrVect,CstrRhs> &cst : csts) {
       simp.load_constraint(cst);
   }
   if (!simp.generate_init_basis(obj,true)) { return Interval::empty_set(); }
   simplex_ret sret = simp.simplex_mat();
   if (sret[INFEASIBLE]) return Interval::empty_set();
   double upper = simp.get_objective_value();

   if (!simp.generate_init_basis(obj,false)) { return Interval::empty_set(); }
   sret = simp.simplex_mat();
   if (sret[INFEASIBLE]) return Interval::empty_set();
   double lower = -simp.get_objective_value();
   return Interval(lower, upper);
}

/** basic start for simplification of polyhedron, using raytracing */
int basic_simplify_polyhedron(int dim,
    IntervalVector &ivbox, std::vector<cstrrhs_status> &clstats,
    CstrVectMap &csts, std::vector<CstrVectMap::iterator> &intmap,
    IntervalVector &inc) {
      Vector start(dim);
      inc &= ivbox;
      if (inc.is_empty()) start=ivbox.mid(); else start=inc.mid();
      /* first we consider all "orthogonal" directions */
      int nbcst = intmap.size();
      std::vector<Interval> alpha(dim+nbcst);
      std::vector<int> alphanum(2*(dim+nbcst));
      for (int j=0;j<nbcst;j++) {
          if (intmap[j]!=csts.end()) {
            alpha[dim+j]=Interval::all_reals();
          } else {
            alpha[dim+j]=Interval::empty_set();
          }
          alphanum[2*(dim+j)]=alphanum[2*(dim+j)+1]=-2;
      }
      for (int i=0;i<dim;i++) { 
          alpha[i] = ivbox[i]-start[i];
	  alphanum[2*i]=2*i; alphanum[2*i+1]=2*i+1;
          for (int j=0;j<nbcst;j++) {
             int u=dim+j;
             if (alpha[u].is_empty()) continue;
             Interval vectstart=alpha[i];
             CstrVectMap::iterator &cs_it = intmap[j];
             if (cs_it->first.vect[i]==0.0) {
                if (vectstart.contains(0.0)) continue;
                alphanum[2*u+1]=-2; alphanum[2*u]=-2;
                alpha[u].set_empty();
                continue;
             }
             Interval inter = vectstart/cs_it->first.vect[i];
             int pos = cs_it->first.vect[i]>0.0 ? 1 : 0;
             if (inter.lb()>alpha[u].lb()) alphanum[2*u+1]=2*i+pos;
                else if (inter.lb()==alpha[u].lb()) alphanum[2*u+1]=-2;
             if (inter.ub()<alpha[u].ub()) alphanum[2*u]=2*i+(1-pos);
                else if (inter.ub()==alpha[u].ub()) alphanum[2*u]=-2;
              alpha[u] &= inter;
          }  
      }
      for (int j=0;j<nbcst;j++) {
          CstrVectMap::iterator &cs_it = intmap[j];
          if (cs_it==csts.end()) continue;
          Interval vectstart=cs_it->second.val-cs_it->first.vect*start;
          for (int i=0;i<dim;i++) {
              if (alpha[i].is_empty()) continue;
              if (cs_it->first.vect[i]==0.0) {
                      if (vectstart.contains(0.0)) continue;
                      alphanum[2*i+1]=-2; alphanum[2*i]=-2;
                      alpha[i].set_empty();
                      continue;
              }
              Interval inter =
		vectstart/(cs_it->first.vect[i]);
              int pos = cs_it->first.vect[i]>0.0 ? 1 : 0;
              if (inter.lb()>alpha[i].lb()) alphanum[2*i+1]=2*(dim+j)+pos;
                else if (inter.lb()==alpha[i].lb()) alphanum[2*i+1]=-2;
              if (inter.ub()<alpha[i].ub()) alphanum[2*i]=2*(dim+j)+(1-pos);
                else if (inter.ub()==alpha[i].ub()) alphanum[2*i]=-2;
              alpha[i] &= inter;
          }
          for (int k=0;k<nbcst;k++) {
             int u=dim+k;
             if (alpha[u].is_empty()) continue;
             CstrVectMap::iterator &cs_it2 = intmap[k];
             double dir = cs_it->first.vect*cs_it2->first.vect;
             if (dir==0.0) {
                if (vectstart.contains(0.0)) continue;
                alphanum[2*u+1]=-2; alphanum[2*u]=-2;
                alpha[u].set_empty();
                continue;
             }
             Interval inter = vectstart/dir;
             int pos = dir>0.0 ? 1 : 0;
             if (inter.lb()>alpha[u].lb()) alphanum[2*u+1]=2*(dim+j)+pos;
                else if (inter.lb()==alpha[u].lb()) alphanum[2*u+1]=-2;
             if (inter.ub()<alpha[u].ub()) alphanum[2*u]=2*(dim+j)+(1-pos);
                else if (inter.ub()==alpha[u].ub()) alphanum[2*u]=-2;
              alpha[u] &= inter;
          }  
      }
      for (int i=0;i<dim+nbcst;i++) {
          if (alpha[i].is_empty()) continue;
          if (inc.is_empty()) {
             inc = start; 
             if (i<dim) inc[i] += alpha[i].mid(); else 
			inc += alpha[i].mid()*intmap[i-dim]->first.vect;
	  }
          int nrcst = alphanum[2*i]/2;
          if (nrcst>=0) {
               if (nrcst<dim) assign_filter(clstats[nrcst],
					rhs_ubnr,rhs_ubfilter);
               else { 
                  cstrrhs_status& st=intmap[nrcst-dim]->second.status;
                  assign_filter(st, rhs_ubnr,rhs_ubfilter);
               }
          }
          nrcst = alphanum[2*i+1]/2;
          if (nrcst>=0) {
               if (nrcst<dim) assign_filter(clstats[nrcst],
					rhs_lbnr,rhs_lbfilter);
               else { 
                  cstrrhs_status& st=intmap[nrcst-dim]->second.status;
                  assign_filter(st, rhs_lbnr,rhs_lbfilter);
               }
          }
      }
      return 0;
}


/** simplify a "polyhedron" defined as an IVbox + vector of (possibly 
    "fuzzy" linear forms + interval... 
    each vector is centered and kept if useful.
    returns -1 if the polyhedron is empty .
    also, gives (if possible) a ivbox included in the polyhedron
    */
int simplify_polyhedron(int dim, IntervalVector &ivbox,
     std::vector<cstrrhs_status> &clstats,
     const std::vector<std::pair<IntervalVector,Interval>> &csts,
     CstrVectMap &csts_rs, IntervalVector &inc) {
   csts_rs.clear();
   clstats.clear();
   inc.set_empty();
   int sz = csts.size(); /* number of constraints, beside the ivbox */
   if (sz==0) { clstats.resize(dim,rhs_nrfilter);
	        inc=ivbox; 
                int nb_non_flat=0;
                for (int i=0;i<dim;i++) 
		    if (!ivbox[i].is_degenerated()) nb_non_flat++;
                return nb_non_flat; }
   clstats.resize(dim, 0);
   Intsimplex simp(ivbox,sz,clstats);
   for (const std::pair<const IntervalVector,Interval> &cst : csts) {
       simp.load_constraint(cst.first,cst.second,0);
   }
//   bool check = simp.constraint_matrix(); /* FIXME! */
//   if (!check) return -1;

   /* we hope that now simp is fully bounded */
   Vector obj(dim);
   double upper,lower;
   int nb_non_flat=0;
   /* 1st, the ivbox */
   for (int i=0;i<dim;i++) {
       if (!simp.generate_init_basis(i,true)) { return -1; }
       simplex_ret sret = simp.simplex_mat();
       if (sret[INFEASIBLE]) { inc.set_empty(); return -1; }
       upper = simp.get_objective_value();
       if (upper<ivbox[i].ub()) {
          clstats[i][UB_R]=true; clstats[i][UB_RT]=true;
       } else
          clstats[i][UB_NR]=true; 
       obj[i]=-1;

       if (!simp.generate_init_basis(i,false)) {  return -1; }
       sret = simp.simplex_mat();
       if (sret[INFEASIBLE]) {  return -1; }
       lower = -simp.get_objective_value();
       if (lower>ivbox[i].lb()) {
          clstats[i][LB_R]=true; clstats[i][LB_RT]=true;
       } else
          clstats[i][LB_NR]=true; 
       ivbox[i] &= Interval(lower,upper);
       if (!ivbox[i].is_degenerated()) nb_non_flat++;
       simp.changeObjRowCol(i,ivbox[i],clstats[i]);
   }
   /* we use the ivbox to "center" the other constraints */
   if (nb_non_flat<=1) { /* mono-dimensionnel : pas de contrainte à garder */
     inc = ivbox;
     return nb_non_flat; 
   }
   for (const std::pair<IntervalVector,Interval> &cst : csts) {
       Interval bounds = cst.second;
       CstrVect cv = traduit_vect(ivbox,cst.first,bounds);
       if (fabs(cv.vdim)<1e-5) { /* almost a dimension */
          cv.vect[cv.bdim]=0.0; /* to get the "rest" */
          bounds+=cv.vect*ivbox;
          /* redundancy of ivbox[i] not completely guaranteed.
             non-redundancy is still ok */
          reset_filter(clstats[cv.bdim],rhs_rfilter);
//          clstats[cv.bdim][UB_R]=false; 
//          clstats[cv.bdim][LB_R]=false; 
          ivbox[cv.bdim]&=bounds;
       } else {
          CstrRhs rhs = { .val = bounds , .status = 0 };
          Interval withiv = cv.vect*ivbox;
          if (bounds.lb()<withiv.lb()) rhs.status[LB_R]=true;
          if (bounds.ub()>withiv.ub()) rhs.status[UB_R]=true;
          if (all_filter(rhs.status,rhs_rfilter)) continue;
          csts_rs.and_constraint(ivbox,cv,rhs);
       }
   }
   return simplify_polyhedron(dim, ivbox, clstats, csts_rs, inc, 0);
}

/* generate an including box from the result of a simplex computation */
/* except for flat dimension, we need the mid of the result to be inside
   the polyhedron and not at the boundary */
IntervalVector generate_including_box(const IntervalVector &ivbox,
		const  Intsimplex &simp,const CstrVectMap &csts,
		const std::vector<CstrVectMap::iterator> &intmap) {
   Vector start = simp.getExtremalPoint(ivbox);
   int dim=ivbox.size();
   IntervalVector direction(ivbox.size(),0.0);
   for (int i=0;i<dim;i++) {
      int col = simp.get_basis()[i];
      bool isneg = simp.is_neg_basis(col);
      if (col<dim) { direction[col]+= (isneg ? -1.0 : 1.0); }
      else { const Vector &v1 = intmap[col]->first.vect;
             direction += (isneg ? -v1 : v1);
           }
   }
   Interval alpha(Interval::all_reals());
   for (int i=0;i<dim;i++) {
      if (ivbox[i].is_degenerated()) continue; /* start _is_ in the flat box */
      Interval k1 = ivbox[i]-start[i]; 
		/* because start is a point ... */
      Interval k2 = direction[i]; 
      if (k2.contains(0.0) && !k1.contains(0.0)) 
		{ return IntervalVector(dim,Interval::empty_set()); }
      if (k2.ub()!=0) { alpha &= k1/(k2.ub()); }
      if (k2.lb()!=0) { alpha &= k1/(k2.lb()); }
      if (alpha.is_degenerated()) {  
  	  return IntervalVector(dim,Interval::empty_set()); 
      }
   }
   for(const CstrVectMap::value_type& csit : csts) {
      Interval sd = csit.first.vect*start;
      Interval k1 = csit.second.val-sd; 
		/* because start is a point ... */
      Interval k2 = csit.first.vect*direction; 
      if (k2.contains(0.0) && !k1.contains(0.0)) 
		{ return IntervalVector(dim,Interval::empty_set()); }
      if (k2.ub()!=0) { alpha &= k1/(k2.ub()); }
      if (k2.lb()!=0) { alpha &= k1/(k2.lb()); }
      if (alpha.is_degenerated()) {  
  	  return IntervalVector(dim,Interval::empty_set()); 
      }
   }
   IntervalVector center= start+alpha.mid()*direction;
   /* check the result by building a small box ? */
   Interval radius = Interval::pos_reals();
   for(const CstrVectMap::value_type& csit : csts) {
       Interval res = csit.second.val - csit.first.vect*center;
       radius &= res/dim;
       if (radius.is_empty() || radius.is_degenerated()) {
          return IntervalVector(dim,Interval::empty_set());
       }
   }
   for (int i=0;i<dim;i++) {
      if (ivbox[i].is_degenerated()) continue; 
      Interval rd(0.0,(ivbox[i]-center[i]).mig());
      rd &= radius;
      if (rd.is_degenerated()) {
          return IntervalVector(dim,Interval::empty_set());
      }
      center[i].inflate(rd.ub());
   }
   return center;
}

/* simplify the polyhedron, returns the number of non_flat variables, or -1
   for empty box. if recheck is true, we don't take any rhs_status as correct
   and check again all results .
   Note : if inc is degenerated, it is "recomputed",
   otherwise it is not modified */
int simplify_polyhedron(int dim, IntervalVector &ivbox,
     std::vector<cstrrhs_status> &clstats, CstrVectMap &csts, 
     IntervalVector &inc,  cstrrhs_status recheck) {
//   std::cerr << "simplify_polyhedron2\n";
   int sz = csts.size(); /* number of constraints, beside the ivbox */
   if (sz==0) { clstats.clear(); 
		clstats.resize(dim,rhs_nrfilter);
	        inc=ivbox; 
                int nb_non_flat=0;
                for (int i=0;i<dim;i++) 
		    if (!ivbox[i].is_degenerated()) nb_non_flat++;
                return nb_non_flat; }
   if (recheck!=0) for (int i=0;i<dim;i++) reset_filter(clstats[i],recheck);
   /* we keep a pointer for all elements of the map (correspondence
      column - element. It is valid except for removed elements */
   std::vector<CstrVectMap::iterator> intmap(sz,csts.end());
   CstrVectMap::iterator csit=csts.begin();
   double cnt=0;
   while (csit!=csts.end()) {
       reset_filter(csit->second.status,recheck);
       intmap[cnt]=csit;
       csit++; cnt++;
   }
//   std::cerr << "mat built : " << simp.get_mat() << "\n" << simp.get_objrow() << "\n";
//   basic_simplify_polyhedron(dim,ivbox,clstats,
//            csts, intmap, inc);
   Intsimplex simp(ivbox,sz,clstats);
   csit=csts.begin();
   while (csit!=csts.end()) {
       intmap[simp.load_constraint(*csit)-dim]=csit;
       csit++;
   }
   double upper,lower;
   int nb_non_flat=0;
   /* 1st, the ivbox */
   for (int i=0;i<dim;i++) {
       simp.extendObjRowCol(i,1.0);
       if (!clstats[i][UB_NR] && !clstats[i][UB_RT]) {
          if (!simp.generate_init_basis(i,true))
				 { inc.set_empty(); return -1; }
          simplex_ret sret = simp.simplex_mat();
//          std::cerr << "result simplex_max1 :  " << simp.get_mat() << "\n" << simp.get_objrow() << "\n";
          if (sret[INFEASIBLE]) { inc.set_empty(); return -1; }
#if 0
          for (int bv=0;bv<dim;bv++) {
               int vbasis=simp.get_basis()[bv];
               int col=vbasis/2;
               if (col<dim) {
                  if (vbasis%2==0) 
			assign_filter(clstats[col],rhs_ubnr,rhs_ubfilter);
                  else assign_filter(clstats[col],rhs_lbnr,rhs_lbfilter);
               } else {
                  cstrrhs_status& st=intmap[col-dim]->second.status;
                  if (vbasis%2==0) 
			assign_filter(st,rhs_ubnr,rhs_ubfilter);
                  else assign_filter(st,rhs_lbnr,rhs_lbfilter);
               }
          }
          if (inc.is_empty()) inc = generate_including_box(ivbox,simp,
			csts,intmap);
#endif
          upper = simp.get_objective_value();
          if (upper<ivbox[i].lb()) { inc.set_empty(); return -1; }
          if (upper>ivbox[i].ub()) {
		assign_filter(clstats[i],rhs_ubnr,rhs_ubfilter);
                upper = ivbox[i].ub();
          }
	  else 
		assign_filter(clstats[i],rhs_ubrt,rhs_ubfilter); 
          /* if we did not need i in the basis, the constraint is redundant */
       } else upper=ivbox[i].ub();
       if (!clstats[i][LB_NR] && !clstats[i][LB_RT]) {
          if (!simp.generate_init_basis(i,false)) { inc.set_empty(); return -1; }
          simplex_ret sret = simp.simplex_mat();
//          std::cerr << "result simplex_max2 :  " << simp.get_mat() << "\n" << simp.get_objrow() << "\n";
          if (sret[INFEASIBLE]) { inc.set_empty(); return -1; }
#if 0
          for (int bv=0;bv<dim;bv++) {
               int vbasis=simp.get_basis()[bv];
               int col=vbasis/2;
               if (col<dim) {
                  if (vbasis%2==0) 
			assign_filter(clstats[col],rhs_ubnr,rhs_ubfilter);
                  else assign_filter(clstats[col],rhs_lbnr,rhs_lbfilter);
               } else {
                  cstrrhs_status& st=intmap[col-dim]->second.status;
                  if (vbasis%2==0) 
			assign_filter(st,rhs_ubnr,rhs_ubfilter);
                  else assign_filter(st,rhs_lbnr,rhs_lbfilter);
               }
          }
          if (inc.is_empty()) inc = generate_including_box(ivbox,simp,
			csts,intmap);
#endif
          lower = -simp.get_objective_value();
          if (lower>ivbox[i].ub()) { inc.set_empty(); return -1; }
          if (lower<ivbox[i].lb()) {
		assign_filter(clstats[i],rhs_lbnr,rhs_lbfilter);
                lower = ivbox[i].lb();
          }
	  else 
		assign_filter(clstats[i],rhs_lbrt,rhs_lbfilter); 
#if 0
          if (!clstats[i][LB_NR]) 
		assign_filter(clstats[i],rhs_lbrt,rhs_lbfilter); 
#endif
       } else lower=ivbox[i].lb();
       if (ivbox[i].is_disjoint(Interval(lower,upper))) {
       std::cerr << "result in ivbox " <<  ivbox[i] << " " << lower << " " << upper << "\n";
        ivbox.set_empty(); inc.set_empty(); return -1;
       }
       ivbox[i] &= Interval(lower,upper);
       if (!ivbox[i].is_degenerated()) nb_non_flat++;
       simp.changeObjRowCol(i,ivbox[i],clstats[i]);
   }
   if (nb_non_flat<=1) { csts.clear(); inc=ivbox; return nb_non_flat; }
   /* now we minimize and keep (or not) the different constraints */
   for (unsigned int i=0;i<intmap.size();i++) {
       CstrVectMap::iterator csts_it=intmap[i];
       if (csts_it == csts.end()) continue;
       simp.extendObjRowCol(i+dim,1.0);
//       std::cerr << "contrainte etudiee..." << csts_it->second.status << "\n";
       if (!csts_it->second.status[UB_NR] && !csts_it->second.status[UB_RT]) {
         if (!simp.generate_init_basis(i+dim,true)) { inc.set_empty(); return -1; }
         simplex_ret sret = simp.simplex_mat();
         if (sret[INFEASIBLE]) { inc.set_empty(); return -1; }
//         std::cerr << "vbasis(" << (i+dim)*2+1 << ") : ";
#if 0
         for (int bv=0;bv<dim;bv++) {
              int vbasis=simp.get_basis()[bv];
//              std::cerr << vbasis << " ";
              int col=vbasis/2;
              if (col<dim) {
                  if (vbasis%2==0) 
			assign_filter(clstats[col],rhs_ubnr,rhs_ubfilter);
                  else assign_filter(clstats[col],rhs_lbnr,rhs_lbfilter);
              } else {
                  cstrrhs_status& st=intmap[col-dim]->second.status;
                  if (vbasis%2==0) 
			assign_filter(st,rhs_ubnr,rhs_ubfilter);
                  else assign_filter(st,rhs_lbnr,rhs_lbfilter);
              }
         }
         if (inc.is_empty()) inc = generate_including_box(ivbox,simp,
		csts,intmap);
#endif
         upper = simp.get_objective_value();
         if (upper<(*csts_it).second.val.lb()) { inc.set_empty(); return -1; }
         if (upper-(*csts_it).second.val.ub()>1e-10) { /* FIXME : threshold */
		assign_filter(csts_it->second.status,rhs_ubnr,rhs_ubfilter);
                upper = (*csts_it).second.val.ub();
         }
	 else 
		assign_filter(csts_it->second.status,rhs_ubrt,rhs_ubfilter); 
          /* if we did not need i in the basis, the constraint is redundant */
#if 0
         if (!csts_it->second.status[UB_NR]) { 
		assign_filter(csts_it->second.status,rhs_ubrt,rhs_ubfilter); 
         }
#endif
       } else upper=(*csts_it).second.val.ub();
       if (!csts_it->second.status[LB_NR] && !csts_it->second.status[LB_RT]) {
         if (!simp.generate_init_basis(i+dim,false)) { inc.set_empty(); return -1; }
         simplex_ret sret = simp.simplex_mat(); 
         if (sret[INFEASIBLE]) { inc.set_empty(); return -1; }
//         std::cerr << "vbasis(" << (i+dim)*2+1 << ") : ";
#if 0
         for (int bv=0;bv<dim;bv++) {
              int vbasis=simp.get_basis()[bv];
//              std::cerr << vbasis << " ";
              int col=vbasis/2;
              if (col<dim) {
                  if (vbasis%2==0) 
			assign_filter(clstats[col],rhs_ubnr,rhs_ubfilter);
                  else assign_filter(clstats[col],rhs_lbnr,rhs_lbfilter);
              } else {
                  cstrrhs_status& st=intmap[col-dim]->second.status;
                  if (vbasis%2==0) 
			assign_filter(st,rhs_ubnr,rhs_ubfilter);
                  else assign_filter(st,rhs_lbnr,rhs_lbfilter);
              }
         }
         if (inc.is_empty()) inc = generate_including_box(ivbox,simp,
		csts,intmap);
#endif
         lower = -simp.get_objective_value();
         if (lower>(*csts_it).second.val.ub()) { inc.set_empty(); return -1; }
         if (lower-(*csts_it).second.val.lb()<-1e-10) { /* FIXME : threshold */
		assign_filter(csts_it->second.status,rhs_lbnr,rhs_lbfilter);
                lower = (*csts_it).second.val.lb();
         }
	 else 
		assign_filter(csts_it->second.status,rhs_lbrt,rhs_lbfilter); 
          /* if we did not need i in the basis, the constraint is redundant */
#if 0
         if (!csts_it->second.status[LB_NR]) { 
		assign_filter(csts_it->second.status,rhs_lbrt,rhs_lbfilter); 
         }
#endif
       } else lower=(*csts_it).second.val.lb();
//       std::cerr << "result " <<  csts_it->second.val << " " << Interval(lower,upper) << "\n";
       csts_it->second.val &= Interval(lower,upper);
       if ((csts_it->second.val).is_empty()) { inc.set_empty(); return -1; }
       simp.changeObjRowCol(i+dim,csts_it->second.val,csts_it->second.status);
       if (all_filter(csts_it->second.status,rhs_rfilter)) {
//          std::cerr << "on elimine " << csts.size() << " " << csts_it->second.status << "\n";
          simp.activeCol(i+dim,false);
          csts.erase(csts_it);
	  intmap[i]=csts.end();
///          std::cerr << csts.size() << "\n";
          continue;
       }
    }
//    std::cerr << "fin simplify poly\n";
    return nb_non_flat;
}

CstrVect::CstrVect(int bdim, double vdim, const Vector &vect)  :
     bdim(bdim), vdim(vdim), vect(vect) {
}

CstrVect traduit_vect(const IntervalVector &box, const Vector &v,
                Interval &bounds) {
     CstrVect ret(-1,0.0,v);
     for (int i=0;i<box.size();i++) {
         if (box[i].is_degenerated()) {
            bounds -= ret.vect[i]*box[i];
            ret.vect[i]=0.0;
         }
     }
     double v1=-2.0, absv1=-1.0, v2=-2.0, absv2=-1.0;
     for (int i=0;i<ret.vect.size();i++) {
        double absvl = fabs(ret.vect[i]);
        if (absvl>absv1) { ret.bdim=i; v2=v1; absv2=absv1;
                           absv1=absvl; v1=ret.vect[i]; }
        else if (absvl>absv2) { absv2=absvl; v2=ret.vect[i]; }
     }
     if (ret.bdim>=0) {
        bounds /= v1;
        for (int i=0;i<ret.vect.size();i++) {
            if (i==ret.bdim) { ret.vect[i]=1.0; continue; }
            if (ret.vect[i]==0.0) continue;
            Interval a = ret.vect[i];
            a /= v1;
            ret.vect[i] = a.mid();
            bounds += (ret.vect[i]-a)*box[i];
        }
        ret.vdim = v2/v1;
     }
     return ret;
}

CstrVect traduit_vect(const IntervalVector &box, const IntervalVector &v,
                Interval &bounds) {
     Vector vmid = v.mid();
     bounds -= (v - vmid)*box;
     return traduit_vect(box,vmid,bounds);
}

static constexpr double threshold=1e-5;

/* ajoute une contrainte déjà "préparée" */
POLYOP_RET CstrVectMap::and_constraint(const IntervalVector &box,
                  const CstrVect &cv, const Interval &bounds,
                  const cstrrhs_status stat) {
    CstrRhs rhs = { .val = bounds, .status = stat };
    return this->and_constraint(box,cv,rhs);
}

POLYOP_RET CstrVectMap::and_constraint(const IntervalVector &box,
                  const CstrVect &cv, const CstrRhs &rhs) {
     const Interval &bounds = rhs.val;
     if (bounds.is_empty()) return POL_EMPTY;
     const cstrrhs_status &stat = rhs.status;
     CstrVectMap::iterator itlw = this->lower_bound(cv);
     if (itlw!=this->end()) {
        if (itlw->first==cv) {
          if (itlw->second.val.is_subset(bounds)) return POL_UNCHANGED;
          if (itlw->second.val.lb()>bounds.lb()) {
             assign_filter(itlw->second.status,stat,rhs_lbfilter);
          }
          if (itlw->second.val.ub()<bounds.ub()) {
             assign_filter(itlw->second.status,stat,rhs_ubfilter);
          }
          itlw->second.val &= bounds;
          if (itlw->second.val.is_empty()) return POL_EMPTY;
          return POL_CHANGED;
        } else if (itlw->first.bdim==cv.bdim) {
          const Vector &v2 = itlw->first.vect;
          double gap = 0.0;
          for (int i=0;i<box.size();i++) {
             gap += fabs(v2[i]-cv.vect[i]);
             if (gap>threshold) break; /* FIXME : give a threshold */
          }
          if (gap<=threshold) {
             Interval u = (bounds + (v2-cv.vect)*box);
             if (itlw->second.val.is_subset(u)) return POL_UNCHANGED;
             if (itlw->second.val.lb()>bounds.lb()) {
               assign_filter(itlw->second.status,stat,rhs_lbfilter);
             }
             if (itlw->second.val.ub()<bounds.ub()) {
               assign_filter(itlw->second.status,stat,rhs_ubfilter);
             }
             itlw->second.val &= u;
             if (itlw->second.val.is_empty()) return POL_EMPTY;
             return POL_CHANGED;
          }
        }
    }
    if (itlw!=this->begin()) {
        itlw--;
        if (itlw->first.bdim==cv.bdim) {
          const Vector &v2 = itlw->first.vect;
          double gap = 0.0;
          for (int i=0;i<box.size();i++) {
             gap += fabs(v2[i]-cv.vect[i]);
             if (gap>threshold) break; /* FIXME : give a threshold */
          }
          if (gap<=threshold) {
             Interval u = (bounds + (v2-cv.vect)*box);
             if (itlw->second.val.is_subset(u)) return POL_UNCHANGED;
             if (itlw->second.val.lb()>bounds.lb()) {
               assign_filter(itlw->second.status,stat,rhs_lbfilter);
             }
             if (itlw->second.val.ub()<bounds.ub()) {
               assign_filter(itlw->second.status,stat,rhs_ubfilter);
             }
             itlw->second.val &= u;
             if (itlw->second.val.is_empty()) return POL_EMPTY;
             return POL_CHANGED;
          }
        }
    }
    itlw = this->insert(itlw,std::make_pair(cv,rhs));
    return POL_CHANGED;
}

/* WARNING : untested code */
double dist_csts_box(const IntervalVector &ivbox,
       const Vector &v1, const Interval& b, const Vector &obj) {
   double lambda=1.0;
   int dim = ivbox.size();
   int actsign = 0; /* unknown sign */
   bool ok=false;
   while (!ok) {
      /* compute : the differential wrt lambda of our sum,
         as well as the value */
      Interval diff(0.0);
      double nextup=lambda, varup=-1.0;
      double nextdown=lambda, vardown=1.0;
      if (lambda>0.0) { 
	 diff = b.ub(); 
         nextdown=0.0; vardown=-b.diam();
      } else {
         diff = b.lb();
         nextup=0.0; varup=b.diam(); /* note: should not be used... */
      }
      for (int i=0;i<dim;i++) {
          double limit=obj[i]-lambda*v1[i];
          if (limit<0.0) {
              diff -= v1[i]*ivbox[i].lb();
              if (v1[i]!=0.0) {
                 double nxt = obj[i]/v1[i];
                 if (nxt>lambda) { /* v1[i]<0 */
                    if (varup==-1.0 || nxt<nextup) {
                        nextup=nxt;
                        varup=-v1[i]*ivbox[i].diam();
		    } else if (nxt==nextup) {
                        varup-=v1[i]*ivbox[i].diam();
                    }
	         } else if (nxt<lambda) { /* v1[i]>0 */
                    if (vardown==-1.0 || nxt>nextdown) {
                        nextdown=nxt;
                        vardown=-v1[i]*ivbox[i].diam();
		    } else if (nxt==nextdown) {
                        vardown-=v1[i]*ivbox[i].diam();
                    }
                 } else  /* should not happen, but if */
		    diff -= v1[i]*(ivbox[i]-ivbox[i].lb());
	      }
          }
          else if (limit>0.0) {
              diff -= v1[i]*ivbox[i].ub();
              if (v1[i]!=0.0) {
                 double nxt = obj[i]/v1[i];
                 if (nxt>lambda) { /* v1[i]>0 */
                    if (varup==-1.0 || nxt<nextup) {
                        nextup=nxt;
                        varup=v1[i]*ivbox[i].diam();
		    } else if (nxt==nextup) {
                        varup+=v1[i]*ivbox[i].diam();
                    }
	         } else if (nxt<lambda) { /* v1[i]<0 */
                    if (vardown==-1.0 || nxt>nextdown) {
                        nextdown=nxt;
                        vardown=v1[i]*ivbox[i].diam();
		    } else if (nxt==nextdown) {
                        vardown+=v1[i]*ivbox[i].diam();
                    }
                 } else  /* should not happen, but if */
		    diff -= v1[i]*(ivbox[i]-ivbox[i].ub());
	      }
          } else diff -= obj[i]*ivbox[i];
      }
      if (diff.contains(0.0)) break;
      if (diff.lb()>0.0) {
          if (actsign==-1) break;
          actsign=1; 
          if (nextdown==lambda) 
		/* we cannot go on. It should not happen,
		so we quit */
	       break;
          lambda=nextdown;
          if (diff.lb()+vardown<=0.0) break; 
      } else if (diff.ub()<0.0) {
          if (actsign==1) break;
          actsign=-1;
          if (nextup==lambda) 
		/* we cannot go on. It should not happen,
		so we quit */
	       break;
          lambda=nextup;
          if (diff.ub()+varup>=0.0) break; 
      }
   } 
   /* compute lambda b + (obj-lambda v1) ivbox */ 
   Interval gap = lambda * b + (obj-lambda*v1) * ivbox;
   return gap.ub();
}


}

