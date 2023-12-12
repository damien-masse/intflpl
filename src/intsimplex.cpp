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


namespace intflpl {

Intsimplex::Intsimplex(int dim, int szalloc, bool useInterval) :
    dim(dim), sz(0), szalloc(szalloc), numobjcol(-1),
    objcolInd(0), status((1<<EMPTY) | 1<<(IDENTLEFT)),
    matInit(dim,dim+szalloc,Interval::zero()),  objcol(dim,Interval::zero()),
    matInitT(dim+szalloc,dim,Interval::zero()),
    objrow_offset(dim+szalloc,Interval::zero()), 
    colstat(dim+szalloc,colstatus_Unused), 
    objrowInit(dim+szalloc,Interval::zero()), useInterval(useInterval),
    LUform(nullptr)
{
    /* copy eye */
      for (int i=0;i<dim;i++) { matInit[i][i]=1.0; matInitT[i][i]=1.0; }
}

Intsimplex::Intsimplex(const IntervalVector &box, int szalloc,
	 std::vector<cstrrhs_status> clstats, bool useInterval) :
    dim(box.size()), sz(0), szalloc(szalloc), numobjcol(-1),
    objcolInd(0), status(1<<IDENTLEFT),
    matInit(dim,dim+szalloc,Interval::zero()), objcol(dim,Interval::zero()),
    matInitT(dim+szalloc,dim,Interval::zero()),
    objrow_offset(dim+szalloc,Interval::zero()), 
    colstat(dim+szalloc,colstatus_Unused),
    objrowInit(dim+szalloc,Interval::zero()), useInterval(useInterval),
    LUform(nullptr)
{
    /* copy eye */
      for (int i=0;i<dim;i++) { 
	matInit[i][i]=1.0;  matInitT[i][i]=1.0;
        objrowInit[i]=box[i];
        colstat[i].reset();
        if (!clstats[i][LB_R]) colstat[i].set(HASLB);
        if (!clstats[i][UB_R]) colstat[i].set(HASUB);
      }     
}

Intsimplex::~Intsimplex() {
      if (LUform!=nullptr) delete(LUform);
}


/* internal functions */

/** change_basis : change the basis of a matrix
 *  row : the row of the pivot, col : new basis column
 *  upper : put the upper bound of the last line at 0
 */
void Intsimplex::change_basis(int oldcol, int ncol) {
    if (oldcol!=ncol) /* sometimes oldcol==ncol if just the sign changes */ {
       if (useInterval) {
		/* FIXME : add useInterval */
       }
       LUform->exchangeColsBasis(oldcol, ncol);
       colstat[oldcol] ^= colstatus_Changebasis;
       colstat[ncol] ^= colstatus_Changebasis;
    } else if (useInterval) {
		/* FIXME : add useInterval */
    }
}

void Intsimplex::recompute_objoffset() {
    LUform->extendXMeqA(objrow_offset);
}


void Intsimplex::load_box(const IntervalVector &box, std::vector<cstrrhs_status> clstats) {
    assert(status[IDENTLEFT]);
    for (int i=0;i<dim;i++) { 
      objrowInit[i]=box[i];
      colstat[i].reset();
      if (!clstats[i][LB_R]) colstat[i].set(HASLB);
      if (!clstats[i][UB_R]) colstat[i].set(HASUB);
    }      
}

int Intsimplex::load_constraint(const Vector &cst, const Interval &val,
                       cstrrhs_status clstat) {
    assert(sz<szalloc);
    assert(status[IDENTLEFT]);
    const int col=dim+sz;
    for (int i=0;i<dim;i++) { 
        matInit[i][col]=cst[i];
    }    
    for (int i=0;i<dim;i++) { 
        matInitT[col][i]=cst[i];
    }    
    objrowInit[col]=val;
    colstat[col].reset();
    if (!clstat[LB_R]) colstat[col].set(HASLB);
    if (!clstat[UB_R]) colstat[col].set(HASUB);
    sz++;
    return col;
}

int Intsimplex::load_constraint(const IntervalVector &cst, const Interval &val,
                       cstrrhs_status clstat) {
    assert(sz<szalloc);
    assert(status[IDENTLEFT]);
    const int col=dim+sz;
    if (useInterval) {
        /* FIXME : useInterval */
    }
    for (int i=0;i<dim;i++) { 
        matInit[i][col]=cst[i];
    }    
    for (int i=0;i<dim;i++) { 
        matInitT[col][i]=cst[i];
    }    
    objrowInit[col]=val;
    colstat[col].reset();
    if (!clstat[LB_R]) colstat[col].set(HASLB);
    if (!clstat[UB_R]) colstat[col].set(HASUB);
    sz++;
    return col;
}

int Intsimplex::load_constraint(const std::pair<CstrVect,CstrRhs> newcst) {
    return this->load_constraint(newcst.first.vect,newcst.second.val,newcst.second.status);
}

void Intsimplex::loadObjRow(const IntervalVector &row, std::vector<cstrrhs_status> clstats) {
    assert(status[IDENTLEFT]);
    for (int i=0;i<dim+sz;i++) {
        objrowInit[i]=row[i];
        colstat[i].reset();
        if (!clstats[i][LB_R]) colstat[i].set(HASLB);
        if (!clstats[i][UB_R]) colstat[i].set(HASUB);
    }
}

void Intsimplex::changeObjRowCol(int col, const Interval &val, cstrrhs_status clstat) {
    assert(status[IDENTLEFT]);
    objrowInit[col]=val;
    colstat[col].reset();
    if (!clstat[LB_R]) colstat[col].set(HASLB);
    if (!clstat[UB_R]) colstat[col].set(HASUB);
}

void Intsimplex::extendObjRowCol(int col, double val) {
    assert(status[IDENTLEFT]);
    objrowInit[col].inflate(val);
}


bool Intsimplex::translateObjRow(const IntervalVector &row) {
    assert(status[CONSBASIS]);
    bool ok=true;
    for (int i=0;i<dim+sz;i++) {
        objrowInit[i]+=row[i];
    }
    /* check the basis */
    for (int i=0;i<dim;i++) {
        int vcol = LUform->getColBasis()[i];
        Interval nobj = objrowInit[vcol]+objrow_offset[vcol];
        if (colstat[vcol][NEGBASIS] && nobj.lb()<0) ok=false;
        if (!colstat[vcol][NEGBASIS] && nobj.ub()>0) ok=false;
    }
    //    if (!objrow[i].contains(0.0)) status[OPTBASIS]=false;
	// FIXME : add a check for OPTBASIS
	// FIXME 2 : change that for parametric LP
    if (!ok) status = 1<<INVALROW;
    return ok;
}


//bool Intsimplex::constraint_matrix() { /* can theoretically be used at anytime */
//    IntervalVector xT(dim, Interval::all_reals());
//    Interval vsave = objrow_offset[lastcol];
//    if (lastcolInd==0) {
//       objrow_offset[lastcol] = Interval::all_reals(); 
//    }
//    bool ret = bwd_mul(objrow,xT,mat,0.1);
//    if (!ret) /* the matrix has been emptied */
//	 { status = (1<<INVALCOL) | (1<<INVALROW) | (1<<UNFEAS); return false; }
//    if (lastcolInd==0) {
//       objrow_offset[lastcol] = vsave;
//    }
//    return true;
//}

Vector Intsimplex::getExtremalPoint(const IntervalVector &ivbox) const {
      assert(status[CONSBASIS]); /* or OPTBASIS ? */
      Vector val(dim,0.0);
      for (int i=0;i<dim;i++)  {
               val[i]=-objrow_offset[i].mid();
               if (!ivbox[i].contains(val[i])) {
		  if ((val[i]-ivbox[i].ub()>1.0) ||
			(val[i]-ivbox[i].lb()< (-1.0)))
		      std::cerr << "erreur vali: " << val[i] << 
		   	    " " << ivbox << "\n" << objrow_offset << "\n" << matInit << objrowInit << "\n" << (*LUform);
	 	  if (val[i]>ivbox[i].ub()) val[i]=ivbox[i].ub();
 		  else if (val[i]<ivbox[i].lb()) val[i]=ivbox[i].lb();
	       }
      }
     return val;
}

/* check the change of basis. R is constraint to enter, c the potentially
 * leaving column. Obj is the objective if needed (if not, obj=R) 
 * ratio is obj[rowC]/R_rowC, and R_rowC is of correct sign wrt asup
 * sign are the ``current'' bounds */
int Intsimplex::check_entry_col_C(int colR, const IntervalVector &R,
            const IntervalVector &obj, const Interval &ratio,
            const Interval &bd, double asup, 
	    const Permut &basis, std::vector<int> &sign,
	    Interval &result, int rowC) const {
    IntervalVector biasedR(dim,0.0);
    if (colR!=numobjcol) {
       biasedR = obj;
       for (int i=0;i<dim;i++) {
           if (i==rowC) continue;
           biasedR[i] -= ratio*obj[i];
           if (biasedR[i].lb()>0) sign[i]=1;
           else 
           if (biasedR[i].ub()<0) sign[i]=-1;
       }
    } /* si colR=lastcol, biasedR = 0 */
    Interval yiRi = (asup ? bd.ub() : bd.lb());
    for (int i=0;i<dim;i++) {
        if (i==rowC) continue;
        yiRi -= R[i]*(sign[i]==1 ? objrowInit[basis[i]].ub() : 
			objrowInit[basis[i]].lb());
    }
    /* yiRi must be in Ri*objrowInit[rowC] */
    Interval goal = R[rowC]*objrowInit[basis[rowC]];
    if (yiRi.ub()>goal.ub()) { /* a bit too high */
       bool ok = yiRi.lb()<=goal.ub(); /* ok=false: clearly too high */
       for (int i=0;i<dim;i++) {
          if (i==rowC) continue;
          if (biasedR[i].contains(0.0)) {
            Interval sr = sign[i]*R[i]; 
		/* positive => change of sign is bad */
            if (sr.ub()<0.0) {
              sign[i] = -sign[i];
              yiRi += sr*objrowInit[basis[i]].diam();
              if (yiRi.ub()<=goal.ub()) break;
              if (yiRi.lb()<=goal.ub()) ok=true;
            }
          } 
       }
       if (yiRi.ub()<goal.lb()) return -2;
			/* we left, we should find something better */
       if (!ok) { /* try again just to "find" something */
         for (int i=0;i<dim;i++) {
            if (i==rowC) continue;
            if (biasedR[i].contains(0.0)) {
              Interval sr = sign[i]*R[i]; 
	  	  /* positive => change of sign is bad */
              if (sr.lb()<0.0) {
                sign[i] = -sign[i];
                yiRi += sr*objrowInit[basis[i]].diam();
                if (yiRi.lb()<=goal.ub()) { ok=true; break; }
              }
            } 
         }
       }
    } else if (yiRi.lb()>bd.ub()) {
       bool ok = yiRi.ub()>=goal.lb(); /* ok=false: clearly too high */
       for (int i=0;i<dim;i++) {
          if (i==rowC) continue;
          if (biasedR[i].contains(0.0)) {
            Interval sr = sign[i]*R[i]; 
		/* negative => change of sign is bad */
            if (sr.lb()>0.0) {
              sign[i] = -sign[i];
              yiRi += sr*objrowInit[basis[i]].diam();
              if (yiRi.lb()>=goal.lb()) break;
              if (yiRi.ub()>=goal.lb()) ok=true;
            }
          } 
       }
       if (yiRi.lb()>goal.ub()) return -2;
			/* we left, we should find something better */
       if (!ok) { /* try again just to "find" something */
         for (int i=0;i<dim;i++) {
            if (i==rowC) continue;
            if (biasedR[i].contains(0.0)) {
              Interval sr = sign[i]*R[i]; 
	  	  /* negative => change of sign is bad */
              if (sr.ub()>0.0) {
                sign[i] = -sign[i];
                yiRi += sr*objrowInit[basis[i]].diam();
                if (yiRi.ub()>=goal.lb()) { ok=true; break; }
              }
            } 
         }
       }

    }
    if (yiRi.is_disjoint(goal)) {
       return -2;
    }
    /* on calcule l'objectif ? */
    result = ratio*(asup ? bd.ub() : bd.lb());
    for (int i=0;i<dim;i++) {
        result = biasedR[i]*
		(sign[i]==1 ? objrowInit[basis[i]].ub() :
				 objrowInit[basis[i]].lb());
    }
    if (yiRi.is_subset(goal)) {
        return 0;
    }
    return -1;
}
	    


/* check the satisfaction/entry of a constraint. */
int Intsimplex::check_entry_col(int col, Interval &res, const Permut &basis,
		std::vector<int> &bsign, bool &asup) const {
    /* first get the constraint  in the current basis */
    IntervalVector R = LUform->MbXeqCol(col);
    const Interval& bd = objrowInit[col];
    /* check satisfiability and emptiness */
    Interval sumInt = Interval::zero();
    Interval sumPt = objrow_offset[col]+bd;
    if (sumPt.contains(0.0)) return -1; /* ok */
    asup = sumPt.ub()<0.0;
    if (asup) { if  (!colstat[col][HASUB]) return -1; }
    else { if (!colstat[col][HASLB]) return -1; }

    std::vector<int> act_sign(dim);
    for (int i=0;i<dim;i++) {
        act_sign[i]=(colstat[basis[i]][NEGBASIS] ? -1 : 1);
        sumInt += R[i]*objrowInit[basis[i]];
    }
    if (sumInt.is_disjoint(bd)) return -2; /* empty */
    /* computation of objective in local basis */
    IntervalVector obj_local(dim);
    if (objcolInd==0) {
       obj_local = LUform->MbXeqA(this->objcol);
    } else if (col==numobjcol) {
       obj_local=R;
    } else {
       obj_local = LUform->MbXeqCol(numobjcol);
       if (objcolInd==-1) obj_local = -obj_local;
    }
    /* now we check the different rows */
    int brow=-1;
    for (int i=0;i<dim;i++) {
        if (colstat[basis[i]][LOCKEDIN]) continue;
        if (asup==(act_sign[i]==1)) { if (R[i].lb()<=0.0) continue; }
                               else { if (R[i].ub()>=0.0) continue; }
        Interval ratio;
        if (col==numobjcol) {
            if (objcolInd==-1) ratio=-1.0; else ratio=1.0;
        } else ratio=obj_local[i]/R[i];
        std::vector<int> sign_new(act_sign);
        Interval result;
        int ret = check_entry_col_C(col, 
		R, obj_local, ratio, objrowInit[basis[i]],
		asup, basis, sign_new, result, i);
        if (ret<0) continue;
        if (brow==-1 || (res.ub()>result.ub() ||
		(res.ub()==result.ub() && basis[i]<basis[brow]))) {
		res=result;
		brow = i;
		bsign.swap(sign_new);
        }
    }
    if (brow==-1) return -3;
    return basis[brow];
}


bool Intsimplex::generate_init_basis_internal(double mult,
	const IntervalVector &obj) {
//    std::cerr << "tagapoum!\n";
    if (LUform==nullptr) {
       LUform = new IntLU(matInit,true,true);
    } else {
       LUform->buildBasis();
    }
    for (int i=0;i<dim;i++) {
       if ((mult*obj[i].mid())>=0.0) {
          if (objrowInit[i].contains(POS_INFINITY)) 
			{ status[UNSAT]=true; return false; }
          objrow_offset[i]=-objrowInit[i].ub();
          colstat[i][NEGBASIS]=false;
       } else {
          if (objrowInit[i].contains(NEG_INFINITY)) 
			{ status[UNSAT]=true; return false; }
          objrow_offset[i]=-objrowInit[i].lb();
          colstat[i][NEGBASIS]=true;
       } 
       colstat[i][INBASIS]=true;
       colstat[i][LOCKEDIN]=false;
    }
    for (int i=dim;i<dim+sz;i++) {
       if  (colstat[i][UNUSED]) continue;
       colstat[i][OUTBASIS]=true;
       colstat[i][LOCKEDIN]=false;
    }
    recompute_objoffset();
    status[IDENTLEFT]=false;
    status[INITBASIS]=true;
    status[CONSBASIS]=true;
//    std::cerr << "tagada!\n";
    return true;

}

bool Intsimplex::generate_init_basis(const IntervalVector &column, bool getub) {
    assert(status[IDENTLEFT]);
    objcolInd=0;
    for (int i=0;i<dim;i++) {
        objcol[i]=(getub ? column[i] : -column[i]);
    }
    return this->generate_init_basis_internal(1.0, objcol);
}

bool Intsimplex::generate_init_basis(int col, bool getub) {
    assert(status[IDENTLEFT]);
    objcolInd=getub ? 1 : -1;
    objcol=(getub ? matInitT[col] : -matInitT[col]);
    numobjcol=col;
    return this->generate_init_basis_internal(getub ? 1.0 : -1.0,objcol);
}


int Intsimplex::check_entry_obj(int colout, const IntervalVector &RowB,
		bool &sign, double &bcoef) const {
    bool isup = !colstat[colout][NEGBASIS];
    int bestin = -1;
    Interval bgap(0.0);
    for (int colin=0;colin<dim+sz;colin++) {
       if (colstat[colin][UNUSED] || colstat[colin][INBASIS])
			 continue; /* cannot use the column */
       if (RowB[colin].contains(0.0)) continue;
      /* contrary to previous case, if the actual basis is positive
         (obj[dim].ub=0 and we want to get a positive new variable in the basis
	 hence decrease obj[col], we need a negative variable */
       bool toupper = (isup == (RowB[colin].ub()<0.0));
       double gap;
       if (toupper) {
   	  if (!colstat[colin][HASUB]) continue;
          gap = objrowInit[colin].ub()+objrow_offset[colin].lb();
       } else {
          if (!colstat[colin][HASLB]) continue;
          gap = objrowInit[colin].lb()+objrow_offset[colin].ub(); 
       } 
       if (gap<0.0) gap=0.0;
       Interval resgap = gap/RowB[colin].mag();
       if (bestin==-1 || resgap.ub()<bgap.ub()) {
          bestin=colin;
          resgap = bgap;
          sign = (toupper ? 1.0 : -1.0);
       }
    }
    bcoef = bgap.ub();
    return bestin;
}


/** utility function simplex_mat :
    apply the dual simplex on a matrix :
      the matrix used for the simplex :
        line 0->dim-1 : rows (basis)
        line dim : objective 
        column 0->dim-1 : start with ivbox
        column dim->lastcol-1 constraints 
        column lastcol : lambdas (either >=0 or <=0)
          [dim][dim+sz] : -objectif
     Eg the matrix
        1     2     -1      0  |  -1
        0     1     -1      1  |   1
     ---------------------------------
      [0,3] [-1,1] [1,3] [-1,0]   -5
     is interpreted with lambda>0 :
      -1 1 |-2 2 |1 -1 | 0  0 |   1
      0  0 |1 -1 |-1 1 | 1 -1 |   1
     ---------------------------------
      3  0 |1  1 |3 -1 | 0  1 |  -5
le pas de calcul donne :
      -1 1 |-1 1 | 0 0 | 1 -1 |   2
      0  0 |1 -1 |-1 1 | 1 -1 |   1
     ---------------------------------
      3  0 |2  0 | 2 0 | 1  0 |  -4
soit 
        1     1      0     -1  |  -2
        0    -1      1     -1  |  -1
     ---------------------------------
      [0,3] [0,2] [0,2]  [0,1]    -4
**/     
simplex_ret Intsimplex::simplex_mat() {
   int nb_iter=0;
   assert(status[CONSBASIS]);
   assert(!status[INVALCOL]);
   assert(!status[INVALROW]);
   simplex_ret retval=0;
//   std::cerr << "simplex_max " << mat << "\n" << objrow << "\n";
//   std::cerr << "lastcol " << lastcol << "\n";
//   for (auto &a : colstat) std::cerr << a << " ";
//   std::cerr << "\n";
   bool debug=false;
   while (nb_iter<maxit) {
       if (nb_iter>3*maxit/4) {
                std::cerr << "nbiter : " << nb_iter << "\n";
                debug=true;
       }
       /* looking for new variable in the basis */
       /* we try to find the best coefficient ...
          if not possible, get the first... */
       int bcolin=-1, bcolout=-1;
       double bres=-1.0; 
       bool b_asup;
       std::vector<int> bestsign(dim,0);
       bool has_uncert=false;
       const Permut &basis = LUform->getColBasis();
       for (int col=0;col<dim+sz;col++) {
     	   if (!colstat[col][OUTBASIS]) continue; 
	   std::vector<int> csign(dim,0);
           Interval res;
	   bool asup;
	   int ret = check_entry_col(col, res, basis, csign, asup);
           if (ret==-1) continue;
           if (ret==-2) {
	      status[UNFEAS]=true;
              return simplexret_Infeasible;  
           }
           if (ret==-3) { has_uncert=true; continue; }
           if (res.is_empty()) {
                std::cerr << "res empty!!!\n";
                std::cerr << matInit << "\n";
                std::cerr << objrowInit << "\n";
                std::cerr << "LU: " << LUform << "\n";
                continue;
           }
           double nobj = res.ub();
           if (debug) std::cerr << "col possible : " << col << " out : " << res << " obj : " << nobj << "\n";
           if (bcolin==-1 || nobj<bres) {
               bcolin=col; bres=nobj; bcolout=ret; bestsign.swap(csign);
	       b_asup = asup;
           }
       }
       if (bcolin==-1) {
            if (!has_uncert) break;  /* optimum reached */
	    retval[APPROXBASIS]=true;
            break;
            /* FIXME : ce n'est pas ça du tout. La proposition serait
			de faire une contraction sur les contraintes
			non respectées */
       } 
			
       if (debug) { std::cerr << "col choisie : " << bcolin << " sortie : " << bcolout << "\n"; }

       for (int row=0;row<dim;row++) {
           colstat[basis[row]][NEGBASIS]=(bestsign[row]==-1);
           objrow_offset[basis[row]]=-(bestsign[row]==1 ?
			objrowInit[basis[row]].ub() :
		        objrowInit[basis[row]].lb());
       }
       colstat[bcolin][NEGBASIS]=!b_asup;
       objrow_offset[bcolin]=-(b_asup ?
			objrowInit[bcolin].ub() :
		        objrowInit[bcolin].lb());
       this->change_basis(bcolin,bcolout);
       this->recompute_objoffset();
       status[INITBASIS]=false;
       nb_iter++;
   }
   if (nb_iter==maxit) {
      std::cerr << "maxiter reached " << maxit << "\n";
      std::cerr << "mat " << matInit << "\n";
      std::cerr << "objrow_offset " << objrow_offset << "\n";
      std::cerr << "objrowInit " << objrowInit << "\n";
      std::cerr << "LU : " << LUform << "\n";
      return simplexret_MaxIter;
   }
//   if (retval[APPROXBASIS]) std::cerr << "sortie approx " << nb_iter << "\n";
   status[OPTBASIS]=true;
   /* no ''real'' check of optimality, could be done as an afterwork */
   return retval;
}

/** use the simple matrix to generate a list of vertices from a 2D-facet
 *  we start with a already optimal matrix and "turn" the basis . 
 *  the ivbox enables to create the coordinates of the vertices ...
 *  locked_basis = not movable basis (all except 2 ?) */
std::list<Vector>
	Intsimplex::generate_vertices_2D(const IntervalVector &ivbox) {
//   std::cout << "ivbox : " << ivbox << "\n";
   assert(status[CONSBASIS]);
   std::list<Vector> lstret;
   const Permut &basis = LUform->getColBasis();
   std::vector<bool> ibasis(2*(dim+sz),false);
		 /* we "turn" until we reach a basis already seen */
   /* here we just "turn" around, hence a basic variable must stay twice */
   int lastin= -1; /* at first we don't have a move */
   /* we consider that the first movable dimension is the «~last in~» */
   for (int i=0;i<dim;i++) {
      if (colstat[basis[i]][LOCKEDIN]) continue;
      lastin=2*basis[i];
      if (colstat[basis[i]][NEGBASIS]) lastin++;
      break; 
   }
   ibasis[lastin] = true;
   bool fini=false;
   double gain=1.0; /* gain>0 : sommet "utile" */
   while (!fini) {
      if (gain>1e-10)  {
          Vector val = this->getExtremalPoint(ivbox);
          lstret.push_back(val);
      }
      /* identifying the exiting basis */
      int basisout=-1;
      int ncol=-1;
      bool bsign;
      for (int i=0;i<dim;i++) {
//          std::cout << "i= " << i << ":";
          /* is this basis movable (not equality) ? */
          if (colstat[basis[i]][LOCKEDIN]) continue;
          if (basis[i]==(lastin/2)) continue;
          /* to modify the place, the last line must be modified,
             __but__ all intervals must stay around 0... hence we try
	     the minimum possible > 0 */
          IntervalVector RowB = LUform->extendXMeqRow(i);
          bool sign;
          double bcoef;
	  
          int ncol = check_entry_obj(i,RowB,sign,bcoef);

          assert (ncol!=-1); /* can't find a move, theoretically not possible
                              as it is bounded 
				FIXME : unbounded modification ? */
//          std::cout << "\n";
          gain=bcoef; basisout=i; bsign=sign; break; /* 2D => only one move ! */
       }
       /* maintenant on se contente de faire entrer/sortir ncol */
       if (ncol==-1) { 
           std::cerr << "basis = " << basis ;
          std::cerr << "\ncolstat = "; for (int q=0;q<dim+sz;q++) std::cerr << colstat[q] << " ";
          std::cerr << "matInit = " << matInit << "\n";
          std::cerr << "LU = " << LUform << "\n";
          status[INVALCOL]=true; return lstret;
       }
       colstat[ncol][NEGBASIS]=(bsign==-1);
       objrow_offset[ncol]=-(bsign==1 ?
			objrowInit[ncol].ub() :
		        objrowInit[ncol].lb());
       this->change_basis(basisout,ncol);
       this->recompute_objoffset();
       lastin = ncol*2; if (bsign==-1) lastin++;
       /* si on est sur la base initiale, on s'arrête... */
       fini=ibasis[lastin];
       ibasis[lastin]=true;
   }
   status[INVALCOL]=true;
   return lstret;
}

}
