/**
 * Matrix manipulations with interval LU (or U L{-1}) decomposition
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
#include "permut.h"
#include "intLU.h"


namespace intflpl {

double IntLU::bestPiv(int &nbR, int &nbC) const {
     nbR=-1;
     nbC=-1;
     double bstVal=0.0;
     for (int r=0;r<nbRows;r++) {
         if (PmRows.is_image(r)) continue;
         double sumRad=0.0, sumMagSq=0.0;
         double minRatio=0.0;
         int bestCol=-1;
         for (int c=0;c<nbCols;c++) {
            if (PmCols.is_image(c)) continue;
            sumRad+=UpM[r][c].rad();
            double mag = UpM[r][c].mag();
            sumMagSq+=(mag*mag);
            if (UpM[r][c].mig()==0.0) continue;
            double ratio = (0.001+UpM[r][c].rad())/UpM[r][c].mig();
            if (bestCol==-1 || ratio<minRatio) {
                minRatio=ratio;
                bestCol=c;
            }
         }
         if (bestCol!=-1) {
            double nvVal = ((sumRad*sumRad+0.001)*minRatio)/sumMagSq;
         std::cout << "row: " << r << " bcol: " << bestCol << " val: " << nvVal << "\n";
            if (nbR==-1 || nvVal<bstVal) {
               nbR=r; nbC=bestCol; bstVal=nvVal;
            }
         }
     }
     return bstVal;
}

void IntLU::computeNewRow(int nbR, int nbC) {
     Interval &piv = UpM[nbR][nbC];
     LowM[nbR][nbR]=1.0;
     PmRows.addNewImg(nbR);
     PmCols.addNewImg(nbC);
     for (int r=0;r<nbRows;r++) {
        if (PmRows.is_image(r)) continue;
        Interval val = UpM[r][nbC]/piv;
        for (int c=0;c<nbRows;c++) {
           if (!PmRows.is_image(c)) continue;
           LowM[r][c] -= LowM[nbR][c]*val;
        }
        for (int c=0;c<nbCols;c++) {
           if (PmCols.is_image(c)) continue;
           UpM[r][c] -= UpM[nbR][c]*val;
        }
        UpM[r][nbC]=0.0;
     }
}

bool IntLU::buildBasisInternal() {
     int nbR, nbC; 
     isLU=true;
     for (int i=0;i<nbRows;i++) { 
       double bv = this->bestPiv(nbR,nbC);
       if (nbR==-1) { isLU=false; break; }
       this->computeNewRow(nbR,nbC);
     }
     return isLU;
}

bool IntLU::buildBasis() {
     /* reset everything */
     UpM=(*Mat);
     LowM.clear();
     PmRows.reset();
     PmCols.reset();
     if (hasId) {
        for (int i=0;i<nbRows;i++) {
           LowM[i][i]=1.0;
           PmRows.addNewImg(i);
           PmCols.addNewImg(i);
        }
        isLU=true;
        return true;
     }
     return this->buildBasisInternal();
}

IntLU::IntLU(const IntervalMatrix& M, bool hasId, bool build) :
     nbRows(M.nb_rows()), nbCols(M.nb_cols()),
     Mat(&M), LowM(nbRows,nbRows,0.0), UpM(M), PmRows(nbRows,nbRows),
     PmCols(nbRows,nbCols), hasId(hasId), isLU(build)
{
     if (!build) return;
     if (hasId) {
        for (int i=0;i<nbRows;i++) {
           LowM[i][i]=1.0;
           PmRows.addNewImg(i);
           PmCols.addNewImg(i);
        }
        return;
     }
     this->buildBasisInternal();
}

const IntervalMatrix &IntLU::getLInv() const {
    return this->LowM;
}
const IntervalMatrix &IntLU::getU() const {
    return this->UpM;
}

Interval IntLU::getDeterminant() const {
    if (!isLU) return Interval::all_reals();
    Interval ret(1.0);
    for (int i=0;i<nbRows;i++) {
       int row= PmRows[i];
       int col= PmCols[i];
       ret *= UpM[row][col];
       ret /= LowM[row][row];
    }
    return ret;
}

IntervalVector IntLU::UinvC(const IntervalVector& A) const {
    IntervalVector res(nbCols,0.0);
    for (int i=nbRows-1;i>=0;i--) {
        int col = PmCols[i];
        int row = PmRows[i];
        Interval &rc=res[col];
        rc=A[row];
        for (int j=nbRows-1;j>i;j--) {
           int col2 = PmCols[j];
           rc-=res[col2]*UpM[row][col2];
        }
        rc/=UpM[row][col];
    }
    return res;
}


IntervalVector IntLU::UBinvC(const IntervalVector& A) const {
    IntervalVector res(nbRows,0.0);
    for (int i=nbRows-1;i>=0;i--) {
        int row = PmRows[i];
        Interval &rc=res[i];
        rc=A[row];
        for (int j=nbRows-1;j>i;j--) {
           int col2 = PmCols[j];
           rc-=res[j]*UpM[row][col2];
        }
        int col = PmCols[i];
        rc/=UpM[row][col];
    }
    return res;
}

IntervalVector IntLU::C_UBinv(const IntervalVector &A) const {
    IntervalVector res(nbRows,0.0);
    for (int c=0;c<=nbRows-1;c++) {
        int row = PmRows[c];
        int col = PmCols[c];
        Interval &rc=res[row];
        rc=A[c];
        for (int c2=0;c2<c;c2++) {
           int row2 = PmRows[c2];
           rc-=res[row2]*UpM[row2][col];
        }
        rc/=UpM[row][col];
    }
    return res;
}

IntervalVector IntLU::C_UBinvR(int row) const {
    IntervalVector res(nbRows,0.0);
    for (int c=row;c<=nbRows-1;c++) {
        int rowC = PmRows[c];
        int colC = PmCols[c];
        Interval &rc=res[rowC];
        rc=(c==row ? 1.0 : 0.0);
        for (int c2=row;c2<c;c2++) {
           int row2 = PmRows[c2];
           rc-=res[row2]*UpM[row2][colC];
        }
        rc/=UpM[rowC][colC];
    }
    return res;
}

IntervalVector IntLU::C_Uinv(const IntervalVector &A) const {
    IntervalVector res(nbRows,0.0);
    for (int c=0;c<=nbRows-1;c++) {
        int row = PmRows[c];
        int col = PmCols[c];
        Interval &rc=res[row];
        rc=A[col];
        for (int c2=0;c2<c;c2++) {
           int row2 = PmRows[c2];
           rc-=res[row2]*UpM[row2][col];
        }
        rc/=UpM[row][col];
    }
    return res;
}


IntervalVector IntLU::UBinvCol(int col) const {
    IntervalVector res(nbRows,0.0);
    int rv = PmCols.rev(col);
    if (rv>=0) {
       res[rv]=1.0;
       return res;
    }
    for (int i=nbRows-1;i>=0;i--) {
        int row = PmRows[i];
        Interval &rc=res[i];
        rc=UpM[row][col];
        for (int j=nbRows-1;j>i;j--) {
           int col2 = PmCols[j];
           rc-=res[j]*UpM[row][col2];
        }
        rc/=UpM[row][PmCols[i]];
    }
    return res;
}

IntervalVector IntLU::MXeqA(const IntervalVector& A) const {
    assert(isLU);
    IntervalVector LA = LowM*A;
    return this->UinvC(LA);
}
IntervalVector IntLU::MbXeqA(const IntervalVector& A) const {
    assert(isLU);
    IntervalVector LA = LowM*A;
    return this->UBinvC(LA);
}
IntervalVector IntLU::MbXeqCol(int col) const {
    assert(isLU);
    return this->UBinvCol(col);
}
IntervalVector IntLU::XMbeqA(const IntervalVector& A) const {
    assert(isLU);
    IntervalVector ret = this->C_UBinv(A);
    return ret*LowM;
}
IntervalVector IntLU::XMeqA(const IntervalVector& A) const {
    assert(isLU);
    IntervalVector ret = this->C_Uinv(A);
    return ret*LowM;
}

IntervalVector IntLU::extendXMeqA(const IntervalVector& A) const {
    assert(isLU);
    IntervalVector X = C_Uinv(A);
    IntervalVector ret(A);
    for (int i=0;i<nbCols;i++) {
      if (PmCols.is_image(i)) continue;
      ret[i] = 0.0;
      for (int j=0;j<nbRows;j++) ret[i] += X[j]*UpM[j][i];
    }
    return ret;
}

IntervalVector IntLU::extendXMeqRow(int row) const {
    assert(isLU);
    IntervalVector X = C_UBinvR(row); /* X[i]!=0 only for i>=row */
    IntervalVector ret(nbCols,0.0);
    for (int i=0;i<nbCols;i++) {
      int a = PmCols.rev(i);
      if (a==row) { ret[a]=1.0; continue;}
      if (a>=0) continue;
      ret[i] = 0.0;
      for (int j=row;j<nbRows;j++) ret[i] += X[j]*UpM[j][i];
    }
    return ret;
}

/* check and adjust a row */
void IntLU::adjustRowUpm(int rowstart, int colstart) {
    double bval=0;
    int bcol=-1;
    int row = PmRows[rowstart];
    IntervalVector &Vrow = UpM[row];
    for (int i=colstart;i<nbRows;i++) {
        Interval &it = Vrow[PmCols[i]];
        if (it.contains(0.0)) continue;
        double val = (1e-3+it.rad())/it.mig();
        if (bcol==-1 || val<bval) {
            bcol=i;
            bval=val/2.0;
        } else bval=bval/2.0; /* we prefer to use the ``first'' available
                                 value */
    }
    if (colstart!=bcol) {
       PmCols.exchangeIdx(colstart,bcol);
    }
    int ncol = PmCols[colstart];
    Interval &piv = UpM[row][ncol];
    for (int i=rowstart+1;i<nbRows;i++) {
       int row2 = PmRows[i];
       if (UpM[row2][ncol]==Interval::zero()) continue;
       LowM[row2] -= (UpM[row2][ncol]/piv)*LowM[row];
       UpM[row2] -= (UpM[row2][ncol]/piv)*UpM[row];
       UpM[row2][ncol] = Interval::zero();
    }
}

void IntLU::exchangeColsBasis(int colIn, int colOut) {
    assert(isLU);
    assert(!PmCols.is_image(colIn));
    int num = PmCols.rev(colOut);
    assert(num>=0);
    PmCols.removeAndAddBack(num,colIn);
    for (int i=num;i<nbRows;i++) {
         this->adjustRowUpm(i,i);
    }
}

const Permut &IntLU::getColBasis() const {
    return this->PmCols;
}

std::ostream &operator<<(std::ostream& os, const IntLU &P) {
     os << "Matrix LU decomposition\n";
     os << "row perm : " << P.PmRows;
     os << "\ncol perm : " << P.PmCols;
     os << "\nLow-1 : " << P.LowM;
     os << "\nUp : " << P.UpM;
     os << "\n";
     return os;
}


}
#undef __TESTLU__
#ifdef __TESTLU__

using namespace intflpl;
int main() {
     IntervalMatrix M(3,5);
     M[0][0]=-1.0; M[0][1]=3.0; M[0][2]=-1.0; M[0][3]=-2.0; M[0][4]=-2.5;
     M[1][0]=1.0; M[1][1]=Interval(-4.0,-4.0); M[1][2]=3.0; M[1][3]=2.0; M[1][4]=-2.5;
     M[2][0]=Interval(1.0,1.0); M[2][1]=-1.0; M[2][2]=0.5; M[2][3]=Interval(-2.0,-2.0); M[2][4]=2.5;
     
     IntLU lu(M,false,true);
     
     std::cout << lu;
     std::cout << "Low-1*Mat: " << (lu.getLInv()*M) << "\n";
     std::cout << "Determinant: " << lu.getDeterminant() << "\n\n";
 
     IntervalVector V(3,0.0);
     V[0]=1.0; V[1]=-0.3; V[2]=-1.5;
     IntervalVector Res = lu.MXeqA(V);
     std::cout << "MX=" << V << "\n";
     std::cout << "Res: " << Res << "\n";
     std::cout << "MRes: " << (M*Res) << "\n\n";

     IntervalVector Res2 = lu.MbXeqA(V);
     std::cout << "Res2: " << Res2 << "\n\n";

     for (int i=0;i<5;i++) {
        std::cout << "col" << i << " " << lu.MbXeqCol(i) << "\n";
     }
     IntervalVector V2(5,0.0);
     V2[0]=1.0; V2[1]=1.5; V2[2]=0.5; V2[3]=3.5; V2[4]=-1.0;
     std::cout << "\nXM'=" << V2 << "\n";
     IntervalVector Res3 = lu.XMeqA(V2);
     std::cout << "Res: " << Res3 << "\n";
     std::cout << "ResM: " << Res3*M << "\n";
     
     std::cout << "ResMbis: " << lu.extendXMeqA(V2) << "\n\n";

     lu.exchangeColsBasis(0,1);
     
     std::cout << lu;
     std::cout << "Low-1*Mat: " << (lu.getLInv()*M) << "\n";
     std::cout << "Determinant: " << lu.getDeterminant() << "\n\n";

     for (int i=0;i<5;i++) {
        std::cout << "col" << i << " " << lu.MbXeqCol(i) << "\n";
     }

     lu.exchangeColsBasis(1,0);
     
     std::cout << lu;
     std::cout << "Low-1*Mat: " << (lu.getLInv()*M) << "\n";
     std::cout << "Determinant: " << lu.getDeterminant() << "\n\n";

     for (int i=0;i<5;i++) {
        std::cout << "col" << i << " " << lu.MbXeqCol(i) << "\n";
     }

     lu.exchangeColsBasis(2,1);
     
     std::cout << lu;
     std::cout << "Low-1*Mat: " << (lu.getLInv()*M) << "\n";
     std::cout << "Determinant: " << lu.getDeterminant() << "\n\n";

     for (int i=0;i<5;i++) {
        std::cout << "col" << i << " " << lu.MbXeqCol(i) << "\n";
     }
     return 0;

}

#endif
