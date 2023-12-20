/**
 * Matrix manipulations with interval LU (or U L{-1}) decomposition
 * ----------------------------------------------------------------------------
 *  \date       2023
 *  \author     Damien Massé
 *  \copyright  Copyright 2023
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#ifndef _FLPL_INTLU_H
#define _FLPL_INTLU_H

#include <vector>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <bitset>

#include "flpl_def.h"
#include "permut.h"


namespace intflpl {

class IntLU {
    public:
      /* reprise d'une matrice existante ;
	 build = true => on crée une décomposition */
      IntLU (const IntervalMatrix& M, bool hasId, bool build);

      /* reset build */
      bool buildBasis();

      /* LU worked */
      bool inLUform() const;

      const IntervalMatrix &getLInv() const;
      const IntervalMatrix &getU() const;
 
      /* get the "inverse"  of the matrix 
         if nbRows=nbCols, the row order is preserved 
         otherwise, the rows are "in built order" */
      IntervalMatrix getInvB() const;

      Interval getDeterminant() const;

      /* solve MX = A. X is a vector column vector of dim nbC, with
         components outside the "basis" eq to 0. dim(A)=nbR, dim(X)=nbC */
      IntervalVector MXeqA(const IntervalVector& A) const;
      /* solve M'X = A, where M' is the restriction of M to the ``basis''
         column, "in built order" */
      IntervalVector MbXeqA(const IntervalVector& A) const;
      /* solve M'X = Mcol, where M' is the restriction of M to the ``basis''
         column, and Mcol the column col of M */
      IntervalVector MbXeqCol(int col) const;
      /* solve XM' = A */
      IntervalVector XMbeqA(const IntervalVector& A) const;
      /* "solve" XM = Ap (limited to basic columns of M) */
      IntervalVector XMeqA(const IntervalVector& A) const;
      /* "extend" a vector : from A, find X such that XM = Ap
         then returns A' so that XM = A' (note: use only UpM and not LowM)  */
      IntervalVector extendXMeqA(const IntervalVector& A) const;
      /* compute a row vector : of A'-1 A (use UpM only) */
      IntervalVector extendXMeqRow(int row) const;

      /* "exchange" two columns in the basis : insert a column in,
         and a column out */
      void exchangeColsBasis(int colIn, int colOut); 

      const Permut &getColBasis() const;

      friend std::ostream &operator<<(std::ostream& os, const IntLU &P);

   private:
      const int nbRows, nbCols;
      const IntervalMatrix *Mat;
      IntervalMatrix LowM;
      IntervalMatrix UpM;
      Permut PmRows, PmCols;
      bool hasId, isLU;

      double bestPiv(int &nbR, int &nbC) const;
      void computeNewRow(int nbR, int nbC);
      /* solve UpM X = A */
      IntervalVector UinvC(const IntervalVector &A) const;
      /* solve UpM_B X = A */
      IntervalVector UBinvC(const IntervalVector &A) const;
      /* solve X UpM_B = A */
      IntervalVector C_UBinv(const IntervalVector &A) const;
      /* "solve" X UpM_B = Row, with Row = 1 on one coef and 0 elsewhere */
      IntervalVector C_UBinvR(int row) const;
      /* "solve" X UpM = Ap, using only A to cols limited to the basis */
      IntervalVector C_Uinv(const IntervalVector &A) const;
      /* solve UpM_B X = Upcol */
      IntervalVector UBinvCol(int col) const;


      bool buildBasisInternal();
      void adjustRowUpm(int rowstart, int colstart);

};

}

#endif
