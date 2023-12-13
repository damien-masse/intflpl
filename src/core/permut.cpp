/**
 * Permutations (or partial permutations) manipulations (for matrices)
 * ----------------------------------------------------------------------------
 *  \date       2023
 *  \author     Damien Massé
 *  \copyright  Copyright 2023
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */


#include <vector>
#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <bitset>

#include "permut.h"


namespace intflpl {

/* representation of partial injective function from [0,nbA-1] to [0,nbB-1],
 * defined on [0,u] */
Permut::Permut(int _nbA, int _nbB) : nbA(_nbA), nbB(_nbB),
		defined(0), direct(_nbA,-1), reverse(_nbB,-1)
{
}

Permut::Permut(int _nbA, int _nbB, const Permut &P):
    nbA(_nbA), nbB(_nbB),
    defined(0), direct(_nbA,-1), reverse(_nbB,-1)
{  size_t i=0;
   int nbC=P.getnbA(); 
   for (;i<(size_t) nbC;i++) {
      if (nbC>=nbA) break;
      if (P.direct[i]<0 || P.direct[i]>=nbB) break;
      this->direct[i]=P.direct[i]; this->reverse[(size_t) P.direct[i]]=i;
   } 
   this->defined=i;
}

bool Permut::is_full() const { return defined==nbA; }

int Permut::nb_defined() const { return defined; }

bool Permut::is_defined(int a) const { return direct[(size_t) a]>=0; }

bool Permut::is_image(int b) const { return reverse[(size_t) b]>=0; }

int Permut::operator[](int a) const { return direct[(size_t) a]; }
int Permut::img(int a) const { return direct[(size_t) a]; }
int Permut::rev(int b) const { return reverse[(size_t) b]; }

void Permut::addNewImg(int b) {
    direct[(size_t) defined]=b;
    reverse[(size_t) b]=defined;
    defined++;
}

void Permut::exchangeIdx(int a1, int a2) {
     int v1 = direct[(size_t) a1];
     int v2 = direct[(size_t) a2];
     direct[(size_t) a2]=v1;
     direct[(size_t) a1]=v2;
     reverse[(size_t) v1]=a2;
     reverse[(size_t) v2]=a1;
}

void Permut::removeAndAddBack(int a, int b2) {
     int i;
     reverse[(size_t) direct[(size_t) a]]=-1;
     for (i=a;i<nbA-1;i++) {
        int v = direct[(size_t)(i+1)];
        if (v<0) break;
        reverse[(size_t) v]=i;
        direct[(size_t) i]=v;
     }
     direct[(size_t) i]=b2;
     if (b2>=0) reverse[(size_t) b2]=i;
}

void Permut::reset() {
     for (int i=0;i<defined;i++) {
        direct[(size_t) i] = -1;
     }
     for (int i=0;i<nbB;i++) {
        reverse[(size_t) i] = -1;
     }
     defined=0;
}

std::ostream &operator<<(std::ostream& os, const Permut &P) {
     os << "[";
     for (int i=0;i<P.nbA-1;i++) { os << P.direct[i] << ","; }
     if (P.nbA>0) os << P.direct[P.nbA-1];
     os << "]";
     return os;
}

Permut &operator>>(Permut &P, int i) {
     P.addNewImg(i);
     return P;
}

}
