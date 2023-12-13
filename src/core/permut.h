/**
 * Permutations (or partial permutations) manipulations (for matrices)
 * ----------------------------------------------------------------------------
 *  \date       2023
 *  \author     Damien Massé
 *  \copyright  Copyright 2023
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#ifndef _FLPL_PERMUT_H
#define _FLPL_PERMUT_H

#include <vector>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <bitset>



namespace intflpl {

/* representation of partial injective function from [0,nbA-1] to [0,nbB-1],
 * defined on [0,u] */
class Permut {
    public:
      /* initialisation à -1 */
      Permut (int _nbA, int _nbB);
      /* reprise d'une contrainte existante */
      Permut (int _nbA, int _nbB, const Permut& P);
      
      bool is_full() const;
      int nb_defined() const;
      bool is_defined(int a) const;
      bool is_image(int b) const;
      int operator[](int a) const;
      int img(int a) const;
      int rev(int b) const;
       
      int getnbA() const;
      int getnbB() const;

      void addNewImg(int b);
      void exchangeIdx(int a1, int a2);
      void removeAndAddBack(int a, int b);

      void reset();

      /* TODO : constuire un iterateur constant */
      
      friend std::ostream &operator<<(std::ostream& os, const Permut &P);
      friend Permut &operator>>(Permut &P, int i);
   
    private:
      const int nbA;
      const int nbB;
      int defined;
      std::vector<int> direct;
      std::vector<int> reverse;
};

inline int Permut::getnbA() const { return this->nbA; }
inline int Permut::getnbB() const { return this->nbB; }

}

#endif
