#ifndef _INTERPOLATION_MATRIX_HPP
#define _INTERPOLATION_MATRIX_HPP

#include "../util/util.hpp"
#include "../common/RNM.hpp"

#include "FESpace.hpp"


template<class Mesh> class GFElement;
template<class Mesh> class GFESpace;

/*
 * Class interpolation Matrix for computation in element
 */
template<class RdHat>
class InterpolationMatrix {
public:
  const int N,np,ncoef;          // dim, nb of point, nb of coef
  bool invariant;                // if same FE for all element
  int k;                         
  KN<RdHat> P;                   // point in ref element
  KNM<R> coef;                    // matrix from value at node to dof - coeff  

  template<class Mesh>
  InterpolationMatrix(const GFESpace<Mesh> &Vh);


private: 
  InterpolationMatrix(const InterpolationMatrix &);
  void operator=(const InterpolationMatrix &);
};


template <class RdHat>
ostream & operator<<(ostream& f,const InterpolationMatrix<RdHat> &M)
{ f<< M.N << " "<< M.k << " "<< M.np << " "<< M.ncoef << endl;
  f<< " = " << M.P ;
  f << "coef=" <<M.coef ;
  return f;
}


/*
 *   Create the interpolation matrix for the FE type
 *
 */
template<class RdHat>
template<class Mesh>
InterpolationMatrix<RdHat>::InterpolationMatrix(const GFESpace<Mesh> &Vh)
  :
  N(Vh.N),
  np(Vh.MaxNbNodePerElement),
  ncoef(Vh.MaxNbNodePerElement),
  invariant(true),
  k(-1),
  P(np),
  coef(ncoef,ncoef)
{ 
  Vh.TFE[0]->init(*this);
}




#endif
