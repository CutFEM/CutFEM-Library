/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/
/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "GTypeOfFE_Sum.hpp"

#include "../common/Mesh3dn.hpp"
#include "../common/Mesh2dn.hpp"
#include "../common/Mesh1dn.hpp"

template <class Mesh>
GTypeOfFESum<Mesh>::GTypeOfFESum(KN<GTypeOfFE<Mesh> const *> &t)
    : GTypeOfFE<Mesh>(t), nbOfFE(t.N()), teb(t) //,
                                                // NN(nbOfFE+1),
                                                // DF(nbOfFE+1) ,
                                                // comp(nbOfFE+1)
{
   // buildData();
}

// template<class Mesh>
// void GTypeOfFESum<Mesh>::buildData() {

// //   map<const GTypeOfFE<Mesh> *,int> m;
//   {
//   int i=nbOfFE, j;
// //   while(i--) { // on va a l'envert pour avoir comp[i] <=i
// //     m[teb[i]]=i;
// //     std::cout << i << "\t" << m[teb[i]] << std::endl;
// //   }

//   // l'ordre comp est important comp est croissant  mais pas de pb.
//   i=nbOfFE;
//   while(i--) {
//     comp[i]= i;//m[teb[i]];
//     std::cout << i << "\t" << comp[i] << "\t" << std::endl;
//   }

//   std::cout << comp << std::endl;
//   std::cout << this->nb_sub_fem << std::endl;
//   // reservatition des intervalles en espaces
//   int n=0,N=0;
//   for ( j=0;j<nbOfFE;j++) { NN[j]=N;N+=teb[j]->N; }
//   NN[nbOfFE] = N;
//   //  reservation des interval en df
//   n=0;
//   for ( j=0;j<nbOfFE;j++) { DF[j]=n;n+=teb[j]->nbDoF; }
//   DF[nbOfFE] = n;

//   int ii=0;
//   for (int i=0;i<nbOfFE;++i) {
//     for (int j=0;j<teb[i]->nb_sub_fem;++j)
//       this->Sub_ToFE[ii++]=teb[i]->Sub_ToFE[j];
//   }

//   int c=0,c0=0;//, fcom=0;
//   for (int i=0;i<this->nb_sub_fem;i++) {
//     int N = this->Sub_ToFE[i]->N;
//     int ndofi = this->Sub_ToFE[i]->nbDoF;
//     //       this->first_comp[i]= fcom;
//     //       this->last_comp[i]= fcom+N;
// //     fcom += N;

//     for(int j=0;j<N;++j) {
//       this->begin_dfcomp[c] = c0 + this->Sub_ToFE[i]->begin_dfcomp[j] ;
//       this->end_dfcomp[c]   = c0 + this->Sub_ToFE[i]->end_dfcomp[j] ;
//       c++;
//     }
//     c0+=ndofi;
//   }

//   cout <<" NbDoF : " << this->nbDoF <<endl;
//   for(int i=0;i<this->N;++i)
//     std::cout << "      comp " << i
// 	      << " ["<<this->begin_dfcomp[i]<<", "<< this->end_dfcomp[i]
// 	      << "[ \n";

//   }
//   // construction de l'interpolation .
// }

// template<class Mesh>
// void GTypeOfFESum<Mesh>::buildInterpolation() {

//   int npi=0;
//   int nci=0;

//   for (int i=0;i<this->nb_sub_fem;i++)
//     {
//       npi +=this->Sub_ToFE[i]->NbPtforInterpolation;
//       nci +=this->Sub_ToFE[i]->NbcoefforInterpolation;
//     }
//   assert(this->NbcoefforInterpolation== nci);

//   this->pInterpolation.init(nci);
//   this->cInterpolation.init(nci);
//   this->dofInterpolation.init(nci);

//   map<RdHat,int,lessRd> mpt;
//   numPtInterpolation.init(npi);
//   int npp=0,kkk=0;
//   KN<RdHat> Ptt(npi);
//   for (int i=0;i<this->nb_sub_fem;i++) {

//     const GTypeOfFE<Mesh> &ti=*this->Sub_ToFE[i];

//     for(int p=0;p<ti.NbPtforInterpolation;++p,++kkk) {
//       Ptt[kkk]=ti.PtInterpolation[p];

//       cout << "    p = "<< p << " [ " << Ptt[kkk]<< "] ,  "
// 	   << kkk<< " "<< npp<< endl;;

//       if( mpt.find(Ptt[kkk]) == mpt.end()) mpt[Ptt[kkk]]=npp++;

//       numPtInterpolation[kkk]=mpt[Ptt[kkk]];
//     }
//   }
//   assert(this->NbPtforInterpolation==0);

//   this->NbPtforInterpolation=npp;
//   this->PtInterpolation.init(npp);
//   for(int i=0;i<npp;++i)
//     this->PtInterpolation[numPtInterpolation[i]]=Ptt[i];

//   int oc=0,odof=0;
//   for (int i=0,k=0;i<this->nb_sub_fem;i++)
//     {
//       const GTypeOfFE<Mesh> &ti=*this->Sub_ToFE[i];
//       for(int j=0;j<ti.NbcoefforInterpolation; ++j,++k)
// 	{
// 	  this->pInterpolation[k]   = numPtInterpolation[ti.pInterpolation[j]];
// 	  this->cInterpolation[k]   = ti.cInterpolation[j]+oc;
// 	  this->dofInterpolation[k] = ti.dofInterpolation[j]+odof;
// 	  this->coefInterpolation[k]= ti.coefInterpolation[j];
// 	}
//       oc += ti.N;
//       odof += ti.nbDoF;
//     }
// }

// template class GTypeOfFESum<Mesh2>;
// template class GTypeOfFESum<Mesh3>;
