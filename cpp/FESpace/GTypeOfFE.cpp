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

#include "GTypeOfFE.hpp"

#include "../common/Mesh3dn.hpp"
#include "../common/Mesh2dn.hpp"
#include "../common/Mesh1dn.hpp"
// #include "interpolationMatrix.hpp"

/*
 *  compute the number of item used
 */
static int nbnode_d(const int nodePerItem[4], const std::vector<int> &nbItem) {
   const int ndf = nbItem[0] * nodePerItem[0] + nbItem[1] * nodePerItem[1] +
                   nbItem[2] * nodePerItem[2] + nbItem[3] * nodePerItem[3];
   return ndf;
}

static int getDFonWhat(const int *data, int nbdf, int on) {
   int kk = 0;
   for (int i = 0; i < nbdf; ++i)
      if (on == data[i])
         ++kk;
   // cout << " on " << on << " k = " << kk << endl;
   return kk;
}

static int startwhat(const std::vector<int> &nitem, int on) {
   int s = 0;
   for (int i = 0; i < on; ++i)
      s += nitem[i];
   return s;
}

dataTypeOfFE::dataTypeOfFE(const std::vector<int> &nitem, const int *Data,
                           int nbdf, int NN)
    : data(Data), dataalloc(0), ndfonVertex(getDFonWhat(data, nbdf, 0)),
      ndfonEdge(getDFonWhat(data, nbdf, startwhat(nitem, 1))),
      ndfonFace(getDFonWhat(data, nbdf, startwhat(nitem, 2))),
      ndfonVolume(getDFonWhat(data, nbdf, startwhat(nitem, 3))), nbDoF(nbdf),
      nbNode(nbnode_d(data + 4 * nbdf, nitem)), nbOfFE(1), N(NN),
      nbNodeOnWhat(data + 4 * nbdf), DFOnWhat(data), DFOfNode(data + nbdf),
      NodeOfDF(data + 2 * nbdf) {}

static int *builddata_d(const std::vector<int> nbItem,
                        const KN<dataTypeOfFE const *> &t) {

   int nbDoF   = 0;
   int dfon[4] = {0, 0, 0, 0};
   int k       = t.N();
   int nbFE    = 0;

   for (int i = 0; i < k; ++i) {
      nbDoF += t[i]->nbDoF;
      dfon[0] += t[i]->ndfonVertex;
      dfon[1] += t[i]->ndfonEdge;
      dfon[2] += t[i]->ndfonFace;
      dfon[3] += t[i]->ndfonVolume;
      nbFE += t[i]->nbOfFE;
   }

   KN<int> NN(k + 1), DF(k + 1), comp(k + 1);
   int n = 0, N = 0;
   for (int j = 0; j < k; j++) {
      N += t[j]->N;
   }

   //  reservation des interval en df
   for (int j = 0; j < k; j++) {
      DF[j] = n;
      n += t[j]->nbDoF;
   }
   DF[k] = n;

   int nwhat = 15;
   KN<int> w(nwhat), nn(nwhat);
   w  = 0;
   nn = 0;

   for (int j = 0; j < k; j++)
      for (int i = 0; i < t[j]->nbDoF; i++) {
         nn[t[j]->DFOnWhat[i]]++;
      }

   int nbn = 0;
   for (int j = 0; j < nwhat; j++) {
      if (nn[j])
         nn[j] = nbn++;
      else
         nn[j] = -1;
   }

   int *data  = new int[nwhat + 7 * nbDoF + N];
   int lgdata = nwhat + 7 * nbDoF + N;

   for (int i = 0; i < 4; ++i)
      data[i] = dfon[i];

   data[4] = nbDoF;
   data[5] = nbn;
   data[6] = N;
   data[7] = nbFE;

   int p = 8;
   for (int i = 0; i <= 3; ++i) { // loop over each kind of item
      int maxN = 0;
      for (int j = 0; j < k; ++j) {
         maxN = std::max(t[j]->nbNodeOnWhat[i], maxN);
      }
      data[p++] = maxN; // save on what item element
   }

   p = 15;
   for (int j = 0; j < k; j++) {
      for (int i = 0; i < t[j]->nbDoF; i++) {
         data[p++] = t[j]->DFOnWhat[i];
      }
   }

   KN<int> dln(nwhat);
   dln = 0;
   for (int j = 0; j < k; j++) {
      int cc = p;
      for (int i = 0; i < t[j]->nbDoF; i++) {
         data[p++] = t[j]->DFOfNode[i] + dln[t[j]->DFOnWhat[i]];
      }

      for (int i = 0; i < t[j]->nbDoF; i++)
         dln[t[j]->DFOnWhat[i]] =
             std::max(dln[t[j]->DFOnWhat[i]], data[cc++] + 1);
   }
   //  Ok si un noeud par what
   for (int j = 0; j < k; j++) {
      for (int i = 0; i < t[j]->nbDoF; i++) {
         data[p++] = nn[t[j]->DFOnWhat[i]];
      }
   }

   for (int i = p + nwhat; i < lgdata; ++i)
      data[i] = 0; // set 0 to the rest
   return data;
}

dataTypeOfFE::dataTypeOfFE(const std::vector<int> &nitemdim,
                           const KN<dataTypeOfFE const *> &t)
    : data(builddata_d(nitemdim, t)), dataalloc(data), ndfonVertex(data[0]),
      ndfonEdge(data[1]), ndfonFace(data[2]), ndfonVolume(data[3]),
      nbDoF(data[4]), nbNode(data[5]), N(data[6]), nbOfFE(data[7]),
      nbNodeOnWhat(data + 8), DFOnWhat(data + 15 + 0 * nbDoF),
      DFOfNode(data + 15 + 1 * nbDoF), NodeOfDF(data + 15 + 2 * nbDoF) {}

//
// static int *builddata_d(const int nbItem[4],
// 			const KN<dataTypeOfFE const*>& t,
// 			const dataTypeOfFE * tt)
// {
//   int nbDoFTime = tt->nbNode;
//   int nbDoF = 0;
//   int dfon[4] = {0,0,0,0};
//   int k = t.N();
//   int kt = nbDoFTime*k;
//   int nbFE = 0;
//
//   for(int i=0;i<k;++i) {
//     nbDoF   += t[i]->nbDoF*nbDoFTime;
//     dfon[0] += t[i]->ndfonVertex*nbDoFTime;
//     dfon[1] += t[i]->ndfonEdge*nbDoFTime;
//     dfon[2] += t[i]->ndfonFace*nbDoFTime;
//     dfon[3] += t[i]->ndfonVolume*nbDoFTime;
//     nbFE += t[i]->nbOfFE;
//   }
//
//   KN<int> NN(k+1), DF(kt+1) , comp(kt+1);
//   int n=0,N=0;
//   for ( int j=0; j<k; j++) { N += t[j]->N;}
//
//   //  reservation des interval en df
//   int ii=0;
//   for(int jt=0;jt<nbDoFTime;++jt)
//     for (int j=0;j<k;j++) {
//       DF[ii] = n ; n += t[j]->nbDoF;}
//   DF[kt] = n;
//
//   int nwhat = 15;
//   KN<int> w(nwhat),nn(nwhat);
//   w=0;
//   nn=0;
//
//   for ( int j=0; j<k; j++)
//     for ( int i=0; i<t[j]->nbDoF; i++) {
//       nn[t[j]->DFOnWhat[i]]++;
//     }
//
//   int nbn=0;
//   for( int j=0; j<nwhat; j++) {
//     if (nn[j]) nn[j] = nbn++;
//     else nn[j] = -1;
//   }
//
//   int * data = new int[nwhat+7*nbDoF*nbDoFTime +N];
//   int lgdata = nwhat+7*nbDoF*nbDoFTime+N;
//
//   for(int i=0;i<4;++i) data[i] = dfon[i];
//
//   data[4] = nbDoF;
//   data[5] = nbn;
//   data[6] = N;
//   data[7] = nbFE;
//
//   int p=8;
//   for(int i=0; i<=3; ++i){            // loop over each kind of item
//     int maxN = 0;
//     for(int j=0; j<k; ++j) {
//       maxN = max(t[j]->nbNodeOnWhat[i], maxN);
//     }
//     data[p++] = maxN;                // save on what item element
//   }
//
//   p = 15;
//   for ( int j=0; j<k; j++) {
//     for(int jt=0;jt<nbDoFTime;++jt){
//       for ( int i=0; i<t[j]->nbDoF; i++){
// 	data[p++] = t[j]->DFOnWhat[i];
//       }
//     }
//   }
//
//   KN<int> dln(nwhat);
//   dln=0;
//   for ( int j=0; j<k; j++) {
//     for(int jt=0;jt<nbDoFTime;++jt){
//       int  cc=p;
//       for ( int i=0; i<t[j]->nbDoF; i++) {
// 	data[p++] = t[j]->DFOfNode[i] + dln[t[j]->DFOnWhat[i]];
// 	// std::cout << data[p-1] << std::endl;
//       }
//
//       for ( int i=0;i<t[j]->nbDoF;i++)
// 	dln[t[j]->DFOnWhat[i]] = Max(dln[t[j]->DFOnWhat[i]],data[cc++]+1);
//     }
//   }
//   // std::cout << " --------------" << std::endl;
//
//   //  Ok si un noeud par what
//   for (int j=0; j<k; j++) {
//     for(int jt=0;jt<nbDoFTime;++jt){
//       for (int  i=0; i<t[j]->nbDoF; i++) {
// 	data[p++] = nn[t[j]->DFOnWhat[i]];
// 	// std::cout << data[p-1] << std::endl;
//       }
//     }
//   }
//   // getchar();
//   for(int i=p+nwhat; i<lgdata; ++i)
//     data[i] = 0;                        // set 0 to the rest
//   return data;
// }
//
//
// dataTypeOfFE::dataTypeOfFE(const int nitemdim[4], const KN<dataTypeOfFE
// const*> &t, 			   const dataTypeOfFE* tt)
//   :
//   data(builddata_d(nitemdim,t, tt)),
//   dataalloc(data),
//   ndfonVertex(data[0]),
//   ndfonEdge(data[1]),
//   ndfonFace(data[2]),
//   ndfonVolume(data[3]),
//   nbDoF(data[4]),
//   nbNode(data[5]),
//   N(data[6]),
//   nbOfFE(data[7]),
//   nbNodeOnWhat(data+8),
//   DFOnWhat(data+15+0*nbDoF),
//   DFOfNode(data+15+1*nbDoF),
//   NodeOfDF(data+15+2*nbDoF)
// {
// }

template class GTypeOfFE<Mesh1>;
template class GTypeOfFE<Mesh2>;
template class GTypeOfFE<Mesh3>;
