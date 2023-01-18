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

#ifndef FESpace_HPP_
#define FESpace_HPP_

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include "../num/util.hpp"
#include "GTypeOfFE_Sum.hpp"

#include "../common/Mesh3dn.hpp"
#include "../common/Mesh2dn.hpp"
#include "../common/Mesh1dn.hpp"

#include "../common/base_interface.hpp"
#include "../common/cut_mesh.hpp"
#include "QuadratureFormular.hpp"

#include "cutFEMConfig.h"
#ifdef USE_MPI
#include "../parallel/cfmpi.hpp"
#endif

template <class Mesh> class GFESpace;
template <class Mesh> class GFElement;
template <class Mesh> class GbaseFElement;

/*
 *  Structure with the FE space information
 *
 */
class DataFENodeDF {
   int *nbref; // pointer  on common ref counter
 public:
   int ndfon[4]; // array with nb of dof per item
   const int nbElement;
   const int nbNode;
   const int nbDoF;
   const int *const NodesOfElement;     // array of indices
   const int *const FirstDfOfNodeData;  // if dof different in elements
   const int *const FirstNodeOfElement; // 0
   const int MaxNbNodePerElement;
   const int MaxNbDFPerElement;
   const int MaxNbDFPerNode;
   bool constDfPerNode = true;
   int ndfonVertex() const { return ndfon[0]; }
   int ndfonEdge() const { return ndfon[1]; }
   int ndfonFace() const { return ndfon[2]; }
   int ndfonTet() const { return ndfon[3]; }

   DataFENodeDF(const DataFENodeDF &m)
       : nbref(m.nbref), nbElement(m.nbElement), nbNode(m.nbNode),
         nbDoF(m.nbDoF), NodesOfElement(m.NodesOfElement),
         FirstDfOfNodeData(m.FirstDfOfNodeData),
         FirstNodeOfElement(m.FirstNodeOfElement),
         MaxNbNodePerElement(m.MaxNbNodePerElement),
         MaxNbDFPerElement(m.MaxNbDFPerElement),
         MaxNbDFPerNode(maxdfon(m.ndfon)), constDfPerNode(m.constDfPerNode) {
      for (int i = 0; i < NbTypeItemElement; ++i)
         ndfon[i] = m.ndfon[i];
      (*nbref)++; // add one to the ref counter
   }

   DataFENodeDF(int andfon[NbTypeItemElement], int anbElement, int anbNode,
                int anbDoF, const int *aNodesOfElement,
                const int *aFirstDfOfNodeData, int aMaxNbNodePerElement,
                int aMaxNbDFPerElement,
                bool cstPerNode = true)
       : nbref(new int(0)), // new ref counter
         nbElement(anbElement), nbNode(anbNode), nbDoF(anbDoF),
         NodesOfElement(aNodesOfElement), FirstDfOfNodeData(aFirstDfOfNodeData),
         FirstNodeOfElement(0), MaxNbNodePerElement(aMaxNbNodePerElement),
         MaxNbDFPerElement(aMaxNbDFPerElement), MaxNbDFPerNode(maxdfon(andfon)),
         constDfPerNode(cstPerNode) {
      for (int i = 0; i < NbTypeItemElement; ++i)
         ndfon[i] = andfon[i];
   }

   DataFENodeDF(int andfon[NbTypeItemElement], int anbElement, int anbNode,
                int anbDoF, const int *aNodesOfElement,
                const int *aFirstDfOfNodeData, int aMaxNbNodePerElement,
                int aMaxNbDFPerElement, int aMaxNbDFPerNode,
                bool cstPerNode = true)
       : nbref(new int(0)), // new ref counter
         nbElement(anbElement), nbNode(anbNode), nbDoF(anbDoF),
         NodesOfElement(aNodesOfElement), FirstDfOfNodeData(aFirstDfOfNodeData),
         FirstNodeOfElement(0), MaxNbNodePerElement(aMaxNbNodePerElement),
         MaxNbDFPerElement(aMaxNbDFPerElement), MaxNbDFPerNode(aMaxNbDFPerNode),
         constDfPerNode(cstPerNode) {
      for (int i = 0; i < NbTypeItemElement; ++i)
         ndfon[i] = andfon[i];
   }

   ~DataFENodeDF() {
      if ((*nbref) == 0) // remove if nbref ==0
      {
         delete nbref;
         delete[] NodesOfElement;
         delete[] FirstDfOfNodeData;
         delete[] FirstNodeOfElement;
      } else
         (*nbref)--;
   }

 private:
   void operator=(const DataFENodeDF &);
};

/*
 *  Basic class for Element of the FESpace
 *
 */
template <class MMesh> class GbaseFElement {
 public:
   typedef MMesh Mesh;
   typedef GFESpace<Mesh> FESpace;
   typedef typename Mesh::Element Element;
   typedef Element E;
   typedef typename E::Rd Rd;
   typedef typename E::RdHat RdHat;
   typedef typename E::RdHatBord RdHatBord;
   typedef GQuadratureFormular<RdHat> QF;
   typedef GQuadratureFormular<RdHatBord> QFB;

   const GFESpace<Mesh> &Vh;
   const Element &T;
   const GTypeOfFE<Mesh> *tfe;
   const int N;      // dim arrived space
   const int number; // my index

   GbaseFElement(const GFESpace<Mesh> &aVh, int k);
   GbaseFElement(const GFESpace<Mesh> &aVh, int k, int kth);

   R EdgeOrientation(int i) const { return T.EdgeOrientation(i); }

   Rd Pt(RdHat Phat) const { return T(Phat); } // Ref to Global
   Rd PtHat(int i) const {
      assert(i < tfe->NbPtforInterpolation);
      RdHat Phat(tfe->Pt_Pi_h(i));
      return Phat;
   }
   Rd Pt(int i) const {
      assert(i < tfe->NbPtforInterpolation);
      RdHat Phat(tfe->Pt_Pi_h(i));
      return T(Phat);
   }

   R getMeasure() const {
      return T.mesure();
   } // mesure felstavat (french spel ;-) )
   R get_measure() const {
      return T.mesure();
   } // mesure felstavat (french spel ;-) )

   Rd map(Rd ip) const { return T(ip); }
   Rd mapToPhysicalElement(Rd ip) const { return T(ip); }

   int degre() const { return tfe->degre(); }
   int index() const { return number; };

   int whichDomain() const { return Vh.whichDomain(number); }
   int get_domain() const { return Vh.get_domain(number); }
   // bool isCut() const {return Vh.isCut(number);}
   // const Interface<Mesh>& getInterface() const { Vh.get_interface(number);}
};

template <class Mesh>
inline GbaseFElement<Mesh>::GbaseFElement(const GFESpace<Mesh> &aVh, int k)
    : Vh(aVh), T(Vh.Th[k]), tfe(aVh.TFE[k]), N(aVh.N), number(k) {}

// for cutFEM
template <class Mesh>
inline GbaseFElement<Mesh>::GbaseFElement(const GFESpace<Mesh> &aVh, int k,
                                          int kth)
    : Vh(aVh), T(Vh.Th[kth]), tfe(aVh.TFE[0]), N(aVh.N), number(k) {}

/*
 *    Class of the finite element of the FESpace
 *
 */
template <class Mesh> class GFElement : public GbaseFElement<Mesh> {
 public:
   typedef typename Mesh::Element Element;
   typedef Element E;
   typedef typename E::Rd Rd;
   typedef typename E::RdHat RdHat;
   typedef typename E::RdHatBord RdHatBord;
   typedef GQuadratureFormular<RdHat> QF;
   typedef GQuadratureFormular<RdHatBord> QFB;

   friend class GFESpace<Mesh>;
   const int *p; // indices nodes
   const int nb; // nb of nodes

   GFElement(const GFESpace<Mesh> *VVh, int k);
   GFElement(const GFESpace<Mesh> *VVh, int k, int kth);

   int NbDoF() const { return this->tfe->nbDoF; }
   int NbNode() const { return nb; }
   int operator[](int i) const;         //{ return  p[i]  Numdu noeud
   int operator()(int i, int df) const; // Num du DoF du noeud i de df local df
   int operator()(int df) const;
   int loc2glb(int i) const { return (*this)(i); }
   int loc2glb(int i, int it) const {
      return (*this)(i) + it * this->Vh.NbDoF();
   }

   void BF(const Rd &P, RNMK_ &val) const;
   void BF(const What_d whatd, const Rd &P, RNMK_ &val) const;
   void BF(const What_d whatd, const Rd &P, RNMK_ &val, const RNM_ &J) const;

   void set(InterpolationMatrix<RdHat> &M) const {
      this->tfe->set(this->Vh.Th, this->T, M, 0, 0, 0);
   }

   KN_<R> &Pi_h(KNM_<R> vpt, RN_ &vdf) const;

   int DFOnWhat(int i) const { return this->tfe->DFOnWhat[i]; }

   // df is the df in element
   int NodeOfDF(int df) const { return this->tfe->NodeOfDF[df]; } // a node
   int DFOfNode(int df) const {
      return this->tfe->DFOfNode[df];
   } // the df of the node

   int dfcbegin(int ic) const { return this->tfe->begin_dfcomp[ic]; }
   int dfcend(int ic) const { return this->tfe->end_dfcomp[ic]; }

   // GFElement(GFElement&&) = default;
   // GFElement& operator=(GFElement&& v) = default;
   // GFElement(const GFElement&) = delete;
   // GFElement& operator=(const GFElement& v) = delete;
};

template <class Mesh>
GFElement<Mesh>::GFElement(const GFESpace<Mesh> *VVh, int k)
    : GbaseFElement<Mesh>(*VVh, k), p(this->Vh.PtrFirstNodeOfElement(k)),
      nb(this->Vh.NbOfNodesInElement(k)) {}
template <class Mesh>
GFElement<Mesh>::GFElement(const GFESpace<Mesh> *VVh, int k, int kth)
    : GbaseFElement<Mesh>(*VVh, k, kth), p(this->Vh.PtrFirstNodeOfElement(k)),
      nb(this->Vh.NbOfNodesInElement(k)) {}

template <class Mesh> inline int GFElement<Mesh>::operator[](int i) const {
   return p ? p[i] : ((&this->T[i]) - this->Vh.Th.vertices);
}

template <class Mesh>
inline int GFElement<Mesh>::operator()(int i, int df) const {
   return this->Vh.FirstDFOfNode(p ? p[i]
                                   : ((&this->T[i]) - this->Vh.Th.vertices)) +
          df;
}

template <class Mesh> inline int GFElement<Mesh>::operator()(int df) const {

   return operator()(NodeOfDF(df), DFOfNode(df));
}

template <class Mesh>
inline void GFElement<Mesh>::BF(const Rd &P, RNMK_ &val) const {
   this->tfe->FB(Fop_D0 | Fop_D1, this->T, P, val);
}

template <class Mesh>
inline void GFElement<Mesh>::BF(const What_d whatd, const Rd &P,
                                RNMK_ &val) const {
   this->tfe->FB(whatd, this->T, P, val);
}

template <class Mesh>
inline void GFElement<Mesh>::BF(const What_d whatd, const Rd &P, RNMK_ &val,
                                const RNM_ &J) const {
   this->tfe->FB(whatd, this->T, P, val, J);
}

template <class Mesh>
KN_<R> &GFElement<Mesh>::Pi_h(KNM_<R> vpt, KN_<R> &vdf) const {
   // compute  the interpolation
   // in : vpt  value of componant j at point p : vpt(p,j)
   // out: vdf  value du the degre of freedom
   const KN<IPJ> ipj(this->tfe->ipj_Pi_h);
   const KN<Rd> PtHat(this->tfe->Pt_Pi_h);
   KN<R> Aipj(ipj.N());
   this->tfe->get_Coef_Pi_h(*this, Aipj);

   vdf = R();
   for (int i = 0; i < Aipj.N(); i++) {
      const IPJ &ipj_i(ipj[i]);
      vdf[ipj_i.i] += Aipj[i] * vpt(ipj_i.p, ipj_i.j);
   }
   return vdf;
}

/*
 *   Class for the FESpace
 *
 */
template <class MMesh> class GFESpace : public DataFENodeDF {
 public:
   typedef MMesh Mesh;
   typedef GFElement<Mesh> FElement;
   typedef typename Mesh::Element Element;
   typedef typename Element::Face Face;
   typedef typename Mesh::BorderElement BorderElement;
   typedef typename Mesh::Rd Rd;
   typedef GTypeOfFE<Mesh> TypeOfFE;
   typedef GQuadratureFormular<typename Element::RdHat> QFElement;
   typedef GQuadratureFormular<typename BorderElement::RdHat> QFBorderElement;

   const Mesh &Th;
   KN<const GTypeOfFE<Mesh> *> TFE;
   const int N;        // dim espace d'arrive
   const int Nproduit; // 1 if non constant Max number df par node. else Max
                       // number df par node..
   const BasisFctType basisFctType;
   const int polynomialOrder;

   GFESpace(const Mesh &TTh, const GTypeOfFE<Mesh> &tfe)
       : DataFENodeDF(this->BuildDFNumbering(
             TTh, tfe.ndfonVertex, tfe.ndfonEdge, tfe.ndfonFace,
             tfe.ndfonVolume, tfe.nbNodeOnWhat[0], tfe.nbNodeOnWhat[1],
             tfe.nbNodeOnWhat[2], tfe.nbNodeOnWhat[3], tfe.N)),
         Th(TTh), TFE(1, 0, &tfe), N(tfe.N),
         Nproduit(FirstDfOfNodeData ? 1 : MaxNbDFPerNode),
         basisFctType(tfe.basisFctType), polynomialOrder(tfe.polynomialOrder) {}

   GFESpace(const ActiveMesh<Mesh> &TTh, const GFESpace &vh)
       : DataFENodeDF(vh.BuildDFNumbering(TTh)), Th(TTh.Th),
         TFE(1, 0, vh.TFE(0)), N(TFE[0]->N),
         Nproduit(FirstDfOfNodeData ? 1 : MaxNbDFPerNode),
         basisFctType(vh.basisFctType), polynomialOrder(vh.polynomialOrder) {}
   DataFENodeDF BuildDFNumbering(const ActiveMesh<Mesh> &TTh) const;

   DataFENodeDF BuildDFNumbering(const Mesh &TTh, int dfon[NbTypeItemElement],
                                 int nndon[NbTypeItemElement], int N = 1);

   DataFENodeDF BuildDFNumbering(const Mesh &TTh, int ndfv, int ndfe, int ndff,
                                 int ndft, int nndv, int nnde, int nndf,
                                 int nndt, int N = 1) {
      int dfon[NbTypeItemElement] = {ndfv, ndfe, ndff, ndft};
      int ndon[NbTypeItemElement] = {nndv, nnde, nndf, nndt};
      return BuildDFNumbering(TTh, dfon, ndon, N);
   }

   const int *PtrFirstNodeOfElement(int k) const {
      return NodesOfElement ? NodesOfElement + (FirstNodeOfElement
                                                    ? FirstNodeOfElement[k]
                                                    : k * MaxNbNodePerElement)
                            : 0;
   }

   int NbOfNodesInElement(int k) const {
      return FirstNodeOfElement // 0 if const dof per element
                 ? FirstNodeOfElement[k + 1] - FirstNodeOfElement[k]
                 : MaxNbNodePerElement;
   }

   int FirstDFOfNode(int i) const {
      return (FirstDfOfNodeData ? FirstDfOfNodeData[i] : i * Nproduit);
   }

   ~GFESpace() {}

   virtual FElement operator[](int k) const { return FElement(this, k); }
   FElement operator[](const Element &K) const { return FElement(this, Th(K)); }

   int operator()(int k) const { return NbOfNodesInElement(k); }
   int operator()(int k, int i) const { //  the node i of element k
      return NodesOfElement ? *(PtrFirstNodeOfElement(k) + i) : Th(k, i);
   }

   virtual int idxElementFromBackMesh(int k) const { return k; }
   virtual int idxElementInBackMesh(int k) const { return k; }
   virtual int idxElementFromBackMesh(int k, int i) const { return k; }

   virtual int whichDomain(int k) const { return -1; }
   virtual int get_domain(int k) const { return -1; }
   virtual int getNeighborElement(int k, int &j, int domain = 0) const {
      return Th.ElementAdj(k, j);
   }
   virtual int nbDomain() const { return 1; }
   virtual bool containBackElement(int k) const { return true; }
   // virtual bool isCut(int k) const { return false;}
   // virtual bool isCut() const { return (this->gamma.size()>0);}
   virtual bool isCutSpace() const { return false; }
   virtual std::vector<int> idxAllElementFromBackMesh(int k, int d) const {
      assert(d == -1);
      std::vector<int> v = {idxElementFromBackMesh(k)};
      return v;
   }
   virtual const GFESpace<Mesh> &get_back_space() const { return *this; }
   virtual const ActiveMesh<Mesh> &get_mesh() const {
      assert(0);
      auto m = std::make_shared<ActiveMesh<Mesh>>(Th);
      return *m;
   }

   int NbNode() const { return this->nbNode; }
   int NbDoF() const { return this->nbDoF; }
   int NbElement() const { return this->nbElement; }
   int get_nb_element() const { return this->nbElement; }
   int get_nb_dof() const { return this->nbDoF; }

// int NbInnerFaces() const { return Th.nbInnerFaces();}
#ifdef USE_MPI
   virtual int first_element() const {
      return MPIcf::first_element(this->nbElement);
   }
   virtual int next_element() const {
      return MPIcf::next_element(this->nbElement);
   }
   virtual int last_element() const {
      return MPIcf::last_element(this->nbElement);
   }

   // virtual int first_face() const { return
   // MPIcf::first_element(this->NbInnerFaces());} virtual int next_face() const
   // {  return MPIcf::next_element(this->NbInnerFaces());} virtual int
   // last_face() const {  return MPIcf::last_element(this->NbInnerFaces());}

   virtual int first_boundary_element() const { return MPIcf::my_rank(); }
   virtual int next_boundary_element() const { return MPIcf::size(); }
   virtual int last_boundary_element() const { return this->Th.nbBrdElmts(); }
#else
   virtual int first_element() const { return 0; }
   virtual int next_element() const { return 1; }
   virtual int last_element() const { return this->nbElement; }

   // virtual int first_face() const { return 0;}
   // virtual int next_face() const {return 1;}
   // virtual int last_face() const { return this->NbInnerFaces();}

   virtual int first_boundary_element() const { return 0; }
   virtual int next_boundary_element() const { return 1; }
   virtual int last_boundary_element() const { return this->Th.nbBrdElmts(); }
#endif

   virtual void info() const {
      // std::cout << "FESpace \t" << this << std::endl;
      // std::cout << "----------------------------------" << std::endl;
      std::cout << " nb node    \t" << NbNode() << std::endl;
      std::cout << " nb dof     \t" << NbDoF() << std::endl;
      std::cout << " nb element \t" << NbElement() << std::endl;
   }
};

template <typename Mesh>
DataFENodeDF
GFESpace<Mesh>::BuildDFNumbering(const ActiveMesh<Mesh> &TTh) const {

   int ndfon[4] = {this->ndfonVertex(), this->ndfonEdge(), this->ndfonFace(),
                   this->ndfonTet()};
   int nSub     = TTh.get_nb_domain();

   int nt = 0, nv = 0, ndof = 0;
   int maxNodePerElement = this->MaxNbNodePerElement;
   int maxDFPerElement   = this->MaxNbDFPerElement;

   // to know where are the nodes
   KN<int> NodeOnWhat(maxNodePerElement);
   {
      int c          = 0;
      const int nk[] = {Element::nv, Element::ne, Element::nf, Element::nt};
      for (int i = 0; i < 4; ++i) {
         for (int k = 0; k < nk[i]; ++k) {
            for (int j = 0; j < this->TFE(0)->nbNodeOnWhat[i]; ++j) {
               NodeOnWhat(c++) = i;
            }
         }
      }
   }

   // p => idx nodes of elements : Vh.nbElement*NodesPerElement
   // node are not the vertex of the mesh but node for the finite element
   // they are build in the Space so we can just use this for
   // the cut space
   int *nodeOfElement = new int[maxNodePerElement * TTh.get_nb_element()];
   int *firstDfOfNode = nullptr;
   // nodeOfElement = new int(maxNodePerElement*TTh.get_nb_element());
   // get number of nodes and dof
   {
      // for(int i=0;i<nSub;++i) nts(i) = TTh.get_nb_element(i);
      nt = TTh.get_nb_element();
      std::vector<std::vector<int>> idxP_globToLoc(nSub);
      for (int i = 0; i < nSub; ++i)
         idxP_globToLoc[i].assign(this->nbNode, -1);
      int c = 0;
      for (int k = 0; k < TTh.get_nb_element(); ++k) {

         int kb     = TTh.idxElementInBackMesh(k);
         int dom    = TTh.get_domain_element(k);
         int nbNode = (*this)[0].NbNode();

         for (int j = 0; j < nbNode; ++j) {
            int iglb = (*this)(kb, j);

            int val            = idxP_globToLoc[dom][iglb];
            nodeOfElement[c++] = (val == -1) ? nv : val;
            if (val == -1) {
               idxP_globToLoc[dom][iglb] = nv++;
               ndof += this->TFE(0)->ndfOn()[NodeOnWhat(j)];
            }
         }
      }
   }

   if (!this->constDfPerNode) {
      int kk = 0, nn = 0;
      int *p        = nodeOfElement;
      firstDfOfNode = new int[nv + 1];
      int *pp       = firstDfOfNode;

      const FElement &FK((*this)[0]);
      const int *ndfon = FK.tfe->ndfOn();

      for (int k = 0; k < nt; ++k) {
         for (int i = 0; i < FK.NbNode(); i++) {
            int onWhat = NodeOnWhat[i];
            int nb     = ndfon[onWhat];

            pp[p[nn++]] = nb;
         }
      }

      for (int n = 0; n < nv; ++n) {
         int ndfn = pp[n];
         pp[n]    = kk;
         kk += ndfn;
      }

      pp[nv] = ndof;
      assert(ndof == kk);
   }

   return DataFENodeDF(ndfon, nt, nv, ndof, nodeOfElement, firstDfOfNode,
                       maxNodePerElement, maxDFPerElement, this->MaxNbDFPerNode,
                       this->constDfPerNode);
}

template <typename Mesh> class CutFESpace : public GFESpace<Mesh> {
 public:
   typedef GFElement<Mesh> FElement;
   typedef typename Mesh::Element Element;
   typedef typename Mesh::BorderElement BorderElement;

   const ActiveMesh<Mesh> &cutTh;
   const GFESpace<Mesh> &backSpace;

   CutFESpace(const ActiveMesh<Mesh> &TTh, const GFESpace<Mesh> &vh)
       : GFESpace<Mesh>(TTh, vh), cutTh(TTh), backSpace(vh) {}

   FElement operator[](int k) const {
      int kb = cutTh.idxElementInBackMesh(k);
      return FElement(this, k, kb);
   }
   const ActiveMesh<Mesh> &get_mesh() const { return cutTh; }
   const GFESpace<Mesh> &get_back_space() const { return backSpace; }
   int get_domain(int k) const { return cutTh.get_domain_element(k); }
   // bool isCut(int k) const { return cutTh.isCut(k);}
   int idxElementInBackMesh(int k) const {
      return cutTh.idxElementInBackMesh(k);
   }
   int idxElementFromBackMesh(int k) const {
      return cutTh.idxElementFromBackMesh(k);
   }
   int idxElementFromBackMesh(int k, int i) const {
      return cutTh.idxElementFromBackMesh(k, i);
   }
   std::vector<int> idxAllElementFromBackMesh(int k, int d) const {
      return cutTh.idxAllElementFromBackMesh(k, d);
   }

   ~CutFESpace() {}
   // const Interface<Mesh>& get_interface(int k) const {
   // cutTh.getInterface(k);}

   // virtual bool isCut(int k) const { return false;}
   // virtual bool isCut() const { return (this->gamma.size()>0);}
   // virtual int whichDomain(int k) const { return -1;}
   // virtual int getNumberOfSubDomain() const { return 1;}
   // virtual int getNeighborElement(int k,int &j, int domain = 0) const {
   // return Th.ElementAdj(k,j);} virtual bool containBackElement(int k)const
   // {return true;} virtual int nbDomain() const {return 1;} virtual bool
   // isCutSpace() const {return false;}
};

template <typename Mesh>
DataFENodeDF GFESpace<Mesh>::BuildDFNumbering(const Mesh &TTh, int dfon[4],
                                              int nndon[4], int N) {

   int *p = 0, *pp = 0;
   bool constndfPerNode  = false;
   int maxNodePerElement = 0;
   int maxDFPerElement   = 0;
   int nbNodes = 0, nbOfDF = 0;
   unsigned int tinfty = -1;
   int nbb             = 0;

   const int nk[] = {Element::nv, Element::ne, Element::nf, Element::nt};
   int nbNodeInK  = nk[0] * nndon[0] + nk[1] * nndon[1] + nk[2] * nndon[2] +
                   nk[3] * nndon[3];
   int keysdim[nbNodeInK];
   int mindf   = 1000;
   int maxdf   = 0;
   int nbnzero = 0;
   for (int dd = 0; dd < 4; ++dd) {
      if (dfon[dd]) {
         nbnzero++;
         mindf = std::min(mindf, dfon[dd]);
         maxdf = std::max(maxdf, dfon[dd]);
         maxDFPerElement += dfon[dd] * nk[dd];
         maxNodePerElement += nndon[dd] * nk[dd];
      }
   }

   if (mindf == maxdf)
      constndfPerNode = true;
   bool nodearevertices = (nbnzero == 1 && dfon[0]);

   if (nodearevertices) {
      nbNodes = TTh.nv;
      nbOfDF  = nbNodes * dfon[0];
   } else {

      p = new int[nbNodeInK * TTh.nt];

      int nodeCounter = 0.;
      std::vector<int> numVertex(TTh.nv, -1);
      std::vector<int> numVol(TTh.nt, -1);

      auto comp = [](const std::array<int, 2> &a, const std::array<int, 2> &b) {
         for (int j = 0; j < 2; j++) {
            if (a[j] > b[j])
               return false;
            if (a[j] < b[j])
               return true;
         }
         return false;
      };
      std::map<std::array<int, 2>, int, decltype(comp)> edge(comp);
      std::map<std::array<int, 2>, int, decltype(comp)> face(comp);

      for (int k = 0; k < TTh.nt; ++k) {
         const Element &K(TTh[k]);
         int ii = 0;

         if (dfon[0] > 0) {
            for (int i = 0; i < Element::nv; ++i) {
               keysdim[ii++] = 0;
               int idx       = TTh(k, i);
               if (numVertex[idx] == -1) {
                  numVertex[idx] = nbNodes;
                  nbOfDF += dfon[0];
                  nbNodes += nndon[0];
               }
               for (int j = 0; j < nndon[0]; ++j) {
                  p[nodeCounter++] = numVertex[idx];
               }
            }
         }

         if (dfon[1] > 0) {
            std::array<int, 2> id_e;
            for (int i = 0; i < Element::ne; ++i) {
               keysdim[ii++] = 1;
               id_e[0]       = TTh(K[Element::nvedge[i][0]]);
               id_e[1]       = TTh(K[Element::nvedge[i][1]]);
               std::sort(id_e.begin(), id_e.end());
               const auto &it = edge.find(id_e);
               int num_node;
               if (it == edge.end()) {
                  num_node   = nbNodes;
                  edge[id_e] = nbNodes;
                  nbOfDF += dfon[1];
                  nbNodes += nndon[1];
               } else {
                  num_node = it->second;
               }
               for (int j = 0; j < nndon[1]; ++j)
                  p[nodeCounter++] = num_node;
            }
         }

         if (dfon[2] > 0) {
            std::array<int, 2> id_e;
            for (int iii, i = 0; i < Element::nf; ++i) {
               keysdim[ii++] = 2;
               if (Element::nf == 1) {
                  id_e[0] = k;
                  id_e[1] = tinfty;
               } else {
                  int kAdj = TTh.ElementAdj(k, iii = i);
                  int kn   = (kAdj == -1) ? --nbb : kAdj;
                  id_e[0]  = k;
                  id_e[1]  = kn;
                  std::sort(id_e.begin(), id_e.end());
               }
               const auto &it = face.find(id_e);
               int num_node;
               if (it == face.end()) {
                  num_node   = nbNodes;
                  face[id_e] = nbNodes;
                  nbOfDF += dfon[2];
                  nbNodes += nndon[2];
               } else {
                  num_node = it->second;
               }
               for (int j = 0; j < nndon[2]; ++j)
                  p[nodeCounter++] = num_node;
            }
         }

         if (dfon[3] > 0) {
            keysdim[ii++] = 3;
            if (numVol[k] == -1) {
               numVol[k] = nbNodes;
               nbOfDF += dfon[3];
               nbNodes += nndon[3];
            }
            for (int j = 0; j < nndon[3]; ++j) {
               p[nodeCounter++] = numVol[k];
            }
         }
      }

      if (!constndfPerNode) {
         pp = new int[nbNodes + 1];

         int kk = 0, nn = 0;
         for (int k = 0; k < TTh.nt; ++k) {
            for (int i = 0; i < nbNodeInK; i++) {
               pp[p[nn++]] = dfon[keysdim[i]];
            }
         }
         for (int n = 0; n < nbNodes; ++n) {
            int ndfn = pp[n];
            pp[n]    = kk;
            kk += ndfn;
         }
         pp[nbNodes] = nbOfDF;
         assert(kk == nbOfDF);
      }
   }
   return DataFENodeDF(dfon, TTh.nt, nbNodes, nbOfDF, p, pp, maxNodePerElement,
                       maxDFPerElement, constndfPerNode);
}

typedef GFESpace<Mesh1> FESpace1;
typedef GFESpace<Mesh2> FESpace2;
typedef GFESpace<MeshQuad2> FESpaceQ2;
typedef GFESpace<MeshHexa> FESpaceQ3;
typedef CutFESpace<Mesh1> CutFESpaceT1;
typedef CutFESpace<Mesh2> CutFESpaceT2;
typedef CutFESpace<MeshQuad2> CutFESpaceQ2;
typedef CutFESpace<MeshHexa> CutFESpaceQ3;
typedef CutFESpace<Mesh3> CutFESpaceT3;
typedef GFESpace<Mesh3> FESpace3;
typedef GFElement<Mesh1> TimeSlab;
typedef GFElement<Mesh2> FElement2;
typedef GFElement<Mesh3> FElement3;

template <int d> struct typeMesh {
   typedef Mesh2 Mesh;
};
template <> struct typeMesh<1> {
   typedef Mesh1 Mesh;
};
template <> struct typeMesh<2> {
   typedef Mesh2 Mesh;
};
template <> struct typeMesh<3> {
   typedef Mesh3 Mesh;
};

#endif
