#ifndef FESpace_HPP_
#define FESpace_HPP_


/*
 *  FESpace.hpp
 *
 */
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
using namespace std;
#include "../util/error.hpp"
#include "../util/util.hpp"
#include "../util/ufunction.hpp"

#include "interpolationMatrix.hpp"
#include "GTypeOfFE_Sum.hpp"
// #include "GTypeOfFE_Rd.hpp"

#include "../common/Mesh3dn.hpp"
#include "../common/Mesh2dn.hpp"
#include "../common/Mesh1dn.hpp"
#include "../common/Interface3dn.hpp"
#include "../common/timeInterface.hpp"

#include "QuadratureFormular.hpp"

#include "cutFEMConfig.h"
#ifdef USE_MPI
#include "../parallel/cfmpi.hpp"
#endif

template<class Mesh> class GFESpace;
template<class Mesh> class GFElement;
template<class Mesh> class GbaseFElement;


/*
 *  Basic class for Element of the FESpace
 *
 */
template<class MMesh>
class GbaseFElement
{
public:
  typedef MMesh  Mesh;
  typedef GFESpace<Mesh>  FESpace;
  typedef typename Mesh::Element Element;
  typedef  Element E;
  typedef typename E::Rd Rd;
  typedef typename E::RdHat RdHat;
  typedef typename E::RdHatBord RdHatBord;
  typedef GQuadratureFormular<RdHat>  QF;
  typedef GQuadratureFormular<RdHatBord>  QFB;

  const GFESpace<Mesh> &Vh;
  const Element &T;
  const GTypeOfFE<Mesh> * tfe;
  const int N;                             // dim arrived space
  const int number;                        // my index

  GbaseFElement(const GFESpace<Mesh> &aVh, int k) ;
  GbaseFElement(const GFESpace<Mesh> &aVh, int k, int kth) ;

  R EdgeOrientation(int i) const {return T.EdgeOrientation(i);}

  Rd Pt(RdHat Phat) const {return T(Phat);}     // Ref to Global
  Rd PtHat(int i) const { assert(i < tfe->NbPtforInterpolation);
    RdHat Phat(tfe->PtInterpolation(i)); return Phat;}
  Rd Pt(int i) const { assert(i < tfe->NbPtforInterpolation);
    RdHat Phat(tfe->Pt_Pi_h(i));    return T(Phat);}

  int whichDomain() const {return Vh.whichDomain(number);}
  bool isCut() const {return Vh.isCut(number);}
  R getMeasure() const {return T.mesure();} // mesure felstavat (french spel ;-) ) 
  Rd map(Rd ip) const {return T(ip);}
  int degre() const {return tfe->degre();}
  int index() const {return number;};

};

template<class Mesh>
inline GbaseFElement<Mesh>::GbaseFElement(  const GFESpace<Mesh> &aVh, int k)
  : Vh(aVh),T(Vh.Th[k]),tfe(aVh.TFE[k]),N(aVh.N),number(k){
}

// for cutFEM
template<class Mesh>
inline GbaseFElement<Mesh>::GbaseFElement(  const GFESpace<Mesh> &aVh, int k, int kth)
  : Vh(aVh),T(Vh.Th[kth]),tfe(aVh.TFE[0]),N(aVh.N),number(k){
}



/*
 *    Class of the finite element of the FESpace
 *
 */
template<class Mesh>
class GFElement : public GbaseFElement<Mesh>
{
public:

  typedef typename Mesh::Element Element;
  typedef  Element E;
  typedef typename E::Rd Rd;
  typedef typename E::RdHat RdHat;
  typedef typename E::RdHatBord RdHatBord;
  typedef GQuadratureFormular<RdHat>  QF;
  typedef GQuadratureFormular<RdHatBord>  QFB;

  friend class GFESpace<Mesh>;
  const int *p;                   // indices nodes
  const int nb;                   // nb of nodes

  GFElement(const GFESpace<Mesh> * VVh,int k) ;
  GFElement(const GFESpace<Mesh> * VVh,int k, int kth) ;


  int NbDoF() const { return this->tfe->nbDoF;}
  int NbNode()const {return nb;}
  int  operator[](int i) const ;       //{ return  p[i]  Numdu noeud
  int  operator()(int i,int df) const ;// Num du DoF du noeud i de df local df
  int  operator()(int df) const ;
  int loc2glb(int i) const {return (*this)(i);}
  int loc2glb(int i, int it) const {return (*this)(i) + it*this->Vh.NbDoF();}


  void BF(const Rd & P,RNMK_ & val) const;
  void BF(const What_d whatd, const Rd & P,RNMK_ & val) const;
  void BF(const What_d whatd, const Rd & P,RNMK_ & val, const RNM_& J) const;

  void set(InterpolationMatrix<RdHat> &M) const {
    this->tfe->set(this->Vh.Th,this->T,M,0,0,0);}

  KN_<R>& Pi_h(KNM_<R> vpt, RN_ & vdf) const ;

  int DFOnWhat(int i) const { return this->tfe->DFOnWhat[i];}

  // df is the df in element
  int NodeOfDF(int df) const { return this->tfe->NodeOfDF[df];} // a node
  int DFOfNode(int df) const { return this->tfe->DFOfNode[df];} // the df of the node

  int dfcbegin(int ic) const { return this->tfe->begin_dfcomp[ic];}
  int dfcend(int ic) const { return this->tfe->end_dfcomp[ic];}

};


template<class Mesh>
GFElement<Mesh>::GFElement(const GFESpace<Mesh> * VVh,int k)
  : GbaseFElement<Mesh>(*VVh,k) ,
    p(this->Vh.PtrFirstNodeOfElement(k)),
    nb(this->Vh.NbOfNodesInElement(k))
{}
template<class Mesh>
GFElement<Mesh>::GFElement(const GFESpace<Mesh> * VVh,int k, int kth)
  : GbaseFElement<Mesh>(*VVh,k,kth) ,
    p(this->Vh.PtrFirstNodeOfElement(k)),
    nb(this->Vh.NbOfNodesInElement(k))
{}


template<class Mesh>
inline int GFElement<Mesh>::operator[](int i) const {
  return  p ? p[i] :  ((&this->T[i])-this->Vh.Th.vertices);}

template<class Mesh>
inline int  GFElement<Mesh>::operator()(int i,int df) const {
  return this->Vh.FirstDFOfNode(p ? p[i] :  ((&this->T[i])-this->Vh.Th.vertices)) + df;
}

template<class Mesh>
  inline int  GFElement<Mesh>::operator()(int df) const {
  return operator()(NodeOfDF(df),DFOfNode(df));
}


template<class Mesh>
inline void GFElement<Mesh>::BF(const Rd & P,RNMK_ & val) const {
  this->tfe->FB(Fop_D0|Fop_D1,this->T,P,val);}

template<class Mesh>
inline void GFElement<Mesh>::BF(const What_d whatd,const Rd & P,RNMK_ & val) const {
  this->tfe->FB(whatd,this->T,P,val);}

template<class Mesh>
inline void GFElement<Mesh>::BF(const What_d whatd,const Rd & P,RNMK_ & val, const RNM_& J) const {
  this->tfe->FB(whatd,this->T,P,val,J);}



template<class Mesh>
KN_<R> & GFElement<Mesh>::Pi_h(KNM_<R> vpt,KN_<R> & vdf)    const
  {
    // compute  the interpolation
    // in : vpt  value of componant j at point p : vpt(p,j)
    // out: vdf  value du the degre of freedom
    const KN<IPJ> ipj(this->tfe->ipj_Pi_h);
    const KN<Rd>  PtHat(this->tfe->Pt_Pi_h);
    KN<R>   Aipj(ipj.N());
    this->tfe->get_Coef_Pi_h(*this, Aipj);

    vdf=R();
    for (int i=0;i<Aipj.N();i++) {
      const IPJ &ipj_i(ipj[i]);
      vdf[ipj_i.i] += Aipj[i] * vpt(ipj_i.p, ipj_i.j);
    }
    return  vdf;
  }



/*
 *   Class for the FESpace
 *
 */
template<class MMesh>
class GFESpace : public DataFENodeDF {
public:
  typedef MMesh Mesh;
  typedef GFElement<Mesh> FElement;
  typedef typename Mesh::Element  Element;
  typedef typename Mesh::BorderElement  BorderElement;
  typedef typename Mesh::Rd  Rd;
  typedef GTypeOfFE<Mesh> TypeOfFE;
  typedef GQuadratureFormular<typename Element::RdHat> QFElement;
  typedef GQuadratureFormular<typename BorderElement::RdHat>  QFBorderElement;
  typedef GenericInterface<Mesh> Interface;


  const Mesh &Th;
  KN<const GTypeOfFE<Mesh> *>  TFE;
  const int N;                     // dim espace d'arrive
  const int Nproduit; // 1 if non constant Max number df par node. else Max number df par node..
  GFESpace const * backSpace = this;
  const PeriodicBC* periodicBC = nullptr;
  KN<const Interface*> gamma;
  KN<const Interface*>* gamma2;


  GFESpace(const Mesh & TTh,
	   const GTypeOfFE<Mesh> & tfe=DataFE<Mesh>::P1,
	   const PeriodicBC* PPeriod = nullptr)
    //	   int nbPeriodicBe = 0, int* periodicBe = 0)
    :
    DataFENodeDF(TTh.BuildDFNumbering(tfe.ndfonVertex,tfe.ndfonEdge,
				      tfe.ndfonFace,tfe.ndfonVolume,
				      tfe.nbNodeOnWhat[0],
				      tfe.nbNodeOnWhat[1],
				      tfe.nbNodeOnWhat[2],
				      tfe.nbNodeOnWhat[3],
				      tfe.N,
				      PPeriod
				      )),
    Th(TTh),
    TFE(1,0,&tfe),
    N(tfe.N),
    Nproduit(FirstDfOfNodeData ? 1 :MaxNbDFPerNode ),
    periodicBC(PPeriod)
  {

  }

    GFESpace(const Mesh & TTh,
      TimeInterface<Mesh>& g,
      const GTypeOfFE<Mesh> & tfe=DataFE<Mesh>::P1,
      const PeriodicBC* PPeriod = nullptr) :
      GFESpace(TTh, tfe, PPeriod)
      {
        this->gamma.resize(g.size());
        for(int i=0;i<g.size();++i) this->gamma[i]= g[i];
      }
  GFESpace(const Mesh & TTh,
    const Interface& g,
    const GTypeOfFE<Mesh> & tfe=DataFE<Mesh>::P1,
    const PeriodicBC* PPeriod = nullptr) :
    GFESpace(TTh, tfe, PPeriod)
     {
      this->gamma.resize(1);
      this->gamma[0]= &g;
    }

      // for CutSpace
  GFESpace( const GFESpace& vh, const DataFENodeDF& data,const PeriodicBC* PPeriod = nullptr) :
    DataFENodeDF(data),
    Th(vh.Th),
    TFE(1,0,vh.TFE(0)),
    N(TFE[0]->N),
    Nproduit(FirstDfOfNodeData ? 1 :MaxNbDFPerNode ),
    periodicBC(PPeriod)
    {
    }

  const int * PtrFirstNodeOfElement(int k) const {
    return NodesOfElement  ?
      NodesOfElement + (FirstNodeOfElement ? FirstNodeOfElement[k] : k*MaxNbNodePerElement)
      : 0;}

  int NbOfNodesInElement(int k)   const {
    return FirstNodeOfElement // 0 if const dof per element
      ?  FirstNodeOfElement[k+1] - FirstNodeOfElement[k] :  MaxNbNodePerElement ;}


  int FirstDFOfNode(int i) const {
    return (FirstDfOfNodeData ? FirstDfOfNodeData[i] : i*Nproduit) ;}
  // int LastDFOfNode(int i)  const {return FirstDfOfNodeData ? FirstDfOfNodeData[i+1] : (i+1)*Nproduit;}
  // int NbDFOfNode(int i)  const {return FirstDfOfNodeData ? FirstDfOfNodeData[i+1]-FirstDfOfNodeData[i] : Nproduit;}

  ~GFESpace(){
  }

  virtual FElement operator[](int k) const { return FElement(this,k);}
  FElement operator[](const Element & K) const { return FElement(this,Th(K));}


  int  operator()(int k)const {return NbOfNodesInElement(k);}
  int  operator()(int k,int i) const { //  the node i of element k
    return NodesOfElement ?  *(PtrFirstNodeOfElement(k) + i)  : Th(k,i)  ;}

  // for mesh build from interface but works for all mesh
  virtual int idxElementFromBackMesh (int k) const { return Th.idxElementFromBackMesh(k) ;}
  virtual int idxElementInBackMesh(int k) const {
    return Th.idxElementInBackMesh(k);}

  virtual int idxElementFromBackMesh(int k,int i) const {
    return idxElementFromBackMesh(k); }

  virtual bool isCut(int k) const { return false;}
  virtual bool isCut() const { return (this->gamma.size()>0);}

  virtual int whichDomain(int k) const { return -1;}

  virtual int getNumberOfSubDomain() const { return 1;}

  virtual int getNeighborElement(int k,int &j, int domain = 0) const {
    return Th.ElementAdj(k,j);
  }
  virtual bool containBackElement(int k)const {
    return true;
  }
  virtual int nbDomain() const {
    return 1;
  }

  // virtual int idxGlob2Loc(int k, int i) const { return k;}
  virtual const GFESpace& getBackSpace() const { return *backSpace;}
  // virtual const  GenericInterface<Mesh>& getInterface(int i) const {assert(0);}
  const Interface& getInterface(int i) const {assert(this->gamma.size() > 0);assert(i<this->gamma.size()); return *this->gamma(i);}

  virtual bool isCutSpace() const {return false;}

  int NbNode() const { return this->nbNode;}
  int NbDoF() const { return this->nbDoF;}
  int NbElement() const { return this->nbElement;}
  #ifdef USE_MPI
  virtual int first_element() const { return MPIcf::first_element(this->nbElement);}
  virtual int next_element() const {  return MPIcf::next_element(this->nbElement);}
  virtual int last_element() const {  return MPIcf::last_element(this->nbElement);}

  virtual int first_boundary_element() const { return MPIcf::my_rank();}
  virtual int next_boundary_element() const { return MPIcf::size();}
  virtual int last_boundary_element() const {return this->Th.nbe;}
  #else
  virtual int first_element() const { return 0;}
  virtual int next_element() const {return 1;}
  virtual int last_element() const { return this->nbElement;}

  virtual int first_boundary_element() const { return 0;}
  virtual int next_boundary_element() const { return 1;}
  virtual int last_boundary_element() const {return this->Th.nbe;}
  #endif




  bool isPeriodic(int lab)const {
    return (periodicBC)? periodicBC->isPeriodic(lab) : false;}


  virtual void info() const {
    std::cout << "FESpace \t" << this << std::endl;
    std::cout << "nb node    \t" << NbNode() << std::endl;
    std::cout << "nb dof     \t" << NbDoF() << std::endl;
    std::cout << "nb element \t" << NbElement() << std::endl;
  }



};







typedef GFESpace<Mesh1> FESpace1;
typedef GFESpace<Mesh2> FESpace2;
typedef GFESpace<Mesh3> FESpace3;
typedef GFElement<Mesh1> TimeSlab;
typedef GFElement<Mesh2> FElement2;
typedef GFElement<Mesh3> FElement3;


template<int d> struct typeMesh {typedef Mesh2 Mesh;};
template<> struct typeMesh<1> {typedef Mesh1 Mesh;};
template<> struct typeMesh<2> {typedef Mesh2 Mesh;};
template<> struct typeMesh<3> {typedef Mesh3 Mesh;};


#endif
