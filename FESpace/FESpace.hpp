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
    RdHat Phat(tfe->Pt_Pi_h(i)); return Phat;
  }
  Rd Pt(int i) const { assert(i < tfe->NbPtforInterpolation);
    RdHat Phat(tfe->Pt_Pi_h(i));
    return T(Phat);
  }

  R getMeasure() const {return T.mesure();} // mesure felstavat (french spel ;-) ) 
  R get_measure() const {return T.mesure();} // mesure felstavat (french spel ;-) ) 

  Rd map(Rd ip) const {return T(ip);}
  Rd mapToPhysicalElement(Rd ip) const {return T(ip);}


  int degre() const {return tfe->degre();}
  int index() const {return number;};

  int whichDomain() const {return Vh.whichDomain(number);}
  int get_domain() const {return Vh.get_domain(number);}
  // bool isCut() const {return Vh.isCut(number);}
  // const Interface<Mesh>& getInterface() const { Vh.get_interface(number);}
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
GFElement<Mesh>::GFElement(const GFESpace<Mesh> * VVh, int k)
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
  typedef typename Element::Face Face;
  typedef typename Mesh::BorderElement  BorderElement;
  typedef typename Mesh::Rd  Rd;
  typedef GTypeOfFE<Mesh> TypeOfFE;
  typedef GQuadratureFormular<typename Element::RdHat> QFElement;
  typedef GQuadratureFormular<typename BorderElement::RdHat>  QFBorderElement;
  // typedef GenericInterface<Mesh> GInterface;


  const Mesh &Th;
  KN<const GTypeOfFE<Mesh> *>  TFE;
  const int N;                     // dim espace d'arrive
  const int Nproduit; // 1 if non constant Max number df par node. else Max number df par node..
  // GFESpace const * backSpace = this;
  const PeriodicBC* periodicBC = nullptr;
  const BasisFctType basisFctType;
  const int polynomialOrder;
  // KN<const GInterface*> gamma;
  // KN<const GInterface*>* gamma2;


  GFESpace(const Mesh & TTh,
    const GTypeOfFE<Mesh> & tfe,//=DataFE<Mesh>::P1,
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
    periodicBC(PPeriod),
    basisFctType(tfe.basisFctType),
    polynomialOrder(tfe.polynomialOrder)
    {
    }




    // GFESpace(const Mesh & TTh,
    //   // TimeInterface<Mesh>& g,
    //   const GTypeOfFE<Mesh> & tfe=DataFE<Mesh>::P1,
    //   const PeriodicBC* PPeriod = nullptr) :
    //   GFESpace(TTh, tfe, PPeriod)
    //   {
    //     // this->gamma.resize(g.size());
    //     // for(int i=0;i<g.size();++i) this->gamma[i]= g[i];
    //   }
    //   GFESpace(const Mesh & TTh,
    //     // const GInterface& g,
    //     const GTypeOfFE<Mesh> & tfe=DataFE<Mesh>::P1,
    //     const PeriodicBC* PPeriod = nullptr) :
    //     GFESpace(TTh, tfe, PPeriod)
    //     {
    //       // this->gamma.resize(1);
    //       // this->gamma[0]= &g;
    //     }

      // for CutSpace
  // GFESpace( const GFESpace& vh, const DataFENodeDF& data,const PeriodicBC* PPeriod = nullptr) :
  //   DataFENodeDF(data),
  //   Th(vh.Th),
  //   TFE(1,0,vh.TFE(0)),
  //   N(TFE[0]->N),
  //   Nproduit(FirstDfOfNodeData ? 1 :MaxNbDFPerNode ),
  //   periodicBC(PPeriod),
  //   basisFctType(vh.basisFctType),
  //   polynomialOrder(vh.polynomialOrder)
  //   {
  //   }

  // GFESpace(const ActiveMesh<Mesh> & TTh, const GTypeOfFE<Mesh> & tfe, const PeriodicBC* PPeriod = nullptr)
  //   :
  //   DataFENodeDF(TTh.Th.BuildDFNumbering(tfe.ndfonVertex,tfe.ndfonEdge,
  //     tfe.ndfonFace,tfe.ndfonVolume,
  //     tfe.nbNodeOnWhat[0],
  //     tfe.nbNodeOnWhat[1],
  //     tfe.nbNodeOnWhat[2],
  //     tfe.nbNodeOnWhat[3],
  //     tfe.N,
  //     PPeriod
  //   )),
  //   Th(TTh.Th),
  //   TFE(1,0,&tfe),
  //   N(tfe.N),
  //   Nproduit(FirstDfOfNodeData ? 1 :MaxNbDFPerNode ),
  //   periodicBC(PPeriod),
  //   basisFctType(tfe.basisFctType),
  //   polynomialOrder(tfe.polynomialOrder)
  //   {
  //   }


    GFESpace(const ActiveMesh<Mesh> & TTh, const GFESpace& vh, const PeriodicBC* PPeriod = nullptr) :
    DataFENodeDF(vh.BuildDFNumbering(TTh)),
    Th(TTh.Th),
    TFE(1,0,vh.TFE(0)),
    N(TFE[0]->N),
    Nproduit(FirstDfOfNodeData ? 1 :MaxNbDFPerNode ),
    periodicBC(PPeriod),
    basisFctType(vh.basisFctType),
    polynomialOrder(vh.polynomialOrder)
    {
    }
    DataFENodeDF BuildDFNumbering(const ActiveMesh<Mesh> & TTh) const;


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
  // virtual int idxElementFromBackMesh (int k)      const {return Th.idxElementFromBackMesh(k) ;}
  // virtual int idxElementInBackMesh(int k)         const {return Th.idxElementInBackMesh(k);}
  // virtual int idxElementFromBackMesh(int k,int i) const {return idxElementFromBackMesh(k); }
  virtual int idxElementFromBackMesh (int k)      const {return k ;}
  virtual int idxElementInBackMesh(int k)         const {return k;}
  virtual int idxElementFromBackMesh(int k,int i) const {return k; }
  // virtual bool faceInDomain(const Face& face, int dom) const {return true;}
  // const Face& get_face(int k) const {return Th.hyper_face(k);}
  virtual int whichDomain(int k) const { return -1;}
  virtual int get_domain(int k) const { return -1;}
  virtual int getNeighborElement(int k,int &j, int domain = 0) const { return Th.ElementAdj(k,j);}
  virtual int nbDomain() const {return 1;}
  virtual bool containBackElement(int k)const {return true;}
  // virtual bool isCut(int k) const { return false;}
  // virtual bool isCut() const { return (this->gamma.size()>0);}
  virtual bool isCutSpace() const {return false;}
  virtual vector<int> idxAllElementFromBackMesh (int k, int d) const { assert(d==-1);vector<int> v = {idxElementFromBackMesh(k)}; return v ;}
  virtual const GFESpace<Mesh>& get_back_space() const {return *this;}
  virtual const ActiveMesh<Mesh>& get_mesh() const { assert(0); return ActiveMesh<Mesh>(Th);}



  int NbNode() const { return this->nbNode;}
  int NbDoF() const { return this->nbDoF;}
  int NbElement() const { return this->nbElement;}
  int get_nb_element() const { return this->nbElement;}
  int get_nb_dof() const { return this->nbDoF;}

  // int NbInnerFaces() const { return Th.nbInnerFaces();}
  #ifdef USE_MPI
  virtual int first_element() const { return MPIcf::first_element(this->nbElement);}
  virtual int next_element() const {  return MPIcf::next_element(this->nbElement);}
  virtual int last_element() const {  return MPIcf::last_element(this->nbElement);}

  // virtual int first_face() const { return MPIcf::first_element(this->NbInnerFaces());}
  // virtual int next_face() const {  return MPIcf::next_element(this->NbInnerFaces());}
  // virtual int last_face() const {  return MPIcf::last_element(this->NbInnerFaces());}

  virtual int first_boundary_element() const { return MPIcf::my_rank();}
  virtual int next_boundary_element() const { return MPIcf::size();}
  virtual int last_boundary_element() const {return this->Th.nbBrdElmts();}
  #else
  virtual int first_element() const { return 0;}
  virtual int next_element() const {return 1;}
  virtual int last_element() const { return this->nbElement;}

  // virtual int first_face() const { return 0;}
  // virtual int next_face() const {return 1;}
  // virtual int last_face() const { return this->NbInnerFaces();}

  virtual int first_boundary_element() const { return 0;}
  virtual int next_boundary_element() const { return 1;}
  virtual int last_boundary_element() const {return this->Th.nbBrdElmts();}
  #endif




  bool isPeriodic(int lab)const {
    return (periodicBC)? periodicBC->isPeriodic(lab) : false;}


  virtual void info() const {
    // std::cout << "FESpace \t" << this << std::endl;
    // std::cout << "----------------------------------" << std::endl;
    std::cout << " nb node    \t" << NbNode() << std::endl;
    std::cout << " nb dof     \t" << NbDoF() << std::endl;
    std::cout << " nb element \t" << NbElement() << std::endl;
  }



};

template<typename Mesh>
DataFENodeDF GFESpace<Mesh>::BuildDFNumbering(const ActiveMesh<Mesh> & TTh) const{

int ndfon[4] = {this->ndfonVertex(), this->ndfonEdge(),this->ndfonFace(), this->ndfonTet()};
int nSub = TTh.get_nb_domain();

int nt = 0, nv = 0, ndof = 0;
int maxNodePerElement = this->MaxNbNodePerElement;
int maxDFPerElement   = this->MaxNbDFPerElement;

// to know where are the nodes
KN<int> NodeOnWhat(maxNodePerElement);
{
  int c= 0;
  const int nk[]={Element::nv,Element::ne,Element::nf,Element::nt};
  for(int i=0;i<4;++i){
    for(int k=0;k<nk[i];++k) {
      for(int j=0;j<this->TFE(0)->nbNodeOnWhat[i];++j) {
        NodeOnWhat(c++) = i;
      }
    }
  }
}

// p => idx nodes of elements : Vh.nbElement*NodesPerElement
// node are not the vertex of the mesh but node for the finite element
// they are build in the Space so we can just use this for
// the cut space
int* nodeOfElement  = new int[maxNodePerElement*TTh.get_nb_element()];
int* firstDfOfNode = nullptr;
// nodeOfElement = new int(maxNodePerElement*TTh.get_nb_element());
// get number of nodes and dof
{
  // for(int i=0;i<nSub;++i) nts(i) = TTh.get_nb_element(i);
  nt = TTh.get_nb_element();
  std::vector<std::vector<int>>  idxP_globToLoc(nSub);
  for(int i=0;i<nSub;++i) idxP_globToLoc[i].assign(this->nbNode, -1);
  int c = 0;
  for(int k=0;k<TTh.get_nb_element();++k){

    int kb = TTh.idxElementInBackMesh(k);
    int dom = TTh.get_domain_element(k);
    int nbNode = (*this)[0].NbNode();

    for( int j = 0; j < nbNode; ++j) {
      int iglb =(*this)(kb, j);

      int val = idxP_globToLoc[dom][iglb];
      nodeOfElement[c++] = (val ==-1)? nv: val;
      if(val ==-1) {
        idxP_globToLoc[dom][iglb] = nv++;
        ndof += this->TFE(0)->ndfOn()[NodeOnWhat(j)];
      }
    }
  }
}

if(!this->constDfPerNode) {
  int kk=0, nn=0;
  int* p = nodeOfElement;
  firstDfOfNode = new int[nv+1];
  int* pp= firstDfOfNode;

  const FElement& FK((*this)[0]);
  const int * ndfon = FK.tfe->ndfOn();

  for(int k=0; k<nt; ++k) {
    for(int i=0; i<FK.NbNode(); i++) {
      int onWhat = NodeOnWhat[i];
      int nb = ndfon[onWhat];

      pp[p[nn++]] = nb;
    }
  }

  for(int n=0; n<nv; ++n) {
    int ndfn=pp[n];
    pp[n] = kk;
    kk += ndfn;
  }

  pp[nv] = ndof;
  assert(ndof == kk);
}

return DataFENodeDF(ndfon,
        nt, nv, ndof,
        nodeOfElement, firstDfOfNode,
        maxNodePerElement, maxDFPerElement,
        this->MaxNbDFPerNode,
        this->constDfPerNode);
}


template<typename Mesh>
class CutFESpace : public GFESpace<Mesh> {
public:
  typedef GFElement<Mesh> FElement;
  typedef typename Mesh::Element  Element;
  typedef typename Mesh::BorderElement  BorderElement;

  const ActiveMesh<Mesh> & cutTh;
  const GFESpace<Mesh>& backSpace;

  CutFESpace(const ActiveMesh<Mesh> & TTh, const GFESpace<Mesh>& vh, const PeriodicBC* PPeriod = nullptr) : GFESpace<Mesh>(TTh, vh, PPeriod), cutTh(TTh), backSpace(vh) {}

  FElement operator[](int k) const {
    int kb = cutTh.idxElementInBackMesh(k);
    return FElement(this,k, kb);
  }
  const ActiveMesh<Mesh>& get_mesh() const { return cutTh;}
  const GFESpace<Mesh>& get_back_space() const {return backSpace;}
  int get_domain(int k) const { return cutTh.get_domain_element(k);}
  // bool isCut(int k) const { return cutTh.isCut(k);}
  int idxElementInBackMesh(int k) const { return cutTh.idxElementInBackMesh(k);}
  int idxElementFromBackMesh (int k) const { return cutTh.idxElementFromBackMesh(k) ;}
  int idxElementFromBackMesh(int k,int i) const { return cutTh.idxElementFromBackMesh(k,i); }
  vector<int> idxAllElementFromBackMesh (int k, int d) const { return cutTh.idxAllElementFromBackMesh(k, d) ;}


  // const Interface<Mesh>& get_interface(int k) const { cutTh.getInterface(k);}

  // virtual bool isCut(int k) const { return false;}
  // virtual bool isCut() const { return (this->gamma.size()>0);}
  // virtual int whichDomain(int k) const { return -1;}
  // virtual int getNumberOfSubDomain() const { return 1;}
  // virtual int getNeighborElement(int k,int &j, int domain = 0) const { return Th.ElementAdj(k,j);}
  // virtual bool containBackElement(int k)const {return true;}
  // virtual int nbDomain() const {return 1;}
  // virtual bool isCutSpace() const {return false;}
};


typedef GFESpace<Mesh1>     FESpace1;
typedef GFESpace<Mesh2>     FESpace2;
typedef GFESpace<MeshQuad2> FESpaceQ2;
typedef GFESpace<MeshHexa>  FESpaceQ3;
typedef CutFESpace<Mesh1>     CutFESpaceT1;
typedef CutFESpace<Mesh2>     CutFESpaceT2;
typedef CutFESpace<MeshQuad2> CutFESpaceQ2;
typedef CutFESpace<MeshHexa>  CutFESpaceQ3;
typedef CutFESpace<Mesh3>     CutFESpaceT3;
typedef GFESpace<Mesh3> FESpace3;
typedef GFElement<Mesh1> TimeSlab;
typedef GFElement<Mesh2> FElement2;
typedef GFElement<Mesh3> FElement3;


template<int d> struct typeMesh {typedef Mesh2 Mesh;};
template<> struct typeMesh<1> {typedef Mesh1 Mesh;};
template<> struct typeMesh<2> {typedef Mesh2 Mesh;};
template<> struct typeMesh<3> {typedef Mesh3 Mesh;};


#endif
