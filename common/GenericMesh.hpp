#ifndef GENERICMESH_HPP_
#define GENERICMESH_HPP_

extern long verbosity;

#include "cassert"
#include "../util/assertion.hpp"
#include "../util/util.hpp"
//#include "../parallel/cfmpi.hpp"
#include <cstdlib>

// #include "RefCounter.hpp"

using namespace ::std;

#include "R3.hpp"
#include "HashTable.hpp"
#include "GenericVertex.hpp"
#include "dataStruct1D.hpp"
#include "dataStruct2D.hpp"
#include "dataStruct3D.hpp"
// #include "cut_method.hpp"
#include "cutFEMConfig.h"
#ifdef USE_MPI
#include "../parallel/cfmpi.hpp"
#endif

inline int maxdfon(const int *dfon){
  return max(max(dfon[0],dfon[1]),max(dfon[2],dfon[3]));}

template<int d> struct typeRd {typedef R0 Rd;};
template<> struct typeRd<1> {typedef R1 Rd;};
template<> struct typeRd<2> {typedef R2 Rd;};
template<> struct typeRd<3> {typedef R3 Rd;};


const int NbTypeItemElement = 4;
const int TypeVertex =0;
const int TypeEdge   =1;
const int TypeFace   =2;
const int TypeVolume =3;



/*
 *  Structure with the FE space information
 *
 */
class DataFENodeDF {
  int  * nbref;                            // pointer  on common ref counter
public:
  int ndfon[4];                            // array with nb of dof per item
  const int nbElement;
  const int nbNode;
  const int nbDoF;
  const int * const NodesOfElement;             // array of indices
  const int * const FirstDfOfNodeData;          // if dof different in elements
  const int * const FirstNodeOfElement;         // 0
  const int MaxNbNodePerElement;
  const int MaxNbDFPerElement;
  const int MaxNbDFPerNode;
  bool constDfPerNode = true;
  int ndfonVertex()const {return ndfon[0];}
  int ndfonEdge()const {return ndfon[1];}
  int ndfonFace()const {return ndfon[2];}
  int ndfonTet()const {return ndfon[3];}

  DataFENodeDF(const DataFENodeDF & m)
    :
    nbref( m.nbref ) ,
    nbElement(m.nbElement),
    nbNode(m.nbNode),
    nbDoF(m.nbDoF),
    NodesOfElement(m.NodesOfElement),
    FirstDfOfNodeData(m.FirstDfOfNodeData),
    FirstNodeOfElement(m.FirstNodeOfElement),
    MaxNbNodePerElement(m.MaxNbNodePerElement),
    MaxNbDFPerElement(m.MaxNbDFPerElement) ,
    MaxNbDFPerNode(maxdfon(m.ndfon)),
    constDfPerNode(m.constDfPerNode)
    {
      for(int i=0;i<NbTypeItemElement;++i)
	ndfon[i]=m.ndfon[i];
      (*nbref)++;                         // add one to the ref counter
    }


  DataFENodeDF(
	       int andfon[NbTypeItemElement],
	       int anbElement,
	       int anbNode,
	       int anbDoF,
	       const int * aNodesOfElement,
	       const int * aFirstDfOfNodeData,
	       int aMaxNbNodePerElement,
	       int aMaxNbDFPerElement,
	       bool cstPerNode = true)
    :
    nbref( new int(0)),// new ref counter
    nbElement(anbElement),
    nbNode(anbNode),
    nbDoF(anbDoF),
    NodesOfElement(aNodesOfElement),
    FirstDfOfNodeData(aFirstDfOfNodeData),
    FirstNodeOfElement(0),
    MaxNbNodePerElement(aMaxNbNodePerElement),
    MaxNbDFPerElement(aMaxNbDFPerElement) ,
    MaxNbDFPerNode(maxdfon(andfon)),
    constDfPerNode(cstPerNode)
  {
    for(int i=0;i<NbTypeItemElement;++i)
      ndfon[i]=andfon[i];
  }

  DataFENodeDF(
	       int andfon[NbTypeItemElement],
	       int anbElement,
	       int anbNode,
	       int anbDoF,
	       const int * aNodesOfElement,
	       const int * aFirstDfOfNodeData,
	       int aMaxNbNodePerElement,
	       int aMaxNbDFPerElement,
	       int aMaxNbDFPerNode,
	       bool cstPerNode = true)
    :
    nbref( new int(0)),// new ref counter
    nbElement(anbElement),
    nbNode(anbNode),
    nbDoF(anbDoF),
    NodesOfElement(aNodesOfElement),
    FirstDfOfNodeData(aFirstDfOfNodeData),
    FirstNodeOfElement(0),
    MaxNbNodePerElement(aMaxNbNodePerElement),
    MaxNbDFPerElement(aMaxNbDFPerElement) ,
    MaxNbDFPerNode(aMaxNbDFPerNode),
    constDfPerNode(cstPerNode)
  {
    for(int i=0;i<NbTypeItemElement;++i)
      ndfon[i]=andfon[i];
  }


 ~DataFENodeDF()
  {
    if ((*nbref)==0) // remove if nbref ==0
       {
	 delete nbref;
         delete [] NodesOfElement;
         delete [] FirstDfOfNodeData;
         delete [] FirstNodeOfElement;
       }
	else  (*nbref)--;
  }
private:
	void operator=(const DataFENodeDF &) ;

};









inline  R1 ExtNormal( GenericVertex<R1> *const v[2],int const f[1])  {
  return f[0]==0 ? R1(-1):R1(1);  }
inline  R2 ExtNormal( GenericVertex<R2> *const v[3],int const f[2])  {
  return R2(*v[f[1]],*v[f[0]]).perp();  }
inline  R3 ExtNormal( GenericVertex<R3> *const v[4],int const f[3])  {
  return R3(*v[f[0]],*v[f[2]])^R3(*v[f[0]],*v[f[1]]) ;  }




template<int N> inline void PermI2J(const void **I,const void **J,int *S)
{
	ffassert(0);
}
template<> inline void PermI2J<1>(const void **I,const void **J,int *S)
{
  S[0]=0;
}
template<> inline void PermI2J<2>(const void **I,const void **J,int *S)
{
  if(I[0]==J[0])
    { assert(I[1]==J[1]);
      S[0]=0;S[1]=1;}
  else
    { assert(I[1]==J[0]&&I[0]==J[1]);
      S[0]=1;S[1]=0;}
}
template<> inline void PermI2J<3>(const void **I,const void **J,int *S)
{
  if(I[0]==J[0]) S[0]=0;
  else if(I[0]==J[1]) S[0]=1;
  else {S[0]=2; assert(I[0]==J[2]) ;}
  if(I[1]==J[0]) S[1]=0;
  else if(I[1]==J[1]) S[1]=1;
  else {S[1]=2; assert(I[1]==J[2]) ; }
  S[2]=3-S[0]-S[1];
  assert(I[2]==J[3-S[0]-S[1]]);
}




// template<typename Mesh>
class PeriodicBC {

public:
  int nbEqui=0;
  KN<int> periodicBE;
  int nbe[6] = {0,0,0,0,0,0};
  bool isP[6] = {false, false, false, false, false, false};

  // lab 2 & 4
  template<typename Mesh>
  PeriodicBC(Mesh& Th) {
    typedef typename Mesh::BorderElement B;

    int idx[6] = {0,0,0,0,0,0};
    periodicBE.init(Th.nbe);

    for(int kb=0;kb<Th.nbe; kb++) {

      const B & BE(Th.be(kb));
      int label = BE.lab;
      nbe[label-1] ++;

      if(label==2) {
	periodicBE(2*idx[1]) = kb;
	idx[1]++;
	isP[1] = true;
      }
      else if(label==4) {
	periodicBE(2*idx[3]+1) = kb;
	idx[3]++;
	isP[3] = true;
      }


    }
    nbEqui = nbe[1];


    // for(int i=0;i<6;++i) std::cout <<nbe[i] << std::endl;
    // std::cout << periodicBE << std::endl;

  };


  bool isPeriodic(int lab)const {assert(lab>0);return isP[lab-1];}

};





template<typename T,typename B,typename V>
class GenericMesh
{
public:
  typedef GenericMesh GMesh;
  typedef T Element;
  typedef typename V::Rd Rd;
  typedef typename Rd::R R;
  typedef V  Vertex;
  typedef B BorderElement;
  typedef typename Element::RdHat RdHat;// for parametrization

  typedef const T* const_element_iterator;
  typedef       T*       element_iterator;
  typedef const V* const_vertex_iterator;
  typedef       V*       vertex_iterator;

  static const int nea=T::nea; //  numbering of adj (4 in Tet,  3 in Tria, 2 in seg)
  static const int nva=T::nva; //  numbering of vertex in Adj hyperface

  int nt,nv,nbe;
  R mes,mesb;
  Rd Pmin,Pmax; // // the bound  of the domain  see BuildBound
  int levelRefinement = 0;
  GMesh const * rootMesh = nullptr;
// protected
  V *vertices;
  T *elements;
  B *borderelements;

  int *TheAdjacencesLink; // to store the adj link  k*nea+i -> k'*nea+i'
  int *BoundaryElementHeadLink; //
  int *ElementConteningVertex;

  // For mesh build from interface
  Uint *ElementIndexInBackMesh;   // from loc to BackMesh
  Uint *VertexIndexInBackMesh;
  Uint *ElementIndexInLocMesh;    // from backMesh to loc

  // For mesh build from coarse mesh
  Uint *ElementIndexInCoarseMesh;

public:
  int nbElmts() const {return nt;}
  int nbBrdElmts() const {return nbe;}
  int nbVertices() const {return nv;}
  const T & operator[](int i) const {return elements[CheckT(i)];}
  const V& operator()(int i) const {return vertices[CheckV(i)];}
  const B& be(int i) const {return borderelements[CheckBE(i)];}
  void  BoundingBox(Rd &pmin,Rd &pmax) const {pmin=Pmin;pmax=Pmax;}
  T & t(int i)  {return elements[CheckT(i)];}
  V & v(int i)  {return vertices[CheckV(i)];}
  B & be(int i) {return borderelements[CheckBE(i)];}


  GenericMesh()
    : nt(0), nv(0), nbe(0), mes(0.), mesb(0.) ,
      levelRefinement(0),
      vertices(0),elements(0),borderelements(0), //bnormalv(0),
      TheAdjacencesLink(0),
      BoundaryElementHeadLink(0),
      ElementConteningVertex(0),
      ElementIndexInBackMesh(0),
      VertexIndexInBackMesh(0),
      ElementIndexInLocMesh(0),
      ElementIndexInCoarseMesh(0)
  {}

  virtual void info() {
    std::cout << "\n ----- Mesh " << this << " info ----- "<< std::endl;
    std::cout << " nb of nodes            : \t" << nv << std::endl;
    std::cout << " nb of elements         : \t" << nt << std::endl;
    std::cout << " nb of border elements  : \t" << nbe << std::endl;

  }

  void set(int mv,int mt,int mbe)
  {
    assert(nt==0 && nv==0 && nbe ==0);
    nt=mt;
    nv=mv;
    nbe=mbe;
    vertices=new V[nv];
    elements= new T[nt];
    borderelements = new B[nbe];

    assert( nt >=0 && elements);
    assert( nv >0 && vertices);
  }
  // #ifdef USE_MPI
  // int first_element() const { return MPIcf::my_rank();}
  // int next_element() const {return MPIcf::size();}
  // #else
  // int first_element() const { return 0;}
  // int next_element() const {return 1;}
  // #endif
  // int last_element() const { return nt;}
  // int last_boundary_element() const { return nbe;}

  int operator()(const T & tt) const {return CheckT(&tt - elements);}
  int operator()(const T * tt) const {return CheckT(tt - elements);}
  int operator()(const V & vv) const {return CheckV(&vv - vertices);}
  int operator()(const V  * vv) const{return CheckV(vv - vertices);}
  int operator()(const B & k) const {return CheckBE(&k - borderelements);}
  int operator()(const B  * k) const{return CheckBE(k - borderelements);}
  int operator()(int it,int j) const {return operator()(elements[it][j]);}// Nu vertex j of triangle it
  int at(int it,int j) const {return operator()(elements[it][j]);}// Nu vertex j of triangle it

  int be(int it,int j) const {return operator()(borderelements[it][j]);}// Nu vertex j of triangle it

  int CheckV (int i) const { ASSERTION(i>=0 && i < nv);  return i;}
  int CheckT (int i) const { ASSERTION(i>=0 && i < nt);  return i;}
  int CheckBE(int i) const { ASSERTION(i>=0 && i < nbe); return i;}

  //for cutMesh
  int idxElementInBackMesh(int k) const {ASSERTION(k>=0 && k < nt);
    return (ElementIndexInBackMesh)?ElementIndexInBackMesh[k] : k;}
  int idxElementFromBackMesh(int k) const {
    return (ElementIndexInLocMesh)?ElementIndexInLocMesh[k] : k;}


  void BuildAdj();
  void Buildbnormalv();
  void BuildBound();
  void BuildjElementConteningVertex();

  virtual DataFENodeDF BuildDFNumbering(int dfon[NbTypeItemElement],
					int nndon[NbTypeItemElement],
					int N=1,
					const PeriodicBC* PPeriod = nullptr)const;
					// int nbPeriodicBe = 0, int* periodicBe = 0 ) const ;

  virtual DataFENodeDF BuildDFNumbering(int dfon[NbTypeItemElement],
					int N=1,
					const PeriodicBC* PPeriod = nullptr)const
					// int nbPeriodicBe = 0, int* periodicBe = 0) const
  {
    return  BuildDFNumbering(dfon, dfon, N, PPeriod);//nbPeriodicBe, periodicBe);
  }

  DataFENodeDF BuildDFNumbering(int ndfv,int ndfe,int ndff,int ndft,
				int nndv,int nnde,int nndf,int nndt,
				int N=1,
				const PeriodicBC* PPeriod = nullptr//)const;
				// int nbPeriodicBe = 0, int* periodicBe = 0
				) const
  { int dfon[NbTypeItemElement]={ndfv,ndfe,ndff,ndft};
    int ndon[NbTypeItemElement]={nndv,nnde,nndf,nndt};
    return  BuildDFNumbering(dfon, ndon, N, PPeriod);//nbPeriodicBe, periodicBe);
  }

  int ElementAdj(int k,int &j) const  {
    int p = TheAdjacencesLink[nea*k+j];
    j=p%nea;
    return p>=0 ? p/nea: -1;
  }

  int GetAllElementAdj(int it,int *tabk) const{ //  get the tab of all adj element (max ne)
    //  and return the size of the tab
    int i=0;
    for(int j=0;j<nea;++j)
      {
	int p = TheAdjacencesLink[nea*it+j];
	if(p >=0 && p/nea!=it) {
	  tabk[i] = p/nea; i++;}
      }
    return i;
  }

  bool isOnBorder(int k) const{
    int i=0;
    for(int j=0;j<nea;++j) {
      int n = TheAdjacencesLink[3*k+j]/3;
      if(n >=0 && n!=k) i++;
    }
    return (i != nea);
  }

  int BoundaryElement(int bbe,int & ItemInK) const {
    int i= BoundaryElementHeadLink[bbe];
    ItemInK = i%nea;
    return i/nea;}

  int BoundaryElement(int bbe) const {
    int ItemInK;
    return BoundaryElement(bbe, ItemInK); }


  template<int N,int M>
  SortArray<int,N> iteme(const int (* const  nu )[N],int k,int i) {
    int nnv[N];
    Element & K(elements[CheckT(k)]);
    ASSERTION(i>=0 && i <M);
    for (int j=0;j<N;++j){
      nnv[j] = operator()(K[nu[i][j]]);
    }

    return SortArray<int,N>(nnv);
  }

  SortArray<int,B::nv> itemadj(int k,int i)
  {
    return iteme<B::nv,T::nea>(T::nvadj,k,i);
  }

  SortArray<int,B::nv> itembe(int k)
  {
    int nnv[B::nv];
    B & K(borderelements[CheckBE(k)]);

    for (int j=0;j<B::nv;++j){
      nnv[j] = operator()(K[j]);
    }

    return SortArray<int,B::nv>(nnv);
  }


  R mesure()const { return mes;}
  R bordermesure()const { return mesb;}
  virtual ~GenericMesh() {
    delete [] ElementConteningVertex;
    delete [] TheAdjacencesLink;
    delete [] BoundaryElementHeadLink;
    delete [] borderelements;
    if(nt>0) delete [] elements;
    delete [] vertices;
    //    delete [] bnormalv;
    if(ElementIndexInBackMesh) delete [] ElementIndexInBackMesh;
    if(VertexIndexInBackMesh)  delete [] VertexIndexInBackMesh;
    if(ElementIndexInLocMesh)  delete [] ElementIndexInLocMesh;
    if(ElementIndexInCoarseMesh) delete [] ElementIndexInCoarseMesh;
  }

private:
  GenericMesh(const GenericMesh &); // pas de construction par copie
  void operator=(const GenericMesh &);// pas affectation par copy
};



template<typename GMesh>
class BuildAdjacencyOfMesh {

  typedef SortArray<int,GMesh::nva> SArray;

  GMesh & mesh;
  HashTable<SArray,int> *h;
  int nk=0, nba=0;

public:
  BuildAdjacencyOfMesh(GMesh &);
  void initializeArray() ;
  void findAdjacencyElement();
  void addFace(const int, const int);
  void addFaceFirstTime(const SArray &);
  void addFaceAlreadySeen(typename HashTable<SArray,int>::iterator, const SArray&);
  void findBoundaryElement();
  void addBoundary(const int);

public :
  ~BuildAdjacencyOfMesh() {
    delete h;
  }
};



template<typename GMesh>
class BuildDofNumberingOfMesh {

  typedef typename GMesh::Element T;
  typedef typename T::Rd Rd;

public :
  const GMesh & Th;

  int nbNodeInK;
  int* nbDofOnItem;
  int* nbNodeOnItem;
  bool constndfPerNode = false;
  bool nodearevertices;
  int maxNbNodePerElement=0;
  int maxNbDFPerElement=0;
  int nbNodes=0;
  int nbOfDF=0;
  int *p = 0, *pp=0;
  int nbnzero=0;
  int N = 1;

  BuildDofNumberingOfMesh(const GMesh& m,
    int ndfon[NbTypeItemElement],
    int nndon[NbTypeItemElement],
    int NN=1)
    : Th(m), nbDofOnItem(ndfon) , nbNodeOnItem(nndon), N(NN){
      for(int i=0;i<4;++i) {
        assert(nndon[i] == 0 || nndon[i] == 1);
      }
    }

  void computeData(const PeriodicBC* PPeriod){// int nbequi) {

    const int nk[]={T::nv,T::ne,T::nf,T::nt};
    nbNodeInK = T::NbNodes(nbNodeOnItem);
    int mindf = 1000;
    int maxdf = 0;
    for (int dd=0;dd<NbTypeItemElement;++dd) {
      if(nbDofOnItem[dd]) {
        nbnzero++;
        mindf = min(mindf, nbDofOnItem[dd]);
        maxdf = max(maxdf, nbDofOnItem[dd]);
        maxNbDFPerElement   += nbDofOnItem[dd]*nk[dd];
        maxNbNodePerElement += nbNodeOnItem[dd]*nk[dd];
      }
    }
    if(mindf == maxdf) constndfPerNode = true;
    nodearevertices = ( nbnzero ==1  && nbDofOnItem[0]) && !PPeriod;//nbequi == 0 ;


  }

  DataFENodeDF buildDataFE();

  public:

  DataFENodeDF createDataFENodeDF(){
    return DataFENodeDF(nbDofOnItem,Th.nt,nbNodes,nbOfDF,p,pp,
			maxNbNodePerElement,maxNbDFPerElement,
			constndfPerNode);
  }


};





class BuildDofNumberingOfMeshHelper {
protected:
  typedef SortArray<unsigned int,2> Key;
  Key* keys;
  int* keysdim, *permArray;
  int ndim[NbTypeItemElement]={0,0,0,0};
  int nbequi = 0;
  int* equibe = 0;

  int of, nbElement, nbb = 0;
  int itemCounterInK, dofCounter = 0;
  HashTable<Key,int>* h;
  HashTable<Key,Key>* equi;
  unsigned int tinfty=-1;

  virtual int nbDofOnItem(const int i) = 0;

  template<typename GMesh>
  void addDofOnVertexOfK(const GMesh& Th, const typename GMesh::Element & K) {
    for(int i=0;i<GMesh::Element::nv;++i) {

      keysdim[itemCounterInK]=0;                     // what item
      keys[itemCounterInK++]=Key(Th(K[i]),tinfty);
    }
  }

  template<typename GMesh>
  void addDofOnVertexOfB(const GMesh& Th, int * v1, int* v2) {
    typedef typename GMesh::BorderElement B;

    for(int i=0;i<B::nv;++i) {

      keysdim[itemCounterInK]=0;                     // what item
      keys[itemCounterInK++]=Key(v1[i], tinfty);
      keys[itemCounterInK++]=Key(v2[i], tinfty);
    }
  }

  template<typename GMesh>
  void addDofOnEdgeOfK(const GMesh& Th, const typename GMesh::Element & K) {
    for(int i=0;i<GMesh::Element::ne;++i) {

      keysdim[itemCounterInK]=1;               // what item (1 = edge)

      keys[itemCounterInK++]=Key(Th(K[GMesh::Element::nvedge[i][0]]),
				 Th(K[GMesh::Element::nvedge[i][1]]));
    }

  }

  template<typename GMesh>
  void addDofOnFaceOfK(const GMesh& Th, const typename GMesh::Element & K) {

    const int k = Th(K);
    for(int ii,i=0;i<GMesh::Element::nf;++i) {

      keysdim[itemCounterInK]=2;               // what item (2 = face)

      if(GMesh::Element::nf == 1)
      keys[itemCounterInK++] = Key(k+of, tinfty);
      else {

        int kAdj = Th.ElementAdj(k,ii=i);
        int iAdj = (kAdj == -1) ? --nbb : kAdj;

        assert( (of + nbb) >= 0);
        keys[itemCounterInK++] = Key(k+of, iAdj+of);
      }
    }
  }

  template<typename GMesh>
  void addDofOnVolumeOfK(const GMesh& Th, const typename GMesh::Element & K) {

    const int k = Th(K);
    keysdim[itemCounterInK]=3;               // what item (3 = volume)
    keys[itemCounterInK++]=Key(k+of, tinfty);


  }

  virtual void fillTheArray() = 0;

public :
  ~BuildDofNumberingOfMeshHelper() {
    if(keys) delete[] keys;
    if(keysdim) delete[] keysdim;
  }

};



template<typename GMesh>
class BuildDofNumberingOfMeshPk : public BuildDofNumberingOfMeshHelper {
public :
  const int nkeys=4+6+4+1;                 // max number of items

  typedef typename GMesh::Element T;
  typedef typename GMesh::BorderElement B;

  BuildDofNumberingOfMesh<GMesh>& builder;

  BuildDofNumberingOfMeshPk(BuildDofNumberingOfMesh<GMesh>& buildDof,
			    const PeriodicBC* PPeriod)
			    // int nbPeriodicBe, int* periodicBe)
    : builder(buildDof)
      //, nbequi(nbPeriodicBe), equibe(periodicBe)
  {
    if(PPeriod) {
      nbequi = PPeriod->nbEqui;//nbPeriodicBe;
      equibe = PPeriod->periodicBE;//periodicBe;
    }

    // if(builder.nodearevertices) {
    //   std::cout << " hello " << std::endl;
    //   builder.nbNodes = builder.Th.nv;
    //   builder.nbOfDF  = builder.nbNodes*builder.nbDofOnItem[0];
    //   std::cout << builder.nbNodes << std::endl;
    //   std::cout << builder.nbOfDF << std::endl;
    // }
    // else {
      setData();
      setPeriodicBC();
      performNodeNumbering();
    // }

  }

  int nbNodeOnItem(const int i) { return builder.nbNodeOnItem[i]; }
  int nbDofOnItem(const int i) { return builder.nbDofOnItem[i]; }


  void setData() {

    nbElement = builder.Th.nt;
    of = builder.Th.nv + 10;
    builder.nbOfDF=0;
    dofCounter=0;

    builder.p =  new int[builder.nbNodeInK*nbElement];
    keys = new Key[nkeys*2];
    keysdim = new int[nkeys*2];
    h = new HashTable<Key,int>(builder.nbNodeInK*nbElement,
			       of+nbElement);

    int  nbmaxeq = 1+builder.nbNodeInK*nbequi;
    int  nbhasheq = nbequi ? of+nbElement : 1;
    equi = new HashTable<Key,Key>(nbmaxeq,nbhasheq);


  }

  void setPeriodicBC() {

    const int nbev= B::nv;

    for( int k=0;k<nbequi;++k) {
      int v1[nbev] = {};             // idx vertex of first boundary
      int v2[nbev] = {};             // idx vertex second boundary

      int be1 = equibe[2*k];
      int be2 = equibe[2*k+1];
      // std::cout << " id boundary  " << be1 << "\t" << be2 << std::endl;

      const B& b1(builder.Th.be(be1));
      const B& b2(builder.Th.be(be2));

      for(int i=0;i<nbev;++i) {
        v1[i] = builder.Th(b1[i]);
        v2[i] = builder.Th(b2[i]);

        // std::cout << v1[i] << "\t" << v2[i] << std::endl;
      }

      itemCounterInK=0;
      if(nbDofOnItem(0)) addDofOnVertexOfB(builder.Th, v1, v2);
      // if(nbDofOnItem(1)) addDofOnEdgeOfB(builder.Th, v1, v2);


      for(int j=0;j<itemCounterInK;) {
        int i0=j++,i1=j++;
        // if(keys[i1]<keys[i0]) swap( keys[i0],keys[i1]);
        if(keys[i0]< keys[i1] )
        equi->add(keys[i0],keys[i1]);
        else
        equi->add(keys[i1],keys[i0]);
      }
    }

    for (typename HashTable<Key,Key>::iterator qe,pe=equi->begin() ; pe != equi->end(); ++pe) {

      // std::cout << pe->k << " \t " << pe->v << std::endl;
      qe=equi->find(pe->v);
      if(qe){
        // std::cout << pe->k << " " << pe->v << " <=> " << qe->k << std::endl;
        pe->v = qe->v;
      }
    }


  }

  void performNodeNumbering() {
    for(int k=0;k<nbElement;++k) {              // loop over element
      const T & K(builder.Th[k]);             // the element
      itemCounterInK=0;

      if(nbDofOnItem(0)) addDofOnVertexOfK(builder.Th, K);
      if(nbDofOnItem(1)) addDofOnEdgeOfK(builder.Th, K);
      if(nbDofOnItem(2)) addDofOnFaceOfK(builder.Th, K);
      if(nbDofOnItem(3)) addDofOnVolumeOfK(builder.Th, K);

      fillTheArray();

    }
    if(!builder.constndfPerNode) {

      builder.pp =  new int[builder.nbNodes+1];

      int kk=0,nn=0;
      for(int k=0; k<nbElement; ++k)
      for(int i=0; i<builder.nbNodeInK; i++)
      builder.pp[builder.p[nn++]]=builder.nbDofOnItem[keysdim[i]];

      for(int n=0; n<builder.nbNodes; ++n) {
        int ndfn=builder.pp[n];
        builder.pp[n]=kk;
        kk += ndfn;
      }
      builder.pp[builder.nbNodes] = builder.nbOfDF;
      assert(kk==builder.nbOfDF);

    }
  }


void fillTheArray() {

  for(int i=0;i<itemCounterInK;i++) {
    Key ki=keys[i];
    // Check if belong to equi Htable pe
    // if yes change ki to pe->v

    typename HashTable<Key,Key>::iterator pe= equi->find(ki);
    if(pe) {   ki= pe->v;
      // std::cout << pe->v << std::endl;
    }


    typename HashTable<Key,int>::iterator pk= h->find(ki); // look for the key
    if(!pk)  {
      pk = h->add(ki,builder.nbNodes);                     // add the key
      builder.nbNodes += nbNodeOnItem(keysdim[i]);
      builder.nbOfDF  += builder.nbDofOnItem[keysdim[i]];  // add the nb of df to add
    }
    for(int j=0;j<nbNodeOnItem(keysdim[i]);++j) {
      builder.p[dofCounter++] = pk->v;
      ndim[keysdim[i]]++;
    }
  }

}

~BuildDofNumberingOfMeshPk() {
  if (h) delete h;
  if (equi) delete equi;
}

};












// enum ElementSignEnum { AllElement, NegElement, PosElement, NoElement};
//
//
// template<typename T>
// class SignElement{
//
//   typedef T Element;
//   static const int nve = Element::nv;
//
//   int sum_;
// public:
//   SignElement() : sum_(0) {}
//   SignElement(const double ls[nve]) : sum_(0) {
//     for (Ubyte i= 0; i < nve; ++i)
//       sum_ += fsign( ls[i]);
//   }
//   SignElement(const byte ls[nve]) : sum_(0) {
//     for (Ubyte i= 0; i < nve; ++i)
//       sum_ += ls[i];
//   }
//
//
//   bool cut() const  {return abs(sum_) != nve;}
//   byte sign() const { return (sum_ == nve) ? 1 : ((sum_ == -nve)? -1 : 0);}
//
// };


//Represents the reference element which is cut by a linear level set function ls.
template<typename T>
class GSignPatternTrait
{
  typedef T Element;
  static const int nve = Element::nv;

  protected:
  Ubyte num_root_vert_;  ///< number of vertices, where the level set function is zero.
  Ubyte num_root_;       ///< number of roots of the level set function; num_root_vert <= num_root
  byte sign_[Element::nv];  //Sign of the level set function of the vertices; \f$\in\{-1,0,1\}\f$
  Ubyte cut_simplex_[Element::nvc]; ///< local number with respect to the reference tetra of the object on the cut: [0,num_root_vert): vertex numbers in (0..2); [num_root_vert, num_root): edge numbers in (0..3). Both parts are sorted in increasing order.
  Ubyte cut_simplex_rep_[Element::nvc]; ///< local number of the object on the cut: (0..5)

  void compute_cuts ();

  public:
  GSignPatternTrait () : num_root_vert_( 3), num_root_( 2) {} ///< Uninitialized default state
  GSignPatternTrait (const byte   ls[nve]) { assign( ls); } ///< Assign the sign pattern on the vertices; throws if ls is identically 0.
  GSignPatternTrait (const double ls[nve]) { assign( ls); } ///< Assign a sign pattern on the vertices; throws if ls is identically 0.
  void assign (const byte   ls[nve]); ///< Assign a sign pattern on the vertices; throws if ls is identically 0.
  void assign (const double ls[nve]); ///< Assign a sign pattern on the vertices; throws if ls is identically 0.

  byte sign (int i) const { return sign_[i]; } ///< -1,0,1; sign of vertex i.

  bool empty () const { return num_root_ == 0; } ///< True, iff there is no intersection.

  bool no_zero_vertex () const { return num_root_vert_ == 0; } ///< True, iff there is no vertex, in which ls vanishes.

  Ubyte num_cut_simplexes () const { return num_root_; } ///< Number of edges and vertices with a root of ls.
  Ubyte num_zero_vertexes () const { return num_root_vert_; } ///< Number of vertices of the tetra that are roots of ls.

  /// Return local number of edges/verts with a root of ls. For edges, [] returns a edge number in 0..5, and () returns an extended vertex number in 4..9.
  ///@{
  Ubyte operator[] (int i) const { return cut_simplex_[i]; }
  Ubyte operator() (int i) const { return cut_simplex_rep_[i]; }
  ///@}


};

template<typename T>
void GSignPatternTrait<T>::assign (const byte ls[Element::nv])
{
  num_root_vert_= num_root_= 0;

  byte sum= 0;
  for (Ubyte i= 0; i < Element::nv; ++i) sum += (sign_[i] = ls[i]);
  if (sum == Element::nv || sum == -Element::nv) return;
  compute_cuts ();
}

template<typename T>
void GSignPatternTrait<T>::assign (const double ls[Element::nv])
{
  num_root_vert_= num_root_= 0;

  byte sum= 0;
  for (Ubyte i= 0; i < Element::nv; ++i)
    sum += (sign_[i] = sign( ls[i]));
  if (sum == Element::nv || sum == -Element::nv) // optimize the case of uncut tetras
    return;

  compute_cuts ();
}

template<typename T>
void GSignPatternTrait<T>::compute_cuts ()
{
  for (Ubyte i= 0; i < Element::nv; ++i){
    if (sign( i) == 0) {
      cut_simplex_[num_root_vert_++]= i;
    }
    num_root_= num_root_vert_;
  }
  for (Ubyte i= 0; i < Element::ne; ++i){
    if (sign( Element::nvedge[i][0])*sign( Element::nvedge[i][1]) == -1) {
      cut_simplex_[num_root_++]= i;
    }
    std::memcpy( cut_simplex_rep_, cut_simplex_, Element::nvc*sizeof(byte));
  }
  for (int i= num_root_vert_; i < num_root_; ++i) cut_simplex_rep_[i]+= Element::nv;

}








#endif
