#ifndef GENERICMESH_HPP_
#define GENERICMESH_HPP_

extern long verbosity;

#include "cassert"
#include "../util/util.hpp"
#include <cstdlib>

// #include "RefCounter.hpp"

using namespace ::std;

#include "R3.hpp"
#include "HashTable.hpp"
#include "GenericVertex.hpp"
#include "dataStruct1D.hpp"
#include "dataStruct2D.hpp"
#include "dataStruct3D.hpp"
#include "cutFEMConfig.h"
#ifdef USE_MPI
#include "../parallel/cfmpi.hpp"
#endif

inline int maxdfon(const int *dfon){
  return max(max(dfon[0],dfon[1]),max(dfon[2],dfon[3]));}




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




template<typename T,typename B,typename V>
class GenericMesh {
public:
  typedef GenericMesh GMesh;
  typedef T Element;
  typedef typename V::Rd Rd;
  typedef typename Rd::R R;
  typedef V  Vertex;
  typedef B BorderElement;
  typedef typename Element::RdHat RdHat;// for parametrization
  typedef typename Element::Face Face;


  static const int nea=T::nea; //  numbering of adj (4 in Tet,  3 in Tria, 2 in seg)
  static const int nva=T::nva; //  numbering of vertex in Adj hyperface

  int nt,nv,nbe;
  R mes,mesb;


// protected
  V *vertices;
  T *elements;
  B *borderelements;

  int *TheAdjacencesLink; // to store the adj link  k*nea+i -> k'*nea+i'
  int *BoundaryElementHeadLink; //

public:
  int nbElmts() const {return nt;}
  int get_nb_element() const {return nt;}
  int nbBrdElmts() const {return nbe;}
  int nbVertices() const {return nv;}
  int nbElements() const {return nt;}
  int NbElement() const {return nt;}

  int nbBorderElements() const {return nbe;}


  const T & operator[](int i) const {return elements[CheckT(i)];}
  const V& operator()(int i) const {return vertices[CheckV(i)];}
  const B& be(int i) const {return borderelements[CheckBE(i)];}

  T & t(int i)  {return elements[CheckT(i)];}
  V & v(int i)  {return vertices[CheckV(i)];}
  B & be(int i) {return borderelements[CheckBE(i)];}


  GenericMesh()
    : nt(0), nv(0), nbe(0), mes(0.), mesb(0.) ,
      vertices(0),elements(0),borderelements(0),
      TheAdjacencesLink(0),
      BoundaryElementHeadLink(0)  {}

  virtual void info() {
    std::cout << " ----- Mesh " << this << " info ----- "<< std::endl;
    std::cout << " nb of nodes            : \t" << nv << std::endl;
    std::cout << " nb of elements         : \t" << nt << std::endl;
    std::cout << " nb of border elements  : \t" << nbe << std::endl;

  }

  void set(int mv,int mt,int mbe) {
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

  #ifdef USE_MPI
  virtual int first_element() const { return MPIcf::first_element(this->nbElements());}
  virtual int next_element() const {  return MPIcf::next_element(this->nbElements());}
  virtual int last_element() const {  return MPIcf::last_element(this->nbElements());}

  virtual int first_boundary_element() const { return MPIcf::my_rank();}
  virtual int next_boundary_element() const { return MPIcf::size();}
  virtual int last_boundary_element() const {return this->nbBrdElmts();}
  #else
  virtual int first_element() const { return 0;}
  virtual int next_element() const {return 1;}
  virtual int last_element() const { return this->nbElements();}

  virtual int first_boundary_element() const { return 0;}
  virtual int next_boundary_element() const { return 1;}
  virtual int last_boundary_element() const {return this->nbBrdElmts();}
  #endif


  int operator()(const T & tt) const {return CheckT(&tt - elements);}
  int operator()(const T * tt) const {return CheckT(tt - elements);}
  int operator()(const V & vv) const {return CheckV(&vv - vertices);}
  int operator()(const V  * vv) const{return CheckV(vv - vertices);}
  int operator()(const B & k) const {return CheckBE(&k - borderelements);}
  int operator()(const B  * k) const{return CheckBE(k - borderelements);}

  int operator()(int it,int j) const {return operator()(elements[it][j]);}// Nu vertex j of triangle it
  int at(int it,int j) const {return operator()(elements[it][j]);}// Nu vertex j of triangle it
  int be(int it,int j) const {return operator()(borderelements[it][j]);}// Nu vertex j of triangle it

  int CheckV (int i) const { assert(i>=0 && i < nv);  return i;}
  int CheckT (int i) const { assert(i>=0 && i < nt);  return i;}
  int CheckBE(int i) const { assert(i>=0 && i < nbe); return i;}


  void BuildAdj();
  void Buildbnormalv();
  void BuildBound();
  void BuildjElementConteningVertex();

  virtual DataFENodeDF BuildDFNumbering(int dfon[NbTypeItemElement], int nndon[NbTypeItemElement], int N=1)const;

  virtual DataFENodeDF BuildDFNumbering(int dfon[NbTypeItemElement], int N=1)const{
    return  BuildDFNumbering(dfon, dfon, N);
  }

  DataFENodeDF BuildDFNumbering(int ndfv,int ndfe,int ndff,int ndft, int nndv,int nnde,int nndf,int nndt, int N=1 ) const {
    int dfon[NbTypeItemElement]={ndfv,ndfe,ndff,ndft};
    int ndon[NbTypeItemElement]={nndv,nnde,nndf,nndt};
    return  BuildDFNumbering(dfon, ndon, N);
  }

  int ElementAdj(int k,int &j) const  {
    int p = TheAdjacencesLink[nea*k+j];
    j=p%nea;
    return p>=0 ? p/nea: -1;
  }

  int GetAllElementAdj(int it,int *tabk) const{ //  get the tab of all adj element (max ne)
    int i=0;
    for(int j=0;j<nea;++j) {
      int p = TheAdjacencesLink[nea*it+j];
      if(p >=0 && p/nea!=it) {
        tabk[i] = p/nea; i++;
      }
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
    return i/nea;
  }

  int BoundaryElement(int bbe) const {
    int ItemInK;
    return BoundaryElement(bbe, ItemInK);
   }


  template<int N,int M>
  SortArray<int,N> iteme(const int (* const  nu )[N],int k,int i) {
    int nnv[N];
    Element & K(elements[CheckT(k)]);
    assert(i>=0 && i <M);
    for (int j=0;j<N;++j){
      nnv[j] = operator()(K[nu[i][j]]);
    }

    return SortArray<int,N>(nnv);
  }
  SortArray<int,B::nv> itemadj(int k,int i)  {
    return iteme<B::nv,T::nea>(T::nvadj,k,i);
  }
  SortArray<int,B::nv> itembe(int k) {
    int nnv[B::nv];
    B & K(borderelements[CheckBE(k)]);
    for (int j=0;j<B::nv;++j){
      nnv[j] = operator()(K[j]);
    }
    return SortArray<int,B::nv>(nnv);
  }

  R mesure()const { return mes;}
  R bordermesure()const { return mesb;}
  double get_mesh_size() const {
    double hh = 1e300;
    for(int k=0;k<nt;++k){
      hh = min((*this)[k].hElement(), hh);
    }
    return hh;
  }

  virtual ~GenericMesh() {
    delete [] TheAdjacencesLink;
    delete [] BoundaryElementHeadLink;
    delete [] borderelements;
    if(nt>0) delete [] elements;
    delete [] vertices;
  }

private:
  GenericMesh(const GenericMesh &); // pas de construction par copie
  void operator=(const GenericMesh &);// pas affectation par copy
};



template<typename GMesh>
class BuildAdjacencyOfMesh {
public:
  typedef SortArray<int,GMesh::nva> SArray;

  GMesh & mesh;
  HashTable<SArray,int> *h;
  int nk=0, nba=0;
  int ne = 0;


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

  void computeData(){

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
    nodearevertices = ( nbnzero ==1  && nbDofOnItem[0]);


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
        // std::cout << of << "\t" << nbb << std::endl;
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

  BuildDofNumberingOfMeshPk(BuildDofNumberingOfMesh<GMesh>& buildDof)
			    // int nbPeriodicBe, int* periodicBe)
    : builder(buildDof)
      //, nbequi(nbPeriodicBe), equibe(periodicBe)
  {
      setData();
      performNodeNumbering();

  }

  int nbNodeOnItem(const int i) { return builder.nbNodeOnItem[i]; }
  int nbDofOnItem(const int i) { return builder.nbDofOnItem[i]; }


  void setData() {

    nbElement = builder.Th.nt;
    of = 2*builder.Th.nv;
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

    typename HashTable<Key,int>::iterator pk= h->find(ki); // look for the key
    if(!pk)  {
      pk = h->add(ki,builder.nbNodes);                     // add the key
      builder.nbNodes += nbNodeOnItem(keysdim[i]);
      builder.nbOfDF  += builder.nbDofOnItem[keysdim[i]];  // add the nb of df to add
    }
    std::cout << keysdim[i] << "\t" << nbNodeOnItem(keysdim[i]) << std::endl;
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








#endif
