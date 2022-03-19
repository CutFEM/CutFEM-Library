#ifndef _CUTFESPACE_HPP
#define _CUTFESPACE_HPP

#include "expression.hpp"
#include "../common/marker.hpp"
#include <unordered_map>

#include "cutFEMConfig.h"
#ifdef USE_MPI
#include "../parallel/cfmpi.hpp"
#endif

template<class M>
class CutGFESpace;
template<class M>
class SubDomainArray;

template<class M>
class GSubDomain {
  typedef M Mesh;
  typedef typename Mesh::SignPattern SignPattern;
  typedef GFESpace<M> FESpace;
  typedef GenericInterface<M> Interface;
  typedef typename FESpace::FElement FElement;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::Rd Rd;
  typedef FunFEM<Mesh> Fun;

public :

  const FESpace& Vh;
  int nt;
  int nv, ndof;
  int NodesPerElement;
  KN<int> NodeOnWhat;
  const int sign;

protected:
  KN<int>   Tri_LocToGlob;

  KN<int>   Tri_GlobToLoc;
  KN<int>   Ver_GlobToLoc;

  KN<int>   NodesOfElement;
  KN<int>   FirstDfOfNodeData;


public:

  GSubDomain(const FESpace& vh) : Vh(vh), nt(0), nv(0){}
  GSubDomain(const FESpace& vh, Fun& ls, int sign_domain) : Vh(vh), sign(sign_domain){
    this->init(Vh);
    this->make(Vh, ls);
  }
  GSubDomain(const FESpace& vh, vector<Fun>& ls, int sign_domain) : Vh(vh), sign(sign_domain){
    this->init(Vh);
    for(int i=0;i<ls.size();++i) this->add(Vh, ls[i]);
    this->finalize(Vh);
  }
  GSubDomain(const FESpace& vh, const Interface& interface, int sign_domain) : Vh(vh), sign(sign_domain){
    this->init(Vh);
    this->make(Vh, interface.ls_sign);
  }
  GSubDomain(const FESpace& vh, const TimeInterface<M>& interface, int sign_domain) : Vh(vh), sign(sign_domain){
    this->init(Vh);
    for(int i=0;i<interface.size();++i) this->add(Vh, interface[i]->ls_sign);
    this->finalize(Vh);
  }

  GSubDomain(const FESpace& vh, int sign_domain) : Vh(vh), sign(sign_domain)  {
    this->init(Vh);
  }
private:
  void init(const FESpace & Vh) {
    NodesPerElement = Vh[0].NbNode();
    Tri_GlobToLoc.resize(Vh.nbElement); Tri_GlobToLoc = -1;
    Tri_LocToGlob.resize(Vh.nbElement);
    Ver_GlobToLoc.resize(Vh.nbNode); Ver_GlobToLoc = -1;
    NodesOfElement.resize(Vh.nbElement*NodesPerElement);
    nt = 0;
    nv = 0; ndof = 0;
    NodeOnWhat.init(NodesPerElement);
    // idQuadTime = 0;

    int c= 0;
    const int nk[]={Element::nv,Element::ne,Element::nf,Element::nt};
    for(int i=0;i<4;++i)
    for(int k=0;k<nk[i];++k)
    for(int j=0;j<Vh.TFE(0)->nbNodeOnWhat[i];++j)
    NodeOnWhat(c++) = i;
  }
  void make(const FESpace& Vh, Fun& levelSet) {
    const int nbNode = Vh[0].NbNode();
    int c=0;
    double loc_ls[levelSet.size(0)];

    for(int k = 0; k < Vh.nbElement; ++k) {
      // levelSet.Fvec(k, loc_ls);
      levelSet.eval(loc_ls, k);
      const SignElement<Element> cut( loc_ls);

      if(cut.sign() == this->sign || cut.cut()) {
        this->addFE(k);

        for( int j = 0; j < nbNode; ++j) {
          NodesOfElement(c++) = this->addVertex( k,j );

        }
      }
    }
    assert(c == NodesPerElement*nt);

    if(!Vh.constDfPerNode) {
      int kk=0, nn=0;
      KN<int>&p(NodesOfElement);
      KN<int>&pp(FirstDfOfNodeData);
      pp.init(nv+1);

      const FElement& FK(Vh[Tri_LocToGlob(0)]);
      const int * ndfon = FK.tfe->ndfOn();
      for(int k=0; k<nt; ++k) {

        for(int i=0; i<nbNode; i++) {
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
    Ver_GlobToLoc.destroy();
    // idQuadTime++;
  }
  void add(const FESpace& Vh, Fun& levelSet) {
    const int nbNode = Vh[0].NbNode();
    int c= NodesPerElement*nt;
    double loc_ls[levelSet.size(0)];

    for(int k = 0; k < Vh.nbElement; ++k) {
      if(FEinSub(k)) continue;
      levelSet.eval(loc_ls, k);
      const SignElement<Element> cut( loc_ls);

      if(cut.sign() == this->sign || cut.cut()) {

        this->addFE(k);

        for( int j = 0; j < nbNode; ++j) {
          NodesOfElement(c++) = this->addVertex( k,j );

        }
      }
    }
    // idQuadTime++;
    assert(c == NodesPerElement*nt);
  }
  void finalize(const FESpace& Vh) {
    if(!Vh.constDfPerNode) {
      const int nbNode = Vh[0].NbNode();
      int kk=0, nn=0;
      KN<int>&p(NodesOfElement);
      KN<int>&pp(FirstDfOfNodeData);
      pp.init(nv+1);

      const FElement& FK(Vh[Tri_LocToGlob(0)]);
      const int * ndfon = FK.tfe->ndfOn();
      for(int k=0; k<nt; ++k) {

        for(int i=0; i<nbNode; i++) {
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
    Ver_GlobToLoc.destroy();
  }

  void make(const FESpace& Vh, const KN<byte>& levelSet) {
    assert(levelSet.size() == Vh.Th.nv);
    const int nbNode = Vh[0].NbNode();
    int c=0;
    double loc_ls[Element::nv];
    for(int k = 0; k < Vh.nbElement; ++k) {

      for(int i=0;i<Element::nv;++i) loc_ls[i] = levelSet(Vh.Th(k,i));
      const SignElement<Element> cut( loc_ls);

      if(cut.sign() == this->sign || cut.cut()) {
        this->addFE(k);

        for( int j = 0; j < nbNode; ++j) {
          NodesOfElement(c++) = this->addVertex( k,j );

        }
      }
    }
    assert(c == NodesPerElement*nt);

    if(!Vh.constDfPerNode) {
      int kk=0, nn=0;
      KN<int>&p(NodesOfElement);
      KN<int>&pp(FirstDfOfNodeData);
      pp.init(nv+1);

      const FElement& FK(Vh[Tri_LocToGlob(0)]);
      const int * ndfon = FK.tfe->ndfOn();
      for(int k=0; k<nt; ++k) {

        for(int i=0; i<nbNode; i++) {
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
    Ver_GlobToLoc.destroy();
   }
  void add (const FESpace& Vh, const KN<byte>& levelSet) {
    const int nbNode = Vh[0].NbNode();
    int c= NodesPerElement*nt;
    double loc_ls[Element::nv];

    for(int k = 0; k < Vh.nbElement; ++k) {
      if(FEinSub(k)) continue;
      for(int i=0;i<Element::nv;++i) loc_ls[i] = levelSet(Vh.Th(k,i));
      const SignElement<Element> cut( loc_ls);

      if(cut.sign() == this->sign || cut.cut()) {

        this->addFE(k);

        for( int j = 0; j < nbNode; ++j) {
          NodesOfElement(c++) = this->addVertex( k,j );

        }
      }
    }
    // idQuadTime++;
    assert(c == NodesPerElement*nt);
  }

  public :
  bool FEinSub(int i) const { return !(Tri_GlobToLoc[i] == -1);}
  // bool FEinSub(int i) const { return Tri_GlobToLoc.find(i) != Tri_GlobToLoc.end();}

  virtual int NbNode() const {return nv;}
  virtual int NbDoF() const {return ndof;}
  virtual int NbElement() const {return nt;}

  FElement operator[](int k) const { return Vh[Tri_LocToGlob(k)];}

private :
  // Access to the finite element
  //------------------------------------------------------------------
  int checkFE(int i) const {assert(i>=0 && i<=nt); return i;}
  int addFE(int i) { if(!FEinSub(i)) {Tri_LocToGlob(nt) = i; Tri_GlobToLoc[i] = nt++;
      return nt-1;}  return Tri_GlobToLoc[i];}


  // Access to the vertices
  //------------------------------------------------------------------
  bool VinSub(int i) const { return (Ver_GlobToLoc[i] != -1);}
  int checkV(int i) const {assert( i>=0 && i<= nv); return i;}

  int addVertex(int k, int j) {
    int i = Vh[k][j];
    if(!VinSub(i)) {
      //Ver_LocToGlob(nv) = i;
      Ver_GlobToLoc[i] = nv++;
      ndof += Vh.TFE(0)->ndfOn()[NodeOnWhat(j)];
      return nv-1;
    }
    return Ver_GlobToLoc[i] ;
  }



public :
  // int idxGlob2Loc(int i) const { return Ver_GlobToLoc[i];}
  int getTriLocToGlob(int i) const { return Tri_LocToGlob(checkFE(i));}
  int getTriGlobToLoc(int i) const { return checkFE(Tri_GlobToLoc[i]);}
  // int getTriGlobToLoc(int i) const { return Tri_GlobToLoc.find(i)->second;}


  int FirstDfOfNode(int i) const {
    return (FirstDfOfNodeData ? FirstDfOfNodeData[i] : i*Vh.Nproduit) ;
  }

  friend class CutGFESpace<M>;
  friend class SubDomainArray<M>;

private:
  GSubDomain(const GSubDomain &); // pas de construction par copie
  void operator=(const GSubDomain &);//
};



class DataCutFENodeDF   {

public :
  struct Indexing {
    int nbNode, nbDoF, nbElement ;
    Indexing(int num[3]) : nbNode(num[0]), nbDoF(num[1]), nbElement(num[2]){}
    Indexing() : nbNode(0), nbDoF(0), nbElement(0){}
  } * sub;

  KN<int> idx0_K;

  DataCutFENodeDF(){}
  DataCutFENodeDF( KNM<int>& numi) {
    int nSub = numi.N();
    idx0_K.init(nSub+1);
    sub = new Indexing[nSub];
    for(int s=0;s<nSub;++s){
      sub[s].nbNode = numi(s,0);
      sub[s].nbDoF = numi(s,1);
      sub[s].nbElement = numi(s,2);
      idx0_K[s+1] = numi(s,2);
    }

  }

  ~DataCutFENodeDF() {
    delete []sub;
  }
};


template<class M>
class SubDomainArray {
  typedef M Mesh;
  typedef GFESpace<M> FESpace;
  typedef GSubDomain<M> SubDomain;
  typedef GenericInterface<M> Interface;
  typedef typename FESpace::Rd Rd;
  typedef typename Mesh::Element Element;
  typedef typename FESpace::FElement FElement;

public:
  KN<SubDomain*> subDomain;
  int nSub;
  bool allocate_memory_subDomain = false;

  SubDomainArray(const GFESpace<M>& vh, const GenericInterface<M>& g, list<int> list_sign)
  : nSub(list_sign.size()) ,subDomain(list_sign.size()) {
    int nn = 0;
    for(auto it = list_sign.begin(); it!= list_sign.end();++it){
      subDomain(nn) = new GSubDomain<M>(vh, g, *it);
      nn++;
    }
    allocate_memory_subDomain = true;
  }
  SubDomainArray(const GFESpace<M>& vh, const TimeInterface<M>& g, list<int> list_sign)
  : nSub(list_sign.size()) ,subDomain(list_sign.size()) {
    int nn = 0;
    for(auto it = list_sign.begin(); it!= list_sign.end();++it){
      subDomain(nn) = new GSubDomain<M>(vh, g, *it);
      nn++;
    }
    allocate_memory_subDomain = true;
  }
  SubDomainArray(const KN<SubDomain*>& sd) : subDomain(sd), nSub(sd.size()){}


  DataFENodeDF BuildDFNumbering(const KN<SubDomain*>& sdomain);
  DataFENodeDF BuildDFNumbering(list<SubDomain*>& sdomain){
    return BuildDFNumbering(KN<SubDomain*>(sdomain));
  }
  DataFENodeDF BuildDFNumbering() {return BuildDFNumbering(subDomain);};

  DataCutFENodeDF BuildDFNumberingCut(const KN<SubDomain*>& sdomain);
  DataCutFENodeDF BuildDFNumberingCut(list<SubDomain*>& sdomain){
    return BuildDFNumberingCut(KN<SubDomain*>(sdomain));
  };
  DataCutFENodeDF BuildDFNumberingCut() {return BuildDFNumberingCut(subDomain);}

  ~SubDomainArray(){
    if(allocate_memory_subDomain){
      for(int i=0;i<nSub;++i) delete subDomain(i);
    }
  }
};


/*
 *   Class for the FESpace
 */
template<class M>
class CutGFESpace : public SubDomainArray<M>,
                    public GFESpace<M> ,
                    public DataCutFENodeDF {
public :
  typedef M Mesh;
  typedef GFESpace<M> FESpace;
  typedef GSubDomain<M> SubDomain;
  typedef GenericInterface<M> Interface;
  typedef typename FESpace::Rd Rd;
  typedef typename Mesh::Element Element;
  typedef typename Element::Face Face;
  typedef typename FESpace::FElement FElement;

  const int idx0 = 0;

  const FESpace& Vh;
  KN<const Marker*> marker;

  int ntCut = 0;
public :


CutGFESpace(const FESpace& vh, const Interface& g, list<int> list_sign) :
  SubDomainArray<M>(vh,g,list_sign),
  FESpace(vh, SubDomainArray<M>::BuildDFNumbering(), vh.periodicBC),
  DataCutFENodeDF(SubDomainArray<M>::BuildDFNumberingCut()),
  Vh(vh)
{
  this->gamma.resize(1); this->gamma[0] = &g;
}

CutGFESpace(const FESpace& vh, const TimeInterface<M>& g, list<int> list_sign) :
  SubDomainArray<M>(vh,g,list_sign),
  FESpace(vh, SubDomainArray<M>::BuildDFNumbering(), vh.periodicBC),
  DataCutFENodeDF(SubDomainArray<M>::BuildDFNumberingCut()),
  Vh(vh)
{
  this->gamma.resize(g.size());
  for(int i=0;i<g.size();++i) this->gamma[i]= g[i];
}


CutGFESpace(std::list<SubDomain*> sdomain, Interface& g) :
SubDomainArray<M>(sdomain),
FESpace((*sdomain.begin())->Vh, SubDomainArray<M>::BuildDFNumbering(sdomain), (*sdomain.begin())->Vh.periodicBC),
DataCutFENodeDF(SubDomainArray<M>::BuildDFNumberingCut()),
Vh((*sdomain.begin())->Vh)
{
  this->gamma.resize(1); this->gamma[0] = &g;
}

  // CutGFESpace(const KN<SubDomain*> sdomain, const Interface& g) :
  //   DataCutFENodeDF(this->BuildDFNumberingCut(sdomain)),
  //   FESpace(sdomain(0)->Vh, this->BuildDFNumbering(sdomain), sdomain(0)->Vh.periodicBC),
  //   Vh(sdomain(0)->Vh),
  //   SubDomainArray<M>(sdomain)
  // {
  //   this->gamma.resize(1); this->gamma[0] = &g;
  // }

  // CutGFESpace(const KN<SubDomain*> sdomain, const TimeInterface<M>& g) :
  //   DataCutFENodeDF(this->BuildDFNumberingCut(sdomain)),
  //   FESpace(sdomain(0)->Vh, this->BuildDFNumbering(sdomain), sdomain(0)->Vh.periodicBC),
  //   Vh(sdomain(0)->Vh),
  //   SubDomainArray<M>(sdomain)
  //   // subDomain(sdomain)
  // {
  //   // nSub = sdomain.size();
  //   this->gamma.resize(g.size());
  //   for(int i=0;i<g.size();++i) this->gamma[i]= g[i];
  // }

  CutGFESpace(std::list<SubDomain*> sdomain, const TimeInterface<M>& g) :
  SubDomainArray<M>(sdomain),
  FESpace((*sdomain.begin())->Vh, SubDomainArray<M>::BuildDFNumbering(sdomain), (*sdomain.begin())->Vh.periodicBC),
  DataCutFENodeDF(SubDomainArray<M>::BuildDFNumberingCut(sdomain)),
  Vh((*sdomain.begin())->Vh)
  {
    this->gamma.resize(g.size());
    for(int i=0;i<g.size();++i) this->gamma[i]= g[i];
  }



  CutGFESpace(const KN<SubDomain*> sdomain, const Marker& mark) :
  SubDomainArray<M>(sdomain),
  FESpace(sdomain(0)->Vh, SubDomainArray<M>::BuildDFNumbering(sdomain), sdomain(0)->Vh.periodicBC),
  DataCutFENodeDF(SubDomainArray<M>::BuildDFNumberingCut(sdomain)),
  Vh(sdomain(0)->Vh)
  {
    // nSub = sdomain.size();
    marker.resize(1); marker[0] = &mark;
  }
  // CutGFESpace(SubDomain& sdomain, const Marker& mark) :
  //   DataCutFENodeDF(this->BuildDFNumberingCut(sdomain)),
  //   FESpace(sdomain.Vh, this->BuildDFNumbering(sdomain), sdomain.Vh.periodicBC),
  //   subDomain(1, &sdomain),
  //   Vh(sdomain.Vh)
  //   {
  //     nSub = 1;//sdomain.size();
  //     marker.resize(1); marker[0] = &mark;
  //   }




  FElement operator[](int k) const {
    int dom_id = (k<sub[0].nbElement)? 0 : 1;
    int kth = this->subDomain(dom_id)->getTriLocToGlob(k - idx0_K[dom_id]);
    return FElement(this,k,kth);
  }

  bool isCutSpace() const {return true;}

  int NbNode() const { return this->nbNode;}
  int NbDoF() const { return this->nbDoF;}
  int NbDoF(int s) const { return this->subDomain(s)->NbDoF();}
  int NbElement() const { return this->nbElement;}
  int NbElement(int s) const { return this->subDomain(s)->NbElement();}

  #ifdef USE_MPI
  virtual int first_element() const { return MPIcf::first_element(this->nbElement);}
  virtual int next_element() const {  return MPIcf::next_element(this->nbElement);}
  virtual int last_element() const {  return MPIcf::last_element(this->nbElement);}
  #else
  virtual int last_element() const { return this->nbElement;}
  #endif

  virtual int idxElementFromBackMesh(int k, int i) const {
    if(i==-1) i = 0;
    assert(i >= 0 && i < 2);
    assert(this->nSub <= 2);
    if(!this->subDomain(i)->FEinSub(k)) return -1;
    int kk = this->subDomain(i)->getTriGlobToLoc(k);
    return idx0_K[i] + kk;
  }
  virtual bool containBackElement(int k)const {
    for(int i=0;i<this->nSub;++i) {
      if(this->subDomain(i)->FEinSub(k)) return true;
    }
    return false;
  }
  // virtual bool faceInDomain(const Face& face, int dom) const {
  //   int k1 = face.get_index_adjacent_element(0);
  //   int k2 = face.get_index_adjacent_element(1);
  //   if(this->subDomain(dom)->FEinSub(k1) && this->subDomain(dom)->FEinSub(k2)) {
  //     return true;
  //   }
  //   return false;
  // }


  virtual int idxElementFromBackMesh(int k) const {
    std::cout << "need to use a domain id to get element" << std::endl;
    assert(0);
    return k;
  }
  virtual int idxElementInBackMesh(int k) const {
    ASSERTION(k>=0 && k < nt);
    return Vh.Th((*this)[k].T);
  }
  virtual int nbDomain() const {
    return this->nSub;
  }
  // virtual int idxGlob2Loc(int k, int i) const {
  //   assert(this->nSub <= 2);
  //   return  this->subDomain(i)->idxGlob2Loc(k);
  // }
  virtual const FESpace& getBackSpace() const {return Vh;}
  virtual int whichDomain(int k) const { assert(this->nSub <= 2);return (k >= sub[0].nbElement);}

  int getNumberOfSubDomain() const { return 2;}
  bool isCut(int k) const {
    // assert(this->nSub == 2);
    int kk = (k<sub[0].nbElement)? k: k-sub[0].nbElement;
    int kback = this->subDomain((k>=sub[0].nbElement))->getTriLocToGlob(kk);
    if(this->nSub == 2){
      return  (this->subDomain(0)->FEinSub(kback) &&
      this->subDomain(1)->FEinSub(kback));
    }
    else{
      for(int i=0;i<this->gamma.size();++i) {
        auto it = this->gamma[i]->face_of_element_.find(kback);
        if(it != this->gamma[i]->face_of_element_.end()) return true;
      }
      return false;
    }

  }

  virtual int getNeighborElement(int k,int &j, int domain = 0) const {
    int k_back = this->idxElementInBackMesh(k);
    int kn_back = this->Th.ElementAdj(k_back,j);
    if(kn_back == -1) return -1;
    else return  this->idxElementFromBackMesh(kn_back, domain);
  }



  virtual void info() const {
    std::cout << "CutFESpace \t" << this << std::endl;
    std::cout << "nb node    \t" << NbNode() << std::endl;
    std::cout << "nb dof     \t" << NbDoF() << std::endl;
    std::cout << "nb element \t" << NbElement() << std::endl;
  }


  ~CutGFESpace(){ }
};

template<class M>
DataCutFENodeDF SubDomainArray<M>::BuildDFNumberingCut(const KN<SubDomain*>& sdomain) {

  this->nSub = sdomain.size();
  KNM<int> numi(this->nSub, 3);
  for(int s=0;s<this->nSub;++s) {
    numi(s,0) = sdomain(s)->NbNode();
    numi(s,1) = sdomain(s)->NbDoF();
    numi(s,2) = sdomain(s)->NbElement();
  }
  // int num1[3] = { sdomain(0)->NbNode(), sdomain(0)->NbDoF(), sdomain(0)->NbElement()};
  // int num2[3] = { sdomain(1)->NbNode(), sdomain(1)->NbDoF(), sdomain(1)->NbElement()};
  return DataCutFENodeDF(numi);
  // return DataCutFENodeDF(num1,num2);
}


template<class M>
DataFENodeDF SubDomainArray<M>::BuildDFNumbering(const KN<SubDomain*>& sdomain) {

  const FESpace& Mh(sdomain(0)->Vh);

  int ndfon[NbTypeItemElement] = {Mh.ndfonVertex(),
				  Mh.ndfonEdge(),
				  Mh.ndfonFace(),
				  Mh.ndfonTet()};
  nSub = sdomain.size();
  KN<int> nts(nSub);
  KN<int> nvs(nSub);
  KN<int> nv_start(nSub);
  KN<int> ndfs(nSub);

  for(int i=0;i<nSub;++i) nts(i) = sdomain(i)->NbElement();
  for(int i=0;i<nSub;++i) nvs(i) = sdomain(i)->NbNode();
  for(int i=0;i<nSub;++i) ndfs(i) = sdomain(i)->NbDoF();

  int nn = 0;
  for(int i=0;i<nSub;++i) {
    nv_start(i) = nn;
    nn += sdomain(i)->NbNode();
  }

  const int nt = nts.sum();
  const int nv = nvs.sum();
  const int ndof = ndfs.sum();

  int *p = 0, *pp=0;

  int maxNode = Mh.MaxNbNodePerElement;
  int maxDF   = Mh.MaxNbDFPerElement;

  p =  new int[maxNode*nt];
  int i = 0;
  for(int s=0;s<nSub;++s)
    for(int j=0; j<maxNode*nts[s];++j,++i)
      p[i] = sdomain(s)->NodesOfElement(j) + nv_start(s);


  // for(    i=0; i<maxNode*nt1;++i )    p[i] = sdomain(0)->NodesOfElement(i);
  // for(int j=0; j<maxNode*nt2;++j,++i) p[i] = sdomain(1)->NodesOfElement(j) + nv1;

  pp =  new int[nv+1];
  // int n0 = sdomain(0)->FirstDfOfNodeData(nv1);
  // for(i=0; i<nv1;++i )         pp[i] = sdomain(0)->FirstDfOfNodeData(i);
  // for(int j=0; j<nv2;++j,++i)  pp[i] = sdomain(1)->FirstDfOfNodeData(j) + n0;
  KN<int> dfBegin(nSub);
  nn = 0;
  for(int s=0;s<nSub;++s) {
    dfBegin(s) = nn;
    nn += sdomain(s)->FirstDfOfNode(nvs[s]);
  }

  i = 0;
  for(int s=0;s<nSub;++s)
    for(int j=0; j<nvs[s];++j,++i)  pp[i] = sdomain(s)->FirstDfOfNode(j) + dfBegin(s);



  pp[nv] = ndof;
  assert(i == nv);

  return DataFENodeDF(ndfon,
		      nt, nv, ndof,
		      p,pp,
		      maxNode,maxDF,
		      Mh.MaxNbDFPerNode,
		      Mh.constDfPerNode);

}


typedef GSubDomain<Mesh2> SubDomain2;
typedef GSubDomain<Mesh3> SubDomain3;
typedef CutGFESpace<Mesh2> CutFESpace2;
typedef CutGFESpace<Mesh3> CutFESpace3;

#endif
