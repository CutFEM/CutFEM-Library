#include "GenericMesh.hpp"


template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::BuildjElementConteningVertex()
{
  const int nkv= T::nv;
  if(!ElementConteningVertex) ElementConteningVertex = new int[nv];

  for(int i=0;i<nv;++i) ElementConteningVertex[i]=-1;

  for (int k=0;k<nt;++k)
    for (int i=0;i<nkv;++i)
      ElementConteningVertex[operator()(elements[k][i])]=k ;

    int kerr=0;
    for(int i=0;i<nv;++i)
      if (ElementConteningVertex[i]<0)  kerr++;
    assert(kerr==0);
}

template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::BuildAdj()
{
  if(TheAdjacencesLink!=0) return ;           //  already build

  BuildAdjacencyOfMesh<GenericMesh<T,B,V>> a(*this);

  ne_ = a.ne - this->nbe;

  // this->BuildInnerFace();

}

// template<typename T,typename B,typename V>
// void GenericMesh<T,B,V>::BuildInnerFace(){
//   inner_faces_ = new Face[ne_];
//   int idx = 0;
//   int iv[T::nva];
//   for(int k=0;k<nt;++k) {
//
//     for(int ie=0;ie<T::nea;++ie) {
//       int je = ie;
//       int kn = this->ElementAdj(k, je);
//
//       // skip boundary edges and only go once
//       if(kn < k) continue;
//       assert(idx < ne_);
//
//       for(int i=0;i<T::nva;++i) iv[i] = this->at(k, T::nvedge[ie][i]);
//
//       inner_faces_[idx].set(vertices, iv, 0);
//       inner_faces_[idx].set_adjacent_element(k, kn);
//
//       idx++;
//     }
//   }
// }
//



template<typename GMesh>
BuildAdjacencyOfMesh<GMesh>::BuildAdjacencyOfMesh(GMesh & m) : mesh(m) {
  initializeArray();

  findAdjacencyElement();

  findBoundaryElement();
}

template<typename GMesh>
void BuildAdjacencyOfMesh<GMesh>::initializeArray() {
  mesh.TheAdjacencesLink = new int[mesh.nea*mesh.nt];
  mesh.BoundaryElementHeadLink = new int[mesh.nbe];
  h = new HashTable<SArray,int>(mesh.nea*mesh.nt,mesh.nv);
}

template<typename GMesh>
void BuildAdjacencyOfMesh<GMesh>::findAdjacencyElement() {
  nk=0, nba=0;ne=0;
  for (int k=0;k<mesh.nt;++k) {
    for (int i=0;i<mesh.nea;++i) {
      addFace(k,i);
    }
  }
}

template<typename GMesh>
void BuildAdjacencyOfMesh<GMesh>::addFace(const int k, const int i) {
  SArray a(mesh.itemadj(k,i));    // sort nodes on the adj item
  typename HashTable<SArray,int>::iterator p = h->find(a);
  if(!p) {
    ne++;
    addFaceFirstTime(a);
  }
  else {
    addFaceAlreadySeen(p,a);
  }
  ++nk;
}

template<typename GMesh>
void BuildAdjacencyOfMesh<GMesh>::addFaceFirstTime(const SArray & a) {
  h->add(a,nk);
  mesh.TheAdjacencesLink[nk]=-1;
  nba++;
}

template<typename GMesh>
void BuildAdjacencyOfMesh<GMesh>::addFaceAlreadySeen(typename HashTable<SArray,int>::iterator p, const SArray& a) {
  assert(p->v>=0);
  mesh.TheAdjacencesLink[nk]=p->v;
  mesh.TheAdjacencesLink[p->v]=nk;

  p->v=-1-nk;
  nba--;
}

template<typename GMesh>
void BuildAdjacencyOfMesh<GMesh>::findBoundaryElement() {
  for (int k=0;k<mesh.nbe;++k) {
    addBoundary(k);
  }
}

template<typename GMesh>
void BuildAdjacencyOfMesh<GMesh>::addBoundary(const int k) {
  int err = 0;
  SArray a(mesh.itembe(k));
  typename HashTable<SArray,int>::iterator p= h->find(a);
  if(!p) {
    err++;
    if(err==1)   cerr << "Err  Border element not in mesh \n";
    if(err<10)   cerr << " \t " << k << " " << a << endl;
  }
  else {
    mesh.BoundaryElementHeadLink[k] = p->v <0 ? -p->v-1 : p->v;
  }
}


template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::BuildBound() {
  mes=0.;
  mesb=0.;

  for (int i=0;i<nt;i++)   mes += this->elements[i].mesure();

  for (int i=0;i<nbe;i++)  mesb += this->be(i).mesure();

  if(vertices && (nv>0)) {
    Pmin=vertices[0];
    Pmax=vertices[0];
    for(int i=1;i<nv;++i) {
      Pmin=Minc(Pmin,vertices[i]);
      Pmax=Maxc(Pmax,vertices[i]);
    }
  }
}





template<typename T,typename B,typename V>
DataFENodeDF GenericMesh<T,B,V>::BuildDFNumbering(
              int nbDofOnItem[NbTypeItemElement],
						  int nbNodeOnItem[NbTypeItemElement],
						  int N,
						  const PeriodicBC* PPeriod
						  //						  int nbPeriodicBe, int* periodicBe
						  ) const
{
  BuildDofNumberingOfMesh<GenericMesh> builderDof(*this, nbDofOnItem, nbNodeOnItem, N);
  builderDof.computeData(PPeriod);//nbPeriodicBe);



  if(!builderDof.nodearevertices){
    BuildDofNumberingOfMeshPk<GenericMesh> ttt(builderDof, PPeriod);//nbPeriodicBe, periodicBe);
  }
  else {

    builderDof.nbNodes = builderDof.Th.nv;
    builderDof.nbOfDF  = builderDof.nbNodes*builderDof.nbDofOnItem[0];
  }


  return builderDof.createDataFENodeDF();

}




template class  GenericMesh<Seg1,BoundaryPoint1,Vertex1>;
template class  GenericMesh<Triangle2,BoundaryEdge2,Vertex2>;
template class  GenericMesh<Tet,Triangle3,Vertex3>;
template class  GenericMesh<Quad2,BoundaryEdge2,Vertex2>;
template class  GenericMesh<Hexa,Quad3,Vertex3>;

template class  BuildAdjacencyOfMesh<GenericMesh<Seg1,BoundaryPoint1,Vertex1>>;
template class  BuildAdjacencyOfMesh<GenericMesh<Triangle2,BoundaryEdge2,Vertex2>>;
template class  BuildAdjacencyOfMesh<GenericMesh<Tet,Triangle3,Vertex3>>;
template class  BuildAdjacencyOfMesh<GenericMesh<Quad2,BoundaryEdge2,Vertex2>>;
template class  BuildAdjacencyOfMesh<GenericMesh<Hexa,Quad3,Vertex3>>;

template class  BuildDofNumberingOfMesh<GenericMesh<Seg1,BoundaryPoint1,Vertex1>>;
template class  BuildDofNumberingOfMesh<GenericMesh<Triangle2,BoundaryEdge2,Vertex2>>;
template class  BuildDofNumberingOfMesh<GenericMesh<Tet,Triangle3,Vertex3>>;
template class  BuildDofNumberingOfMesh<GenericMesh<Quad2,BoundaryEdge2,Vertex2>>;
template class  BuildDofNumberingOfMesh<GenericMesh<Hexa,Quad3,Vertex3>>;
