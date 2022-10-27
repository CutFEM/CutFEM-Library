#include "GenericMesh.hpp"



template<typename T,typename B,typename V>
void GenericMesh<T,B,V>::BuildAdj() {
  if(TheAdjacencesLink!=0) return ;           //  already build

  BuildAdjacencyOfMesh<GenericMesh<T,B,V>> a(*this);

}


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
}



template<typename T,typename B,typename V>
DataFENodeDF GenericMesh<T,B,V>::BuildDFNumbering(
              int nbDofOnItem[NbTypeItemElement],
						  int nbNodeOnItem[NbTypeItemElement],
						  int N ) const
{
  BuildDofNumberingOfMesh<GenericMesh> builderDof(*this, nbDofOnItem, nbNodeOnItem, N);
  builderDof.computeData();



  if(!builderDof.nodearevertices){
    BuildDofNumberingOfMeshPk<GenericMesh> ttt(builderDof);
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
