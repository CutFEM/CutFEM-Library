#include <fstream>
#include <iostream>
#include "../util/ufunction.hpp"
#include "Mesh2dn.hpp"
#include "RNM.hpp"
#include "libmesh5.h"
#include "Interface2dn.hpp"
#include "timeInterface.hpp"



Mesh2::Mesh2(const char * filename) { // read the mesh

  int nt,nv,nbe;
  int ok=1;//load(filename);
  if(ok)
    {
      ifstream f(filename);
      if(!f) {cerr << "Mesh2::Mesh2 Erreur openning " << filename<<endl;exit(1);}
      if(verbosity)
      cout << " Read On file \"" <<filename<<"\""<<  endl;
      f >> nv >> nt >> nbe ;
      this->set(nv,nt,nbe);
      if(verbosity)
      cout << "  -- Nb of Vertex " << nv << " " << " Nb of Triangles " << nt
	   << " , Nb of border edges " << nbe <<  endl;
      assert(f.good() && nt && nv) ;
      for (int i=0;i<nv;i++)
	{
	  f >> this->vertices[i];
	  assert(f.good());
	}
      mes=0;
      for (int i=0;i<nt;i++)
	{
	  this->t(i).Read1(f,this->vertices,nv);
	  mes += t(i).mesure();
	}
      mesb=0.;
      for (int i=0;i<nbe;i++)
	{
	  this->be(i).Read1(f,this->vertices,nv);
	  mesb += be(i).mesure();
	}
    }
  BuildBound();
  BuildAdj();
//   Buildbnormalv();
  BuildjElementConteningVertex();

  if(verbosity)
  cout << "   - mesh mesure = " << mes << " border mesure: " << mesb << endl;
}

Mesh2::Mesh2(int nx, int ny, R orx, R ory, R lx, R ly) {

  // int idQ[2][3] = {{0,1,2},{1,3,2}};
   int idQ2[2][3] = {{0,1,2},{3,2,1}};
   int idQ[2][3] = {{0,1,3},{2,0,3}};

  int mv  = nx * ny;
  int mt  = 2 * (nx - 1)*(ny - 1);
  int mbe = 2 * ((nx - 1) + (ny - 1));
  const R hx = lx / (nx - 1);
  const R hy = ly / (ny - 1);
  this->set(mv,mt,mbe);

  KN<int> iv(4),  indT(3);

  int jt = 0;
  for(int j=0; j<ny-1; j++) {
    for(int i=0; i<nx-1; i++) {

      int id=0;
      for(int jj=j; jj<j+2; ++jj) {
        for(int ii=i; ii<i+2; ++ii) {

          int ivl  =  ii + jj*nx;                          // index
          iv(id++) = ivl;

          vertices[ivl].x = ii*hx + orx;
          vertices[ivl].y = jj*hy + ory;
        }
      }

      for(int l=0;l<2;++l){                                 // create 2 elements
        for(int e=0; e<3; ++e) {
          if(jt==0 || jt == 1 || jt == mt-1 || jt == mt-2) indT(e) = iv(idQ2[l][e]);
          else indT(e) = iv(idQ2[l][e]);
        }

        elements[jt++].set(vertices, indT, 0);
      }
    }
  }



  // create the for borders
  int lab, k=0;
  for(int i=0; i<nx-1;++i) {
    indT(0) = i; indT(1) = i+1;lab = 1;
    for(int j=0;j<2;++j) vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
    borderelements[k++].set(vertices,indT, lab);
  }
  for(int i=0; i<ny-1;++i) {
    indT(0) = (i+1)*nx - 1 ; indT(1) = indT(0) + nx; lab = 2;
    for(int j=0;j<2;++j) vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
    borderelements[k++].set(vertices,indT, lab);
  }
  for(int i=0; i<nx-1;++i) {
    indT(0) = i+nx*(ny-1); indT(1) = indT(0)+1;lab = 3;
    for(int j=0;j<2;++j) vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
    borderelements[k++].set(vertices,indT, lab);
  }
  for(int i=0; i<ny-1;++i) {
    indT(0) = i*nx ; indT(1) = indT(0) + nx; lab = 4;
    for(int j=0;j<2;++j) vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
    borderelements[k++].set(vertices,indT, lab);
  }
  assert(k == nbe);


  BuildBound();
  BuildAdj();
//   Buildbnormalv();
  BuildjElementConteningVertex();
}


void Mesh2::square(int nx, R orx, R lx) {
  assert(!vertices);
  int ny = nx;
  R ory = orx;
  R ly = lx;

  int idQ[2][3] = {{0,1,2},{1,3,2}};
  int mv  = nx * ny;
  int mt  = 2 * (nx - 1)*(ny - 1);
  int mbe = 2 * ((nx - 1) + (ny - 1));
  const R hx = lx / (nx - 1);
  const R hy = ly / (ny - 1);
  this->set(mv,mt,mbe);

  KN<int> iv(4),  indT(3);

  int jt = 0;
  for(int j=0; j<ny-1; j++) {
    for(int i=0; i<nx-1; i++) {

      int id=0;
      for(int jj=j; jj<j+2; ++jj) {
  	for(int ii=i; ii<i+2; ++ii) {

  	  int ivl  =  ii + jj*nx;                          // index
  	  iv(id++) = ivl;

  	  vertices[ivl].x = ii*hx + orx;
  	  vertices[ivl].y = jj*hy + ory;
  	}
      }

      for(int l=0;l<2;++l){                                 // create 2 elements
      	for(int e=0; e<3; ++e) {
      	  indT(e) = iv(idQ[l][e]);
	}

      	elements[jt++].set(vertices, indT, 0);
      }
    }
  }



  // create the for borders
  int lab, k=0;
  for(int i=0; i<nx-1;++i) {
    indT(0) = i; indT(1) = i+1;lab = 1;
    for(int j=0;j<2;++j) vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
    borderelements[k++].set(vertices,indT, lab);
  }
  for(int i=0; i<ny-1;++i) {
    indT(0) = (i+1)*nx - 1 ; indT(1) = indT(0) + nx; lab = 2;
    for(int j=0;j<2;++j) vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
    borderelements[k++].set(vertices,indT, lab);
  }
  for(int i=0; i<nx-1;++i) {
    indT(0) = i+nx*(ny-1); indT(1) = indT(0)+1;lab = 3;
    for(int j=0;j<2;++j) vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
    borderelements[k++].set(vertices,indT, lab);
  }
  for(int i=0; i<ny-1;++i) {
    indT(0) = i*nx ; indT(1) = indT(0) + nx; lab = 4;
    for(int j=0;j<2;++j) vertices[indT(j)].lab = std::max(vertices[indT(j)].lab, lab);
    borderelements[k++].set(vertices,indT, lab);
  }
  assert(k == nbe);


  BuildBound();
  BuildAdj();
//   Buildbnormalv();
  BuildjElementConteningVertex();
}


// Mesh2::Mesh2(const GenericInterface<Mesh2>& gamma) {
Mesh2::Mesh2(const Interface2& gamma) {

  const Mesh2& backMesh(*gamma.backMesh);                   // the back-ground mesh

  const int mt = gamma.nbElement();            // exact approximation (2D)
  const int mv = 3*mt;                         // large approximation
  const int mbe = 0;                           // no boundary element

  this->set(mv,mt,mbe);
  ElementIndexInBackMesh = new Uint[mt];
  ElementIndexInLocMesh = new Uint[backMesh.nt];
  VertexIndexInBackMesh = new Uint[mv];

  nv=0;                                        // reinitialization

  KN<int> foundVertex(backMesh.nv); foundVertex = -1;
  int indT[3], jt=0;

  for(Uint ifac=0; ifac != gamma.nbElement(); ++ifac) {

    const Uint k = gamma.idxElementOfFace(ifac);
    const Element & K = (backMesh)[k];

    for( int i=0; i<3;++i) {
      const int idxG = (backMesh)(K[i]);

      if(foundVertex(idxG) == -1) {

	foundVertex(idxG) = nv;
	vertices[nv].x = K[i].x;
	vertices[nv].y = K[i].y;
	VertexIndexInBackMesh[nv] = idxG;
	indT[i] = nv++;
      }
      else {
	indT[i] = foundVertex(idxG);
      }
    }
    ElementIndexInBackMesh[jt] = k;
    ElementIndexInLocMesh[k] = jt;
    elements[jt++].set(vertices, indT, 0);
  }
  assert(jt == nt);

  BuildBound();
  BuildAdj();
  BuildjElementConteningVertex();

}


// To build mesh for time problem
Mesh2::Mesh2(TimeInterface2& gamma) {
  const Uint nInterface = gamma.size();
  const Mesh2& backMesh(*(gamma[0]->backMesh));                   // the back-ground mesh


  for(Uint i=0;i<nInterface-1;++i) {             // all interface must have the same backMesh
    assert(gamma[i]->backMesh == gamma[i+1]->backMesh);
  }

  int mmt = 0;
  for(Uint i=0;i<nInterface;++i) {
    mmt += gamma[i]->nbElement();            // exact approximation (2D)
  }
  const int mt = 2*mmt;
  const int mv = 3*mt;                         // large approximation
  const int mbe = 0;

  this->set(mv,mt,mbe);
  ElementIndexInBackMesh = new Uint[mt];
  ElementIndexInLocMesh = new Uint[backMesh.nt];
  VertexIndexInBackMesh = new Uint[mv];

  nv=0;                                        // reinitialization
  KN<int> foundVertex(backMesh.nv); foundVertex = -1;
  KN<int> foundElement(backMesh.nt); foundElement = -1;
  int indT[3], jt=0;

  for(Uint iInterface=0; iInterface<nInterface; ++iInterface) {
    for(Uint ifac=0; ifac != gamma[iInterface]->nbElement(); ++ifac) {

      const Uint k = gamma[iInterface]->idxElementOfFace(ifac);
      const Element & K = (backMesh)[k];


      if(foundElement(k) == -1) {

	for( int i=0; i<3;++i) {
	  const int idxG = (backMesh)(K[i]);

	  if(foundVertex(idxG) == -1) {

	    foundVertex(idxG) = nv;
	    vertices[nv].x = K[i].x;
	    vertices[nv].y = K[i].y;
	    VertexIndexInBackMesh[nv] = idxG;
	    indT[i] = nv++;
	  }
	  else {
	    indT[i] = foundVertex(idxG);
	  }
	}
	foundElement(k) = jt;
	ElementIndexInBackMesh[jt] = k;
	ElementIndexInLocMesh[k] = jt;
	elements[jt++].set(vertices, indT, 0);
      }
    }
  }
  nt = jt;

  // we need to get the elements that change domain during the time slab
  for(int k=0; k<backMesh.nt;++k) {
    if(foundElement(k) ==-1) {

      const Element & K = (backMesh)[k];

      // We only need to look at the sign of one of the node
      bool changeSign = false;    // excluded by default
      for(Uint i=0; i<nInterface-1;++i) {
	int idx1 = (backMesh)(K[0]);
	//	if( (gamma[i].levelSet)[idx1] * (gamma[i+1].levelSet)[idx1] <= 0 )
	if( (gamma[i]->ls_sign)[idx1] * (gamma[i+1]->ls_sign)[idx1] <= 0 )

	  changeSign = true;
      }
      if(changeSign) {
      	for( int i=0; i<3;++i) {
      	  const int idxG = (backMesh)(K[i]);

      	  if(foundVertex(idxG) == -1) {

      	    foundVertex(idxG) = nv;
      	    vertices[nv].x = K[i].x;
      	    vertices[nv].y = K[i].y;
      	    VertexIndexInBackMesh[nv] = idxG;
      	    indT[i] = nv++;
      	  }
      	  else {
      	    indT[i] = foundVertex(idxG);
      	  }
      	}
	if(jt>=mt) assert(0);
	foundElement(k) = jt;
      	ElementIndexInBackMesh[jt] = k;
      	ElementIndexInLocMesh[k] = jt;
      	elements[jt++].set(vertices, indT, 0);
      }
      nt = jt;
    }

  }


  BuildBound();
  BuildAdj();
  BuildjElementConteningVertex();

}




Mesh2::Mesh2(const Mesh2& Th, std::string whatToDo) {
  assert(whatToDo == "refine");
  this->rootMesh = &Th;
  this->levelRefinement = Th.levelRefinement+1;
  int of = Th.nv+10;
  unsigned int tinfty=-1;
  const int nkeys=6+4;                 // max number of items
  const int nkv= Element::nv;
  const int nke= Element::ne;
  const int nnodeK = nkv + nke;
  int nbNodes=0;

  int mt = Element::nref * Th.nt;
  int mv = Element::nref * Th.nv;
  int mbe = Th.nbe * BorderElement::nref;  this->set(mv,mt,mbe);
  ElementIndexInCoarseMesh = new Uint[mt];

  typedef SortArray<unsigned int,2> Key;
  Key keys[nkeys];
  int keysdim[nkeys];
  HashTable<Key,int> h(nnodeK*Th.nt,of+Th.nt);
  int idx[nkeys];

  int kt = 0;
  for(int k=0;k<Th.nt;++k) {
    int m = 0;
    const Element & K(Th[k]);                            // the element
    for(int i=0;i<nkv;++i) {      // loop over the nodes
      keysdim[m]=0;               // what item
      keys[m++]=Key(Th(K[i]),tinfty);  // index saved
    }
    for(int i=0;i<nke;++i) {      // loop over edges
      keysdim[m]=1;               // what item
      keys[m++]=Key(Th(K[Element::nvedge[i][0]]),
      Th(K[Element::nvedge[i][1]])); // sorted id
    }

    for(int i=0;i<m;i++) {
      Key ki=keys[i];//,kio=ki;
      typename HashTable<Key,int>::iterator pk= h.find(ki); // look for the key
      if(!pk) {                          // if not found

        if(keysdim[i] == 0) { // on node
          for(int d=0;d < Rd::d;++d) vertices[nbNodes][d] = K[i][d];
          // M::vertices[nbNodes].lab = K[i].lab;
          // std::cout << vertices[nbNodes]  << std::endl;
        }
        else {
          int j = i - nkv;
          const int i0 = Element::nvedge[j][0], i1 = Element::nvedge[j][1];
          for(int d=0;d < Rd::d;++d) vertices[nbNodes][d] = 0.5*(K[i0][d] + K[i1][d]);
          // M::vertices[nbNodes].lab = K[i].lab;
          // std::cout << vertices[nbNodes]  << std::endl;
        }

        pk = h.add(ki,nbNodes++);        // add the key
      }
      idx[i] = pk->v;
    }
    int idxK[nkv];

    for(int i=0;i<Element::nref;++i) {
      for(int j=0;j<nkv;++j) idxK[j] = idx[Element::refElement[i][j]];
      ElementIndexInCoarseMesh[kt] = k;
      elements[kt++].set(vertices, idxK, 0);
      // std::cout << elements[kt-1] << std::endl;
    }
  }



  nv = nbNodes;
  assert(kt == nt);  //  take care of the boundary elements
  int n = 0;
  for(int ke=0; ke<Th.nbe;++ke) {                    // loop over boundary element
    int kf, k = Th.BoundaryElement(ke,kf);           // k is the element index
    const Element & K(Th[k]);
    const BorderElement & Be(Th.be(ke));


    if(Rd::d == 2) {
      int idxBE[3];

      Key k0(Th(K[Element::nvedge[kf][0]]),tinfty);
      typename HashTable<Key,int>::iterator pk = h.find(k0);
      idxBE[0] = pk->v;

      Key k1(Th(K[Element::nvedge[kf][1]]),tinfty);
      pk =  h.find(k1);
      idxBE[2] = pk->v;

      Key ki(Th(Be[0]), Th(Be[1]));
      pk =  h.find(ki);
      if(!pk) assert(0);

      idxBE[1] = pk->v;

      borderelements[n++].set(vertices,idxBE  , Be.lab);
      borderelements[n++].set(vertices,idxBE+1, Be.lab);

    }
    else {
      int idxK[6];
      for(int i=0;i<BorderElement::nv;++i) {
        Key k0(Th(Be[i]),tinfty);
        typename HashTable<Key,int>::iterator pk = h.find(k0);
        idxK[i] = pk->v;
      }
      for(int i=0;i<BorderElement::ne;++i) {

        Key k0(Th(Be[BorderElement::nvedge[i][0]]),
        Th(Be[BorderElement::nvedge[i][1]]));

        typename HashTable<Key,int>::iterator pk = h.find(k0);
        idxK[i+3] = pk->v;
      }
      int idxBE[3];
      for(int j=0;j<BorderElement::nref;++j) {
        for(int i=0;i<3;++i) idxBE[i] = idxK[BorderElement::refElement[j][i]];

        borderelements[n++].set(vertices,idxBE, Be.lab);

      }
    }
  }


  BuildBound();
  BuildAdj();
  //   Buildbnormalv();
  BuildjElementConteningVertex();

  // std::cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << std::endl;
  // std::cout << "  -- Mesh3 : "<< " d "<< 3  << ", n Tet " << nt << ", n Vtx "
  // 	    << nv << " n Bord " << nbe << std::endl;
  ffassert(mes>=0);

}
