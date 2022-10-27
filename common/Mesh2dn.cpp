#include <fstream>
#include <iostream>
#include "../util/ufunction.hpp"
#include "Mesh2dn.hpp"
#include "RNM.hpp"
#include "libmesh5.h"



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

}
MeshQuad2::MeshQuad2(int nx, int ny, R orx, R ory, R lx, R ly) {

  // int idQ[2][3] = {{0,1,2},{1,3,2}};
   // int idQ2[2][3] = {{0,1,2},{3,2,1}};
   int indQ[4] = {0,1,3,2};

  int mv  = nx * ny;
  int mt  = (nx - 1)*(ny - 1);
  int mbe = 2 * ((nx - 1) + (ny - 1));
  const R hx = lx / (nx - 1);
  const R hy = ly / (ny - 1);
  this->set(mv,mt,mbe);

  KN<int> iv(4),  indT(4);

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
      for(int e=0; e<4; ++e) {
        indT(e) = iv(indQ[e]);
      }
      elements[jt++].set(vertices, indT, 0);
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
}
