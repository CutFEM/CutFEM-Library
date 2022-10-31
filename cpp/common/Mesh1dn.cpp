#include <fstream>
#include <iostream>
#include "../util/ufunction.hpp"
#include "../util/error.hpp"
#include "RNM.hpp"
#include "libmesh5.h"
#include "Mesh1dn.hpp"

long verbosity=2;

Mesh1::Mesh1(const char * filename) { // read the mesh

  int nt,nv,nbe;
  int ok=0;//load(filename);
  if(ok)
    {
      ifstream f(filename);
      if(!f) {cerr << "Mesh1::Mesh1 Erreur openning " << filename<<endl;exit(1);}
      if(verbosity)
	cout << " Read On file \"" <<filename<<"\""<<  endl;
      f >> nv >> nt >> nbe ;
      this->set(nv,nt,nbe);
      if(verbosity)
	  cout << "  -- Nb of Vertex " << nv << " " << " Nb of Seg " << nt
	       << " , Nb of border Vertex " << nbe <<  endl;
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






Mesh1::Mesh1(int nx, R orx, R lx) {

  int mv  = nx;
  int mt  = (nx - 1);
  int mbe = 2;
  const R hx = lx / (nx - 1);

  this->set(mv,mt,mbe);

  for(int i=0; i<nx; i++) {
    vertices[i].x = i*hx + orx;
  }

  for(int k=0; k<nx-1; k++) {
    int iv[2] = {k,k+1};
    elements[k].set(vertices, iv, 0);
  }

  int iv[1] = {0};
  be(0).set(vertices,iv, 0);
  iv[0] = nx-1;
  be(1).set(vertices,iv, 1);




  BuildBound();
  BuildAdj();
}
