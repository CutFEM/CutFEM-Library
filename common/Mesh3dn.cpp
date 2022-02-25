#include <fstream>
#include <iostream>
#include <cstring>
#include "libmesh5.h"
#include "../util/ufunction.hpp"
#include "../util/error.hpp"
#include "RNM.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "Interface3dn.hpp"
#include "timeInterface.hpp"




Mesh3::Mesh3(const string  filename)
{
  int ok=load(filename);
  if(ok) {
    std::cout << " could not load the file " << std::endl;
    //   ifstream f(filename.c_str());
    //   if(!f) {
    //     cerr << "  --  Mesh3::Mesh3 Erreur openning " << filename<<endl;
    //     ffassert(0);exit(1);}
    //   if(verbosity>2)
    //     cout << "  -- Mesh3:  Read On file \"" <<filename<<"\""<<  endl;
    //   if(filename.rfind(".msh")==filename.length()-4)
    //     readmsh(f);
    //   else {
    //     std::cout << " not a good format" << std::endl;
    //   }
  }

  BuildBound();
  if(nt > 0){
    BuildAdj();
    BuildjElementConteningVertex();
  }

  verbosity = 3;
  if(verbosity>2)
  cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;
  if(verbosity)
  cout << "  -- Mesh3 : "<<filename  << ", d "<< 3
  << ", n Tet " << nt << ", n Vtx "
  << nv << " n Bord " << nbe << endl;
  ffassert(mes>=0); // add F. Hecht sep 2009.
}



int Mesh3::load(const string & filename) {

  int ver,inm,dim;
  int lf=filename.size()+20;

  KN<char>  fileb(lf),filef(lf);
  char * pfile;
  strcpy(filef,filename.c_str());
  strcpy(fileb,filef);
  strcat(filef,".mesh");
  strcat(fileb,".meshb");
  //    int bin;
  if( (inm=GmfOpenMesh(pfile=fileb, GmfRead,&ver,&dim)) ) {
    //bin=true;
    std::cout << " opening " << (char *) fileb << std::endl;
  }
  else if( (inm=GmfOpenMesh(pfile=filef, GmfRead,&ver,&dim)) ) {
    //  bin=false;
    std::cout << " opening " << (char *) fileb << std::endl;
  }
  else {
    cerr << " Erreur ouverture file " << (char *) fileb
    << " " << (char *) filef << endl;
    return   1;
  }
  int nv,nt,neb;
  nv = GmfStatKwd(inm,GmfVertices);
  nt = GmfStatKwd(inm,GmfTetrahedra);
  neb= GmfStatKwd(inm,GmfTriangles);
  this->set(nv,nt,neb);
  if(verbosity>1)
  cout << "  -- Mesh3(load): "<<pfile <<", ver " << ver << ", d "<< dim
  << ", nt " << nt << ", nv " << nv << " nbe:  = " << nbe << endl;
  if(dim  != 3) {
    cerr << "Err dim == " << dim << " !=3 " <<endl;
    return 2; }
    if( nv<=0 && (nt <0 || nbe <=0)  ) {
      cerr << " missing data "<< endl;
      return 3;
    }

    int iv[4],lab;
    float cr[3];
    // read vertices
    GmfGotoKwd(inm,GmfVertices);
    int mxlab=0;
    int mnlab=0;
    for(int i=0;i<nv;++i) {
      if(ver<2) {
        GmfGetLin(inm,GmfVertices,&cr[0],&cr[1],&cr[2],&lab);
        vertices[i].x=cr[0];
        vertices[i].y=cr[1];
        vertices[i].z=cr[2];}
        else
        GmfGetLin(inm,GmfVertices,&vertices[i].x,&vertices[i].y,&vertices[i].z,&lab);
        vertices[i].lab=lab;

        mxlab= max(mxlab,lab);
        mnlab= min(mnlab,lab);
      }

      //    /* read mesh triangles */
      if(nbe > 0) {
        if(mnlab==0 && mxlab==0 ) {
          int kmv=0;
          mesb=0;
          GmfGotoKwd(inm,GmfTriangles);
          for(int i=0;i<nbe;++i) {
            GmfGetLin(inm,GmfTriangles,&iv[0],&iv[1],&iv[2],&lab);
            for(int j=0;j<3;++j)
            if(!vertices[iv[j]-1].lab) {
              vertices[iv[j]-1].lab=1;
              kmv++;
            }
            for (int j=0;j<3;++j) iv[j]--;
            this->be(i).set(this->vertices,iv,lab);
            mesb += this->be(i).mesure();
          }

          if(kmv&& verbosity>1)
          cout << " Aucun label Hack (FH)  ??? => 1 sur les triangle frontiere "<<endl;
        }
        else {
          mesb=0;
          GmfGotoKwd(inm,GmfTriangles);
          for(int i=0;i<nbe;++i) {
            GmfGetLin(inm,GmfTriangles,&iv[0],&iv[1],&iv[2],&lab);
            for (int j=0;j<3;++j) iv[j]--;
            this->be(i).set(this->vertices,iv,lab);
            mesb += this->be(i).mesure();
          }
        }
      }

      if(nt>0) {
        /* read mesh tetrahedra */
        GmfGotoKwd(inm,GmfTetrahedra);
        for(int i=0;i<nt;++i) {
          GmfGetLin(inm,GmfTetrahedra,&iv[0],&iv[1],&iv[2],&iv[3],&lab);
          assert( iv[0]>0 && iv[0]<=nv && iv[1]>0 && iv[1]<=nv && iv[2]>0 && iv[2]<=nv && iv[3]>0 && iv[3]<=nv);
          for (int j=0;j<4;j++) iv[j]--;
          this->elements[i].set(vertices,iv,lab);
          mes += this->elements[i].mesure();
        }
      }
      GmfCloseMesh(inm);
      return(0); // OK

    }


    void Mesh3::readmsh(ifstream & f) {
      f >> nv >> nt >> nbe;
      if(verbosity>2)
      cout << " GRead : nv " << nv << " " << nt << " " << nbe << endl;
      this->vertices = new Vertex[nv];
      this->elements = new Element [nt];
      this->borderelements = new BorderElement[nbe];
      for (int k=0; k<nv; k++) {
        Vertex & P = this->vertices[k];
        f >> P.x >>P.y >> P.z >> P.lab ;
      }
      mes=0.;
      mesb=0.;

      if(nt != 0) {
        for (int k=0; k<nt; k++) {
          int i[4],lab;
          Element & K(this->elements[k]);
          f >> i[0] >> i[1] >> i[2] >> i[3] >> lab;
          K.set(this->vertices,i,lab);
          mes += K.mesure();

        }
      }
      for (int k=0; k<nbe; k++) {
        int i[4],lab;
        BorderElement & K(this->borderelements[k]);
        f >> i[0] >> i[1] >> i[2]  >> lab;
        K.set(this->vertices,i,lab);
        mesb += K.mesure();
      }
    }

    Mesh3::Mesh3(int nnv, int nnt, int nnbe, Vertex3 *vv, Tet *tt, Triangle3 *bb) {

      nv = nnv;
      nt = nnt;
      nbe =nnbe;

      vertices = vv;
      elements = tt;
      borderelements = bb;

      mes=0.;
      mesb=0.;

      for (int i=0;i<nt;i++)
      mes += this->elements[i].mesure();

      for (int i=0;i<nbe;i++)
      mesb += this->be(i).mesure();

      //  Add FH to be consitant we all constructor ...  July 09
      BuildBound();
      if(nt > 0){
        BuildAdj();
        // 	  Buildbnormalv();
        BuildjElementConteningVertex();
      }
      //  end add
      // if(verbosity>1)
      // 	  cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;

      //    if(verbosity>1)
      //      cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;
      assert(mes>=0.);
    }

    Mesh3::Mesh3(int nnv, int nnbe, Vertex3 *vv, Triangle3 *bb) {

      nv = nnv;
      nbe =nnbe;

      vertices = vv;
      borderelements = bb;

      mes=0.;
      mesb=0.;

      for (int i=0;i<nbe;i++)
      mesb += this->be(i).mesure();

      //  Add FH to be consitant we all constructor ...  July 09
      BuildBound();
      if(nt > 0){
        BuildAdj();
        //     Buildbnormalv();
        BuildjElementConteningVertex();
      }
      //  end add

      if(verbosity>1)
      cout << "  -- End of Construct  mesh3: mesure = " << mes << " border mesure " << mesb <<  endl;
      ffassert(mes>=0); // add F. Hecht sep 2009.
    }


    int  signe_permutation(int i0,int i1,int i2,int i3){
      int p=1;
      if(i0>i1) Exchange(i0,i1), p = -p;
      if(i0>i2) Exchange(i0,i2), p = -p;
      if(i0>i3) Exchange(i0,i3), p = -p;
      if(i1>i2) Exchange(i1,i2), p = -p;
      if(i1>i3) Exchange(i1,i3), p = -p;
      if(i2>i3) Exchange(i2,i3), p = -p;
      return p;
    }

    Mesh3::Mesh3(int nx, int ny, int nz, R orx, R ory, R orz, R lx, R ly, R lz) {

      const int idQ[6][4] = {{1,3,2,5},{1,2,4,5},{0,1,2,4},
      {6,3,5,2},{6,2,5,4},{7,3,5,6}};
      int mv  = nx * ny * nz;
      int mt  = 6 * (nx - 1)*(ny - 1)*(nz - 1);
      int mbe = 4 * ((nx - 1) * (ny - 1) + (nx - 1) * (nz - 1) +(nz - 1) * (ny - 1) );
      const R hx = lx / (nx - 1);
      const R hy = ly / (ny - 1);
      const R hz = lz / (nz - 1);

      this->set(mv,mt,mbe);

      KN<int> iv(8),  indT(4);

      int jt = 0;
      for(int k=0; k<nz-1; k++) {
        for(int j=0; j<ny-1; j++) {
          for(int i=0; i<nx-1; i++) {

            int id=0;
            for(int kk=k; kk<k+2;++kk) {
              for(int jj=j; jj<j+2; ++jj) {
                for(int ii=i; ii<i+2; ++ii) {

                  int ivl  =  ii + jj*nx + kk*nx*ny;                    // index
                  iv(id++) = ivl;

                  vertices[ivl].x = ii*hx + orx;
                  vertices[ivl].y = jj*hy + ory;
                  vertices[ivl].z = kk*hz + orz;
                }
              }
            }

            for(int l=0;l<6;++l){                                        // create 2 elements
              for(int e=0; e<4; ++e) {
                indT(e) = iv(idQ[l][e]);
              }
              elements[jt++].set(vertices, indT, 0);
            }
          }
        }
      }

      assert(jt == nt);


      const R xmin = orx, ymin = ory, zmin = orz;
      const R xmax = orx + lx, ymax = ory + ly, zmax = orz + lz;

      // create the for borders
      int kf=0;
      int indV[3],indKV[3];
      const R step = 1e-12;
      for( int jt=0; jt<nt; ++jt) {                     // loop over element
        const Element & K((*this)[jt]);
        for (int k=0;k<4;++k) {                         // loop over faces
          for( int i=0;i<3;++i){
            indV[i] = Element::nvface[k][i];
            indKV[i] = (*this)(K[indV[i]]);
          }

          const int iv1 = Element::nvface[k][0];
          const int iv2 = Element::nvface[k][1];
          const int iv3 = Element::nvface[k][2];

          //  Faces on the plane x=orx
          //-----------------------------------------------------
          if( (fabs(K[iv1].x - xmin) < step)   && (fabs(K[iv2].x - xmin) < step)
          && (fabs(K[iv3].x - xmin) < step) ) {

            for(int j=0;j<3;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 1);
            borderelements[kf++].set(vertices, indKV, 1);
          }

          //  Faces on the plane x=orx + glx
          //-----------------------------------------------------
          else if( (fabs(K[iv1].x - xmax) < step)   && (fabs(K[iv2].x - xmax) < step)
          && (fabs(K[iv3].x - xmax) < step) ) {

            for(int j=0;j<3;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 2);
            borderelements[kf++].set(vertices, indKV, 2);
          }

          //  Faces on the plane y=ory
          //-----------------------------------------------------
          else if( (fabs(K[iv1].y - ymin) < step)   && (fabs(K[iv2].y - ymin) < step)
          && (fabs(K[iv3].y - ymin) < step) ) {

            for(int j=0;j<3;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 3);
            borderelements[kf++].set(vertices, indKV,3);
          }

          //  Faces on the plane y=ory + gly
          //-----------------------------------------------------
          else if( (fabs(K[iv1].y - ymax) < step)   && (fabs(K[iv2].y - ymax) < step)
          && (fabs(K[iv3].y - ymax) < step) ) {

            for(int j=0;j<3;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 4);
            borderelements[kf++].set(vertices, indKV,4);
          }
          //  Faces on the plane z=orz
          //-----------------------------------------------------
          else if( (fabs(K[iv1].z - zmin) < step)   && (fabs(K[iv2].z - zmin) < step)
          && (fabs(K[iv3].z - zmin) < step) ) {

            for(int j=0;j<3;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 5);
            borderelements[kf++].set(vertices, indKV,5);
          }

          //  Faces on the plane z=orz + glz
          //-----------------------------------------------------
          else if( (fabs(K[iv1].z - zmax) < step)   && (fabs(K[iv2].z - zmax) < step)
          && (fabs(K[iv3].z - zmax) < step) ) {

            for(int j=0;j<3;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 6);
            borderelements[kf++].set(vertices, indKV,6);
          }
        } // end loop over faces
      }

      assert(kf == nbe);


      BuildBound();
      if(nt > 0){
        BuildAdj();
        //     Buildbnormalv();
        BuildjElementConteningVertex();
      }



      // std::cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << std::endl;
      // std::cout << "  -- Mesh3 : "<< " d "<< 3  << ", n Tet " << nt << ", n Vtx "
      // 	    << nv << " n Bord " << nbe << std::endl;
      ffassert(mes>=0);

    }

    // To build mesh for time problem
    Mesh3::Mesh3( TimeInterface3 &gamma) {

      const Uint nInterface = gamma.size();
      const Mesh3& backMesh(*(gamma[0]->backMesh));                   // the back-ground mesh


      for(Uint i=0;i<nInterface-1;++i) {             // all interface must have the same backMesh
        assert(gamma[i]->backMesh == gamma[i+1]->backMesh);
      }

      int mmt = 0;
      for(Uint i=0;i<nInterface;++i) {
        mmt += gamma[i]->nbElement();            // exact approximation (2D)
      }
      const int mt = backMesh.nt;//2*mmt;
      const int mv = 4*mt;                         // large approximation
      const int mbe = 0;

      if (mt==0) { nt=0;nv=0;nbe=0; return;}

      this->set(mv,mt,mbe);
      nv=0;  nt=0;
      ElementIndexInBackMesh = new Uint[mt];
      VertexIndexInBackMesh = new Uint[mv];
      ElementIndexInLocMesh = new Uint[backMesh.nt];


      KN<int> foundVertex(backMesh.nv); foundVertex = -1;
      KN<byte> createdTet(backMesh.nt); createdTet = static_cast<byte>(0);
      int indT[4];


      for(Uint iInterface=0; iInterface<nInterface; ++iInterface) {
        for(Uint ifac=0; ifac != gamma[iInterface]->nbElement(); ++ifac) {

          const int k = gamma[iInterface]->idxElementOfFace(ifac);
          const Element & K = (backMesh)[k];

          if(createdTet(k)) continue;                 // already created

          for( int i=0; i<4;++i) {
            const int idxG = (backMesh)(K[i]);

            if(foundVertex(idxG) == -1) {

              foundVertex(idxG) = nv;
              vertices[nv].x = K[i].x;
              vertices[nv].y = K[i].y;
              vertices[nv].z = K[i].z;
              VertexIndexInBackMesh[nv] = idxG;
              indT[i] = nv++;
            }
            else {
              indT[i] = foundVertex(idxG);
            }
          }
          ElementIndexInBackMesh[nt] = k;
          ElementIndexInLocMesh[k] = nt;

          elements[nt++].set(vertices, indT, 0);
          createdTet(k) = static_cast<byte>(1);
        }
      }

      // we need to get the elements that change domain during the time slab
      for(int k=0; k<backMesh.nt;++k) {
        if(createdTet(k) == static_cast<byte>(0)) {

          const Element & K = (backMesh)[k];

          // We only need to look at the sign of one of the node
          bool changeSign = false;    // excluded by default
          for(Uint i=0; i<nInterface-1;++i) {
            int idx1 = (backMesh)(K[i]);
            // if( (gamma[i].levelSet)[idx1] * (gamma[i+1].levelSet)[idx1] <= 0 )
            if( (gamma[i]->ls_sign)[idx1] * (gamma[i+1]->ls_sign)[idx1] <= 0 )
            changeSign = true;
          }
          if(changeSign) {
            for( int i=0; i<4;++i) {
              const int idxG = (backMesh)(K[i]);

              if(foundVertex(idxG) == -1) {

                foundVertex(idxG) = nv;
                vertices[nv].x = K[i].x;
                vertices[nv].y = K[i].y;
                vertices[nv].z = K[i].z;
                VertexIndexInBackMesh[nv] = idxG;
                indT[i] = nv++;
              }
              else {
                indT[i] = foundVertex(idxG);
              }
            }
            if(nt>=mt) assert(0);//ElementIndexInBackMesh.resize(2*mt);//assert(0);

            createdTet(k) = static_cast<byte>(1);
            ElementIndexInBackMesh[nt] = k;
            ElementIndexInLocMesh[k] = nt;
            elements[nt++].set(vertices, indT, 0);
          }
        }
      }
      BuildBound();
      if(nt > 0){
        BuildAdj();
        BuildjElementConteningVertex();
      }


    }

    Mesh3::Mesh3(const Interface& gamma) {

      const Mesh3& backMesh(*gamma.backMesh);                  // the back-ground mesh

      const int mt = gamma.nbElement();            // large approximation
      const int mv = 4*mt;                         // very large approximation
      const int mbe = 0;

      if (mt==0) { nt=0;nv=0;nbe=0; return;}

      this->set(mv,mt,mbe);
      nv=0;  nt=0;
      ElementIndexInBackMesh = new Uint[mt];
      VertexIndexInBackMesh = new Uint[mv];
      ElementIndexInLocMesh = new Uint[backMesh.nt];


      KN<int> foundVertex(backMesh.nv); foundVertex = -1;
      KN<byte> createdTet(backMesh.nt); createdTet = static_cast<byte>(0);
      int indT[4];



      for(Uint ifac=0; ifac != gamma.nbElement(); ++ifac) {

        const int k = gamma.idxElementOfFace(ifac);
        const Element & K = (backMesh)[k];

        if(createdTet(k)) continue;                 // already created

        for( int i=0; i<4;++i) {
          const int idxG = (backMesh)(K[i]);

          if(foundVertex(idxG) == -1) {

            foundVertex(idxG) = nv;
            vertices[nv].x = K[i].x;
            vertices[nv].y = K[i].y;
            vertices[nv].z = K[i].z;
            VertexIndexInBackMesh[nv] = idxG;
            indT[i] = nv++;
          }
          else {
            indT[i] = foundVertex(idxG);
          }
        }
        ElementIndexInBackMesh[nt] = k;
        ElementIndexInLocMesh[k] = nt;

        elements[nt++].set(vertices, indT, 0);
        createdTet(k) = static_cast<byte>(1);
      }

      BuildBound();
      if(nt > 0){
        BuildAdj();
        BuildjElementConteningVertex();
      }

    }




MeshHexa::MeshHexa(int nx, int ny, int nz, R orx, R ory, R orz, R lx, R ly, R lz) {

  const int idQ[8] = {0,1,3,2,4,5,7,6};
  int mv  = nx * ny * nz;
  int mt  = (nx - 1)*(ny - 1)*(nz - 1);
  int mbe = 2 * ((nx - 1) * (ny - 1) + (nx - 1) * (nz - 1) +(nz - 1) * (ny - 1) );
  const R hx = lx / (nx - 1);
  const R hy = ly / (ny - 1);
  const R hz = lz / (nz - 1);

  this->set(mv,mt,mbe);

  KN<int> iv(8),  indT(8);

  int jt = 0;
  for(int k=0; k<nz-1; k++) {
    for(int j=0; j<ny-1; j++) {
      for(int i=0; i<nx-1; i++) {

        int id=0;
        for(int kk=k; kk<k+2;++kk) {
          for(int jj=j; jj<j+2; ++jj) {
            for(int ii=i; ii<i+2; ++ii) {

              int ivl  =  ii + jj*nx + kk*nx*ny;                    // index
              iv(id++) = ivl;

              vertices[ivl].x = ii*hx + orx;
              vertices[ivl].y = jj*hy + ory;
              vertices[ivl].z = kk*hz + orz;
            }
          }
        }
        for(int e=0; e<8; ++e) {
          indT(e) = iv(idQ[e]);
        }
        elements[jt++].set(vertices, indT, 0);
      }
    }
  }

  assert(jt == nt);


  const R xmin = orx, ymin = ory, zmin = orz;
  const R xmax = orx + lx, ymax = ory + ly, zmax = orz + lz;

  // create the for borders
  int kf=0;
  int indV[4],indKV[4];
  const R step = 1e-12;
  for( int jt=0; jt<nt; ++jt) {                     // loop over element
    const Element & K((*this)[jt]);
    for (int k=0;k<6;++k) {                         // loop over faces
      for( int i=0;i<4;++i){
        indV[i] = Element::nvface[k][i];
        indKV[i] = (*this)(K[indV[i]]);
      }

      const int iv1 = Element::nvface[k][0];
      const int iv2 = Element::nvface[k][1];
      const int iv3 = Element::nvface[k][2];
      const int iv4 = Element::nvface[k][3];

      //  Faces on the plane x=orx
      //-----------------------------------------------------
      if( (fabs(K[iv1].x - xmin) < step) && (fabs(K[iv2].x - xmin) < step)
      &&  (fabs(K[iv3].x - xmin) < step) && (fabs(K[iv4].x - xmin) < step) ) {

        for(int j=0;j<4;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 1);
        borderelements[kf++].set(vertices, indKV, 1);
      }

      //  Faces on the plane x=orx + glx
      //-----------------------------------------------------
      else if( (fabs(K[iv1].x - xmax) < step)   && (fabs(K[iv2].x - xmax) < step)
      && (fabs(K[iv3].x - xmax) < step) && (fabs(K[iv4].x - xmax) < step) ) {

        for(int j=0;j<4;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 2);
        borderelements[kf++].set(vertices, indKV, 2);
      }

      //  Faces on the plane y=ory
      //-----------------------------------------------------
      else if( (fabs(K[iv1].y - ymin) < step)   && (fabs(K[iv2].y - ymin) < step)
      && (fabs(K[iv3].y - ymin) < step) && (fabs(K[iv4].y - ymin) < step)) {

        for(int j=0;j<4;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 3);
        borderelements[kf++].set(vertices, indKV,3);
      }

      //  Faces on the plane y=ory + gly
      //-----------------------------------------------------
      else if( (fabs(K[iv1].y - ymax) < step)   && (fabs(K[iv2].y - ymax) < step)
      && (fabs(K[iv3].y - ymax) < step) && (fabs(K[iv4].y - ymax) < step)) {

        for(int j=0;j<4;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 4);
        borderelements[kf++].set(vertices, indKV,4);
      }
      //  Faces on the plane z=orz
      //-----------------------------------------------------
      else if( (fabs(K[iv1].z - zmin) < step)   && (fabs(K[iv2].z - zmin) < step)
      && (fabs(K[iv3].z - zmin) < step) && (fabs(K[iv4].z - zmin) < step) ) {

        for(int j=0;j<4;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 5);
        borderelements[kf++].set(vertices, indKV,5);
      }

      //  Faces on the plane z=orz + glz
      //-----------------------------------------------------
      else if( (fabs(K[iv1].z - zmax) < step)   && (fabs(K[iv2].z - zmax) < step)
      && (fabs(K[iv3].z - zmax) < step) && (fabs(K[iv4].z - zmax) < step)) {

        for(int j=0;j<4;++j) vertices[indKV[j]].lab = std::max(vertices[indKV[j]].lab, 6);
        borderelements[kf++].set(vertices, indKV,6);
      }
    } // end loop over faces
  }

  assert(kf == nbe);


  BuildBound();
  if(nt > 0){
    BuildAdj();
    //     Buildbnormalv();
    BuildjElementConteningVertex();
  }

  ffassert(mes>=0);

}




    //
    // Mesh3::Mesh3(const Mesh3& Th, std::string whatToDo) {
    //   this->levelRefinement = Th.levelRefinement+1;
    //   assert(whatToDo == "refine");
    //   int of = Th.nv+10;
    //   unsigned int tinfty=-1;
    //   const int nkeys=6+4;                 // max number of items
    //   const int nkv= Element::nv;
    //   const int nke= Element::ne;
    //   const int nnodeK = nkv + nke;
    //   int nbNodes=0;
    //
    //   int mt = Element::nref * Th.nt;
    //   int mv = Element::nref * Th.nv;
    //   int mbe = Th.nbe * BorderElement::nref;
    //
    //   this->set(mv,mt,mbe);
    //   ElementIndexInCoarseMesh = new Uint[mt];
    //
    //   typedef SortArray<unsigned int,2> Key;
    //   Key keys[nkeys];
    //   int keysdim[nkeys];
    //   HashTable<Key,int> h(nnodeK*Th.nt,of+Th.nt);
    //   int idx[nkeys];
    //
    //   int kt = 0;
    //   for(int k=0;k<Th.nt;++k) {
    //     int m = 0;
    //     const Element & K(Th[k]);                            // the element
    //     for(int i=0;i<nkv;++i) {      // loop over the nodes
    //       keysdim[m]=0;               // what item
    //       keys[m++]=Key(Th(K[i]),tinfty);  // index saved
    //     }
    //     for(int i=0;i<nke;++i) {      // loop over edges
    //       keysdim[m]=1;               // what item
    //       keys[m++]=Key(Th(K[Element::nvedge[i][0]]),
    // 		    Th(K[Element::nvedge[i][1]])); // sorted id
    //     }
    //
    //     for(int i=0;i<m;i++) {
    //       Key ki=keys[i];//,kio=ki;
    //       typename HashTable<Key,int>::iterator pk= h.find(ki); // look for the key
    //       if(!pk) {                          // if not found
    //
    // 	if(keysdim[i] == 0) { // on node
    // 	  for(int d=0;d < Rd::d;++d) vertices[nbNodes][d] = K[i][d];
    // 	  // M::vertices[nbNodes].lab = K[i].lab;
    // 	  // std::cout << vertices[nbNodes]  << std::endl;
    // 	}
    // 	else {
    // 	  int j = i - nkv;
    // 	  const int i0 = Element::nvedge[j][0], i1 = Element::nvedge[j][1];
    // 	  for(int d=0;d < Rd::d;++d) vertices[nbNodes][d] = 0.5*(K[i0][d] + K[i1][d]);
    // 	  // M::vertices[nbNodes].lab = K[i].lab;
    // 	  // std::cout << vertices[nbNodes]  << std::endl;
    // 	}
    //
    // 	pk = h.add(ki,nbNodes++);        // add the key
    //       }
    //       idx[i] = pk->v;
    //     }
    //     int idxK[nkv];
    //
    //     for(int i=0;i<Element::nref;++i) {
    //       for(int j=0;j<nkv;++j) idxK[j] = idx[Element::refElement[i][j]];
    //       ElementIndexInCoarseMesh[kt] = k;
    //       elements[kt++].set(vertices, idxK, 0);
    //       // std::cout << elements[kt-1] << std::endl;
    //     }
    //   }
    //   nv = nbNodes;
    //   assert(kt == nt);
    //
    //
    //
    //   //  take care of the boundary elements
    //   int n = 0;
    //   for(int ke=0; ke<Th.nbe;++ke) {                    // loop over boundary element
    //     int kf, k = Th.BoundaryElement(ke,kf);           // k is the element index
    //     const Element & K(Th[k]);
    //     const BorderElement & Be(Th.be(ke));
    //
    //
    //     if(Rd::d == 2) {
    //       int idxBE[3];
    //
    //       Key k0(Th(K[Element::nvedge[kf][0]]),tinfty);
    //       typename HashTable<Key,int>::iterator pk = h.find(k0);
    //       idxBE[0] = pk->v;
    //
    //       Key k1(Th(K[Element::nvedge[kf][1]]),tinfty);
    //       pk =  h.find(k1);
    //       idxBE[2] = pk->v;
    //
    //       Key ki(Th(Be[0]), Th(Be[1]));
    //       pk =  h.find(ki);
    //       if(!pk) assert(0);
    //
    //       idxBE[1] = pk->v;
    //
    //       borderelements[n++].set(vertices,idxBE  , Be.lab);
    //       borderelements[n++].set(vertices,idxBE+1, Be.lab);
    //
    //     }
    //     else {
    //       int idxK[6];
    //       for(int i=0;i<BorderElement::nv;++i) {
    //       	Key k0(Th(Be[i]),tinfty);
    //       	typename HashTable<Key,int>::iterator pk = h.find(k0);
    //       	idxK[i] = pk->v;
    //       }
    //       for(int i=0;i<BorderElement::ne;++i) {
    //
    //       	Key k0(Th(Be[BorderElement::nvedge[i][0]]),
    //     	       Th(Be[BorderElement::nvedge[i][1]]));
    //
    //     	typename HashTable<Key,int>::iterator pk = h.find(k0);
    //     	idxK[i+3] = pk->v;
    //       }
    //       int idxBE[3];
    //       for(int j=0;j<BorderElement::nref;++j) {
    //       	for(int i=0;i<3;++i) idxBE[i] = idxK[BorderElement::refElement[j][i]];
    //
    //       	borderelements[n++].set(vertices,idxBE, Be.lab);
    //
    //       }
    //     }
    //   }
    //
    //
    //   BuildBound();
    //   BuildAdj();
    // //   Buildbnormalv();
    //   BuildjElementConteningVertex();
    //
    //   std::cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << std::endl;
    //   std::cout << "  -- Mesh3 : "<< " d "<< 3  << ", n Tet " << nt << ", n Vtx "
    //   	    << nv << " n Bord " << nbe << std::endl;
    //   ffassert(mes>=0);
    //
    // }
    //
