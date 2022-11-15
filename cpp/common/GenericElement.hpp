#ifndef _GENERIC_ELEMENT_HPP
#define _GENERIC_ELEMENT_HPP

#include "Label.hpp"
#include "R3.hpp"
#include "GenericVertex.hpp"

const double UnSetMesure = -1e+200;

inline R1 ExtNormal(GenericVertex<R1> *const v[2], int const f[1]) {
   return f[0] == 0 ? R1(-1) : R1(1);
}
inline R2 ExtNormal(GenericVertex<R2> *const v[3], int const f[2]) {
   return R2(*v[f[1]], *v[f[0]]).perp();
}
inline R3 ExtNormal(GenericVertex<R3> *const v[4], int const f[3]) {
   return R3(*v[f[0]], *v[f[2]]) ^ R3(*v[f[0]], *v[f[1]]);
}

template <typename Data> class GenericElement : public Label {
 public:
   typedef typename Data::V Vertex;
   typedef typename Data::Face Face;
   typedef typename Data::V::Rd Rd;
   typedef typename Data::RdHat RdHat;         // for parametrization
   typedef typename Data::RdHatBord RdHatBord; // for parametrization
   typedef typename Rd::R R;

   static const int nv       = Data::NbOfVertices; // nb  of vertices
   static const int ne       = Data::NbOfEdges;    // nb  of edges
   static const int nf       = Data::NbOfFaces;    // nb of faces
   static const int nt       = Data::NbOfTet;      // nb of tets
   static const int nitem    = nv + ne + nf + nt;
   static const int nva      = Data::NbOfVertexOnHyperFace;
   static const int nea      = Data::NbOfAdjElem;
   static const int d        = Rd::d;
   static const int nref     = Data::NbOfRef;
   static const int nvOnFace = Data::NvOnFace;

   static const int (*const nvedge)[2];        // idx nodes on edges
   static const int (*const nvface)[nvOnFace]; // idx nodes on faces

   static const int (*const nvhyperFace)[nva];

   static const int (*const onWhatBorder)[nitem];   //
   static const int (*const commonVertOfEdges)[ne]; //
   static const int (*const refElement)[nv];        //
   static const int (*const faceOfEdge)[2];
   static const int (*const edgeOfFace)[nvOnFace];

   static const int (*const nvadj)[nva]; //
   static const int nitemdim[4];
   static const int ParaviewNumCell = Data::ParaviewNumCell;

   // for cut
   static const int nb_ntcut        = Data::NbNtCut;
   static const int nvc             = Data::NbOfVerticesCut;
   static const int nb_sign_pattern = Data::NbSignPattern;

   static int oppVertOfEdge(int edge, int vert) {
      return vert == nvedge[edge][0] ? nvedge[edge][1] : nvedge[edge][0];
   }

 public:
   Vertex *vertices[nv]; // array 3 pointer to vertex
   R mes;

 public:
   GenericElement() {}
   // GenericElement(Vertex * v0,int * iv,int r, double mss=UnSetMesure)
   // {
   //   for(int i=0;i<nv;++i)
   //     vertices[i]=v0+iv[i];
   //   mes=(mss!=UnSetMesure) ? mss : Data::mesure(vertices);
   //   lab=r;
   //   assert(mss==UnSetMesure && mes>0);
   // }

   const Vertex &operator[](int i) const {
      assert(i >= 0 && i < nv);
      return *vertices[i];
   } // to see triangle as a array of vertex

   Vertex &operator[](int i) {
      assert(i >= 0 && i < nv);
      return *vertices[i];
   } // to see triangle as a array of vertex

   const Vertex &at(int i) const {
      return *vertices[i];
   } // to see triangle as a array of vert

   Vertex &at(int i) {
      return *vertices[i];
   } // to see triangle as a array of vert

   GenericElement &set(Vertex *v0, int *iv, int r, double mss = UnSetMesure) {
      for (int i = 0; i < nv; ++i)
         vertices[i] = v0 + iv[i];
      mes = (mss != UnSetMesure) ? mss : Data::mesure(vertices);
      lab = r;
      if (mss != UnSetMesure || mes <= 0) {
         for (int i = 0; i < nv; ++i)
            std::cout << *vertices[i] << std::endl;
         std::cout << mss << "\t" << mes << std::endl;
         getchar();
      }
      assert(mss == UnSetMesure && mes > 0);
      return *this;
   }

   void set_face(int ifac, Face &face) const {
      // int iv[nva];
      // for(int i=0;i<nva;++i) {iv[i]= nvface[ifac][i];}
      for (int i = 0; i < nva; ++i) {
         face.vertices[i] = vertices[nvface[ifac][i]];
      }
      face.mes = Data::mesure(face.vertices);
      face.lab = 0;
   }

   std::istream &Read1(std::istream &f, Vertex *v0, int n) {
      int iv[nv], ir, err = 0;
      for (int i = 0; i < nv; ++i) {
         f >> iv[i];
         iv[i]--;
         if (!(iv[i] >= 0 && iv[i] < n))
            err++;
      }
      f >> ir;
      if (err || !f.good()) {
         std::cerr << " Erreur GenericElement::Read1 " << nv << " " << n
                   << "  : ";
         for (int j = 0; j < nv; ++j)
            std::cerr << iv[j] << " ";
         std::cerr << " , " << ir << std::endl;
         std::abort();
      }

      set(v0, iv, ir);
      return f;
   }

   Rd Edge(int i) const {
      assert(i >= 0 && i < ne);
      return Rd(at(nvedge[i][0]), at(nvedge[i][1]));
   } // opposite edge vertex i

   Rd N(int i) const {
      return ExtNormal(vertices, nvadj[i]) /
             (ExtNormal(vertices, nvadj[i]).norm());
   }
   Rd N_notNormalized(int i) const { return ExtNormal(vertices, nvadj[i]); }

   Rd PBord(int i, RdHatBord P) const { return Data::PBord(nvadj[i], P); }

   // THIS IS DIFFERENT FOR RECTANGLES
   virtual Rd operator()(const RdHat &Phat) const = 0;
   //{
   //   Rd r= (1.-Phat.sum())*(*(Rd*) vertices[0]);
   //   for (int i=1;i<nv;++i)
   //     r+=  Phat[i-1]*(*(Rd*) vertices[i]);
   //   return r;
   // }
   Rd barycenter() const {
      Rd Q;
      for (int i = 1; i < nv; ++i)
         Q += *(Rd *)vertices[i];
      return 1. / nv * Q;
   }

   int faceOrient(int i) const { // def the permutatution of orient the face
      int fo             = 1;
      const Vertex *f[3] = {&at(nvface[i][0]), &at(nvface[i][1]),
                            &at(nvface[i][2])};
      if (f[0] > f[1])
         fo = -fo, Exchange(f[0], f[1]);
      if (f[1] > f[2]) {
         fo = -fo, Exchange(f[1], f[2]);
         if (f[0] > f[1])
            fo = -fo, Exchange(f[0], f[1]);
      }
      return fo;
   }

   // THIS IS DIFFERENT FOR RECTANGLES

   // Permutation
   //  static const int PERM_FACE[6][3] =
   //  {{0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}};
   // virtual int facePermutation(int i) const = 0;
   // // {//def the permutatution of orient the face
   // //   const Vertex * f[3]={&at(nvface[i][0]), &at(nvface[i][1]),
   // &at(nvface[i][2])};
   // //   if(f[0] < f[1]) {
   // //     if(f[1] < f[2]) return 0;      // 0,1,2
   // //     if(f[0] < f[2]) return 1;      // 0,2,1
   // //     else return 4;                 // 2,0,1
   // //   }
   // //   if(f[0] < f[2]) return 2;        // 1 0 2
   // //   if(f[1] < f[2]) return 3;        // 1,2,0
   // //   else return 5;                   // 2,1,0
   // //
   // // }

   int EdgeOrientation(int i) const {
      return 2 * (&at(nvedge[i][0]) < &at(nvedge[i][1])) - 1;
   }

   int edgePermutation(int i) const {
      return &at(nvedge[i][1]) < &at(nvedge[i][0]);
   } // 0 : no permutation

   R lenEdge(int i) const {
      assert(i >= 0 && i < 3);
      Rd E = Edge(i);
      return sqrt((E, E));
   }

   R hElement() const {
      double h = 0;
      for (int i = 0; i < ne; ++i)
         h += lenEdge(i);
      return h / ne;
   }
   R hMax() const {
      double h = 0;
      for (int i = 0; i < ne; ++i)
         h = max(h, lenEdge(i));
      return h;
   }
   R mesure() const { return mes; }
   R measure() const { return mes; }
   R get_h() const {
      double h = 0;
      for (int i = 0; i < ne; ++i)
         h += lenEdge(i);
      return h / ne;
   }
   Rd map(const RdHat &Phat) const { return (*this)(Phat); }
   Rd mapToPhysicalElement(const RdHat &Phat) const { return (*this)(Phat); }

   //   static  int NbNodes(int c)  // use the bit i of c to say if node in
   //   objet of dim  i existe { int c0=(c&1)!=0, c1=(c&2)!=0, c2=(c&4)!=0,
   //   c3=(c&8)!=0;
   //     return nv*c0 +ne*c1 +nf*c2 + nt*c3 ;}

   //   static  int NbNodes(const int c[4])  // use the bit i of c to say if
   //   node in objet of dim  i existe { int c0=(c[0])!=0, c1=(c[1])!=0,
   //   c2=(c[2])!=0, c3=(c[3])!=0;
   //     return nv*c0 +ne*c1 +nf*c2 + nt*c3 ;}

   static int NbNodes(const int c[4]) {
      return nv * c[0] + ne * c[1] + nf * c[2] + nt * c[3];
   }

   void Renum(Vertex *v0, int *r) {
      for (int i = 0; i < nv; i++)
         vertices[i] = v0 + r[vertices[i] - v0];
   }

   void Change(Vertex *vold, Vertex *vnew) {
      for (int i = 0; i < nv; i++)
         vertices[i] = vnew + vertices[i] - vold;
   }

   // Rd n(int i) const //  unit exterior normal
   //   {Rd E=Edge(i);return Rd(E.y,-E.x)/Norme2(E);}
   //   friend std::ostream& operator <<(std::ostream& f, const GenericElement &
   //   K )
   //   {
   //     for(int i=0;i<nv;++i)
   //       for(int j=i+1;j<nv;++j)
   // 	f << (R3) K[i] << "\n" << (R3) K[j] << std::endl << std::endl;
   //     return f;
   //   }

 private:
   // pas de copie
   GenericElement(const GenericElement &);
   GenericElement &operator=(const GenericElement &);
};

#endif
