#ifndef _GENERIC_ELEMENT_HPP
#define _GENERIC_ELEMENT_HPP
#include <cassert>
#include "global.hpp"
#include "GenericVertex.hpp"

inline R1 ExtNormal(GenericVertex<R1> *const v[2], const std::vector<int> &f) {
   return f[0] == 0 ? R1(-1) : R1(1);
}
inline R2 ExtNormal(GenericVertex<R2> *const v[3], const std::vector<int> &f) {
   return R2(*v[f[1]], *v[f[0]]).perp();
}
inline R3 ExtNormal(GenericVertex<R3> *const v[4], const std::vector<int> &f) {
   return R3(*v[f[0]], *v[f[2]]) ^ R3(*v[f[0]], *v[f[1]]);
}

template <typename Data> class GenericElement : public Label {
 public:
   typedef typename Data::V Vertex;
   typedef typename Data::Face Face;
   typedef typename Data::V::Rd Rd;
   typedef typename Data::RdHat RdHat;
   typedef typename Data::RdHatBord RdHatBord;
   typedef typename Rd::R R;

   static const int nv              = Data::NbOfVertices;
   static const int ne              = Data::NbOfEdges;
   static const int nf              = Data::NbOfFaces;
   static const int nt              = Data::NbOfTet;
   static const int nitem           = nv + ne + nf + nt;
   static const int nva             = Data::NbOfVertexOnHyperFace;
   static const int nea             = Data::NbOfAdjElem;
   static const int d               = Rd::d;
   static const int nvOnFace        = Data::NvOnFace;
   static const int ParaviewNumCell = Data::ParaviewNumCell;
   static const int nb_ntcut        = Data::NbNtCut;
   static const int nvc             = Data::NbOfVerticesCut;
   static const int nb_sign_pattern = Data::NbSignPattern;

   static const std::vector<std::vector<int>> nvedge;
   static const std::vector<std::vector<int>> nvface;
   static const std::vector<std::vector<int>> edgeOfFace;
   static const std::vector<std::vector<int>> faceOfEdge;
   static const std::vector<std::vector<int>> nvhyperFace;
   static const std::vector<std::vector<int>> commonVertOfEdges;
   static std::vector<int> itemTopology() { return {nv, ne, nf, nt}; }
   static int oppVertOfEdge(int edge, int vert) {
      return vert == nvedge[edge][0] ? nvedge[edge][1] : nvedge[edge][0];
   }

 public:
   Vertex *vertices[nv];
   R mes;

 public:
   GenericElement() {}

   const Vertex &operator[](int i) const {
      assert(i >= 0 && i < nv);
      return *vertices[i];
   }

   Vertex &operator[](int i) {
      assert(i >= 0 && i < nv);
      return *vertices[i];
   }

   const Vertex &at(int i) const { return *vertices[i]; }

   Vertex &at(int i) { return *vertices[i]; }

   GenericElement &set(Vertex *v0, int *iv, int r,
                       double mss = globalVariable::UnSetMesure) {
      for (int i = 0; i < nv; ++i)
         vertices[i] = v0 + iv[i];
      mes = (mss != globalVariable::UnSetMesure) ? mss : Data::mesure(vertices);
      lab = r;
      if (mss != globalVariable::UnSetMesure || mes <= 0) {
         for (int i = 0; i < nv; ++i)
            std::cout << *vertices[i] << std::endl;
         std::cout << mss << "\t" << mes << std::endl;
         getchar();
      }
      assert(mss == globalVariable::UnSetMesure && mes > 0);
      return *this;
   }

   void set_face(int ifac, Face &face) const {
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
      return ExtNormal(vertices, nvhyperFace[i]) /
             (ExtNormal(vertices, nvhyperFace[i]).norm());
   }
   Rd N_notNormalized(int i) const {
      return ExtNormal(vertices, nvhyperFace[i]);
   }

   Rd PBord(int i, RdHatBord P) const { return Data::PBord(nvhyperFace[i], P); }

   // THIS IS DIFFERENT FOR RECTANGLES
   virtual Rd operator()(const RdHat &Phat) const = 0;

   Rd barycenter() const {
      Rd Q;
      for (int i = 1; i < nv; ++i)
         Q += *(Rd *)vertices[i];
      return 1. / nv * Q;
   }

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

 private:
   // pas de copie
   GenericElement(const GenericElement &);
   GenericElement &operator=(const GenericElement &);
};

#endif
