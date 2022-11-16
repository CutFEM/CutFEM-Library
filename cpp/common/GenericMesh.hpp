#ifndef GENERICMESH_HPP_
#define GENERICMESH_HPP_

extern long verbosity;

#include "cassert"
#include "../util/util.hpp"
#include <cstdlib>

#include "dataStruct1D.hpp"
#include "dataStruct2D.hpp"
#include "dataStruct3D.hpp"
#include "sort_array.hpp"
#include "cutFEMConfig.h"
#ifdef USE_MPI
#include "../parallel/cfmpi.hpp"
#endif

inline int maxdfon(const int *dfon) {
   return std::max(std::max(dfon[0], dfon[1]), std::max(dfon[2], dfon[3]));
}

const int NbTypeItemElement = 4;
const int TypeVertex        = 0;
const int TypeEdge          = 1;
const int TypeFace          = 2;
const int TypeVolume        = 3;

template <typename T, typename B, typename V> class GenericMesh {
 public:
   typedef GenericMesh GMesh;
   typedef T Element;
   typedef typename V::Rd Rd;
   typedef typename Rd::R R;
   typedef V Vertex;
   typedef B BorderElement;
   typedef typename Element::RdHat RdHat; // for parametrization
   typedef typename Element::Face Face;

   static const int nea =
       T::nea; //  numbering of adj (4 in Tet,  3 in Tria, 2 in seg)
   static const int nva = T::nva; //  numbering of vertex in Adj hyperface

   int nt, nv, nbe;
   R mes, mesb;

   // protected
   V *vertices;
   T *elements;
   B *borderelements;

   int *TheAdjacencesLink;       // to store the adj link  k*nea+i -> k'*nea+i'
   int *BoundaryElementHeadLink; //

 public:
   int nbElmts() const { return nt; }
   int get_nb_element() const { return nt; }
   int nbBrdElmts() const { return nbe; }
   int nbVertices() const { return nv; }
   int nbElements() const { return nt; }
   int NbElement() const { return nt; }

   int nbBorderElements() const { return nbe; }

   const T &operator[](int i) const { return elements[CheckT(i)]; }
   const V &operator()(int i) const { return vertices[CheckV(i)]; }
   const B &be(int i) const { return borderelements[CheckBE(i)]; }

   T &t(int i) { return elements[CheckT(i)]; }
   V &v(int i) { return vertices[CheckV(i)]; }
   B &be(int i) { return borderelements[CheckBE(i)]; }

   GenericMesh()
       : nt(0), nv(0), nbe(0), mes(0.), mesb(0.), vertices(0), elements(0),
         borderelements(0), TheAdjacencesLink(0), BoundaryElementHeadLink(0) {}

   virtual void info() {
      std::cout << " ----- Mesh " << this << " info ----- " << std::endl;
      std::cout << " nb of nodes            : \t" << nv << std::endl;
      std::cout << " nb of elements         : \t" << nt << std::endl;
      std::cout << " nb of border elements  : \t" << nbe << std::endl;
   }

   void set(int mv, int mt, int mbe) {
      assert(nt == 0 && nv == 0 && nbe == 0);
      nt             = mt;
      nv             = mv;
      nbe            = mbe;
      vertices       = new V[nv];
      elements       = new T[nt];
      borderelements = new B[nbe];

      assert(nt >= 0 && elements);
      assert(nv > 0 && vertices);
   }

#ifdef USE_MPI
   virtual int first_element() const {
      return MPIcf::first_element(this->nbElements());
   }
   virtual int next_element() const {
      return MPIcf::next_element(this->nbElements());
   }
   virtual int last_element() const {
      return MPIcf::last_element(this->nbElements());
   }

   virtual int first_boundary_element() const { return MPIcf::my_rank(); }
   virtual int next_boundary_element() const { return MPIcf::size(); }
   virtual int last_boundary_element() const { return this->nbBrdElmts(); }
#else
   virtual int first_element() const { return 0; }
   virtual int next_element() const { return 1; }
   virtual int last_element() const { return this->nbElements(); }

   virtual int first_boundary_element() const { return 0; }
   virtual int next_boundary_element() const { return 1; }
   virtual int last_boundary_element() const { return this->nbBrdElmts(); }
#endif

   int operator()(const T &tt) const { return CheckT(&tt - elements); }
   int operator()(const T *tt) const { return CheckT(tt - elements); }
   int operator()(const V &vv) const { return CheckV(&vv - vertices); }
   int operator()(const V *vv) const { return CheckV(vv - vertices); }
   int operator()(const B &k) const { return CheckBE(&k - borderelements); }
   int operator()(const B *k) const { return CheckBE(k - borderelements); }

   int operator()(int it, int j) const {
      return operator()(elements[it][j]);
   } // Nu vertex j of triangle it
   int at(int it, int j) const {
      return operator()(elements[it][j]);
   } // Nu vertex j of triangle it
   int be(int it, int j) const {
      return operator()(borderelements[it][j]);
   } // Nu vertex j of triangle it

   int CheckV(int i) const {
      assert(i >= 0 && i < nv);
      return i;
   }
   int CheckT(int i) const {
      assert(i >= 0 && i < nt);
      return i;
   }
   int CheckBE(int i) const {
      assert(i >= 0 && i < nbe);
      return i;
   }

   void BuildAdj();
   void Buildbnormalv();
   void BuildBound();

   int ElementAdj(int k, int &j) const {
      int p = TheAdjacencesLink[nea * k + j];
      j     = p % nea;
      return p >= 0 ? p / nea : -1;
   }

   int GetAllElementAdj(int it, int *tabk)
       const { //  get the tab of all adj element (std::max ne)
      int i = 0;
      for (int j = 0; j < nea; ++j) {
         int p = TheAdjacencesLink[nea * it + j];
         if (p >= 0 && p / nea != it) {
            tabk[i] = p / nea;
            i++;
         }
      }
      return i;
   }

   bool isOnBorder(int k) const {
      int i = 0;
      for (int j = 0; j < nea; ++j) {
         int n = TheAdjacencesLink[3 * k + j] / 3;
         if (n >= 0 && n != k)
            i++;
      }
      return (i != nea);
   }

   int BoundaryElement(int bbe, int &ItemInK) const {
      int i   = BoundaryElementHeadLink[bbe];
      ItemInK = i % nea;
      return i / nea;
   }

   int BoundaryElement(int bbe) const {
      int ItemInK;
      return BoundaryElement(bbe, ItemInK);
   }

   template <int N, int M>
   SortArray<int, N> iteme(const std::vector<std::vector<int>> &nu, int k,
                           int i) {
      int nnv[N];
      Element &K(elements[CheckT(k)]);
      assert(i >= 0 && i < M);
      for (int j = 0; j < N; ++j) {
         nnv[j] = operator()(K[nu[i][j]]);
      }

      return SortArray<int, N>(nnv);
   }
   SortArray<int, B::nv> itemadj(int k, int i) {
      return iteme<B::nv, T::nea>(T::nvhyperFace, k, i);
   }
   SortArray<int, B::nv> itembe(int k) {
      int nnv[B::nv];
      B &K(borderelements[CheckBE(k)]);
      for (int j = 0; j < B::nv; ++j) {
         nnv[j] = operator()(K[j]);
      }
      return SortArray<int, B::nv>(nnv);
   }

   R mesure() const { return mes; }
   R bordermesure() const { return mesb; }
   double get_mesh_size() const {
      double hh = 1e300;
      for (int k = 0; k < nt; ++k) {
         hh = std::min((*this)[k].hElement(), hh);
      }
      return hh;
   }

   ~GenericMesh() {
      delete[] TheAdjacencesLink;
      delete[] BoundaryElementHeadLink;
      delete[] borderelements;
      if (nt > 0)
         delete[] elements;
      delete[] vertices;
   }

 private:
   GenericMesh(const GenericMesh &);    // pas de construction par copie
   void operator=(const GenericMesh &); // pas affectation par copy
};

template <typename GMesh> class BuildAdjacencyOfMesh {
 public:
   typedef SortArray<int, GMesh::nva> SArray;

   GMesh &mesh;
   std::map<SArray, int> h;
   int nk = 0, nba = 0;
   int ne = 0;

   BuildAdjacencyOfMesh(GMesh &);
   void initializeArray();
   void findAdjacencyElement();
   void addFace(const int, const int);
   void addFaceFirstTime(const SArray &);
   void addFaceAlreadySeen(typename std::map<SArray, int>::iterator,
                           const SArray &);
   void findBoundaryElement();
   void addBoundary(const int);

 public:
};

#endif
