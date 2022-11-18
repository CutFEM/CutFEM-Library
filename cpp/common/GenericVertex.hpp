#ifndef COMMON_GENERIC_VERTEX_HPP
#define COMMON_GENERIC_VERTEX_HPP

#include "Label.hpp"
#include "point.hpp"
/*
 *    Class of the generic vertex
 *
 */
template <typename Rn> class GenericVertex : public Rn, public Label {
   template <typename T, typename B, typename V> friend class GenericMesh;

   friend inline std::ostream &operator<<(std::ostream &f,
                                          const GenericVertex &v) {
      f << (const Rn &)v << ' ' << (const Label &)v;
      return f;
   }
   friend inline std::istream &operator>>(std::istream &f, GenericVertex &v) {
      f >> (Rn &)v >> (Label &)v;
      return f;
   }

 public:
   typedef Rn Rd;
   static const int d = Rd::d;

   GenericVertex() : Rd(), Label(){};                         //,normal(0) {};
   GenericVertex(const Rd &P, int r = 0) : Rd(P), Label(r){}; //,normal(0){}

 private: // pas de copie pour ne pas prendre l'adresse
   GenericVertex(const GenericVertex &);
   void operator=(const GenericVertex &);
};

#endif
