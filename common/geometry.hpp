#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "R3.hpp"
#include "Mesh2dn.hpp"


namespace geometry{

  template<typename Rd>
  struct Segment{
    Rd A;
    Rd B;

    Segment(Rd a, Rd b) : A(a) , B(b) {}
    bool is_between(const Rd C) const ;
  };

  // typedef Segment<R2> Segment2;
  // equation cartesienne d'une droite
  typedef struct coefDroite{
    double pente;	//taux d'accroissement de la droite
    double ord_or;	//ordonnée à l'origine de la droite
  } Droite;

  struct Interval {
    double v_min, v_max;
    Interval(double a, double b) :v_min(a), v_max(b){}
    Interval operator*(double x) const {return Interval(x*v_min, x*v_max);}
    Interval operator+(double x) const {return Interval(x+v_min, x+v_max);}
  };

  Droite equation(const R2 a,const R2 b) ;
  bool p_boncote(const R2 a, const R2 b, const R2 c, const R2 p) ;
  bool p_dans_triangle(const typename Mesh2::Element& K, const R2 P);
  bool p_dans_triangle(const R2 P, R2 a, R2 b, R2 c);
  int find_triangle_contenant_p(const Mesh2& Th, const R2 P, int k_init = 0);


  // bool check_intersect(R2 u, R2 v, R2 w, R& t) {
  //
  //   R det = -u.x*v.y + u.y*v.x;
  //   if(fabs(det) < 1e-14) return false;
  //
  //   t = -w[0]*v[1] + w[1]*v[0];
  //   t /= det;
  //   R t2 = -w[0]*u[1] + w[1]*u[0];
  //   t2 /= det;
  //   if(t >= 0 && t <= 1 && t2 >=0 && t2 <= 1) return true;
  //   return false;
  // }


  // template argument is RdHat::d

  R3 map_point_to_simplex(const R3 N[4], const R3 P);
  R2 map_point_to_simplex(const R2 N[3], const R2 P);
  // maps to face simplex
  R3 map_point_to_simplex(const R3 N[3], const R2 P);
  R2 map_point_to_simplex(const R2 N[2], const R1 P);


  R measure_hyper_simplex(R2 N[2]);
  R measure_hyper_simplex(R3 N[3]);

  //change to template maybe
  template<int d>
  R mesure_simplex(R2 N[d+1]);
  template<int d>
  R mesure_simplex(R3 N[d+1]);
  // R mesure_simplex(R2 N[2]);
  // R mesure_simplex(R3 N[3]);
}

#endif
