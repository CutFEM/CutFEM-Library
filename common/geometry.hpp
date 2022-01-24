#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "R3.hpp"
#include "Mesh2dn.hpp"

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



Droite equation(const R2 a,const R2 b) ;
bool p_boncote(const R2 a, const R2 b, const R2 c, const R2 p) ;
bool p_dans_triangle(const typename Mesh2::Element& K, const R2 P);
int find_triangle_contenant_p(const Mesh2& Th, const R2 P, int k_init = 0);





#endif
