#include "geometry.hpp"


template<> bool Segment<R2>::is_between(const R2 C) const {
  assert(0);
  return true;
}


Droite equation(const R2 a,const R2 b) {
	Droite d;
	if (!(b.x-a.x)) d.pente=b.x;
	else d.pente=(b.y-a.y) / (b.x-a.x);
	d.ord_or=a.y - a.x * d.pente;
	return d;
}

bool p_boncote(const R2 a, const R2 b, const R2 c, const R2 p) {
  //teste si la droite est perpendiculaire à l'axe des abcisses
	//et retourne le resultat en conséquence

	if(a.x==b.x){
    return ((p.x-a.x) * (c.x-a.x) >= 0);
  }

  Droite d;
	d=equation(a,b);
	return ((d.pente*c.x+d.ord_or-c.y)*(d.pente*p.x+d.ord_or-p.y) >= 0);
}

bool p_dans_triangle(const typename Mesh2::Element& K, const R2 P){
  return ( (p_boncote( K[0], K[1], K[2], P))
        && (p_boncote( K[2], K[1], K[0], P))
        && (p_boncote( K[0], K[2], K[1], P)) );
}

int find_triangle_contenant_p(const Mesh2& Th, const R2 P, int k_init){

  int idx_elt = k_init;
  int count = 0;

  // std::ofstream plot;
  // plot.open("searchElement.dat", std::ofstream::out);

  while (!p_dans_triangle(Th[idx_elt],P)) {


    const typename Mesh2::Element& K(Th[idx_elt]);

    // for(int i=0;i<3;++i) {
    //   plot << K[i] << std::endl;
    // }
    // plot << K[0] << std::endl;
    // plot << std::endl;
    // plot << std::endl;
    //passe au triangle voisin qui est dans la "direction" de p
    //première arête
    int iface;
    if (! p_boncote(K[0],K[1],K[2],P)) {
      iface = 2;
    }
    else {
      if (! p_boncote(K[1],K[2],K[0],P)) {
        iface = 0;
      }
      else{
        iface = 1;
      }
    }
    idx_elt = Th.ElementAdj(idx_elt, iface);
  }
  const typename Mesh2::Element& K(Th[idx_elt]);

  // for(int i=0;i<3;++i) {
  //   plot << K[i] << std::endl;
  // }
  // plot << K[0] << std::endl;
  // plot << std::endl;
  // plot << std::endl;
  // plot.close();

  return idx_elt;
}
