#ifndef _FIXE_LABEL_HPP
#define _FIXE_LABEL_HPP

#include "Mesh3dn.hpp"


template<typename Mesh>
void setLabelBorder(Mesh& Th, R3 normal, int newLabel) {
  typedef typename Mesh::Rd Rd;
  for( int ifac = 0; ifac < Th.nbe; ifac+=1) {
     // for( int ifac = Th.first_element(); ifac < Th.last_boundary_element(); ifac+=Th.next_element()) {
    typename Mesh::BorderElement & face(Th.be(ifac));
    int ifaceK; // index of face of triangle corresp to edge (0,1,2)
    const int k = Th.BoundaryElement(ifac, ifaceK); // index of element (triangle), ifaceK gets modified inside

    int ib = ifaceK;
    if(Th.ElementAdj(k,ib) != -1) continue; // not on the boundary. Because mesh built with buildlayers in freefem++ 

    const typename Mesh::Element& T(Th[k]);
    Rd normal_ext = T.N(ifaceK);

    Rd diff_normal(normal - normal_ext);
    if(diff_normal.norm() < 1e-15){
      face.lab = newLabel;
    }
  }
}



#endif
