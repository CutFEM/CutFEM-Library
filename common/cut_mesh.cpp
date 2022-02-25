#include "cut_mesh.hpp"
#include "Interface2dn.hpp"
#include "Interface3dn.hpp"

template<>
void Cut_Mesh<Triangle2,BoundaryEdge2,Vertex2>::add(const Rn& level_set, int sign_domain, int dom_id){

  const int nve = Th[0].nv;
  double loc_ls[nve];

  for(int k = 0; k < Th.nt; ++k) {
    if(check_exist(k, dom_id)) continue;

    for(int i=0;i<nve;++i) {
      loc_ls[i] = level_set(Th(k,i));
    }
    const SignElement<Element> signK( loc_ls);

    if(signK.sign() == sign_domain || signK.cut()) {

      // const Partition2& cut =  Partition2(Th[k], loc_ls);
      // ElementSignEnum part = cut.what_part(dom_id);
      // nt_paraview[dom_id] += cut.element_end(part) - cut.element_begin(part);

      idx_in_background_mesh[dom_id][nt] = k;
      idx_from_background_mesh[dom_id][k] = nt;
      this->nt+= 1;
    }
  }

}

template<>
void Cut_Mesh<Tet,Triangle3,Vertex3>::add(const Rn& level_set, int sign_domain, int dom_id){

  const int nve = Th[0].nv;
  double loc_ls[nve];

  for(int k = 0; k < Th.nt; ++k) {
    if(check_exist(k, dom_id)) continue;

    for(int i=0;i<nve;++i) loc_ls[i] = level_set(Th(k,i));
    const SignElement<Element> signK( loc_ls);

    if(signK.sign() == sign_domain || signK.cut()) {

      idx_in_background_mesh[dom_id][nt] = k;
      idx_from_background_mesh[dom_id][k] = nt;
      this->nt+= 1;

    }
  }
}
