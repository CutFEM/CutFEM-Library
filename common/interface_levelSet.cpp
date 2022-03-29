#include "interface_levelSet.hpp"

template<>
void InterfaceLevelSet<Mesh2>::cut_partition(Physical_Partition<Element>& local_partition, vector<ElementIdx>& new_element_idx, std::list<int>& erased_element, int sign_part)const {
  new_element_idx.resize(0);
  erased_element.resize(0);
  byte ls[3];
  const Element& T(local_partition.T);
  int kb = (*this->backMesh)(T);
  for(int k=0; k<local_partition.nb_element(0);++k){ // 0 is just useless. Not needed by this class

    // BUILD ELEMENT
    const CutElement<Element> K = local_partition.get_element(k);
    for(int i=0;i<3;++i) ls[i] = util::sign(fun.eval(kb, K[i]));

    //COMPUTE THE PARTITION
    const RefPartition<Triangle2>& patch(RefPartition<Triangle2>::instance(ls));

    // interface i CUT THIS LOCAL ELEMENT k
    if(patch.is_cut()) {
      erased_element.push_back(k);

      // LOOP OVER ELEMENTS IN PATCH
      // need to get only the part corresponding to the sign
      for(auto it = patch.element_begin(); it != patch.element_end(); ++it) {

        // if(patch.whatSign(it) != sign_part &&  patch.whatSign(it) != 0) continue;
        if(patch.whatSign(it) != sign_part ) continue;

        // create the Nodes
        // std::cout << " index node to create " << std::endl;
        int idx_begin = local_partition.nb_node();
        ElementIdx idx_array(idx_begin, idx_begin+1, idx_begin+2);
        for(int i=0; i<3;++i) {
          Uint idx = (*it)[i];
          // std::cout << idx << std::endl;
//
          if(idx < 3) {
            local_partition.add_node(K[idx]);
            // std::cout << K[idx] << std::endl;

          }
          else{
            int i0 = Triangle2::nvedge[idx - 3][0];
            int i1 = Triangle2::nvedge[idx - 3][1];
            local_partition.add_node(get_intersection_node(kb,K[i0], K[i1]));
            // std::cout << get_intersection_node(kb,K[i0], K[i1]) << std::endl;
          }
        }
        // ADD THE INDICES
        new_element_idx.push_back(idx_array);
      }

      // std::cout << " local element " << k << " is cut" << std::endl;
    }

    else {
      // std::cout << " local element " << k << " is not cut" << std::endl;
      // has to be removed if not in domain
      if(patch.whatSign() != sign_part) {
        erased_element.push_back(k);
      };

    }
  }
}



template<>
bool InterfaceLevelSet<Mesh3>::isCutFace(int k, int ifac) const {

  for(int e=0;e<Element::Face::ne;++e){
    int id_edge = Element::edgeOfFace[ifac][e];
    int i1 = Element::nvedge[id_edge][0];
    int i2 = Element::nvedge[id_edge][1];
    if(ls_[this->backMesh->at(k,i1)] * ls_[this->backMesh->at(k,i2)] < 0) return true;
  }
  return false;
}

template<>
bool InterfaceLevelSet<Mesh2>::isCutFace(int k, int ifac) const {
  int i1 = Element::nvedge[ifac][0];
  int i2 = Element::nvedge[ifac][1];
  if(ls_[this->backMesh->at(k,i1)] * ls_[this->backMesh->at(k,i2)] < 0) return true;
  return false;

}

template<>
bool InterfaceLevelSet<MeshHexa>::isCutFace(int k, int ifac) const {
  for(int e=0;e<Element::Face::ne;++e){
    int id_edge = Element::edgeOfFace[ifac][e];
    int i1 = Element::nvedge[id_edge][0];
    int i2 = Element::nvedge[id_edge][1];
    if(ls_[this->backMesh->at(k,i1)] * ls_[this->backMesh->at(k,i2)] < 0) return true;
  }
  return false;

}

template<>
bool InterfaceLevelSet<MeshQuad2>::isCutFace(int k, int ifac) const {
  int i1 = Element::nvedge[ifac][0];
  int i2 = Element::nvedge[ifac][1];
  if(ls_[this->backMesh->at(k,i1)] * ls_[this->backMesh->at(k,i2)] < 0) return true;
return false;

}
