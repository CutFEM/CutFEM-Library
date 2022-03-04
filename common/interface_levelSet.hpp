#ifndef LEVELSET_INTERFACE_HPP_
#define LEVELSET_INTERFACE_HPP_

#include "base_interface.hpp"
// #include "../FESpace/expression.hpp"


template<typename M>
class Interface_LevelSet  : public Interface<M> {

  typedef M Mesh;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::Rd Rd;
  static const int nve = Rd::d;
  typedef FaceInterface<nve> Face;
  typedef SortArray<int, Element::Rd::d+1> ElementIdx;

  KN<byte> ls_sign;
  KN<double> ls_;
  // const FunFEM<Mesh>& fun;
  const FunFEMVirtual& fun;

public:
  // Interface_LevelSet(const Mesh & MM, const KN<double>& lss, int label = 0)
  // : Interface<M>(MM) , ls_(lss)
  // {
  //   make_patch(ls_, label);
  // }
  // Interface_LevelSet(const Mesh & MM, const FunFEM<Mesh>& lss, int label = 0)
  Interface_LevelSet(const Mesh & MM, const FunFEMVirtual& lss, int label = 0)

  : Interface<M>(MM) , ls_(lss.v), fun(lss)
  {
    make_patch(ls_, label);
  }

  SignElement<Element> get_SignElement(int k) const {
    byte loc_ls[Element::nv];
    for(int i=0;i<Element::nv;++i) {
      int iglb = this->backMesh->at(k,i);
      loc_ls[i] = ls_sign[iglb];
    }
    return SignElement<Element>(loc_ls);
  }

  Partition<Element> get_partition(int k) const {
    double loc_ls[Element::nv];
    for(int i=0;i<Element::nv;++i) {
      int iglb = this->backMesh->at(k,i);
      loc_ls[i] = ls_[iglb];
    }

    return Partition<Element>((*this->backMesh)[k], loc_ls);
  }

  void cut_partition(Physical_Partition<Element>& local_partition, vector<ElementIdx>& new_element_idx, std::list<int>& erased_element, int sign_part)const ;

private:
  void make_patch (const KN<double>& ls, int label);
  const Face  make_face(const typename RefPatch<Element>::FaceIdx& ref_tri,
			    const typename Mesh::Element& K,
			    const double lset[Element::nv],
			    int label);
  Rd make_normal(const typename Mesh::Element& K, const double lset[Element::nv]);

  Rd get_intersection_node(int k, const Rd A, const Rd B) const {
    double fA = fun.eval(k, A);
    double fB = fun.eval(k, B);
    double t = -fA/(fB-fA);
    return (1-t) * A + t * B;
  }

};


template<typename M>
void Interface_LevelSet<M>::make_patch(const KN<double>& ls, int label) {

  assert(this->backMesh);
  this->faces_.resize( 0);                          // reinitialize arrays
  this->vertices_.resize(0);
  this->element_of_face_.resize(0);
  this->outward_normal_.resize(0);
  this->face_of_element_.clear();

  const Mesh& Th = *(this->backMesh) ;
  util::copy_levelset_sign( ls, ls_sign);

  const Uint nb_vertex_K = Element::nv;
  double loc_ls[nb_vertex_K];
  byte   loc_ls_sign[nb_vertex_K];

  for (int k=0; k<this->backMesh->nbElmts(); k++) {                      // loop over elements

    const typename Mesh::Element & K(Th[k]);

    for (Uint i= 0; i < K.nv; ++i) {
      loc_ls_sign[i] = ls_sign[Th(K[i])];
      loc_ls     [i] = ls     [Th(K[i])];
      // std::cout << loc_ls[i] << std::endl;
    }

    const RefPatch<Element>& cut =  RefPatch<Element>::instance( loc_ls_sign);
    if (cut.empty()) continue;
    // int iii=0;
    for (typename RefPatch<Element>::const_face_iterator it= cut.face_begin(), end= cut.face_end();
    it != end; ++it) {
      // std::cout << " face numero " << iii++ << std::endl;
      this->face_of_element_[k] = this->element_of_face_.size();
      this->faces_.push_back( make_face(*it, K, loc_ls, label));
      this->element_of_face_.push_back(k);
      this->outward_normal_.push_back(make_normal(K, loc_ls));
    }
  }
}

template<typename M>
const typename Interface_LevelSet<M>::Face
Interface_LevelSet<M>::make_face (const typename RefPatch<Element>::FaceIdx& ref_tri,
				  const typename Mesh::Element& K,
				  const double lset[Element::nv],
				  int label)
{

  Uint loc_vert_num;
  Uint triIdx[nve];


  for (Uint j= 0; j < nve; ++j) {
    loc_vert_num = ref_tri[j];

    if (loc_vert_num < K.nv) {                            // zero vertex

      const Uint idx = (*this->backMesh)(K[loc_vert_num]);
      Rd Q = (*this->backMesh)(idx);
      this->vertices_.push_back(Q);
      triIdx[j] = this->vertices_.size() - 1;
      assert(0);
    }
    else { // genuine edge vertex

      const Ubyte i0 = Mesh::Element::nvedge[loc_vert_num - K.nv][0],
      	i1 = Mesh::Element::nvedge[loc_vert_num - K.nv][1];

      const double t = lset[i0]/(lset[i0] - lset[i1]);
      Rd Q = (1.0 - t) * ((Rd) K[i0]) + t * ((Rd) K[i1]); // linear interpolation
      this->vertices_.push_back(Q);
      triIdx[j] = this->vertices_.size() - 1;

      this->edge_of_node_.push_back(loc_vert_num - K.nv);
    }
  }
  return Face(triIdx, label);
}



template<typename M>
typename Interface_LevelSet<M>::Rd
Interface_LevelSet<M>::make_normal (const typename Mesh::Element& K, const double lset[Element::nv]) {

  Rd grad[Element::nv];
  // K.Gradlambda(grad);
  Rd normal_ls;
  // for(int i=0; i<Element::nv;++i) {
  //   normal_ls += grad[i]*lset[i];
  // }
  // normal_ls /= normal_ls.norm();
  return normal_ls;
}

template<typename M>
void Interface_LevelSet<M>::cut_partition(Physical_Partition<Element>& local_partition, vector<ElementIdx>& new_element_idx, std::list<int>& erased_element, int sign_part)const {
  new_element_idx.resize(0);
  erased_element.resize(0);
  byte ls[3];
  const Element& T(local_partition.T);
  int kb = (*this->backMesh)(T);
  for(int k=0; k<local_partition.nb_element();++k){

    // BUILD ELEMENT
    const CutElement2 K = local_partition.get_element(k);
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




#endif
