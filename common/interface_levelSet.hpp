#include "base_interface.hpp"



template<typename M>
class Interface_LevelSet  : public Interface<M> {

  typedef M Mesh;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::Rd Rd;
  static const int nve = Rd::d;
  typedef FaceInterface<nve> Face;

  KN<byte> ls_sign;

public:
  Interface_LevelSet(const Mesh & MM, const KN<double>& ls, int label = 0)
  : Interface<M>(MM)
  {
    make_patch(ls, label);
  }


private:
  void make_patch (const KN<double>& ls, int label);
  const Face  make_face(const typename RefPatch<Element>::FaceIdx& ref_tri,
			    const typename Mesh::Element& K,
			    const double lset[Element::nv],
			    int label);
  Rd make_normal(const typename Mesh::Element& K, const double lset[Element::nv]);

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
