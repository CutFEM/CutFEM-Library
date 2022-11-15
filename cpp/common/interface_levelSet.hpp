#ifndef LEVELSET_INTERFACE_HPP_
#define LEVELSET_INTERFACE_HPP_

#include "base_interface.hpp"
// #include "../FESpace/expression.hpp"


template<typename M>
class InterfaceLevelSet  : public Interface<M> {

  typedef M Mesh;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::Rd Rd;
  static const int nve = Rd::d;
  typedef FaceInterface<nve> Face;
  typedef SortArray<Ubyte, Element::Rd::d+1> ElementIdx;

  KN<byte> ls_sign;
  KN<double> ls_;
  // const FunFEMVirtual& fun;

public:
  // InterfaceLevelSet(const Mesh & MM, const KN<double>& lss, int label = 0)
  // : Interface<M>(MM) , ls_(lss)
  // {
  //   make_patch(ls_, label);
  // }
  // InterfaceLevelSet(const Mesh & MM, const FunFEM<Mesh>& lss, int label = 0)
  InterfaceLevelSet(const Mesh & MM, const FunFEMVirtual& lss, int label = 0)

  : Interface<M>(MM) , ls_(lss.getArray())//, fun(lss)
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
  Partition<typename Element::Face> get_partition_face(const typename Element::Face& face, int k, int ifac) const {
    double loc_ls[Element::Face::nv];
    for(int i=0;i<Element::Face::nv;++i) {
      int j = Element::nvhyperFace[ifac][i];
      int iglb = this->backMesh->at(k,j);
      loc_ls[i] = ls_[iglb];
    }
    return Partition<typename Element::Face>(face, loc_ls);
  }
  bool isCutFace(int k, int ifac) const;

  void cut_partition(Physical_Partition<Element>& local_partition, vector<ElementIdx>& new_element_idx, std::list<int>& erased_element, int sign_part)const {assert(0);} ;

private:
  void make_patch (const KN<double>& ls, int label);
  const Face  make_face(const typename RefPatch<Element>::FaceIdx& ref_tri,
			    const typename Mesh::Element& K,
			    const double lset[Element::nv],
			    int label);
  Rd make_normal(const typename Mesh::Element& K, const double lset[Element::nv]);

  // Rd get_intersection_node(int k, const Rd A, const Rd B) const {
  //   double fA = fun.eval(k, A);
  //   double fB = fun.eval(k, B);
  //   double t = -fA/(fB-fA);
  //   return (1-t) * A + t * B;
  // }
  R measure(const Face& f) const {
    Rd l[nve];
    for(int i=0;i<nve;++i) l[i] = this->vertices_[f[i]];
    return geometry::measure_hyper_simplex(l);
  };

  Rd mapToPhysicalFace(int ifac, const typename Element::RdHatBord x ) const {
    Rd N[nve];
    for(int i=0;i<nve;++i) N[i] = this->vertices_[this->faces_[ifac][i]];
    return geometry::map_point_to_simplex(N, x);
  }

};


template<typename M>
void InterfaceLevelSet<M>::make_patch(const KN<double>& ls, int label) {

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
const typename InterfaceLevelSet<M>::Face
InterfaceLevelSet<M>::make_face (const typename RefPatch<Element>::FaceIdx& ref_tri,
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
typename InterfaceLevelSet<M>::Rd
InterfaceLevelSet<M>::make_normal (const typename Mesh::Element& K, const double lset[Element::nv]) {

  Rd grad[Element::nv];
  K.Gradlambda(grad);
  Rd normal_ls;
  for(int i=0; i<Mesh::Element::nv;++i) {
    normal_ls += grad[i]*lset[i];
  }
  normal_ls /= normal_ls.norm();
  return normal_ls;
}

template<typename Mesh>
void TimeInterface<Mesh>::init(int i, const Mesh & Th, const FunFEMVirtual& ls) {
  assert(0 <= i && i < n_);
  if(interface_[i]) {
    delete interface_[i];
  }
  interface_[i] = new InterfaceLevelSet<Mesh>(Th,ls);
}
template<typename Mesh>
void TimeInterface<Mesh>::init(const Mesh & Th, const KN<FunFEMVirtual>& ls) {
  for(int i=0;i<n_;++i){
    if(interface_[i]) { delete interface_[i];}
  }
  if(n_ != ls.size()) {
    n_ = ls.size();
    interface_.resize(n_);
  }
  for(int i=0;i<n_;++i){
    interface_[i] = new InterfaceLevelSet<Mesh>(Th,ls[i]);
  }
}



#endif
