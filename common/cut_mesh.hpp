#ifndef _CUT_MESH_HPP
#define _CUT_MESH_HPP

#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "base_interface.hpp"


struct CBorder {
  CBorder() {}
};
const  CBorder INTEGRAL_BOUNDARY;
struct CFacet {
  CFacet() {}
};
const CFacet INTEGRAL_INNER_FACET;
const CFacet INTEGRAL_INNER_EDGE_2D;
const CFacet INTEGRAL_INNER_FACE_3D;

struct CRidge {
  CRidge() {}
};
const CRidge INTEGRAL_INNER_RIDGE;
const CRidge INTEGRAL_INNER_NODE_2D;
const CRidge INTEGRAL_INNER_EDGE_3D;

struct CExtension {
  CExtension() {}
};
const CExtension INTEGRAL_EXTENSION;

// class enum{Explicit_Partition}


template<typename E >
struct Cut_Part {
  static const int dim = E::RdHat::d;
  typedef SortArray<Ubyte, dim+1> ElementIdx;           // the vertices of a triangle of the cut:
  typedef const ElementIdx* const_element_iterator;
  typedef       ElementIdx*       element_iterator;
  typedef typename E::Rd Rd;
  typedef typename E::RdHat RdHat;


  const Virtual_Partition<E>* partition_;
  const Partition<E> ip;
  const Physical_Partition<E> pp;
  int sign_cut_;

  Cut_Part(const Partition<E> p, int s) : ip(p), sign_cut_(s),pp(p.T) {
    partition_ = &ip;
  }
  Cut_Part(const Physical_Partition<E> p, int s) : pp(p), sign_cut_(s),ip(p.T) {
    partition_ = &pp;
  }
  Cut_Part(const Cut_Part<E>& p) :pp(p.pp), ip(p.ip) {
    if(p.partition_ == &p.pp) partition_ = &pp;
    else partition_ = &ip;
  }

  // GETTERS
  int get_sign() const {return sign_cut_;}
  int get_sign_node(int i) const {return partition_->get_sign_node(i);}
  void get_list_node (vector<typename E::Rd>& node) const{ partition_->get_list_node(node, sign_cut_); }
  CutElement<E> get_element(int k) const {return partition_->get_element(k);}
  Rd get_vertex(const_element_iterator it, const int i) const {return partition_->get_vertex(it,i);}
  int get_nb_element() const {return partition_->nb_element(sign_cut_);}
  int get_local_domain_id() const {
    if(sign_cut_==0) return -1;  // not cout
    else return (sign_cut_==-1);
  }

  // OTHER METHODS
  // GIVE THE MEASURE OF THE CUT PART IN Rd
  double measure() const {return partition_->measure(sign_cut_);}
  double measure(const_element_iterator it) const {return partition_->measure(it);}

  // //GIVE THE MEASURE OF CUT PART OF A FACE IN RdBord
  // double measureBord(int ifac) const {return partition_->measureBord(sign_cut_, ifac);}

  bool multi_interface() const {return partition_ == &pp;}

  Rd mapToPhysicalElement(const_element_iterator it, const RdHat Phat) const {
    return partition_->mapToPhysicalElement(it, Phat);
  }

  // ITERATORS
  const_element_iterator element_begin () const{
    return partition_->element_begin(sign_cut_);
  }
  const_element_iterator element_end () const{
    return partition_->element_end(sign_cut_);
  }
  const_element_iterator other_side_element_begin () const{
    assert(sign_cut_ != 0);
    int other_side_cut = (sign_cut_==1);
    return partition_->element_begin(other_side_cut);
  }
  const_element_iterator other_side_element_end () const{
    assert(sign_cut_ != 0);
    int other_side_cut = (sign_cut_==1);
    return partition_->element_end(other_side_cut);
  }


};


template<typename Mesh>
class ActiveMesh {

public:
  typedef typename Mesh::Element Element;
  typedef typename Element::Face Face;
  typedef typename Mesh::Rd Rd;
  typedef typename Mesh::BorderElement BorderElement;
  typedef typename Element::RdHat RdHat;// for parametrization
  typedef SortArray<int, 2> pairIndex;
  static const int nea=Element::nea; //  numbering of adj (4 in Tet,  3 in Tria, 2 in seg)
  static const int nva=Element::nva; //  numbering of vertex in Adj hyperface

  const Mesh& Th;

  std::vector<std::vector<int>> idx_in_background_mesh_;   // [domain][idxK_cutMesh] -> idxK_backMesh
  std::vector<std::map<int,int>>  idx_from_background_mesh_; // [domain](idxK_backMesh) -> idxK_cutMesh
  std::vector<std::map<std::pair<int,int>, std::vector<std::pair<const Interface<Mesh>*, int>>>> interface_id_; // [time_quad](domain_id, idx_k) -> [n_interface](interface, sign)
  std::vector<int> idx_element_domain;


  // For time problem
  int nb_quadrature_time_;
  // map des elements that are not always in the active mesh
  std::vector<std::vector<std::map<int, bool>>> in_active_mesh_; // [dom][itq][idx_element] -> true/false

public:
  // Create a CutMesh without cut on the backMesh
  // Usefull if wanna add sub domains
  ActiveMesh(const Mesh& th) : Th(th) {
    idx_in_background_mesh_.reserve(10);
    idx_from_background_mesh_.reserve(10);
    idx_in_background_mesh_.resize(1);
    idx_from_background_mesh_.resize(1);
    idx_in_background_mesh_[0].resize(Th.nt);
    interface_id_.resize(5);
    nb_quadrature_time_ = 1;
    for(int k = 0; k < Th.nt; ++k) {
      idx_in_background_mesh_[0][k] = k;
      idx_from_background_mesh_[0][k] = k;
    }
    idx_element_domain.push_back(0);
    idx_element_domain.push_back(Th.nt);
    in_active_mesh_.resize(10);
    for(int i=0;i<10;++i) in_active_mesh_[i].resize(nb_quadrature_time_);
  }
  // Give the background mesh and a sign Function defined on the mesh nodes
  // Will create 2 subdomains
  ActiveMesh(const Mesh& th, const Interface<Mesh>& interface) :  Th(th) {
    idx_in_background_mesh_.reserve(10);
    idx_from_background_mesh_.reserve(10);
    idx_in_background_mesh_.resize(2);
    idx_from_background_mesh_.resize(2);
    interface_id_.resize(1);
    nb_quadrature_time_ = 1;
    this->init(interface);
  }

  ActiveMesh(const Mesh& th, const TimeInterface<Mesh>& interface) :  Th(th) {
    idx_in_background_mesh_.reserve(10);
    idx_from_background_mesh_.reserve(10);
    idx_in_background_mesh_.resize(2);
    idx_from_background_mesh_.resize(2);
    nb_quadrature_time_ = interface.size();
    interface_id_.resize(nb_quadrature_time_);
    in_active_mesh_.resize(10);   // Added to avoid stack overflow error
    for(int i=0;i<10;++i) in_active_mesh_[i].resize(nb_quadrature_time_);
    this->init(interface);
  }

  void truncate(const Interface<Mesh>& interface, int sign_domain);
  void truncate(const TimeInterface<Mesh>& interface, int sign_domain);
  void add(const Interface<Mesh>& interface, int sign_domain);
  void createSurfaceMesh(const Interface<Mesh>& interface);
  void createSurfaceMesh(const TimeInterface<Mesh>& interface);
private:
  void init(const Interface<Mesh>& interface);
  void init(const TimeInterface<Mesh>& interface);
  bool check_exist(int k, int dom) const {
    const auto it = idx_from_background_mesh_[dom].find(k);
    if(it == idx_from_background_mesh_[dom].end()) return false;
    else return true;
  }
  int idxK_begin(int i) const {return this->idx_element_domain[i];}
  int idxK_in_domain(int k, int i) const {
    return k - this->idx_element_domain[i];
  }
  Physical_Partition<Element> build_local_partition(const int k, int t=0) const ;
  Physical_Partition<Face> build_local_partition(Face& face, const int k, int ifac, int t=0) const ;
public:
  DataFENodeDF BuildDFNumbering(int ndfv,int ndfe,int ndff,int ndft, int nndv,int nnde,int nndf,int nndt, int N=1, const PeriodicBC* PPeriod = nullptr) const {
        assert(0);
        return Th.BuildDFNumbering(ndfv,ndfe,ndff,ndft,nndv,nnde,nndf,nndt,N,PPeriod);
      }

  int nbElmts() const {return get_nb_element();}
  int NbElement() const {return get_nb_element();}

  int nbBrdElmts() const {return Th.nbe;}
  int nbVertices() const {return Th.nv;}

  const Element & operator[](int i) const {
    int k = idxElementInBackMesh(i);
    return Th[k];
  }
  const BorderElement& be(int i) const {return Th.be(i);}

  int get_nb_domain() const {return idx_in_background_mesh_.size();}
  int get_nb_element(int i) const {
    assert(i>=0 && i<this->get_nb_domain());
    return idx_in_background_mesh_[i].size();
  }
  int get_nb_element() const {
    int s = 0;
    for(int i=0;i<this->get_nb_domain();++i) {
      s+= get_nb_element(i);
    }
    return s;
  }
  int get_domain_element(const int k) const {
    for(int i=0;i<this->get_nb_domain();++i) {
      if(k<idx_element_domain[i+1]) {return i;}
    }
    assert(0);
  }


  bool isCut(int k, int t) const {
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_[t].find(std::make_pair(domain, kloc));
    return (it != interface_id_[t].end());
  }
  bool isCutFace(int k, int ifac, int t) const {
    if(!isCut(k,t)) return false;
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_[t].find(std::make_pair(domain, kloc));
    assert(it != interface_id_[t].end());
    int s = it->second.at(0).second;  // 0 because no multi cut
    if(s == 0) return false; // means surface mesh
    int kb = this->idxElementInBackMesh(k);
    return it->second.at(0).first->isCutFace(kb,ifac);
  }
  bool isStabilizeElement(int k) const {
    // is cut or not always active
    for(int i=0;i<nb_quadrature_time_;++i){
      if(this->isCut(k,i) || this->isInactive(k,i)) return true;
    }
    return false;
  }
  bool isInactive(int k, int t) const {
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = in_active_mesh_[domain][t].find(kloc);
    if(it == in_active_mesh_[domain][t].end()) return false;
    return true;
  }
  const Interface<Mesh>& get_interface(int k, int t) const {
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_[t].find(std::make_pair(domain, kloc));
    assert(it != interface_id_[t].end());
    return *(it->second.at(0).first);
  }
  Partition<Element> get_partition(int k, int t) const {
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_[t].find(std::make_pair(domain, kloc));
    assert(it != interface_id_[t].end());
    int kb = this->idxElementInBackMesh(k);
    return it->second.at(0).first->get_partition(kb);
  }
  int get_sign_cut(int k, int t) const {
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_[t].find(std::make_pair(domain, kloc));
    assert(it != interface_id_[t].end());
    return it->second.at(0).second;
  }
  Cut_Part<Element> get_cut_part(int k, int t) const {
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_[t].find(std::make_pair(domain, kloc));
    // if not cut build a partition that consider full element
    if(it == interface_id_[t].end()) {
      return Cut_Part<Element>(Partition<Element>((*this)[k]), -1);
    }
    int kb = this->idxElementInBackMesh(k);
    if(it->second.size() == 1)
    return Cut_Part<Element>(it->second.at(0).first->get_partition(kb), it->second.at(0).second);
    else
    return Cut_Part<Element>(this->build_local_partition(k), 0);
  }
  Cut_Part<typename Element::Face> get_cut_face(Face& face, int k, int ifac, int t) const {

    // BUILD THE FACE
    // In the class mesh the inner faces are not built
    int kb = this->idxElementInBackMesh(k);
    int iv[Face::nv];
    for(int i=0;i<Face::nv;++i) iv[i] = Th(kb, Element::nvhyperFace[ifac][i]);
    face.set(Th.vertices, iv, 0);

    // GET THE INTERFACE
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_[t].find(std::make_pair(domain, kloc));
    assert(it != interface_id_[t].end());

    if(it->second.size() == 1)
    return Cut_Part<Face>(it->second.at(0).first->get_partition_face(face,kb,ifac), it->second.at(0).second);
    else
    return Cut_Part<Face>(this->build_local_partition(face, k, ifac), 0);
  }

  vector<int> idxAllElementFromBackMesh(int k, int d) const {
    std::vector<int> idx(0);
    if(d!=-1) {
      int ret = idxElementFromBackMesh(k,  d);
      assert(ret != -1);idx.push_back(ret);
      return idx;
      }
    for(int i=0;i<get_nb_domain();++i) {
      int ret = idxElementFromBackMesh(k,  i);
      if(ret != -1) idx.push_back(ret);
    }
    assert(idx.size()>0 && idx.size() < 3);
    return idx;
  }

  vector<int> getAllDomainId(int k) const {
    std::vector<int> idx(0);
    for(int i=0;i<get_nb_domain();++i) {
      int ret = idxElementFromBackMesh(k,  i);
      if(ret != -1) idx.push_back(i);
    }
    assert(idx.size()>0 && idx.size() < 3);
    return idx;
  }

  int idxElementFromBackMesh(int k) const {
    assert(0);
    return -1;
  }
  int idxElementFromBackMesh(int k, int i) const {
    if(i==-1) assert(0);
    auto it = idx_from_background_mesh_[i].find(k);
    if(it ==  idx_from_background_mesh_[i].end()) return -1;
    return idxK_begin(i) + it->second;
  }

  int idxElementInBackMesh(const int k) const {
    int i = this->get_domain_element(k);
    int l = idxK_in_domain(k,i);
    return idx_in_background_mesh_[i][l];
  }
  int idxElementInBackMesh(const int k,int i) const {
    int l = idxK_in_domain(k,i);
    return idx_in_background_mesh_[i][l];
  }
  int ElementAdj(const int k,int &j) const {
    int domain = get_domain_element(k);
    int kb  = this->idxElementInBackMesh(k);
    int kbn = this->Th.ElementAdj(kb,j);
    if(kbn == -1) return -1;

    return this->idxElementFromBackMesh(kbn, domain);
  }

  void info() const {
    std::cout << " ------------------------------- " << std::endl;
    std::cout << " Cut Mesh has  \t" << get_nb_domain() << " domains" << std::endl;
    for(int i=0;i<get_nb_domain();++i) {
      std::cout << " nb elements in \t" << i << " => " << this->get_nb_element(i) << std::endl;
    }
    std::cout << " nb elements in total => \t" << this->get_nb_element() << std::endl;
  }

  #ifdef USE_MPI
  virtual int first_element() const { return MPIcf::first_element(this->get_nb_element());}
  virtual int next_element() const {  return MPIcf::next_element(this->get_nb_element());}
  virtual int last_element() const {  return MPIcf::last_element(this->get_nb_element());}

  virtual int first_boundary_element() const { return MPIcf::my_rank();}
  virtual int next_boundary_element() const { return MPIcf::size();}
  virtual int last_boundary_element() const {return this->Th.nbBrdElmts();}
  #else
  virtual int first_element() const { return 0;}
  virtual int next_element() const {return 1;}
  virtual int last_element() const { return this->get_nb_element();}

  virtual int first_boundary_element() const { return 0;}
  virtual int next_boundary_element() const { return 1;}
  virtual int last_boundary_element() const {return this->Th.nbBrdElmts();}
  #endif



};



//  constructor for basic 2 subdomains problem {1, -1}
template<typename Mesh>
void ActiveMesh<Mesh>::init(const Interface<Mesh>& interface){

  idx_in_background_mesh_[0].reserve(Th.nt);
  idx_in_background_mesh_[1].reserve(Th.nt);
  idx_element_domain.push_back(0);
  int nt0=0, nt1=0;
  for(int k = 0; k < Th.nt; ++k) {

    const SignElement<Element> signK = interface.get_SignElement(k);

    if(signK.cut()) {
      idx_in_background_mesh_[0].push_back(k);
      idx_from_background_mesh_[0][k] = nt0;
      idx_in_background_mesh_[1].push_back(k);
      idx_from_background_mesh_[1][k] = nt1;
      interface_id_[0][std::make_pair(0,nt0)].push_back(std::make_pair(&interface,  1));
      interface_id_[0][std::make_pair(1,nt1)].push_back(std::make_pair(&interface, -1));

      nt0++; nt1++;

    }
    else {
      int s = signK.sign();
      int& nnt = (s > 0)? nt0 : nt1;
      idx_in_background_mesh_[(s<0)].push_back(k);
      idx_from_background_mesh_[(s<0)][k] = nnt;

      nnt++;
    }
  }
  idx_in_background_mesh_[0].resize(nt0);
  idx_in_background_mesh_[1].resize(nt1);
  idx_in_background_mesh_[0].shrink_to_fit();
  idx_in_background_mesh_[1].shrink_to_fit();

  idx_element_domain.push_back(nt0);
  idx_element_domain.push_back(nt0+nt1);
  in_active_mesh_.resize(10);
  for(int i=0;i<10;++i) in_active_mesh_[i].resize(nb_quadrature_time_);

}

//  constructor a subdopmain corresponding to the positive sign
template<typename Mesh>
void ActiveMesh<Mesh>::add(const Interface<Mesh>& interface, int sign_domain){

  int dom_size = this->get_nb_domain();
  idx_element_domain.resize(0);
  // int sign_domain = -1;
  // Initialize the first new subdomain domain
  // and clear old array with element indices.
  {
    idx_in_background_mesh_.resize(dom_size+1);
    idx_from_background_mesh_.resize(dom_size+1);
    for( int d=0; d<dom_size+1;++d){
      idx_in_background_mesh_[d].resize(0);
    }
    int nt_max = idx_from_background_mesh_[0].size();
    idx_in_background_mesh_[dom_size].reserve(nt_max);
  }
  std::vector<int> nt(2*dom_size, 0.);
  int new_dom_id = dom_size;
  for( int d=0; d<dom_size;++d){
    bool sub_is_cut = false;
    for(auto it_k = idx_from_background_mesh_[d].begin(); it_k != idx_from_background_mesh_[d].end() ;) {

      int kb  = it_k->first;
      int k   = it_k->second;

      auto it_gamma = interface_id_[0].find(std::make_pair(d, k));
      const SignElement<Element> signK = interface.get_SignElement(kb);

      int nb_interface = (it_gamma == interface_id_[0].end())? 0 :it_gamma->second.size();
      std::vector<const Interface<Mesh>*> old_interface(nb_interface);
      std::vector<int> ss(nb_interface);
      for(int i=0;i<nb_interface;++i) old_interface[i] = it_gamma->second[i].first;
      for(int i=0;i<nb_interface;++i) ss[i] = it_gamma->second[i].second;

      if(it_gamma != interface_id_[0].end()) {
        auto ittt =   interface_id_[0].erase(it_gamma);
      }

      if(signK.sign() == sign_domain || signK.cut()) {

        // Initialize first time we find a cut element
        if(!sub_is_cut && new_dom_id !=dom_size){
           new_dom_id++;
           idx_in_background_mesh_.resize(new_dom_id+1);
           idx_from_background_mesh_.resize(new_dom_id+1);
           idx_in_background_mesh_[new_dom_id].resize(0);
           int nt_max = idx_from_background_mesh_[d+1].size();
           idx_in_background_mesh_[new_dom_id].reserve(nt_max);
         }


        sub_is_cut = true;
        idx_in_background_mesh_[new_dom_id].push_back(kb);
        idx_from_background_mesh_[new_dom_id][kb] = nt[new_dom_id];

        for(int i=0; i<nb_interface;++i) {
          interface_id_[0][std::make_pair(new_dom_id,nt[new_dom_id])].push_back(std::make_pair(old_interface[i],  ss[i]));
        }


        if(!signK.cut()) {
          it_k = idx_from_background_mesh_[d].erase(it_k);
        }
        else {
          idx_in_background_mesh_[d].push_back(kb);
          it_k->second = nt[d];

          // need to change k and add new interface
          // attach all interface to new element
          for(int i=0; i<nb_interface;++i) {
            interface_id_[0][std::make_pair(d,nt[d])].push_back(std::make_pair(old_interface[i],  ss[i]));
          }
          interface_id_[0][std::make_pair(new_dom_id,nt[new_dom_id])].push_back(std::make_pair(&interface,  sign_domain));
          interface_id_[0][std::make_pair(d,nt[d])].push_back(std::make_pair(&interface,  -sign_domain));
          nt[d]++;
          it_k++;
        }

        nt[new_dom_id]++;
      }
      else {
        // std::cout << " in old domain " << std::endl;
        idx_in_background_mesh_[d].push_back(kb);
        it_k->second = nt[d];

        // change the key of the interface_id map
        for(int i=0;i<nb_interface;++i){
          interface_id_[0][std::make_pair(d,nt[d])].push_back(std::make_pair(old_interface[i],  ss[i]));
        }
        it_k++;
        nt[d]++;
      }
    }

    // if(sub_is_cut && d+1 !=dom_size){
    //    new_dom_id++;
    //    idx_in_background_mesh_.resize(new_dom_id+1);
    //    idx_from_background_mesh_.resize(new_dom_id+1);
    //    idx_in_background_mesh_[new_dom_id].resize(0);
    //    int nt_max = idx_from_background_mesh_[d+1].size();
    //    idx_in_background_mesh_[new_dom_id].reserve(nt_max);
    //  }
  }


  idx_element_domain.push_back(0);
  for( int d=0; d<new_dom_id+1;++d){
    idx_in_background_mesh_[d].resize(nt[d]);
    idx_in_background_mesh_[d].shrink_to_fit();
    int sum_nt = idx_element_domain[d] + nt[d];
    idx_element_domain.push_back(sum_nt);
  }

}

template<typename Mesh>
void ActiveMesh<Mesh>::truncate(const Interface<Mesh>& interface,int sign_domain_remove){

  int dom_size = this->get_nb_domain();
  idx_element_domain.resize(0);

  {
    for( int d=0; d<dom_size;++d){
      idx_in_background_mesh_[d].resize(0);
      int nt_max = idx_from_background_mesh_[d].size();
      idx_in_background_mesh_[d].reserve(nt_max);
    }
  }

  std::vector<int> nt(dom_size, 0.);
  for( int d=0; d<dom_size;++d){
    for(auto it_k = idx_from_background_mesh_[d].begin(); it_k != idx_from_background_mesh_[d].end() ;) {

      int kb  = it_k->first;
      int k   = it_k->second;

      // std::cout << "domain \t" << d << " element back " << kb << "\t => loc id " << k << std::endl;
      auto it_gamma = interface_id_[0].find(std::make_pair(d, k));
      const SignElement<Element> signK = interface.get_SignElement(kb);

      // REMOVE THE ELEMENT IN THE INPUT DOMAIN
      if(signK.sign() == sign_domain_remove) {

        it_k = idx_from_background_mesh_[d].erase(it_k);

        continue;
      }


      // SAVE AND ERASE OLD INTERFACES
      int nb_interface = (it_gamma == interface_id_[0].end())? 0 :it_gamma->second.size();
      std::vector<const Interface<Mesh>*> old_interface(nb_interface);
      std::vector<int> ss(nb_interface);
      for(int i=0;i<nb_interface;++i) old_interface[i] = it_gamma->second[i].first;
      for(int i=0;i<nb_interface;++i) ss[i] = it_gamma->second[i].second;
      if(it_gamma != interface_id_[0].end()) {
        auto ittt =   interface_id_[0].erase(it_gamma);
      }

      // SET NEW INDICES AND PUT BACK INTERFACES
      idx_in_background_mesh_[d].push_back(kb);
      it_k->second = nt[d];
      for(int i=0; i<nb_interface;++i) {
        interface_id_[0][std::make_pair(d,nt[d])].push_back(std::make_pair(old_interface[i],  ss[i]));
      }
      // IS CUT SO NEED TO ADD INTERFACE AND SIGN
      if (signK.cut()){
        interface_id_[0][std::make_pair(d,nt[d])].push_back(std::make_pair(&interface,  -sign_domain_remove));
      }
      nt[d]++;
      it_k++;
    }
  }

  idx_element_domain.push_back(0);
  for( int d=0; d<dom_size;++d){
    idx_in_background_mesh_[d].resize(nt[d]);
    idx_in_background_mesh_[d].shrink_to_fit();
    int sum_nt = idx_element_domain[d] + nt[d];
    idx_element_domain.push_back(sum_nt);
  }


}


template<typename Mesh>
void ActiveMesh<Mesh>::createSurfaceMesh(const Interface<Mesh>& interface){

  int dom_size = this->get_nb_domain();
  idx_element_domain.resize(0);

  {
    for( int d=0; d<dom_size;++d){
      idx_in_background_mesh_[d].resize(0);
      int nt_max = idx_from_background_mesh_[d].size();
      idx_in_background_mesh_[d].reserve(nt_max);
    }
  }

  std::vector<int> nt(dom_size, 0.);
  for( int d=0; d<dom_size;++d){
    for(auto it_k = idx_from_background_mesh_[d].begin(); it_k != idx_from_background_mesh_[d].end() ;) {

      int kb  = it_k->first;
      int k   = it_k->second;

      // std::cout << "domain \t" << d << " element back " << kb << "\t => loc id " << k << std::endl;
      auto it_gamma = interface_id_[0].find(std::make_pair(d, k));
      const SignElement<Element> signK = interface.get_SignElement(kb);

      // REMOVE THE ELEMENT IN THE INPUT DOMAIN
      if(!signK.cut()) {

        it_k = idx_from_background_mesh_[d].erase(it_k);
        continue;
      }


      // SAVE AND ERASE OLD INTERFACES
      int nb_interface = (it_gamma == interface_id_[0].end())? 0 :it_gamma->second.size();
      std::vector<const Interface<Mesh>*> old_interface(nb_interface);
      std::vector<int> ss(nb_interface);
      for(int i=0;i<nb_interface;++i) old_interface[i] = it_gamma->second[i].first;
      for(int i=0;i<nb_interface;++i) ss[i] = it_gamma->second[i].second;
      if(it_gamma != interface_id_[0].end()) {
        auto ittt =   interface_id_[0].erase(it_gamma);
      }

      // SET NEW INDICES AND PUT BACK INTERFACES
      idx_in_background_mesh_[d].push_back(kb);
      it_k->second = nt[d];
      for(int i=0; i<nb_interface;++i) {
        interface_id_[0][std::make_pair(d,nt[d])].push_back(std::make_pair(old_interface[i],  ss[i]));
      }
      // IS CUT SO NEED TO ADD INTERFACE AND SIGN
      interface_id_[0][std::make_pair(d,nt[d])].push_back(std::make_pair(&interface,  0));
      nt[d]++;
      it_k++;
    }
  }

  idx_element_domain.push_back(0);
  for( int d=0; d<dom_size;++d){
    idx_in_background_mesh_[d].resize(nt[d]);
    idx_in_background_mesh_[d].shrink_to_fit();
    int sum_nt = idx_element_domain[d] + nt[d];
    idx_element_domain.push_back(sum_nt);
  }

}

template<typename Mesh>
void ActiveMesh<Mesh>::createSurfaceMesh(const TimeInterface<Mesh>& interface){

  int n_tid = interface.size();
  nb_quadrature_time_ = n_tid;
  in_active_mesh_.resize(10);
  for(int i=0;i<10;++i) in_active_mesh_[i].resize(nb_quadrature_time_);

  int dom_size = this->get_nb_domain();
  idx_element_domain.resize(0);

  {
    for( int d=0; d<dom_size;++d){
      idx_in_background_mesh_[d].resize(0);
      int nt_max = idx_from_background_mesh_[d].size();
      idx_in_background_mesh_[d].reserve(nt_max);
    }
  }

  std::vector<int> nt(dom_size, 0.);
  for( int d=0; d<dom_size;++d){
    for(auto it_k = idx_from_background_mesh_[d].begin(); it_k != idx_from_background_mesh_[d].end() ;) {

      int kb  = it_k->first;
      int k   = it_k->second;


      bool active_element = false;
      for(int t=0;t<interface.size()-1;++t){
        const SignElement<Element> signKi  = interface(t)->get_SignElement(kb);
        const SignElement<Element> signKii = interface(t+1)->get_SignElement(kb);

        if(signKi.cut() || signKii.cut() || signKi.sign()*signKii.sign()<=0) {
          active_element = true;
          break;
        }
      }

      if(!active_element){
        it_k = idx_from_background_mesh_[d].erase(it_k);
        continue;
      }


      // SAVE AND ERASE OLD INTERFACES EACH TIME STEP
      for( int it=0; it< n_tid; ++it) {
        auto it_gamma = interface_id_[it].find(std::make_pair(d, k));
        int nb_interface = (it_gamma == interface_id_[it].end())? 0 :it_gamma->second.size();
        std::vector<const Interface<Mesh>*> old_interface(nb_interface);
        std::vector<int> ss(nb_interface);
        for(int i=0;i<nb_interface;++i) old_interface[i] = it_gamma->second[i].first;
        for(int i=0;i<nb_interface;++i) ss[i] = it_gamma->second[i].second;
        if(it_gamma != interface_id_[it].end()) {
          auto ittt =   interface_id_[it].erase(it_gamma);
        }
        // PUT BACK INTERFACES
        for(int i=0; i<nb_interface;++i) {
          interface_id_[it][std::make_pair(d,nt[d])].push_back(std::make_pair(old_interface[i],  ss[i]));
        }
        // IS CUT SO NEED TO ADD INTERFACE AND SIGN
        const SignElement<Element> signK  = interface(it)->get_SignElement(kb);

        if (signK.cut()) {
          interface_id_[it][std::make_pair(d,nt[d])].push_back(std::make_pair(interface[it],  0));
        } else {
          in_active_mesh_[d][it][nt[d]] = false;
        }
      }
      it_k->second = nt[d];
      idx_in_background_mesh_[d].push_back(kb);
      nt[d]++;
      it_k++;
    }
  }

  idx_element_domain.push_back(0);
  for( int d=0; d<dom_size;++d){
    idx_in_background_mesh_[d].resize(nt[d]);
    idx_in_background_mesh_[d].shrink_to_fit();
    int sum_nt = idx_element_domain[d] + nt[d];
    idx_element_domain.push_back(sum_nt);
  }

}


//  constructor for basic 2 subdomains problem {1, -1}
template<typename Mesh>
void ActiveMesh<Mesh>::init(const TimeInterface<Mesh>& interface){

  int n_tid = interface.size();
  nb_quadrature_time_ = n_tid;
  in_active_mesh_.resize(10);
  for(int i=0;i<10;++i) in_active_mesh_[i].resize(nb_quadrature_time_);

  idx_in_background_mesh_[0].reserve(Th.nt);
  idx_in_background_mesh_[1].reserve(Th.nt);
  idx_element_domain.push_back(0);
  int nt0=0, nt1=0;
  for(int k = 0; k < Th.nt; ++k) {

    bool active_element = false;
    int s;
    for(int t=0;t<interface.size()-1;++t){
      const SignElement<Element> signKi  = interface(t)->get_SignElement(k);
      const SignElement<Element> signKii = interface(t+1)->get_SignElement(k);
      s = signKi.sign();
      if(signKi.cut() || signKii.cut() || signKi.sign()*signKii.sign()<=0) {
        active_element = true;
        break;
      }
    }

    if(active_element) {
      for( int it=0; it< n_tid; ++it) {
        const SignElement<Element> signK  = interface(it)->get_SignElement(k);
        int st = signK.sign();
        if(!signK.cut() && st == -1) {in_active_mesh_[0][it][nt0] = false; }
        if(!signK.cut() && st == 1)  {in_active_mesh_[1][it][nt1] = false; }

        if (signK.cut()){
          interface_id_[it][std::make_pair(0,nt0)].push_back(std::make_pair(interface[it],  1));
          interface_id_[it][std::make_pair(1,nt1)].push_back(std::make_pair(interface[it], -1));
        }
      }
      idx_in_background_mesh_[0].push_back(k);
      idx_from_background_mesh_[0][k] = nt0;
      idx_in_background_mesh_[1].push_back(k);
      idx_from_background_mesh_[1][k] = nt1;
      nt0++; nt1++;
    }
    else {
      int dom_add = (s<0);
      int dom_rm  = (s>0);
      int& nnt = (s > 0)? nt0 : nt1;
      idx_in_background_mesh_[dom_add].push_back(k);
      idx_from_background_mesh_[dom_add][k] = nnt;
      nnt++;
    }
  }
  idx_in_background_mesh_[0].resize(nt0);
  idx_in_background_mesh_[1].resize(nt1);
  idx_in_background_mesh_[0].shrink_to_fit();
  idx_in_background_mesh_[1].shrink_to_fit();

  idx_element_domain.push_back(nt0);
  idx_element_domain.push_back(nt0+nt1);

}

template<typename Mesh>
void ActiveMesh<Mesh>::truncate(const TimeInterface<Mesh>& interface,int sign_domain_remove){

  int n_tid = interface.size();
  nb_quadrature_time_ = n_tid;
  in_active_mesh_.resize(10);
  for(int i=0;i<10;++i) in_active_mesh_[i].resize(nb_quadrature_time_);

  int dom_size = this->get_nb_domain();
  idx_element_domain.resize(0);

  {
    for( int d=0; d<dom_size;++d){
      idx_in_background_mesh_[d].resize(0);
      int nt_max = idx_from_background_mesh_[d].size();
      idx_in_background_mesh_[d].reserve(nt_max);
    }
  }

  std::vector<int> nt(dom_size, 0.);
  for( int d=0; d<dom_size;++d){
    for(auto it_k = idx_from_background_mesh_[d].begin(); it_k != idx_from_background_mesh_[d].end() ;) {

      int kb  = it_k->first;
      int k   = it_k->second;

      bool active_element = false;
      int s;
      for(int t=0;t<interface.size()-1;++t){
        const SignElement<Element> signKi  = interface(t)->get_SignElement(kb);
        const SignElement<Element> signKii = interface(t+1)->get_SignElement(kb);
        s = signKi.sign();
        if(signKi.cut() || signKii.cut() || signKi.sign()*signKii.sign()<=0) {
          active_element = true;
          break;
        }
      }
      // REMOVE THE ELEMENT IN THE INPUT DOMAIN
      if(s == sign_domain_remove && !active_element) {
        it_k = idx_from_background_mesh_[d].erase(it_k);
        continue;
      }



      // SAVE AND ERASE OLD INTERFACES
      for( int it=0; it< n_tid; ++it) {
        auto it_gamma = interface_id_[it].find(std::make_pair(d, k));
        int nb_interface = (it_gamma == interface_id_[it].end())? 0 :it_gamma->second.size();
        std::vector<const Interface<Mesh>*> old_interface(nb_interface);
        std::vector<int> ss(nb_interface);
        for(int i=0;i<nb_interface;++i) old_interface[i] = it_gamma->second[i].first;
        for(int i=0;i<nb_interface;++i) ss[i] = it_gamma->second[i].second;
        if(it_gamma != interface_id_[it].end()) {
          auto ittt =   interface_id_[it].erase(it_gamma);
        }
        // PUT BACK INTERFACES
        for(int i=0; i<nb_interface;++i) {
          interface_id_[it][std::make_pair(d,nt[d])].push_back(std::make_pair(old_interface[i],  ss[i]));
        }
        // IS CUT SO NEED TO ADD INTERFACE AND SIGN
        const SignElement<Element> signK  = interface(it)->get_SignElement(kb);
        if (signK.cut()){
          interface_id_[it][std::make_pair(d,nt[d])].push_back(std::make_pair(interface[it],  -sign_domain_remove));
        }
        else if(signK.sign() == sign_domain_remove) {
          in_active_mesh_[d][it][nt[d]] = false;
        }
      }

      // SET NEW INDICES AND PUT BACK INTERFACES
      idx_in_background_mesh_[d].push_back(kb);
      it_k->second = nt[d];
      nt[d]++;
      it_k++;
    }
  }

  idx_element_domain.push_back(0);
  for( int d=0; d<dom_size;++d){
    idx_in_background_mesh_[d].resize(nt[d]);
    idx_in_background_mesh_[d].shrink_to_fit();
    int sum_nt = idx_element_domain[d] + nt[d];
    idx_element_domain.push_back(sum_nt);
  }

}



template<typename Mesh>
Physical_Partition<typename Mesh::Element> ActiveMesh<Mesh>::build_local_partition(const int k, int t) const {

  int nvc = Element::Rd::d+1;
  typedef SortArray<Ubyte, Element::Rd::d+1> ElementIdx;
  // GET THE ELEMENT TO BE CUT
  const Element& K((*this)[k]);
  Physical_Partition<Element> partition(K);
  std::vector<ElementIdx>  elements_idx;

  Ubyte iv[Element::nb_ntcut][nvc];

  // INITIALIZE THE LOCAL_PARTITION WITH K
  for( int e=0;e<Element::nb_ntcut ;e++){
    for(int i=0; i<nvc; ++i) {
      iv[e][i] = (i+2*e)%Element::nv;
    }
  }
  for(int i=0; i<Element::nv; ++i) {
    partition.add_node(K[i]);
  }
  for( int e=0;e<Element::nb_ntcut ;e++){
    // partition.add_element(ElementIdx(iv[e]));
    elements_idx.push_back(ElementIdx(iv[e]));
  }

  if(!isCut(k, t)) return partition;

  // START CUTTING PROCEDURE
  std::list<int> erased_element;
  vector<ElementIdx> new_element_idx;

  int domain = get_domain_element(k);
  int kloc = idxK_in_domain(k, domain);
  auto it = interface_id_[t].find(std::make_pair(domain, kloc));

  // std::cout << " ELEMENT \t" << k << "in domain " << domain << std::endl;
  // std::cout << it->second.size() << " interfaces " << std::endl;
  for(int i=0; i<it->second.size(); ++i) {

    // FRACTURE i DOES NOT CUT THIS ELEMENT => NEXT
  //   if(!cut_lines_[i]->is_cut_element(k)) continue;
  //   std::cout << " is cut by fracture \t" << i << std::endl;
  //   // WE PERFORM THE CUT
     int s = it->second[i].second;
     const Interface<Mesh>& interface(*(it->second[i].first));
     interface.cut_partition(partition, new_element_idx, erased_element,s);

    // SORT THE LIST FROM HIGHER TO LOWER INDEX
    // TO NOT DESTROYED INDICES WHEN ERASING
    erased_element.sort(std::greater<int>());
    // std::cout << " nb K to erase " << erased_element.size() << std::endl;
    // ERASING THE CUT ELEMENTS
    for(auto it = erased_element.begin(); it != erased_element.end(); ++it){
      // std::cout << " erased element " << *it << std::endl;
      // partition.erase_element(*it);
      elements_idx.erase(elements_idx.begin()+(*it));
    }
    // CREATING THE NEW ELEMENTS
     for(auto it = new_element_idx.begin(); it != new_element_idx.end(); ++it){
       // partition.add_element(*it);
       elements_idx.push_back(ElementIdx(*it));

     }
     // std::cout << " iteration " << i << "nb K " << partition.nb_element() << std::endl;
  }
  if(elements_idx.size() >= partition.max_size()){
    std::cout << " Need to increase maximum size of array element_idx in Physical partition" << std::endl;
    assert(0);
  }
  for(int i=0;i<elements_idx.size();++i){
    partition.set_element(i, elements_idx[i]);
  }
  return partition;
}

template<typename Mesh>
Physical_Partition<typename Mesh::Element::Face> ActiveMesh<Mesh>::build_local_partition(Face& face, const int k, int ifac, int t) const {

  Physical_Partition<typename Element::Face> partition(face);

  return partition;

}
typedef ActiveMesh<Mesh2> ActiveMeshT2;
typedef ActiveMesh<Mesh3> ActiveMeshT3;
typedef ActiveMesh<MeshQuad2> ActiveMeshQ2;
typedef ActiveMesh<MeshHexa> ActiveMeshQ3;










#endif
