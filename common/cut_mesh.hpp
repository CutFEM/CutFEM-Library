#ifndef _CUT_MESH_HPP
#define _CUT_MESH_HPP

#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "base_interface.hpp"

// class enum{Explicit_Partition}


template<typename E >
struct Cut_Part {

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

  int get_sign() const {return sign_cut_;}
  // const Partition<E>& get_partition() const {return partition_;}
  void get_list_node (vector<typename E::Rd>& node) const{ partition_->get_list_node(node, sign_cut_); }
  CutElement2 get_element(int k) const {return partition_->get_element(k);}
  int nb_element() const {return partition_->nb_element();}


  Cut_Part(const Cut_Part<E>& p) :pp(p.pp), ip(p.ip) {
    if(p.partition_ == &p.pp) partition_ = &pp;
    else partition_ = &ip;
  }

  bool multi_interface() const {return partition_ == &pp;}

};


template<typename Mesh>
class Cut_Mesh {

public:
  typedef typename Mesh::Element Element;
  typedef typename Mesh::Rd Rd;
  typedef typename Mesh::BorderElement BorderElement;
  typedef typename Element::RdHat RdHat;// for parametrization
  static const int nea=Element::nea; //  numbering of adj (4 in Tet,  3 in Tria, 2 in seg)
  static const int nva=Element::nva; //  numbering of vertex in Adj hyperface


  typedef SignPatternTrait2 SignPattern;
  typedef RefPatch2 RefPatch;
  typedef RefPartition2 RefPartition;
  // typedef Partition2 Partition;

  const Mesh& Th;

  // std::vector<std::vector<int>> idx_in_background_mesh_;   // [domain][idxK_cutMesh] -> idxK_backMesh
  std::vector<std::vector<int>> idx_in_background_mesh_;   // [domain][idxK_cutMesh] -> idxK_backMesh
  std::vector<std::map<int,int>>  idx_from_background_mesh_; // [domain](idxK_backMesh) -> idxK_cutMesh

  std::map<std::pair<int,int>, std::vector<std::pair<const Interface<Mesh>*, int>>> interface_id_; // (domain_id, idx_k) -> [time_quad][n_interface](interface, sign)

  std::vector<int> idx_element_domain;

public:
  // Create a CutMesh without cut on the backMesh
  // Usefull if wanna add sub domains
  Cut_Mesh(const Mesh& th) : Th(th) {
    idx_in_background_mesh_.reserve(10);
    idx_from_background_mesh_.reserve(10);
    idx_in_background_mesh_.resize(1);
    idx_from_background_mesh_.resize(1);
    idx_in_background_mesh_[0].resize(Th.nt);

    for(int k = 0; k < Th.nt; ++k) {
      idx_in_background_mesh_[0][k] = k;
      idx_from_background_mesh_[0][k] = k;
    }
    idx_element_domain.push_back(0);
    idx_element_domain.push_back(Th.nt);
  }
  // Give the background mesh and a sign Function defined on the mesh nodes
  // Will create 2 subdomains
  Cut_Mesh(const Mesh& th, const Interface<Mesh>& interface) :  Th(th) {
    idx_in_background_mesh_.reserve(10);
    idx_from_background_mesh_.reserve(10);
    idx_in_background_mesh_.resize(2);
    idx_from_background_mesh_.resize(2);
    this->init(interface);
  }

  void truncate(const Interface<Mesh>& interface, int sign_domain);
  void add(const Interface<Mesh>& interface);
  void create_surface_mesh(const Interface<Mesh>& interface);
  void create_surface_mesh(const Time_Interface<Mesh>& interface);
private:
  void init(const Interface<Mesh>& interface);

  bool check_exist(int k, int dom) const {
    const auto it = idx_from_background_mesh_[dom].find(k);
    if(it == idx_from_background_mesh_[dom].end()) return false;
    else return true;
  }
  int idxK_begin(int i) const {return this->idx_element_domain[i];}
  int idxK_in_domain(int k, int i) const {
    return k - this->idx_element_domain[i];
  }
  Physical_Partition<Element> build_local_partition(const int k) const ;

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
  int get_domain_element(int k) const {
    for(int i=0;i<this->get_nb_domain();++i) {
      if(k<idx_element_domain[i+1]) return i;
    }
    assert(0);
  }

  bool isCut(int k) const {
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_.find(std::make_pair(domain, kloc));
    return (it != interface_id_.end());
  }
  bool isCut(int k, int domain) const {
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_.find(std::make_pair(domain, kloc));
    return (it != interface_id_.end());
  }


  const Interface<Mesh>& get_interface(int k) const {
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_.find(std::make_pair(domain, kloc));
    assert(it != interface_id_.end());
    return *(it->second.at(0).first);
  }
  Partition<Element> get_partition(int k) const {
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_.find(std::make_pair(domain, kloc));
    assert(it != interface_id_.end());
    int kb = this->idxElementInBackMesh(k);
    return it->second.at(0).first->get_partition(kb);
  }
  int get_sign_cut(int k) const {
    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_.find(std::make_pair(domain, kloc));
    assert(it != interface_id_.end());
    return it->second.at(0).second;
  }
  Cut_Part<Element> get_cut_part(int k) const {

    int domain = get_domain_element(k);
    int kloc = idxK_in_domain(k, domain);
    auto it = interface_id_.find(std::make_pair(domain, kloc));
    assert(it != interface_id_.end());
    int kb = this->idxElementInBackMesh(k);

    if(it->second.size() == 1)
    return Cut_Part<Element>(it->second.at(0).first->get_partition(kb), it->second.at(0).second);
    else
    return Cut_Part<Element>(this->build_local_partition(k), 0);
  }

  int idxElementFromBackMesh(int k) const {
    assert(0);
    return Th.idxElementFromBackMesh(k);
  }
  int idxElementFromBackMesh(int k, int i) const {
    if(i==-1) assert(0);
    auto it = idx_from_background_mesh_[i].find(k);
    if(it ==  idx_from_background_mesh_[i].end()) return -1;
    return idxK_begin(i) + it->second;
  }

  int idxElementInBackMesh(int k) const {
    int i = this->get_domain_element(k);
    int l = idxK_in_domain(k,i);
    return idx_in_background_mesh_[i][l];
  }
  int ElementAdj(int k,int &j, int domain = 0) const {
    int kb  = this->idxElementInBackMesh(k);
    int kbn = this->Th.ElementAdj(kb,j);
    if(kbn == -1) return -1;
    this->idxElementFromBackMesh(kbn, domain);
  }

};



//  constructor for basic 2 subdomains problem {1, -1}
template<typename Mesh>
void Cut_Mesh<Mesh>::init(const Interface<Mesh>& interface){

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
      interface_id_[std::make_pair(0,nt0)].push_back(std::make_pair(&interface,  1));
      interface_id_[std::make_pair(1,nt1)].push_back(std::make_pair(&interface, -1));

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

  std::cout << " nb elements in 0 \t" << this->get_nb_element(0) << std::endl;
  std::cout << " nb elements in 1 \t" << this->get_nb_element(1) << std::endl;
  std::cout << " nb elements  \t" << this->get_nb_element() << std::endl;
}

template<typename Mesh>
void Cut_Mesh<Mesh>::add(const Interface<Mesh>& interface){

  int dom_size = this->get_nb_domain();
  idx_element_domain.resize(0);
  int sign_domain = 1;
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


      // std::cout << "domain \t" << d << " element back " << kb << "\t => loc id " << k << std::endl;
      auto it_gamma = interface_id_.find(std::make_pair(d, k));


      const SignElement<Element> signK = interface.get_SignElement(kb);

      int nb_interface = (it_gamma == interface_id_.end())? 0 :it_gamma->second.size();
      std::vector<const Interface<Mesh>*> old_interface(nb_interface);
      std::vector<int> ss(nb_interface);
      for(int i=0;i<nb_interface;++i) old_interface[i] = it_gamma->second[i].first;
      for(int i=0;i<nb_interface;++i) ss[i] = it_gamma->second[i].second;

      if(it_gamma != interface_id_.end()) {
        // std::cout << " found old interface " << std::endl;
        auto ittt =   interface_id_.erase(it_gamma);
      }

      if(signK.sign() == sign_domain || signK.cut()) {
        sub_is_cut = true;

        // if(new_dom_id == 3) {
        //   std::cout << nt[new_dom_id] << "\t" << kb << " nb interface found " << nb_interface << std::endl;
        // }
        idx_in_background_mesh_[new_dom_id].push_back(kb);
        idx_from_background_mesh_[new_dom_id][kb] = nt[new_dom_id];

        for(int i=0; i<nb_interface;++i) {
          interface_id_[std::make_pair(new_dom_id,nt[new_dom_id])].push_back(std::make_pair(old_interface[i],  ss[i]));
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
            // interface_id_[std::make_pair(new_dom_id,nt[new_dom_id])].push_back(std::make_pair(old_interface[i],  ss[i]));
            interface_id_[std::make_pair(d,nt[d])].push_back(std::make_pair(old_interface[i],  ss[i]));
          }
          interface_id_[std::make_pair(new_dom_id,nt[new_dom_id])].push_back(std::make_pair(&interface,  sign_domain));
          interface_id_[std::make_pair(d,nt[d])].push_back(std::make_pair(&interface,  -sign_domain));
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
          interface_id_[std::make_pair(d,nt[d])].push_back(std::make_pair(old_interface[i],  ss[i]));
        }
        it_k++;
        nt[d]++;
      }
    }

    if(sub_is_cut && d+1 !=dom_size){
       new_dom_id++;
       idx_in_background_mesh_.resize(new_dom_id+1);
       idx_from_background_mesh_.resize(new_dom_id+1);
       idx_in_background_mesh_[new_dom_id].resize(0);
       int nt_max = idx_from_background_mesh_[d+1].size();
       idx_in_background_mesh_[new_dom_id].reserve(nt_max);

     }
  }

  idx_element_domain.push_back(0);
  for( int d=0; d<new_dom_id+1;++d){
    idx_in_background_mesh_[d].resize(nt[d]);
    idx_in_background_mesh_[d].shrink_to_fit();
    int sum_nt = idx_element_domain[d] + nt[d];
    idx_element_domain.push_back(sum_nt);
    // std::cout << sum_nt << "\t" ;
  }
// std::cout << std::endl;
  // for( int d=0; d<dom_size +1;++d){
  //   std::cout << " domain " << d << "\t -------------" << std::endl;
  //   for(int i=0;i<nt[d];++i) {
  //     std::cout << i << "\t" << idx_in_background_mesh_[d][i] << std::endl;
  //   }
  //
  // }

//   for( int d=0; d<idx_element_domain.size();++d){
//   std::cout << idx_element_domain[d] << std::endl;
// }
  // for(auto it=interface_id_.begin(); it != interface_id_.end();++it) {
  //   std::cout << it->first.first << "\t" << it->first.second << "\t" << it->second[0].second << std::endl;
  // }

  for( int d=0; d<new_dom_id+1;++d){
  std::cout << " nb elements in  \t" << d << "\t" << this->get_nb_element(d) << std::endl;
  }
  // std::cout << " nb elements in 0 \t" << this->get_nb_element(0) << std::endl;
  // std::cout << " nb elements in 1 \t" << this->get_nb_element(1) << std::endl;

  std::cout << " nb elements  \t" << this->get_nb_element() << std::endl;

}

template<typename Mesh>
void Cut_Mesh<Mesh>::truncate(const Interface<Mesh>& interface,int sign_domain_remove){

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
      auto it_gamma = interface_id_.find(std::make_pair(d, k));
      const SignElement<Element> signK = interface.get_SignElement(kb);

      // REMOVE THE ELEMENT IN THE INPUT DOMAIN
      if(signK.sign() == sign_domain_remove) {

        it_k = idx_from_background_mesh_[d].erase(it_k);

        continue;
      }


      // SAVE AND ERASE OLD INTERFACES
      int nb_interface = (it_gamma == interface_id_.end())? 0 :it_gamma->second.size();
      std::vector<const Interface<Mesh>*> old_interface(nb_interface);
      std::vector<int> ss(nb_interface);
      for(int i=0;i<nb_interface;++i) old_interface[i] = it_gamma->second[i].first;
      for(int i=0;i<nb_interface;++i) ss[i] = it_gamma->second[i].second;
      if(it_gamma != interface_id_.end()) {
        auto ittt =   interface_id_.erase(it_gamma);
      }

      // SET NEW INDICES AND PUT BACK INTERFACES
      idx_in_background_mesh_[d].push_back(kb);
      it_k->second = nt[d];
      for(int i=0; i<nb_interface;++i) {
        interface_id_[std::make_pair(d,nt[d])].push_back(std::make_pair(old_interface[i],  ss[i]));
      }
      // IS CUT SO NEED TO ADD INTERFACE AND SIGN
      if (signK.cut()){
        interface_id_[std::make_pair(d,nt[d])].push_back(std::make_pair(&interface,  -sign_domain_remove));
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
// std::cout << std::endl;
  // for( int d=0; d<dom_size +1;++d){
  //   std::cout << " domain " << d << "\t -------------" << std::endl;
  //   for(int i=0;i<nt[d];++i) {
  //     std::cout << i << "\t" << idx_in_background_mesh_[d][i] << std::endl;
  //   }
  //
  // }

//   for( int d=0; d<idx_element_domain.size();++d){
//   std::cout << idx_element_domain[d] << std::endl;
// }
  // for(auto it=interface_id_.begin(); it != interface_id_.end();++it) {
  //   std::cout << it->first.first << "\t" << it->first.second << "\t" << it->second[0].second << std::endl;
  // }

  for( int d=0; d<dom_size;++d){
  std::cout << " nb elements in  \t" << d << "\t" << this->get_nb_element(d) << std::endl;
  }
  // std::cout << " nb elements in 0 \t" << this->get_nb_element(0) << std::endl;
  // std::cout << " nb elements in 1 \t" << this->get_nb_element(1) << std::endl;

  std::cout << " nb elements  \t" << this->get_nb_element() << std::endl;

}


template<typename Mesh>
void Cut_Mesh<Mesh>::create_surface_mesh(const Interface<Mesh>& interface){

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
      auto it_gamma = interface_id_.find(std::make_pair(d, k));
      const SignElement<Element> signK = interface.get_SignElement(kb);

      // REMOVE THE ELEMENT IN THE INPUT DOMAIN
      if(!signK.cut()) {

        it_k = idx_from_background_mesh_[d].erase(it_k);
        continue;
      }


      // SAVE AND ERASE OLD INTERFACES
      int nb_interface = (it_gamma == interface_id_.end())? 0 :it_gamma->second.size();
      std::vector<const Interface<Mesh>*> old_interface(nb_interface);
      std::vector<int> ss(nb_interface);
      for(int i=0;i<nb_interface;++i) old_interface[i] = it_gamma->second[i].first;
      for(int i=0;i<nb_interface;++i) ss[i] = it_gamma->second[i].second;
      if(it_gamma != interface_id_.end()) {
        auto ittt =   interface_id_.erase(it_gamma);
      }

      // SET NEW INDICES AND PUT BACK INTERFACES
      idx_in_background_mesh_[d].push_back(kb);
      it_k->second = nt[d];
      for(int i=0; i<nb_interface;++i) {
        interface_id_[std::make_pair(d,nt[d])].push_back(std::make_pair(old_interface[i],  ss[i]));
      }
      // IS CUT SO NEED TO ADD INTERFACE AND SIGN
      interface_id_[std::make_pair(d,nt[d])].push_back(std::make_pair(&interface,  0));
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

  for( int d=0; d<dom_size;++d){
    std::cout << " nb elements in  \t" << d << "\t" << this->get_nb_element(d) << std::endl;
  }
  std::cout << " nb elements  \t" << this->get_nb_element() << std::endl;

}

template<typename Mesh>
void Cut_Mesh<Mesh>::create_surface_mesh(const Time_Interface<Mesh>& interface){

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
      auto it_gamma = interface_id_.find(std::make_pair(d, k));

      bool active_element = false;
      for(int t=0;t<interface.size()-1;++t){
        const SignElement<Element> signKi  = interface(t)->get_SignElement(kb);
        const SignElement<Element> signKii = interface(t+1)->get_SignElement(kb);
        if(signKi.cut() || signKii.cut() || signKi.sign()*signKii.sign()<=0) {
          active_element = true;
          break;
        }
      }

      // REMOVE THE ELEMENT IN THE INPUT DOMAIN
      if(!active_element ) {
        it_k = idx_from_background_mesh_[d].erase(it_k);
        continue;
      }


      // SAVE AND ERASE OLD INTERFACES
      int nb_interface = (it_gamma == interface_id_.end())? 0 :it_gamma->second.size();
      std::vector<const Interface<Mesh>*> old_interface(nb_interface);
      std::vector<int> ss(nb_interface);
      for(int i=0;i<nb_interface;++i) old_interface[i] = it_gamma->second[i].first;
      for(int i=0;i<nb_interface;++i) ss[i] = it_gamma->second[i].second;
      if(it_gamma != interface_id_.end()) {
        auto ittt =   interface_id_.erase(it_gamma);
      }

      // SET NEW INDICES AND PUT BACK INTERFACES
      idx_in_background_mesh_[d].push_back(kb);
      it_k->second = nt[d];
      for(int i=0; i<nb_interface;++i) {
        interface_id_[std::make_pair(d,nt[d])].push_back(std::make_pair(old_interface[i],  ss[i]));
      }
      // IS CUT SO NEED TO ADD INTERFACE AND SIGN
      interface_id_[std::make_pair(d,nt[d])].push_back(std::make_pair(interface[0],  0));
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

  for( int d=0; d<dom_size;++d){
    std::cout << " nb elements in  \t" << d << "\t" << this->get_nb_element(d) << std::endl;
  }
  std::cout << " nb elements  \t" << this->get_nb_element() << std::endl;

}




template<typename Mesh>
Physical_Partition<typename Mesh::Element> Cut_Mesh<Mesh>::build_local_partition(const int k) const {

  int nvc = Element::Rd::d+1;
  typedef SortArray<int, Element::Rd::d+1> ElementIdx;
  // GET THE ELEMENT TO BE CUT
  const Element& K((*this)[k]);
  Physical_Partition<Element> partition(K);

  int iv[Element::nb_ntcut][nvc];

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
    partition.add_element(ElementIdx(iv[e]));
    // std::cout << " initial element " << e << std::endl;
    // for(int i=0; i<nvc; ++i) {
    //   std::cout << K[iv[e][i]] << "\t";
    // }
    // std::cout << std::endl;
  }

  if(!isCut(k)) return partition;

  // START CUTTING PROCEDURE
  std::list<int> erased_element;
  vector<ElementIdx> new_element_idx;

  int domain = get_domain_element(k);
  int kloc = idxK_in_domain(k, domain);
  auto it = interface_id_.find(std::make_pair(domain, kloc));

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
      partition.erase_element(*it);
    }
    // CREATING THE NEW ELEMENTS
     for(auto it = new_element_idx.begin(); it != new_element_idx.end(); ++it){
       partition.add_element(*it);
     }
     // std::cout << " iteration " << i << "nb K " << partition.nb_element() << std::endl;
  }
  return partition;
}

typedef Cut_Mesh<Mesh2> Cut_MeshT2;
typedef Cut_Mesh<Mesh3> Cut_MeshT3;
typedef Cut_Mesh<MeshQuad2> Cut_MeshQ2;
typedef Cut_Mesh<MeshHexa> Cut_MeshQ3;










#endif
