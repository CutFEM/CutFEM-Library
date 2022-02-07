#include "macroElement.hpp"
#include "../common/SparseMatMap.hpp"

MacroElement::MacroElement(const FESpace2& vh, const double C) : GMacro() , Vh(vh) {
  double h = Vh[0].T.lenEdge(0);
  double meas = Vh[0].T.mesure();
  nb_element_0 = 0;
  nb_element_1 = 0;
  tol = C * h*h;//meas;

  std::cout << " tolerance \t" << tol << std::endl;

  find_small_element();
  std::cout << nb_element_0 << " \t in Omega 1 " << std::endl;
  std::cout << nb_element_1 << " \t in Omega 2 " << std::endl;
  find_root_element();
  std::cout << " root found " << std::endl;


}

void MacroElement::find_small_element() {
  for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
    if(!Vh.isCut(k)) continue;

    const FElement2& FK(Vh[k]);
    const int kb = Vh.Th(FK.T);
    const int domain = FK.whichDomain();

    CutData2 cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
    const Partition2& cutK =  Partition2(FK.T, cutData);  // build the cut
    ElementSignEnum the_part = cutK.what_part(domain);
    double areaCut = cutK.mesure(domain);

    if(areaCut < tol) {
      // small_or_fat_K[k] = small;   // -1 to small elements
      // idx_small_K_temp.push_back(k);
      // chain_position[k] = 0;
      // idxElement2idxSmallElement[k] = small_element.size();
      // small_element.push_back(SmallElement(k));

      small_element[k] = SmallElement(k);

      if(domain == 0) { nb_element_0++;}else{ nb_element_1++;}
    }
  }
}
void MacroElement::find_root_element(){
  vector<std::pair<int,int>> idx_small_K_temp(small_element.size());
  map<int,int> chain_position;
  vector<int> small_or_fat_K(Vh.nbElement);


  for(int i=0;i<small_or_fat_K.size();++i) small_or_fat_K[i] = i;
  int ii = 0;
  for(auto it=small_element.begin(); it!= small_element.end();++it) {
    idx_small_K_temp[ii++] = std::make_pair(it->second.index, it->first);;
    chain_position[it->second.index] = 0;
    small_or_fat_K[it->second.index] = small;
  }

  ii = 0;
  while (idx_small_K_temp.size() > 0) {
    int nb_of_small_K_left = idx_small_K_temp.size();
    // std::cout << nb_of_small_K_left << std::endl;
    for (int i=nb_of_small_K_left-1;i>=0;--i) {

      int k = idx_small_K_temp[i].first;
      int idx_Ks = idx_small_K_temp[i].second;
      SmallElement& Ks(small_element[idx_Ks]);
      const FElement2& FK(Vh[k]);
      int k_back = Vh.Th(FK.T);
      int the_domain = FK.whichDomain();

      int kn, ie;
      bool found_fat_neigh = false;
      for(int ifac = 0; ifac < 3; ++ifac) {    //loop over the edges / faces

        int ifacn = ifac;
        int kn_back = Vh.Th.ElementAdj(k_back,ifacn);
        if(kn_back == -1) continue;
        int kn_tmp = Vh.idxElementFromBackMesh(kn_back, the_domain);   // not in the domain
        if(kn_tmp ==-1) continue;

        if(small_or_fat_K[kn_tmp] == kn_tmp) {   // means kn_tmp is a root fat
          kn = kn_tmp;
          ie = (k < kn)? ifac : ifacn;
          int kk = (k < kn)?k: kn;
          chain_position[k] = 1;
          Ks.setChainPosition(1);
          found_fat_neigh = true;
          auto it = macro_element.find(kn);
          if(it != macro_element.end()){ // already exist

            it->second.add(k, std::make_pair(kk,ie));
          }
          else{
            macro_element[kn] = MElement(kn);
            macro_element[kn].add(k, std::make_pair(kk, ie));
          }
          break;
        }
        else if((small_or_fat_K[kn_tmp] != small)) {  // means it is a new fat

          if(chain_position[k] == 0 ||  chain_position[k] > chain_position[kn_tmp]){
            kn = kn_tmp;
            ie = (k < kn)? ifac : ifacn;;
            chain_position[k] = chain_position[kn_tmp]+1;
            Ks.setChainPosition(chain_position[kn_tmp]+1);
            found_fat_neigh = true;
          }
        }

      }

      // if(isFat(kn)) {
      if(found_fat_neigh) {
        int kk = (k<kn)? k:kn;
        small_or_fat_K[k] = small_or_fat_K[kn];
        idx_small_K_temp.erase(idx_small_K_temp.begin()+i); // remove the element.
        Ks.setRoot(small_or_fat_K[kn]);

        if(chain_position[k] != 1){
          int idx_root = small_or_fat_K[kn];
          auto it = macro_element.find(idx_root);
          it->second.add(k, std::make_pair(kk, ie));

        }
      }
    }
    ii++;
  }
}

MacroElementSurface::MacroElementSurface(const Interface2& gh, const double C) : GMacro() , interface(gh) {
  double h = (*interface.backMesh)[0].lenEdge(0);
  tol = C * h;

  std::cout << " tolerance macro surface\t" << tol << std::endl;
  find_small_element();
  std::cout << " Found " << small_element.size() << " small elements " << std::endl;
  find_root_element();


}

void MacroElementSurface::find_small_element() {
  for(int iface=interface.first_element(); iface<interface.last_element(); iface+= interface.next_element()) {

    const Face& face = interface[iface];
    const int kb = interface.idxElementOfFace(iface);
    const R meas = interface.computeDx(face).norm();

    if(meas > tol) continue;

    small_element[iface] = SmallElement(iface);
  }
}
void MacroElementSurface::find_root_element(){
  const Mesh2& Th(*interface.backMesh);

  for( auto it=small_element.begin(); it!=small_element.end(); ++it) {

    int iface = it->first;
    const Face& face = interface[iface];
    const int kb = interface.idxElementOfFace(iface);
    int root=-1, chain_position, ie;

    for(int i=0;i<2;++i) {
      chain_position = 0;
      int idx_node = face[i];
      ie = interface.edge_of_node_[idx_node];

      int root_K = check_direction(kb, ie, chain_position);
      assert(interface.face_of_element_.find(root_K) != interface.face_of_element_.end());
      root = interface.face_of_element_.find(root_K)->second;
      if(i==0) {
        it->second.setRoot(root);
        it->second.setChainPosition(chain_position);
        it->second.setEdgeDirection(ie);
      }
    }

    if(chain_position < it->second.chain_position) {
      it->second.setRoot(root);
      it->second.setChainPosition(chain_position);
      it->second.setEdgeDirection(ie);
    }

  }

  for( auto it=small_element.begin(); it!=small_element.end(); ++it) {

    int idx_root = it->second.index_root;
    auto pRoot = macro_element.find(idx_root);
    int k_loc  = it->first;
    int k  = interface.idxElementOfFace(k_loc);
    int ie = it->second.idx_edge_to_root;
    int je = ie;
    int kn = Th.ElementAdj(k, je);
    int kn_loc = interface.idxFaceOfElement(kn);
    assert(kn != -1);
    int kk = (k < kn)? k_loc : kn_loc;
    int ke = (k < kn)? ie : je;
    if(pRoot != macro_element.end()) {   // already exist
      pRoot->second.add(k_loc, std::make_pair(kk, ke));
    }
    else{
      macro_element[idx_root] = MElement(idx_root);
      macro_element[idx_root].add(k_loc, std::make_pair(kk, ke));
    }
  }

}
int MacroElementSurface::check_direction(const int k, const int ie, int& chain_position){
  assert(chain_position < interface.nbElement());
  const Mesh2& Th(*interface.backMesh);
  int e1 = ie;
  int kn = Th.ElementAdj(k,e1);
  chain_position++;
  int kn_loc = interface.idxFaceOfElement(kn);
  if(small_element.find(kn_loc) == small_element.end()) { // fat element
    return kn;
  }
  // find which edge to look at
  int iface = interface.face_of_element_.find(kn)->second;
  const Face& face = interface[iface];
  int ie_next = (interface.edge_of_node_[face[0]]==e1)?
  interface.edge_of_node_[face[1]] : interface.edge_of_node_[face[0]];

  return check_direction(kn, ie_next, chain_position);
}
