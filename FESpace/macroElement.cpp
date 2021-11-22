#include "MacroElement.hpp"
#include "../common/SparseMatMap.hpp"

MacroElement::MacroElement(const FESpace2& vh, const double C) : Vh(vh) {
  double h = Vh[0].T.lenEdge(0);
  double meas = Vh[0].T.mesure();
  nb_element_0 = 0;
  nb_element_1 = 0;
  tol = C * meas;

  std::cout << " tolerance \t" << tol << std::endl;

  find_small_K();
  std::cout << nb_element_0 << " \t in Omega 1 " << std::endl;
  std::cout << nb_element_1 << " \t in Omega 2 " << std::endl;
  find_root_element();
  std::cout << " root found " << std::endl;


}

void MacroElement::find_small_K() {
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
      idxElement2idxSmallElement[k] = small_element.size();
      small_element.push_back(SmallElement(k));

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
    idx_small_K_temp[ii++] = std::make_pair(it->index, ii);
    chain_position[it->index] = 0;
    small_or_fat_K[it->index] = small;
  }


  ii = 0;
  while (idx_small_K_temp.size() > 0) {
    int nb_of_small_K_left = idx_small_K_temp.size();
    std::cout << nb_of_small_K_left << std::endl;
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
          chain_position[k] = 1;
          Ks.setChainPosition(1);
          found_fat_neigh = true;
          idxRoot2idxMElement[kn] = macro_element.size();

          macro_element.push_back(MElement(kn, k, ie));
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
          int idx_ME = idxRoot2idxMElement[idx_root];
          macro_element[idx_ME].add(kk, ie);
        }
      }
    }
    ii++;
  }
}

void MacroElement::do_extension(const std::map<std::pair<int,int>,int>::const_iterator& it){
  int idx_Ks = it->first.first;
  int id_e = it->first.second;
  int handle = it->second;
  int idx_Kf = getIndexRootElement(idx_Ks);
  const FElement2& Ks(Vh[idx_Ks]);
  const FElement2& Kf(Vh[idx_Kf]);
  int ic0 = 0;
  for(int ic=0;ic< Kf.tfe->Sub_ToFE.size();++ic){
    if(Kf.tfe->Sub_ToFE[ic] == &DataFE<Mesh2>::RT0){
      do_extension_RT0(Ks,Kf,id_e,ic);
    }
    else if (Kf.tfe->Sub_ToFE[ic] == &DataFE<Mesh2>::BDM1){
      do_extension_BDM1(Ks,Kf,id_e,ic);
    }
    else if (Kf.tfe->Sub_ToFE[ic] == &DataFE<Mesh2>::P0){
      do_extension_P0(Ks,Kf,id_e,ic0);
    }
    else {
      assert(0);
    }
    ic0 += Kf.tfe->Sub_ToFE[ic]->N;
  }
}
void MacroElement::do_extension_P0  (const FElement2& Ks, const FElement2& Kf, int id_e, int ic){

  int idx0 = Kf.dfcbegin(ic);   // the pressure index
  if (S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] == 0) return;
  S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
  St[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
  S[std::make_pair(Ks.loc2glb(idx0),Kf.loc2glb(idx0))] = 1;
  St[std::make_pair(Kf.loc2glb(idx0),Ks.loc2glb(idx0))] = 1;

  dof2rm.insert(Ks.loc2glb(idx0));

}
void MacroElement::do_extension_RT0 (const FElement2& Ks, const FElement2& Kf, int id_e, int ic){
  int ndof = Kf.NbDoF();
  KNM<double> val(1, ndof);
  evaluate_dofRT0(Ks, id_e, Kf, val);
  int id_df = id_e;
  S [std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))]  = 0.;        // 1st component small element
  St[std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
  for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
    S[std::make_pair(Ks.loc2glb(id_df),Kf.loc2glb(i))]  = val(0,i);        // 1st component small element
    St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_df))] = val(0,i);        // 1st component small element
  }

  dof2rm.insert(Ks.loc2glb(id_df));

}
void MacroElement::do_extension_BDM1(const FElement2& Ks, const FElement2& Kf, int id_e, int ic){
  int ndofOnEdge = Kf.tfe->ndfonEdge;
  int ndof = Kf.NbDoF();
  KNM<double> val(ndofOnEdge, ndof);
  evaluate_dofBDM1(Ks, id_e, Kf, val);

  for(int df = 0; df < ndofOnEdge; ++df) {
    // if(  df==1) {
      int id_df = ndofOnEdge*id_e + df;
      S [std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))]  = 0.;        // 1st component small element
      St[std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
      for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
        // for(int i = 0; i < 6; ++i) {
        // if( (df==0 && i%2==0) || (df==1 && i%2==1) ) {
        // if(df==1 && i%2==1 ) {
        // if(  df==1) {
        S[std::make_pair(Ks.loc2glb(id_df),Kf.loc2glb(i))]  = val(df,i);        // 1st component small element
        St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_df))] = val(df,i);        // 1st component small element
        // }

      }
      // getchar();
      dof2rm.insert(Ks.loc2glb(id_df));
    // }
  }

}

void MacroElement::make_S(){
  typedef typename Mesh2::Partition Partition;
  typedef CutData2 CutData;

  for(int i=0; i<Vh.NbDoF(); ++i){
    S[std::make_pair(i ,i)] = 1;
    St[std::make_pair(i ,i)] = 1;
  }
  dof2rm.clear();

  for(auto it=element_edge_handle.begin(); it !=element_edge_handle.end();++it) {

    int idx_Ks = it->first.first;
    int id_e = it->first.second;
    int handle = it->second;
    if(handle == good) {
      continue;         // do nothing
    }

    if(handle == extension) {
      do_extension(it);
    }
  }
  std::cout << " removing \t" << dof2rm.size() << "  dof " << std::endl;
}

void MacroElement::evaluate_dofRT0(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val) {

  val = 0.;
  KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
  What_d Fop = Fwhatd(0);
  double meas = FKs.T.lenEdge(e);
  R2 normal = FKs.T.EdgeOrientation(e)*FKs.T.N(e);
  for(int iq=0;iq<QF.getNbrOfQuads();++iq) {

    QuadraturePoint1d ip_1d(QF[iq]);
    R2 ip_KsHat = FKs.T.toKref(ip_1d, e);
    R2 mip_Ks   = FKs.map(ip_KsHat);
    R2 ip_KfHat = FKf.T.toKref(mip_Ks);

    FKf.BF(Fop_D0, ip_KfHat, bf);

    for(int i=0;i<3;++i) {
      val(0,i) += meas*ip_1d.getWeight()*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
    }
  }
}
void MacroElement::evaluate_dofBDM1(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val) {

  val = 0.;
  KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
  What_d Fop = Fwhatd(0);
  double meas = FKs.T.lenEdge(e);
  // R2 normal = FKs.T.EdgeOrientation(e)*FKs.T.N(e);
  R2 normal = -FKs.T.Edge(e).perp(); // contain the mesure of edge

  for(int iq=0;iq<QF.getNbrOfQuads();++iq) {

    QuadraturePoint1d ip_1d(QF[iq]);
    R2 ip_KsHat = FKs.T.toKref(ip_1d, e);
    R2 mip_Ks   = FKs.map(ip_KsHat);
    R2 ip_KfHat = FKf.T.toKref(mip_Ks);

    FKf.BF(Fop_D0, ip_KfHat, bf);

    for(int i=0;i<6;++i) {
      val(0,i) += FKs.T.EdgeOrientation(e)*ip_1d.getWeight()*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
      val(1,i) += (-6*ip_1d.x+3)          *ip_1d.getWeight()*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
    }
  }
}

void MacroElement::precond(std::map<std::pair<int,int>,double>& A, Rn & rhs) {
  int t0 = MPIcf::Wtime();
  int N = Vh.NbDoF();
  std::map<std::pair<int,int>,double> R;
  multiply(N, St, A, R);
  A.clear();
  multiply(N, R, S, A);
  Rn b(N);
  multiply(N, N, St, rhs, b);
  rhs = b;

  std::cout << " initial size \t" << rhs.size() << std::endl;

  // need to remove the bad dof
  removeDF(N, A, rhs);

  std::cout << " size after \t" << rhs.size() << std::endl;

  std::cout << " time precond \t" << MPIcf::Wtime() - t0 << std::endl;
}
void MacroElement::removeDF( int N, std::map<std::pair<int,int>,double>& A, Rn& b){

  std::map<std::pair<int,int>,double>  C;

  int i0=0, j=0, i=0;
  for( auto it = dof2rm.begin();it != dof2rm.end();++it) {
    int iend = *it;
    while(j < iend) {
      P [make_pair(i,j)] = 1;
      Pt[make_pair(j,i)] = 1;
      ++i, ++j;
    }
    j += 1;
  }

  while(j<N){
    P[make_pair(i,j)] = 1;
    Pt[make_pair(j,i)] = 1;
    ++i;
    ++j;
  }

  int ndf = dof2rm.size();
  int nline = N - ndf;
  int ncol  = N;

  SparseMatrixRC<double> AA (N    ,N   ,A );
  SparseMatrixRC<double> PP (nline,ncol,P );
  SparseMatrixRC<double> PPt(ncol,nline,Pt);
  multiply(PP, AA, C);
  SparseMatrixRC<double> CC(nline,ncol,C);
  multiply(CC, PPt, A);

  Rn x(nline, 0.);
  multiply(nline, ncol, P, b, x);
  b.resize(nline);
  b = x;
}


void MacroElement::tag_extension_edges() {

  // get all the dof of small Elements
  std::map<int, std::pair<int, int>> df2fix;
  element_edge_handle.clear();
  // First set the trivial dof (good and trivial exhaust)
  for(auto it = small_element.begin(); it!=small_element.end();++it) {
    int k = it->index;
    const FElement2& FK(Vh[k]);
    int the_domain = FK.whichDomain();
    int pos_k = it->chain_position;

    for(int e=0;e<3;++e) {
      int df = FK.loc2glb(e);
      int je = e;
      int kn = Vh.getNeighborElement(k, je, the_domain);
      auto it_df = df2fix.find(df);
      bool df_seen = (it_df != df2fix.end());

      if(isRootFat(kn)) {  // triavial good
        df2fix[df] = std::make_pair(k, 3*e+good);
        continue;
      }
      else{
        int pos_kn = (kn == -1) ? pos_k + 1 : small_element[idxElement2idxSmallElement[kn]].chain_position;
        if (pos_kn > pos_k || (pos_kn == pos_k && kn > k) ){
          df2fix[df] = std::make_pair(k, 3*e+extension);  // default value
        }
      }
    }
  }

  // create data structure element -> whatToDoOnEdges
  for(auto it = df2fix.begin(); it!=df2fix.end();++it) {
    int k = it->second.first;
    int handle = it->second.second%3;
    int id_e = it->second.second/3;
    element_edge_handle[std::make_pair(k, id_e)] = handle;
  }

  // for(auto it = element_edge_handle.begin(); it!=element_edge_handle.end();++it) {
  //   std::cout << it->first.first << "\t"
  //             << it->first.second << ") = \t"
  //             << it->second << std::endl;
  //   }
}
void MacroElement::tag_exhaust_edges() {

  element_edge_handle.clear();
  // get all the dof of small Elements
  std::map<int, std::pair<int, int>> df2fix;
  std::set<int> exhaust_element;
  // vector<std::pair<int,int>> small_K_temp = idx_small_K;

  vector<int> small_K_temp(small_element.size());
  int ii = 0;
  for(auto it=small_element.begin(); it!= small_element.end();++it) {
    small_K_temp[ii++] = it->index;
  }



  // First set the trivial dof (good and trivial exhaust)
  // for(auto it = idx_small_K.begin(); it!=idx_small_K.end();++it) {
  //   int k = it->first;
  // for (int i=idx_small_K.size()-1;i>=0;--i) {
  for (int i=small_element.size()-1;i>=0;--i) {
    // int k = idx_small_K[i].first;
    int k = small_element[i].index;//.first;
    const FElement2& FK(Vh[k]);
    int the_domain = FK.whichDomain();
    bool is_exhaust = false;
    int count = 0;
    int notGoodEdge= -1;
    for(int e=0;e<3;++e) {
      int df = FK.loc2glb(e);
      int je = e;
      int kn = Vh.getNeighborElement(k, je, the_domain);

      auto it = df2fix.find(df);
      bool df_seen = (it != df2fix.end());
      if(!df_seen) {
        df2fix[df] = std::make_pair(k, 3*e+extension);  // default value
      }
      if(kn == -1) {  // trivial exhaust
        df2fix[df] = std::make_pair(k, 3*e+exhaust);
        notGoodEdge = e;
        is_exhaust = true;
        continue;
      }
      if(isRootFat(kn)) {  // triavial good
        df2fix[df] = std::make_pair(k, 3*e+good);
        count++;
        continue;
      }
    }
    // set element as exhaust
    if(is_exhaust) {
       if(count != 2) exhaust_element.insert(k);
       if(count != 2) small_K_temp.erase(small_K_temp.begin()+i);
    }
    if(count == 2) {
      int df = FK.loc2glb(notGoodEdge);
      // df2fix.erase(df);
      df2fix[df] = std::make_pair(k, 3*notGoodEdge+good);
      small_K_temp.erase(small_K_temp.begin()+i);
    }
  }

  // eliminating element when found exhaust edge
  while (small_K_temp.size() > 0) {

    int nb_of_small_K_left = small_K_temp.size();
    std::cout << nb_of_small_K_left << std::endl;
    for (int i=nb_of_small_K_left-1;i>=0;--i) {
      int k = small_K_temp[i];//first;

      const FElement2& FK(Vh[k]);
      int the_domain = FK.whichDomain();
      bool is_exhaust = false;
      for(int e=0;e<3;++e) {
        int df = FK.loc2glb(e);

        int je = e;
        int kn = Vh.getNeighborElement(k, je, the_domain);

        auto it_exhaust = exhaust_element.find(kn);
        bool neighIsExhaust = (it_exhaust != exhaust_element.end() ); // is neighbor exhaust?
        if(neighIsExhaust) {  // if my neighbor is exhaust
          df2fix[df] = std::make_pair(k, 3*e+exhaust);
          is_exhaust = true;
          break;  // can stop now, only one exhaust edge
        }

      }
      if(is_exhaust) {
        // then I can be remove from the search list
        exhaust_element.insert(k);
        small_K_temp.erase(small_K_temp.begin()+i);
      }
    }
  }

  // create data structure element -> whatToDoOnEdges
  for(auto it = df2fix.begin(); it!=df2fix.end();++it) {
    int k = it->second.first;
    int handle = it->second.second%3;
    int id_e = it->second.second/3;
    element_edge_handle[std::make_pair(k, id_e)] = handle;
  }

  // for(auto it = element_edge_handle.begin(); it!=element_edge_handle.end();++it) {
  //   std::cout << it->first.first << "\t"
  //             << it->first.second << ") = \t"
  //             << it->second << std::endl;
  //   }
    // getchar();
}


MacroElementSurface::MacroElementSurface(const Interface2& gh, const double C) : interface(gh) {
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

    small_element[kb] = SmallElement(iface);
  }
}

void MacroElementSurface::find_root_element(){

  for( auto it=small_element.begin(); it!=small_element.end(); ++it) {

    const int kb = it->first;
    int iface = interface.face_of_element_.find(kb)->second;
    const Face& face = interface[iface];

    int root=-1, chain_position;
    // check where is the closest fat element, in the two directions
    // need to set those two variables
    for(int i=0;i<2;++i) {
      chain_position = 0;
      int idx_node = face[i];
      int ie = interface.edge_of_node_[idx_node];

      int root_K = check_direction(kb, ie, chain_position);
      assert(interface.face_of_element_.find(root_K) != interface.face_of_element_.end());

      root = interface.face_of_element_.find(root_K)->second;
      if(i==0) {
        it->second.setRoot(root);
        it->second.setChainPosition(chain_position);
        it->second.setEdgeDirection(i);
      }
    }
    if(chain_position < it->second.chain_position) {
      it->second.setRoot(root);
      it->second.setChainPosition(chain_position);
      it->second.setEdgeDirection(1);
    }

    // std::cout << iface << "\t" << chain_position << std::endl;

  }


  // create all root elements
  for( auto it=small_element.begin(); it!=small_element.end(); ++it) {

    int idx_root = it->second.index_root;
    auto pRoot = macro_element.find(idx_root);
    if(pRoot != macro_element.end()) {   // already exist
      pRoot->second.add(it->second.index, it->second.idx_edge_to_root);
    }
    else{
      macro_element[idx_root] = MElement(idx_root, it->second.index, it->second.idx_edge_to_root);
    }
  }

  // for( auto it=macro_element.begin(); it!=macro_element.end(); ++it) {
  //   std::cout << " Root : " << it->second.idx_root_element << std::endl;
  //   std::cout << " Element : " << std::endl;
  //
  //   for(int i = 0 ; i< it->second.idx_element.size(); ++i ){
  //     std::cout << it->second.idx_element[i] << std::endl;
  //   }
  //   std::cout << std::endl;
  // }

}

int MacroElementSurface::check_direction(const int k, const int ie, int& chain_position){
  assert(chain_position < interface.nbElement());
  const Mesh2& Th(*interface.backMesh);
  int e1 = ie;
  int kn = Th.ElementAdj(k,e1);
  chain_position++;
  if(small_element.find(kn) == small_element.end()) { // fat element
    return kn;
  }
  // find which edge to look at
  int iface = interface.face_of_element_.find(kn)->second;
  const Face& face = interface[iface];
  int ie_next = (interface.edge_of_node_[face[0]]==e1)?
  interface.edge_of_node_[face[1]] : interface.edge_of_node_[face[0]];

  return check_direction(kn, ie_next, chain_position);
}

//   ii = 0;
//   while (idx_small_K_temp.size() > 0) {
//     int nb_of_small_K_left = idx_small_K_temp.size();
//     std::cout << nb_of_small_K_left << std::endl;
//     for (int i=nb_of_small_K_left-1;i>=0;--i) {
//
//       int k = idx_small_K_temp[i].first;
//       int idx_Ks = idx_small_K_temp[i].second;
//       SmallElement& Ks(small_element[idx_Ks]);
//       const FElement2& FK(Vh[k]);
//       int k_back = Vh.Th(FK.T);
//       int the_domain = FK.whichDomain();
//
//       int kn, ie;
//       bool found_fat_neigh = false;
//       for(int ifac = 0; ifac < 3; ++ifac) {    //loop over the edges / faces
//
//         int ifacn = ifac;
//         int kn_back = Vh.Th.ElementAdj(k_back,ifacn);
//         if(kn_back == -1) continue;
//         int kn_tmp = Vh.idxElementFromBackMesh(kn_back, the_domain);   // not in the domain
//         if(kn_tmp ==-1) continue;
//
//         if(small_or_fat_K[kn_tmp] == kn_tmp) {
//           kn = kn_tmp;
//           ie = ifac;
//           chain_position[k] = 1;
//           Ks.setChainPosition(1);
//           found_fat_neigh = true;
//           idxRoot2idxMElement[kn] = macro_element.size();
//
//           macro_element.push_back(MElement(kn, k, ie));
//           break;
//         }
//         else if((small_or_fat_K[kn_tmp] != small)) {
//
//           if(chain_position[k] == 0 ||  chain_position[k] > chain_position[kn_tmp]){
//             kn = kn_tmp;
//             ie = ifac;
//             chain_position[k] = chain_position[kn_tmp]+1;
//             Ks.setChainPosition(chain_position[kn_tmp]+1);
//             found_fat_neigh = true;
//           }
//         }
//
//       }
//
//       // if(isFat(kn)) {
//       if(found_fat_neigh) {
//         edges_to_stabilize.push_back(std::make_pair(k, ie));
//         small_or_fat_K[k] = small_or_fat_K[kn];
//         idx_small_K_temp.erase(idx_small_K_temp.begin()+i); // remove the element.
//         Ks.setRoot(small_or_fat_K[kn]);
//
//         if(chain_position[k] != 1){
//           int idx_root = small_or_fat_K[kn];
//           int idx_ME = idxRoot2idxMElement[idx_root];
//           macro_element[idx_ME].add(k, ie);
//         }
//       }
//   }
//   ii++;
// }
// }



// void find_face_to_stabilize() {
//   int ii = 0;
//   // vector<int> idx_small_K_temp = idx_small_K;
//   while (idx_small_K_temp.size() > 0) {
//
//     int nb_of_small_K_left = idx_small_K_temp.size();
//     for (int i=nb_of_small_K_left-1;i>=0;--i) {
//
//       int k = idx_small_K_temp[i];
//       const FElement2& FK(Vh[k]);
//       int k_back = Vh.Th(FK.T);
//       int the_domain = FK.whichDomain();
//
//       int kn, ie;
//       bool found_fat_neigh = false;
//       for(int ifac = 0; ifac < 3; ++ifac) {    //loop over the edges / faces
//
//         int ifacn = ifac;
//         int kn_back = Vh.Th.ElementAdj(k_back,ifacn);
//         if(kn_back == -1) continue;
//         int kn_tmp = Vh.idxElementFromBackMesh(kn_back, the_domain);   // not in the domain
//         if(kn_tmp ==-1) continue;
//
//
//         if(isRootFat(kn_tmp)) {
//           kn = kn_tmp;
//           ie = ifac;
//           chain_position[k] = 1;
//           found_fat_neigh = true;
//
//           idxRoot2idxMElement[kn] = macro_element.size();
//           macro_element.push_back(MElement(kn, k, ie));
//           break;
//         }
//         else if(isFat(kn_tmp)) {
//           if(chain_position[k] == 0 ||  chain_position[k] > chain_position[kn_tmp]){
//             kn = kn_tmp;
//             ie = ifac;
//             chain_position[k] = chain_position[kn_tmp]+1;
//             found_fat_neigh = true;
//           }
//         }
//
//       }
//
//       // if(isFat(kn)) {
//       if(found_fat_neigh) {
//         edges_to_stabilize.push_back(std::make_pair(k, ie));
//         idx_small_K.push_back(std::make_pair(k, small_or_fat_K[kn]));
//         small_or_fat_K[k] = small_or_fat_K[kn];
//         idx_small_K_temp.erase(idx_small_K_temp.begin()+i); // remove the element.
//
//
//         if(chain_position[k] != 1){
//           int idx_root = small_or_fat_K[kn];
//           int idx_ME = idxRoot2idxMElement[idx_root];
//           macro_element[idx_ME].add(k, ie);
//         }
//       }
//
//     }
//     ii++;
//   }
//   small_or_fat_K.clear();
// }
//

//
// void print(std::string filename = "EdgeToStabilize.dat") {
//   std::ofstream plot;
//   plot.open(filename.c_str(), std::ofstream::out);
//
//   for(std::vector<std::pair<int,int>>::iterator it = edges_to_stabilize.begin(); it != edges_to_stabilize.end(); ++it){
//
//     int k = it->first;
//     int ifac = it->second;
//
//     R2 p1 = Vh[k].T[Mesh2::Element::nvedge[ifac][0]];
//     R2 p2 = Vh[k].T[Mesh2::Element::nvedge[ifac][1]];
//
//     plot << p1 << "\n" << p2 << "\n \n";
//   }
//   plot.close();
// }



// void make_S(){
//   typedef typename Mesh2::Partition Partition;
//   typedef CutData2 CutData;
//
//   for(int i=0; i<Vh.NbDoF(); ++i){
//     S[std::make_pair(i ,i)] = 1;
//     St[std::make_pair(i ,i)] = 1;
//   }
//   dof2rm.clear();
//   for(auto it=small_element.begin(); it !=small_element.end();++it) {
//
//     int idx_Ks = it->first;
//     int idx_Kf = it->second;
//
//     const FElement2& Ks(Vh[idx_Ks]);
//     const FElement2& Kf(Vh[idx_Kf]);
//
//     int ks_back = Vh.Th(Ks.T);
//     int kf_back = Vh.Th(Kf.T);
//     int the_domain = Ks.whichDomain();
//     int ndofOnEdge = Kf.tfe->ndfonEdge;
//     int ndof = Kf.NbDoF();
//     KNM<double> val(ndofOnEdge, ndof);
//     int k = idx_Ks;
//     int pos_k = chain_position[k];
//
//     for(int e=0;e<3;++e) {     // loop over edges
//
//       int je = e;
//       int kn_back = Vh.Th.ElementAdj(ks_back, je);
//       if(kn_back == -1) continue;
//       int kn = Vh.idxElementFromBackMesh(kn_back, the_domain);
//       if(isRootFat(kn)) continue; // good dof, common with a Fat root element
//
//       int pos_kn = (kn == -1) ? pos_k + 1 : chain_position[kn];
//       if (pos_kn < pos_k) continue;
//       if (pos_kn == pos_k && kn < k) continue;
//
//       evaluate_dofBDM1(Ks, e, Kf, val);
//
//       for(int df = 0; df < ndofOnEdge; ++df) {
//         int id_e = ndofOnEdge*e + df;
//         S [std::make_pair(Ks.loc2glb(id_e),Ks.loc2glb(id_e))]  = 0.;        // 1st component small element
//         St[std::make_pair(Ks.loc2glb(id_e),Ks.loc2glb(id_e))] = 0.;        // 1st component small element
//         // for(int ic=0; ic<1;++ic) {
//           // for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
//           for(int i = 0; i < 6; ++i) {
//             if( (df==0 && i%2==0) || (df==1 && i%2==1) ) {
//             // if(  df==1) {
//             S[std::make_pair(Ks.loc2glb(id_e),Kf.loc2glb(i))]  = val(df,i);        // 1st component small element
//             St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_e))] = val(df,i);        // 1st component small element
//           }
//           // }
//         }
//         // getchar();
//         dof2rm.insert(Ks.loc2glb(id_e));
//       }
//     }
//     // int idx0 = Kf.dfcbegin(2);   // the pressure index
//     // if (S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] == 0) continue;
//     // S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
//     // St[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
//     // S[std::make_pair(Ks.loc2glb(idx0),Kf.loc2glb(idx0))] = 1;
//     // St[std::make_pair(Kf.loc2glb(idx0),Ks.loc2glb(idx0))] = 1;
//     //
//     // dof2rm.insert(Ks.loc2glb(idx0));
//
//   }
//   std::cout << " removing \t" << dof2rm.size() << "  dof " << std::endl;
// }


// void precondDiag(int N,std::map<std::pair<int,int>,double>& NL,std::map<std::pair<int,int>,double>& DF, Rn & rhs) {
//
//   SparseMatrixRC<double> B (N    ,N   ,NL );
//   NL.clear();
//   std::map<std::pair<int,int>,double> Pl;
//
//   // create the diagonal Matrix
//   for(int i=0;i<B.n;++i){
//     for(int k=B.p[i];k<B.p[i+1];++k){
//       Pl[std::make_pair(i,i)] += B.a[k];
//     }
//   }
//   for(int i=0;i<B.n;++i){
//     Pl[std::make_pair(i,i)] = 1./Pl[std::make_pair(i,i)];
//   }
//
//   SparseMatrixRC<double> PPl(N,N,Pl);
//   SparseMatrixRC<double> DDF(N,N,DF);
//   multiply(PPl, DDF, DF);
//
//   Rn x(N, 0.);
//   multiply(N, N, Pl, rhs, x);
//   rhs.resize(N);
//   rhs = x;
// }
//
// void modifyA(std::map<std::pair<int,int>,double>& A, Rn & rhs){
//   typedef typename Mesh2::Partition Partition;
//   typedef CutData2 CutData;
//
//   // for(int i=0; i<Vh.NbDoF(); ++i){
//   //   S[std::make_pair(i ,i)] = 1;
//   //   St[std::make_pair(i ,i)] = 1;
//   // }
//   // dof2rm.clear();
//   for(auto it=idx_small_K.begin(); it !=idx_small_K.end();++it) {
//
//     int idx_Ks = it->first;
//     int idx_Kf = it->second;
//
//     const FElement2& Ks(Vh[idx_Ks]);
//     const FElement2& Kf(Vh[idx_Kf]);
//
//     int ks_back = Vh.Th(Ks.T);
//     int kf_back = Vh.Th(Kf.T);
//     int the_domain = Ks.whichDomain();
//     int ndofOnEdge = Kf.tfe->ndfonEdge;
//     int ndof = Kf.NbDoF();
//     KNM<double> val(ndofOnEdge, ndof);
//     int k = idx_Ks;
//     int pos_k = chain_position[k];
//
//     for(int e=0;e<3;++e) {     // loop over edges
//
//       int je = e;
//       int kn_back = Vh.Th.ElementAdj(ks_back, je);
//       if(kn_back == -1) continue;
//       int kn = Vh.idxElementFromBackMesh(kn_back, the_domain);
//       if(isRootFat(kn)) continue; // good dof, common with a Fat root element
//
//       int pos_kn = (kn == -1) ? pos_k + 1 : chain_position[kn];
//       if (pos_kn < pos_k) continue;
//       if (pos_kn == pos_k && kn < k) continue;
//
//       evaluate_dofBDM1(Ks, e, Kf, val);
//       // for(int df = 0; df < ndofOnEdge; ++df) {
//       for(int df = 2; df < 2; ++df) {
//         int id_e = ndofOnEdge*e + df;
//         S[std::make_pair(Ks.loc2glb(id_e),Ks.loc2glb(id_e))]  = 0.;        // 1st component small element
//         St[std::make_pair(Ks.loc2glb(id_e),Ks.loc2glb(id_e))] = 0.;        // 1st component small element
//         for(int ic=0; ic<2;++ic) {
//           for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
//             S[std::make_pair(Ks.loc2glb(id_e),Kf.loc2glb(i))]  = val(df,i);        // 1st component small element
//             St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_e))] = val(df,i);        // 1st component small element
//           }
//         }
//         dof2rm.insert(Ks.loc2glb(id_e));
//       }
//     }
//     int idx0 = Kf.dfcbegin(2);   // the pressure index
//     if (S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] == 0) continue;
//     S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
//     St[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
//     S[std::make_pair(Ks.loc2glb(idx0),Kf.loc2glb(idx0))] = 1;
//     St[std::make_pair(Kf.loc2glb(idx0),Ks.loc2glb(idx0))] = 1;
//
//     dof2rm.insert(Ks.loc2glb(idx0));
//
//   }
//   std::cout << " removing \t" << dof2rm.size() << "  dof " << std::endl;
// }
//
//

/*
void exhaust_algo(const std::map<std::pair<int,int>,double>& A, const Rn& b, R (*divfun)(const R2, const int, const int), bool grabs1Edge){
   typedef typename Mesh2::Partition Partition;
   typedef CutData2 CutData;
   typedef GFESpace<Mesh2> FESpace;
   typedef typename FESpace::FElement FElement;
   typedef typename FElement::QF QF;
   typedef typename QF::QuadraturePoint QuadraturePoint;
   const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));


   std::set<std::pair<std::pair<int,int>,float>> e_pair_glb; // ((idx_small_edge,idx_fat_edge),coef)
   std::set<std::pair<int,float>> inner_exhaust_pair; // (idx_small_edge,rhs_value)

   for(int loop = 0; loop<2; ++loop) {

     for(auto it=idx_small_K.begin(); it !=idx_small_K.end();++it) {

     int idx_Ks = it->first;
     int idx_Kf = it->second;

     const FElement2& Ks(Vh[idx_Ks]);
     const FElement2& Kf(Vh[idx_Kf]);

     int ks_back = Vh.Th(Ks.T);
     int kf_back = Vh.Th(Kf.T);
     int the_domain = Ks.whichDomain(); // 0 or 1

     // if (chain_position[idx_Ks] == 2) continue;

     // [Sets e_exhaust if it is found, otherwise it is pathological case =-1]
     int e_exhaust = -1;
     for(int e=0;e<3;++e) { // loop over every edge
       int je = e;
       int kn_back = Vh.Th.ElementAdj(ks_back, je);
       if(kn_back == -1) continue;
       int kn = Vh.idxElementFromBackMesh(kn_back, the_domain);
       if(isRootFat(kn)) continue; // good dof, common with a Fat root element

       if(kn == -1) e_exhaust = Ks.loc2glb(e); // [sets e_exhaust]
     }
     if(e_exhaust != -1 && loop == 0) continue; // [skips all elements with exhaust first time]
     if(e_exhaust == -1 && loop == 1) continue; // [skips all elements without exhaust 2nd time]

     std::set<std::pair<std::pair<int,int>,float>> e_pair_loc;   // size 1 or 2, gets filled with changed edge and replacement edge
     // [Priority order: skip root fat edge, set e_exhaust (happens once), checks if edge is already replaced]
     for(int e=0;e<3;++e) {
       int je = e;
       int kn_back = Vh.Th.ElementAdj(ks_back, je);
       if(kn_back == -1) continue;
       int kn = Vh.idxElementFromBackMesh(kn_back, the_domain);
       if(isRootFat(kn)) { // val = 1 since it is not replaced by a combination
         e_pair_loc.insert(std::make_pair(std::make_pair(e,Ks.loc2glb(e)),1.));
         // e_pair_glb.insert(std::make_pair(std::make_pair(Ks.loc2glb(e),Ks.loc2glb(e)),1.));
         continue; // good dof, common with a Fat root element
       }

       if(kn == -1 || e_exhaust==-1) { // exhaust pipe edge (OR element has no exhaust)

         int check_exhaust = e_exhaust;
         e_exhaust = Ks.loc2glb(e);
         A[std::make_pair(e_exhaust,e_exhaust)] = 1*Ks.T.EdgeOrientation(e);


         R val = 0;
         const R meas = Ks.getMeasure();
         for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
           QuadraturePoint ip(qf[ipq]); // integration point
           R2 mip = Ks.map(ip);
           const R Cint = meas * ip.getWeight();
           val += Cint * divfun(mip,0,the_domain);
         } b[e_exhaust] = val; //std::cout << val << std::endl;

         if(check_exhaust==-1) {
           std::cout << ">>>>>> INNER EXHAUST: " << e_exhaust << std::endl;
           e_pair_glb.insert(std::make_pair(std::make_pair(e_exhaust,e_exhaust),1));
           // inner_exhaust_pair.insert(std::make_pair(e_exhaust,val));
         } else {
           // for (auto it = inner_exhaust_pair.begin(); it != inner_exhaust_pair.end(); ++it) {
           //   int loc_inner_exhaust = -1;
           //   for (int e2=0;e2<3;e2++) {
           //     if(it->first == Ks.loc2glb(e2)) {
           //       loc_inner_exhaust = e2;
           //       break;
           //     }
           //   } if(loc_inner_exhaust != -1) {
           //     // b[e_exhaust] -= it->second;
           //     // std::cout << b[e_exhaust] << std::endl;
           //     e_pair_loc.insert(std::make_pair(std::make_pair(e_exhaust,it->first),Ks.T.EdgeOrientation(loc_inner_exhaust)));
           //   }
           // }
         }
         continue; // skip rest of code to the next edge
       }
       bool alreadyReplaced = false;
       for (auto it = e_pair_glb.begin(); it != e_pair_glb.end(); ++it) {
         if (it->first.first == Ks.loc2glb(e)) { //  && it->first.second != Ks.loc2glb(e) // has been added before
           e_pair_loc.insert(std::make_pair(std::make_pair(e,it->first.second),it->second));
           alreadyReplaced = true;
           if(grabs1Edge) break;
         }
       } if(alreadyReplaced) continue;

       if(grabs1Edge) {
         // find the parallel edge of the fat connected element
         R2 nS = Ks.T.EdgeOrientation(e)*Ks.T.N(e);
         int idx_eF = -1;
         for(int eF=0;eF<3;++eF){
           R2 u = Kf.T.Edge(eF);
           if( fabs((u , nS)) < 1e-13) {
             idx_eF = eF;
             break;
           }
         }
         assert(idx_eF != -1);

         // double sign = (nS, Kf.T.EdgeOrientation(idx_eF)*Kf.T.N(idx_eF));
         A[std::make_pair(Ks.loc2glb(e),Ks.loc2glb(e))] = 1;
         A[std::make_pair(Ks.loc2glb(e),Kf.loc2glb(idx_eF))] = -1; // sign==1 always in regular mesh
         b[Ks.loc2glb(e)] = 0;

         e_pair_loc.insert(std::make_pair(std::make_pair(e,Kf.loc2glb(idx_eF)),1));
         e_pair_glb.insert(std::make_pair(std::make_pair(Ks.loc2glb(e),Kf.loc2glb(idx_eF)),1));
       } else {
         A[std::make_pair(Ks.loc2glb(e),Ks.loc2glb(e))] = 1;
         b[Ks.loc2glb(e)] = 0;

         double val;
         for(int i = Kf.dfcbegin(0); i < Kf.dfcend(0); ++i) {
           // [takes poly of bf on the small side - doesnt work]
           // int ei = i;
           // int kfn_back = Vh.Th.ElementAdj(kf_back, ei);
           // if(kfn_back == ks_back) continue; // [skip the shared edge, val here is +-1..?!]

           RT0_dof_eval(Kf, i, Ks, e, &val);

           A[std::make_pair(Ks.loc2glb(e),Kf.loc2glb(i))] = -val;

           e_pair_loc.insert(std::make_pair(std::make_pair(e,Kf.loc2glb(i)),val));
           e_pair_glb.insert(std::make_pair(std::make_pair(Ks.loc2glb(e),Kf.loc2glb(i)),val));
         }
       }
     }

     if(e_pair_loc.size()>0) {
       std::cout << "On " << Ks.loc2glb(3) << " changed: " << std::endl;
       for (auto e_ch = e_pair_loc.begin(); e_ch != e_pair_loc.end(); ++e_ch) {
         std::cout << " from " << Ks.loc2glb(e_ch->first.first) << " to " << e_ch->second << "* " << e_ch->first.second << "\n";
       } // [if from=to, it means that this edge is unchanged (had a rootfat)]
     }

     std::cout << " out exhaust: " << e_exhaust << std::endl;
     for (auto e_ch = e_pair_loc.begin(); e_ch != e_pair_loc.end(); ++e_ch) {
       A[std::make_pair(e_exhaust,e_ch->first.second)] += e_ch->second*Ks.T.EdgeOrientation(e_ch->first.first);
       if (e_exhaust==1589 || e_exhaust==1586) {
         std::cout << e_ch->first.second << " : " << A[std::make_pair(e_exhaust,e_ch->first.second)] << " " << e_ch->second << std::endl;
       }
     }

     // if(loop==0) {
     //   for (auto e_ch = e_pair_loc.begin(); e_ch != e_pair_loc.end(); ++e_ch) {
     //     e_pair_glb.insert(std::make_pair(std::make_pair(e_exhaust,e_ch->first.second),A[std::make_pair(e_exhaust,e_ch->first.second)]));
     //   }
     // }

     // pressure
     int idx0 = 3;//Kf.dfcbegin(2);   // the pressure index (note Kf.dfcbegin(1)=0 since vector element)
     A[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 1;
     A[std::make_pair(Ks.loc2glb(idx0),Kf.loc2glb(idx0))] = -1;
     b[Ks.loc2glb(idx0)] = 0;

     } // end of for

   } // end of outer for loop < 2
 }
*/








//
//
//
//
//
//
//
//
//
//
//
// #ifndef _MACRO_ELEMENT_HPP
// #define _MACRO_ELEMENT_HPP
//
//
// #include "CutFESpace.hpp"
//
// struct MElement {
// public:
//
//   int idx_root_element;
//   vector<int> idx_element;
//   vector<std::pair<int,int>> inner_edge;
//
//   MElement(int idx_root, int idx_K, int idx_edge){
//     idx_root_element = idx_root;
//     idx_element.push_back(idx_root);
//     this->add(idx_K, idx_edge);
//   }
//
//   void add(int idx_K, int idx_edge) {
//     idx_element.push_back(idx_K);
//     inner_edge.push_back(std::make_pair(idx_K, idx_edge));
//   }
//
//   void print() const {
//     std::cout << " Root idx   : \t " << idx_root_element << std::endl;
//     std::cout << " small element : \n";
//     for(auto it = inner_edge.begin(); it != inner_edge.end();++it) {
//       std::cout << it->first << " with edge " << it->second << std::endl;
//     }
//   }
// };
//
// struct SmallElement {
//   int index;
//   int index_root;
//   int chain_position;
//
//   SmallElement(int idx = -1, int idx_root = -1) :index(idx), index_root(idx_root), chain_position(0) {}
//   void setRoot(int i) { index_root = i;}
//   void setChainPosition(int i) { chain_position = i;}
//
// };
// class MacroElement {
//
// const QuadratureFormular1d& QF = QF_GaussLegendre2;
//
// public:
//
//   enum{small=-1};
//   const int good = 0, extension = 1, exhaust = 2;
//
//   vector<int> small_or_fat_K; // all the element [-1 or idx_root]
//   vector<int> idx_small_K_temp;   // list of just the small element
//
//   vector<std::pair<int,int>> idx_small_K;   // list of just the small element
//   vector<std::pair<int,int>> edges_to_stabilize;
//   map<int,int> chain_position;
//
//   set<int> dof2rm;
//   // maybe only need one map <root, MElement>
//   vector<MElement> macro_element;
//   map<int, int> idxRoot2idxMElement;
//
//   vector<SmallElement> small_element;
//   map<int, int> idxElement2idxSmallElement;
//   // for exhaust algo
//   std::map<std::pair<int ,int>, int> element_edge_handle;
//
//   std::map<std::pair<int,int>,double> S, St;
//   std::map<std::pair<int,int>,double> P, Pt;
//
//   const FESpace2& Vh;
//
//   R tol;
//   int nb_element_0, nb_element_1;
//
//   MacroElement(const FESpace2& vh, const double C) : Vh(vh) , small_or_fat_K(vh.nbElement){
//     double h = Vh[0].T.lenEdge(0);
//     double meas = Vh[0].T.mesure();
//     nb_element_0 = 0;
//     nb_element_1 = 0;
//     tol = C * h*h;//meas;
//
//     for(int i=0;i<small_or_fat_K.size();++i) small_or_fat_K[i] = i;
//
//     // std::cout << small_or_fat_K.size()<< std::endl;
//     std::cout << " tolerance \t" << tol << std::endl;
//
//     find_small_K();
//     // std::cout << " Found " << idx_small_K_temp.size() << " small elements" << std::endl;
//     std::cout << nb_element_0 << " \t in Omega 1 " << std::endl;
//     std::cout << nb_element_1 << " \t in Omega 2 " << std::endl;
//     find_face_to_stabilize() ;
//     // std::cout << " Number of edge found \t " << edges_to_stabilize.size() << std::endl;
//
//     // for(int k=0; k<idx_small_K.size(); ++k){
//     //
//     //   int idx = idx_small_K[k].first;
//     //   std::cout << idx << "\t" <<   Vh.whichDomain(idx) << "\t"
//     //             << chain_position[idx] << std::endl;
//     //
//     // }
//
//
//   }
//
//   int getRootElement(int k) const {
//     auto it = idxElement2idxSmallElement.find(k);
//     if(it == idxElement2idxSmallElement.end()){
//       return k;
//     }else {
//       return small_element[it->second].index_root;
//     }
//   }
//
//   bool isFat(int k) const {
//     return (small_or_fat_K[k] != small);
//   }
//
//   bool isRootFat(int k) const {
//     return (small_or_fat_K[k] == k);
//     // return (idxRoot2idxMElement.find(k) != idxRoot2idxMElement.end());
//   }
//
//   void find_small_K() {
//     for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
//       if(!Vh.isCut(k)) continue;
//
//       const FElement2& FK(Vh[k]);
//       const int kb = Vh.Th(FK.T);
//       const int domain = FK.whichDomain();
//
//       CutData2 cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
//       const Partition2& cutK =  Partition2(FK.T, cutData);  // build the cut
//       ElementSignEnum the_part = cutK.what_part(domain);
//       double areaCut = cutK.mesure(domain);
//
//       if(areaCut < tol) {
//         small_or_fat_K[k] = small;   // -1 to small elements
//         idx_small_K_temp.push_back(k);
//         chain_position[k] = 0;
//         idxElement2idxSmallElement[k] = small_element.size();
//         small_element.push_back(SmallElement(k));
//
//         if(domain == 0) { nb_element_0++;}else{ nb_element_1++;}
//       }
//     }
//   }
//
// //   void find_root_element(){
// //     vector<std::pair<int,int>> idx_small_K_temp(small_element.size());
// //     map<int,int> chain_position;
// //     vector<int> small_or_fat_K(Vh.nbElement);
// //
// //
// //     for(int i=0;i<small_or_fat_K.size();++i) small_or_fat_K[i] = i;
// //     int ii = 0;
// //     for(auto it=small_element.begin(); it!= small_element.end();++it) {
// //       idx_small_K_temp[ii++] = std::make_pair(it->index, ii);
// //       chain_position[it->index] = 0;
// //       small_or_fat_K[it->index] = small;
// //     }
// //
// //
// //     ii = 0;
// //     while (idx_small_K_temp.size() > 0) {
// //       int nb_of_small_K_left = idx_small_K_temp.size();
// //       for (int i=nb_of_small_K_left-1;i>=0;--i) {
// //
// //         int k = idx_small_K_temp[i].first;
// //         int idx_Ks = idx_small_K_temp[i].second;
// //         SmallElement& Ks(small_element[idx_Ks]);
// //         const FElement2& FK(Vh[k]);
// //         int k_back = Vh.Th(FK.T);
// //         int the_domain = FK.whichDomain();
// //
// //         int kn, ie;
// //         bool found_fat_neigh = false;
// //         for(int ifac = 0; ifac < 3; ++ifac) {    //loop over the edges / faces
// //
// //           int ifacn = ifac;
// //           int kn_back = Vh.Th.ElementAdj(k_back,ifacn);
// //           if(kn_back == -1) continue;
// //           int kn_tmp = Vh.idxElementFromBackMesh(kn_back, the_domain);   // not in the domain
// //           if(kn_tmp ==-1) continue;
// //
// //           if(small_or_fat_K[kn_tmp] == kn_tmp) {
// //             kn = kn_tmp;
// //             ie = ifac;
// //             chain_position[k] = 1;
// //             Ks.setChainPosition(1);
// //             found_fat_neigh = true;
// //             idxRoot2idxMElement[kn] = macro_element.size();
// //
// //             macro_element.push_back(MElement(kn, k, ie));
// //             break;
// //           }
// //           else if((small_or_fat_K[kn_tmp] != small)) {
// //
// //             if(chain_position[k] == 0 ||  chain_position[k] > chain_position[kn_tmp]){
// //               kn = kn_tmp;
// //               ie = ifac;
// //               chain_position[k] = chain_position[kn_tmp]+1;
// //               Ks.setChainPosition(chain_position[kn_tmp]+1);
// //               found_fat_neigh = true;
// //             }
// //           }
// //
// //         }
// //
// //         // if(isFat(kn)) {
// //         if(found_fat_neigh) {
// //           small_or_fat_K[k] = small_or_fat_K[kn];
// //           idx_small_K_temp.erase(idx_small_K_temp.begin()+i); // remove the element.
// //           Ks.setRoot(small_or_fat_K[kn]);
// //
// //           if(chain_position[k] != 1){
// //             int idx_root = small_or_fat_K[kn];
// //             int idx_ME = idxRoot2idxMElement[idx_root];
// //             macro_element[idx_ME].add(k, ie);
// //           }
// //         }
// //     }
// //     ii++;
// //   }
// // }
//
//
//   void find_face_to_stabilize() {
//     int ii = 0;
//     // vector<int> idx_small_K_temp = idx_small_K;
//     while (idx_small_K_temp.size() > 0) {
//
//       int nb_of_small_K_left = idx_small_K_temp.size();
//       for (int i=nb_of_small_K_left-1;i>=0;--i) {
//
//         int k = idx_small_K_temp[i];
//         const FElement2& FK(Vh[k]);
//         int k_back = Vh.Th(FK.T);
//         int the_domain = FK.whichDomain();
//
//         int kn, ie;
//         bool found_fat_neigh = false;
//         for(int ifac = 0; ifac < 3; ++ifac) {    //loop over the edges / faces
//
//           int ifacn = ifac;
//           int kn_back = Vh.Th.ElementAdj(k_back,ifacn);
//           if(kn_back == -1) continue;
//           int kn_tmp = Vh.idxElementFromBackMesh(kn_back, the_domain);   // not in the domain
//           if(kn_tmp ==-1) continue;
//
//
//           if(isRootFat(kn_tmp)) {
//             kn = kn_tmp;
//             ie = ifac;
//             chain_position[k] = 1;
//             found_fat_neigh = true;
//
//             idxRoot2idxMElement[kn] = macro_element.size();
//             macro_element.push_back(MElement(kn, k, ie));
//             break;
//           }
//           else if(isFat(kn_tmp)) {
//             if(chain_position[k] == 0 ||  chain_position[k] > chain_position[kn_tmp]){
//               kn = kn_tmp;
//               ie = ifac;
//               chain_position[k] = chain_position[kn_tmp]+1;
//               found_fat_neigh = true;
//             }
//           }
//
//         }
//
//         // if(isFat(kn)) {
//         if(found_fat_neigh) {
//           edges_to_stabilize.push_back(std::make_pair(k, ie));
//           idx_small_K.push_back(std::make_pair(k, small_or_fat_K[kn]));
//           small_or_fat_K[k] = small_or_fat_K[kn];
//           idx_small_K_temp.erase(idx_small_K_temp.begin()+i); // remove the element.
//
//
//           if(chain_position[k] != 1){
//             int idx_root = small_or_fat_K[kn];
//             int idx_ME = idxRoot2idxMElement[idx_root];
//             macro_element[idx_ME].add(k, ie);
//           }
//         }
//
//       }
//       ii++;
//     }
//     small_or_fat_K.clear();
//   }
//   //
//   // void print(std::string filename = "EdgeToStabilize.dat") {
//   //   std::ofstream plot;
//   //   plot.open(filename.c_str(), std::ofstream::out);
//   //
//   //   for(std::vector<std::pair<int,int>>::iterator it = edges_to_stabilize.begin(); it != edges_to_stabilize.end(); ++it){
//   //
//   //     int k = it->first;
//   //     int ifac = it->second;
//   //
//   //     R2 p1 = Vh[k].T[Mesh2::Element::nvedge[ifac][0]];
//   //     R2 p2 = Vh[k].T[Mesh2::Element::nvedge[ifac][1]];
//   //
//   //     plot << p1 << "\n" << p2 << "\n \n";
//   //   }
//   //   plot.close();
//   // }
//   // void make_S2(){
//   //   typedef typename Mesh2::Partition Partition;
//   //   typedef CutData2 CutData;
//   //
//   //   for(int i=0; i<Vh.NbDoF(); ++i){
//   //     S[std::make_pair(i ,i)] = 1;
//   //     St[std::make_pair(i ,i)] = 1;
//   //   }
//   //   dof2rm.clear();
//   //   for(auto it=small_element.begin(); it !=small_element.end();++it) {
//   //
//   //     int idx_Ks = it->first;
//   //     int idx_Kf = it->second;
//   //
//   //     const FElement2& Ks(Vh[idx_Ks]);
//   //     const FElement2& Kf(Vh[idx_Kf]);
//   //
//   //     int ks_back = Vh.Th(Ks.T);
//   //     int kf_back = Vh.Th(Kf.T);
//   //     int the_domain = Ks.whichDomain();
//   //     int ndofOnEdge = Kf.tfe->ndfonEdge;
//   //     int ndof = Kf.NbDoF();
//   //     KNM<double> val(ndofOnEdge, ndof);
//   //     int k = idx_Ks;
//   //     int pos_k = chain_position[k];
//   //
//   //     for(int e=0;e<3;++e) {     // loop over edges
//   //
//   //       int je = e;
//   //       int kn_back = Vh.Th.ElementAdj(ks_back, je);
//   //       if(kn_back == -1) continue;
//   //       int kn = Vh.idxElementFromBackMesh(kn_back, the_domain);
//   //       if(isRootFat(kn)) continue; // good dof, common with a Fat root element
//   //
//   //       int pos_kn = (kn == -1) ? pos_k + 1 : chain_position[kn];
//   //       if (pos_kn < pos_k) continue;
//   //       if (pos_kn == pos_k && kn < k) continue;
//   //
//   //       evaluate_dofBDM1(Ks, e, Kf, val);
//   //
//   //       for(int df = 0; df < ndofOnEdge; ++df) {
//   //         int id_e = ndofOnEdge*e + df;
//   //         S [std::make_pair(Ks.loc2glb(id_e),Ks.loc2glb(id_e))]  = 0.;        // 1st component small element
//   //         St[std::make_pair(Ks.loc2glb(id_e),Ks.loc2glb(id_e))] = 0.;        // 1st component small element
//   //         // for(int ic=0; ic<1;++ic) {
//   //           // for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
//   //           for(int i = 0; i < 6; ++i) {
//   //             if( (df==0 && i%2==0) || (df==1 && i%2==1) ) {
//   //             // if(  df==1) {
//   //             S[std::make_pair(Ks.loc2glb(id_e),Kf.loc2glb(i))]  = val(df,i);        // 1st component small element
//   //             St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_e))] = val(df,i);        // 1st component small element
//   //           }
//   //           // }
//   //         }
//   //         // getchar();
//   //         dof2rm.insert(Ks.loc2glb(id_e));
//   //       }
//   //     }
//   //     // int idx0 = Kf.dfcbegin(2);   // the pressure index
//   //     // if (S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] == 0) continue;
//   //     // S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
//   //     // St[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
//   //     // S[std::make_pair(Ks.loc2glb(idx0),Kf.loc2glb(idx0))] = 1;
//   //     // St[std::make_pair(Kf.loc2glb(idx0),Ks.loc2glb(idx0))] = 1;
//   //     //
//   //     // dof2rm.insert(Ks.loc2glb(idx0));
//   //
//   //   }
//   //   std::cout << " removing \t" << dof2rm.size() << "  dof " << std::endl;
//   // }
//   // void make_S(){
//   //   typedef typename Mesh2::Partition Partition;
//   //   typedef CutData2 CutData;
//   //
//   //   for(int i=0; i<Vh.NbDoF(); ++i){
//   //     S[std::make_pair(i ,i)] = 1;
//   //     St[std::make_pair(i ,i)] = 1;
//   //   }
//   //   dof2rm.clear();
//   //   for(auto it=small_element.begin(); it !=small_element.end();++it) {
//   //
//   //     int idx_Ks = it->first;
//   //     int idx_Kf = it->second;
//   //
//   //     const FElement2& Ks(Vh[idx_Ks]);
//   //     const FElement2& Kf(Vh[idx_Kf]);
//   //
//   //     int ks_back = Vh.Th(Ks.T);
//   //     int kf_back = Vh.Th(Kf.T);
//   //     int the_domain = Ks.whichDomain();
//   //     int ndofOnEdge = Kf.tfe->ndfonEdge;
//   //     int ndof = Kf.NbDoF();
//   //     KNM<double> val(ndofOnEdge, ndof);
//   //     int k = idx_Ks;
//   //     int pos_k = chain_position[k];
//   //
//   //     for(int e=0;e<3;++e) {     // loop over edges
//   //
//   //       int je = e;
//   //       int kn_back = Vh.Th.ElementAdj(ks_back, je);
//   //       if(kn_back == -1) continue;
//   //       int kn = Vh.idxElementFromBackMesh(kn_back, the_domain);
//   //       if(isRootFat(kn)) continue; // good dof, common with a Fat root element
//   //
//   //       int pos_kn = (kn == -1) ? pos_k + 1 : chain_position[kn];
//   //       if (pos_kn < pos_k) continue;
//   //       if (pos_kn == pos_k && kn < k) continue;
//   //
//   //       evaluate_dofBDM1(Ks, e, Kf, val);
//   //
//   //       for(int df = 0; df < ndofOnEdge; ++df) {
//   //         int id_e = ndofOnEdge*e + df;
//   //         S [std::make_pair(Ks.loc2glb(id_e),Ks.loc2glb(id_e))]  = 0.;        // 1st component small element
//   //         St[std::make_pair(Ks.loc2glb(id_e),Ks.loc2glb(id_e))] = 0.;        // 1st component small element
//   //         // for(int ic=0; ic<1;++ic) {
//   //           // for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
//   //           for(int i = 0; i < 6; ++i) {
//   //             if( (df==0 && i%2==0) || (df==1 && i%2==1) ) {
//   //             // if(  df==1) {
//   //             S[std::make_pair(Ks.loc2glb(id_e),Kf.loc2glb(i))]  = val(df,i);        // 1st component small element
//   //             St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_e))] = val(df,i);        // 1st component small element
//   //           }
//   //           // }
//   //         }
//   //         // getchar();
//   //         dof2rm.insert(Ks.loc2glb(id_e));
//   //       }
//   //     }
//   //     // int idx0 = Kf.dfcbegin(2);   // the pressure index
//   //     // if (S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] == 0) continue;
//   //     // S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
//   //     // St[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
//   //     // S[std::make_pair(Ks.loc2glb(idx0),Kf.loc2glb(idx0))] = 1;
//   //     // St[std::make_pair(Kf.loc2glb(idx0),Ks.loc2glb(idx0))] = 1;
//   //     //
//   //     // dof2rm.insert(Ks.loc2glb(idx0));
//   //
//   //   }
//   //   std::cout << " removing \t" << dof2rm.size() << "  dof " << std::endl;
//   // }
//   //
//   // void evaluate_dofRT0(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val) {
//   //
//   //   val = 0.;
//   //   KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
//   //   What_d Fop = Fwhatd(0);
//   //   double meas = FKs.T.lenEdge(e);
//   //   R2 normal = FKs.T.EdgeOrientation(e)*FKs.T.N(e);
//   //   for(int iq=0;iq<QF.getNbrOfQuads();++iq) {
//   //
//   //     QuadraturePoint1d ip_1d(QF[iq]);
//   //     R2 ip_KsHat = FKs.T.toKref(ip_1d, e);
//   //     R2 mip_Ks   = FKs.map(ip_KsHat);
//   //     R2 ip_KfHat = FKf.T.toKref(mip_Ks);
//   //
//   //     FKf.BF(Fop_D0, ip_KfHat, bf);
//   //
//   //     for(int i=0;i<3;++i) {
//   //       val(0,i) += meas*ip_1d.getWeight()*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
//   //     }
//   //   }
//   // }
//   //
//   // void evaluate_dofBDM1(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val) {
//   //
//   //   val = 0.;
//   //   KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
//   //   What_d Fop = Fwhatd(0);
//   //   double meas = FKs.T.lenEdge(e);
//   //   // R2 normal = FKs.T.EdgeOrientation(e)*FKs.T.N(e);
//   //   R2 normal = -FKs.T.Edge(e).perp(); // contain the mesure of edge
//   //
//   //   for(int iq=0;iq<QF.getNbrOfQuads();++iq) {
//   //
//   //     QuadraturePoint1d ip_1d(QF[iq]);
//   //     R2 ip_KsHat = FKs.T.toKref(ip_1d, e);
//   //     R2 mip_Ks   = FKs.map(ip_KsHat);
//   //     R2 ip_KfHat = FKf.T.toKref(mip_Ks);
//   //
//   //     FKf.BF(Fop_D0, ip_KfHat, bf);
//   //
//   //     for(int i=0;i<6;++i) {
//   //       val(0,i) += FKs.T.EdgeOrientation(e)*ip_1d.getWeight()*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
//   //       val(1,i) += (-6*ip_1d.x+3)          *ip_1d.getWeight()*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
//   //     }
//   //   }
//   // }
//   //
//   // void precond(std::map<std::pair<int,int>,double>& A, Rn & rhs) {
//   //   int t0 = MPIcf::Wtime();
//   //   int N = Vh.NbDoF();
//   //   std::map<std::pair<int,int>,double> R;
//   //   multiply(N, St, A, R);
//   //   A.clear();
//   //   multiply(N, R, S, A);
//   //   Rn b(N);
//   //   multiply(N, N, St, rhs, b);
//   //   rhs = b;
//   //
//   //   std::cout << " initial size \t" << rhs.size() << std::endl;
//   //
//   //   // need to remove the bad dof
//   //   removeDF(N, A, rhs);
//   //
//   //   std::cout << " size after \t" << rhs.size() << std::endl;
//   //
//   //   std::cout << " time precond \t" << MPIcf::Wtime() - t0 << std::endl;
//   // }
//   //
//   // void removeDF( int N, std::map<std::pair<int,int>,double>& A, Rn& b){
//   //
//   //   std::map<std::pair<int,int>,double>  C;
//   //
//   //   int i0=0, j=0, i=0;
//   //   for( auto it = dof2rm.begin();it != dof2rm.end();++it) {
//   //     int iend = *it;
//   //     while(j < iend) {
//   //       P [make_pair(i,j)] = 1;
//   //       Pt[make_pair(j,i)] = 1;
//   //       ++i, ++j;
//   //     }
//   //     j += 1;
//   //   }
//   //
//   //   while(j<N){
//   //     P[make_pair(i,j)] = 1;
//   //     Pt[make_pair(j,i)] = 1;
//   //     ++i;
//   //     ++j;
//   //   }
//   //
//   //   int ndf = dof2rm.size();
//   //   int nline = N - ndf;
//   //   int ncol  = N;
//   //
//   //   SparseMatrixRC<double> AA (N    ,N   ,A );
//   //   SparseMatrixRC<double> PP (nline,ncol,P );
//   //   SparseMatrixRC<double> PPt(ncol,nline,Pt);
//   //   multiply(PP, AA, C);
//   //   SparseMatrixRC<double> CC(nline,ncol,C);
//   //   multiply(CC, PPt, A);
//   //
//   //   Rn x(nline, 0.);
//   //   multiply(nline, ncol, P, b, x);
//   //   b.resize(nline);
//   //   b = x;
//   // }
//   //
//   // void precondDiag(int N,std::map<std::pair<int,int>,double>& NL,std::map<std::pair<int,int>,double>& DF, Rn & rhs) {
//   //
//   //   SparseMatrixRC<double> B (N    ,N   ,NL );
//   //   NL.clear();
//   //   std::map<std::pair<int,int>,double> Pl;
//   //
//   //   // create the diagonal Matrix
//   //   for(int i=0;i<B.n;++i){
//   //     for(int k=B.p[i];k<B.p[i+1];++k){
//   //       Pl[std::make_pair(i,i)] += B.a[k];
//   //     }
//   //   }
//   //   for(int i=0;i<B.n;++i){
//   //     Pl[std::make_pair(i,i)] = 1./Pl[std::make_pair(i,i)];
//   //   }
//   //
//   //   SparseMatrixRC<double> PPl(N,N,Pl);
//   //   SparseMatrixRC<double> DDF(N,N,DF);
//   //   multiply(PPl, DDF, DF);
//   //
//   //   Rn x(N, 0.);
//   //   multiply(N, N, Pl, rhs, x);
//   //   rhs.resize(N);
//   //   rhs = x;
//   // }
//   //
//   // void modifyA(std::map<std::pair<int,int>,double>& A, Rn & rhs){
//   //   typedef typename Mesh2::Partition Partition;
//   //   typedef CutData2 CutData;
//   //
//   //   // for(int i=0; i<Vh.NbDoF(); ++i){
//   //   //   S[std::make_pair(i ,i)] = 1;
//   //   //   St[std::make_pair(i ,i)] = 1;
//   //   // }
//   //   // dof2rm.clear();
//   //   for(auto it=idx_small_K.begin(); it !=idx_small_K.end();++it) {
//   //
//   //     int idx_Ks = it->first;
//   //     int idx_Kf = it->second;
//   //
//   //     const FElement2& Ks(Vh[idx_Ks]);
//   //     const FElement2& Kf(Vh[idx_Kf]);
//   //
//   //     int ks_back = Vh.Th(Ks.T);
//   //     int kf_back = Vh.Th(Kf.T);
//   //     int the_domain = Ks.whichDomain();
//   //     int ndofOnEdge = Kf.tfe->ndfonEdge;
//   //     int ndof = Kf.NbDoF();
//   //     KNM<double> val(ndofOnEdge, ndof);
//   //     int k = idx_Ks;
//   //     int pos_k = chain_position[k];
//   //
//   //     for(int e=0;e<3;++e) {     // loop over edges
//   //
//   //       int je = e;
//   //       int kn_back = Vh.Th.ElementAdj(ks_back, je);
//   //       if(kn_back == -1) continue;
//   //       int kn = Vh.idxElementFromBackMesh(kn_back, the_domain);
//   //       if(isRootFat(kn)) continue; // good dof, common with a Fat root element
//   //
//   //       int pos_kn = (kn == -1) ? pos_k + 1 : chain_position[kn];
//   //       if (pos_kn < pos_k) continue;
//   //       if (pos_kn == pos_k && kn < k) continue;
//   //
//   //       evaluate_dofBDM1(Ks, e, Kf, val);
//   //       // for(int df = 0; df < ndofOnEdge; ++df) {
//   //       for(int df = 2; df < 2; ++df) {
//   //         int id_e = ndofOnEdge*e + df;
//   //         S[std::make_pair(Ks.loc2glb(id_e),Ks.loc2glb(id_e))]  = 0.;        // 1st component small element
//   //         St[std::make_pair(Ks.loc2glb(id_e),Ks.loc2glb(id_e))] = 0.;        // 1st component small element
//   //         for(int ic=0; ic<2;++ic) {
//   //           for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
//   //             S[std::make_pair(Ks.loc2glb(id_e),Kf.loc2glb(i))]  = val(df,i);        // 1st component small element
//   //             St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_e))] = val(df,i);        // 1st component small element
//   //           }
//   //         }
//   //         dof2rm.insert(Ks.loc2glb(id_e));
//   //       }
//   //     }
//   //     int idx0 = Kf.dfcbegin(2);   // the pressure index
//   //     if (S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] == 0) continue;
//   //     S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
//   //     St[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
//   //     S[std::make_pair(Ks.loc2glb(idx0),Kf.loc2glb(idx0))] = 1;
//   //     St[std::make_pair(Kf.loc2glb(idx0),Ks.loc2glb(idx0))] = 1;
//   //
//   //     dof2rm.insert(Ks.loc2glb(idx0));
//   //
//   //   }
//   //   std::cout << " removing \t" << dof2rm.size() << "  dof " << std::endl;
//   // }
//   //
//   void tag_extension_edges() {
//
//     // get all the dof of small Elements
//     std::map<int, std::pair<int, int>> df2fix;
//     element_edge_handle.clear();
//     // First set the trivial dof (good and trivial exhaust)
//     for(auto it = small_element.begin(); it!=small_element.end();++it) {
//       int k = it->index;
//       const FElement2& FK(Vh[k]);
//       int the_domain = FK.whichDomain();
//       int pos_k = it->chain_position;
//
//       for(int e=0;e<3;++e) {
//         int df = FK.loc2glb(e);
//         int je = e;
//         int kn = Vh.getNeighborElement(k, je, the_domain);
//         auto it_df = df2fix.find(df);
//         bool df_seen = (it_df != df2fix.end());
//
//         if(isRootFat(kn)) {  // triavial good
//           df2fix[df] = std::make_pair(k, 3*e+good);
//           continue;
//         }
//         else{
//           // int pos_kn = (kn == -1) ? pos_k + 1 : chain_position[kn];
//           int pos_kn = (kn == -1) ? pos_k + 1 : small_element[idxElement2idxSmallElement[kn]].chain_position;
//           if (pos_kn > pos_k || (pos_kn == pos_k && kn > k) ){
//             df2fix[df] = std::make_pair(k, 3*e+extension);  // default value
//           }
//         }
//       }
//     }
//
//     // create data structure element -> whatToDoOnEdges
//     for(auto it = df2fix.begin(); it!=df2fix.end();++it) {
//       int k = it->second.first;
//       int handle = it->second.second%3;
//       int id_e = it->second.second/3;
//       element_edge_handle[std::make_pair(k, id_e)] = handle;
//     }
//
//     for(auto it = element_edge_handle.begin(); it!=element_edge_handle.end();++it) {
//       std::cout << it->first.first << "\t"
//                 << it->first.second << ") = \t"
//                 << it->second << std::endl;
//       }
//   }
//   void tag_exhaust_edges() {
//
//     element_edge_handle.clear();
//     // get all the dof of small Elements
//     std::map<int, std::pair<int, int>> df2fix;
//     std::set<int> exhaust_element;
//     vector<std::pair<int,int>> small_K_temp = idx_small_K;
//
//     // vector<int> small_K_temp(small_element.size());
//     // int ii = 0;
//     // for(auto it=small_element.begin(); it!= small_element.end();++it) {
//     //   small_K_temp[ii++] = it->index;
//     // }
//
//
//
//     // First set the trivial dof (good and trivial exhaust)
//     // for(auto it = idx_small_K.begin(); it!=idx_small_K.end();++it) {
//     //   int k = it->first;
//     for (int i=idx_small_K.size()-1;i>=0;--i) {
//     // for (int i=small_element.size()-1;i>=0;--i) {
//       int k = idx_small_K[i].first;
//       // int k = small_element[i].index;//.first;
//       const FElement2& FK(Vh[k]);
//       int the_domain = FK.whichDomain();
//       bool is_exhaust = false;
//       int count = 0;
//       int notGoodEdge= -1;
//       for(int e=0;e<3;++e) {
//         int df = FK.loc2glb(e);
//         int je = e;
//         int kn = Vh.getNeighborElement(k, je, the_domain);
//
//         auto it = df2fix.find(df);
//         bool df_seen = (it != df2fix.end());
//         if(!df_seen) {
//           df2fix[df] = std::make_pair(k, 3*e+extension);  // default value
//         }
//         if(kn == -1) {  // trivial exhaust
//           df2fix[df] = std::make_pair(k, 3*e+exhaust);
//           notGoodEdge = e;
//           is_exhaust = true;
//           continue;
//         }
//         if(isRootFat(kn)) {  // triavial good
//           df2fix[df] = std::make_pair(k, 3*e+good);
//           count++;
//           continue;
//         }
//       }
//       // set element as exhaust
//       if(is_exhaust) {
//          if(count != 2) exhaust_element.insert(k);
//          if(count != 2) small_K_temp.erase(small_K_temp.begin()+i);
//       }
//       if(count == 2) {
//         int df = FK.loc2glb(notGoodEdge);
//         // df2fix.erase(df);
//         df2fix[df] = std::make_pair(k, 3*notGoodEdge+good);
//         small_K_temp.erase(small_K_temp.begin()+i);
//       }
//     }
//
//     // eliminating element when found exhaust edge
//     while (small_K_temp.size() > 0) {
//
//       int nb_of_small_K_left = small_K_temp.size();
//       std::cout << nb_of_small_K_left << std::endl;
//       for (int i=nb_of_small_K_left-1;i>=0;--i) {
//         int k = small_K_temp[i].first;
//
//         const FElement2& FK(Vh[k]);
//         int the_domain = FK.whichDomain();
//         bool is_exhaust = false;
//         for(int e=0;e<3;++e) {
//           int df = FK.loc2glb(e);
//
//           int je = e;
//           int kn = Vh.getNeighborElement(k, je, the_domain);
//
//           auto it_exhaust = exhaust_element.find(kn);
//           bool neighIsExhaust = (it_exhaust != exhaust_element.end() ); // is neighbor exhaust?
//           if(neighIsExhaust) {  // if my neighbor is exhaust
//             df2fix[df] = std::make_pair(k, 3*e+exhaust);
//             is_exhaust = true;
//             break;  // can stop now, only one exhaust edge
//           }
//
//         }
//         if(is_exhaust) {
//           // then I can be remove from the search list
//           exhaust_element.insert(k);
//           small_K_temp.erase(small_K_temp.begin()+i);
//         }
//       }
//     }
//
//     // create data structure element -> whatToDoOnEdges
//     for(auto it = df2fix.begin(); it!=df2fix.end();++it) {
//       int k = it->second.first;
//       int handle = it->second.second%3;
//       int id_e = it->second.second/3;
//       element_edge_handle[std::make_pair(k, id_e)] = handle;
//     }
//
//     for(auto it = element_edge_handle.begin(); it!=element_edge_handle.end();++it) {
//       std::cout << it->first.first << "\t"
//                 << it->first.second << ") = \t"
//                 << it->second << std::endl;
//       }
//   }
//
