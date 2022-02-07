#include "extension.hpp"


void Extension::tag_extension_edges(const MacroElement& macro) {

  // get all the dof of small Elements
  std::vector<std::pair<int, int>> df2fix;
  const FESpace2& Vh(macro.Vh);
  mapMacro[&Vh] = &macro;
  int k0 = k_begin(Vh);

  // First set the trivial dof (good and trivial exhaust)
  for(auto it = macro.small_element.begin(); it!=macro.small_element.end();++it) {
    int k = it->second.index;
    const FElement2& FK(Vh[k]);
    int the_domain = FK.whichDomain();
    int pos_k = it->second.chain_position;

    for(int e=0;e<3;++e) {
      int je = e;
      int kn = Vh.getNeighborElement(k, je, the_domain);
      if(macro.isRootFat(kn) ) {  // triavial good
        df2fix.push_back(std::make_pair(k, 3*e+good));
        continue;
      }
      else{
        auto it_neigh =  macro.small_element.find(kn);
        // int pos_kn = (kn == -1) ? pos_k + 1 : macro.small_element[kn].chain_position;
        int pos_kn = (kn == -1) ? pos_k + 1 : it->second.chain_position;

        if (pos_kn > pos_k || (pos_kn == pos_k && kn > k) ){
            df2fix.push_back(std::make_pair(k, 3*e+extension));  // default value
        }
      }
    }
  }
  // create data structure element -> whatToDoOnEdges
  for(auto it = df2fix.begin(); it!=df2fix.end();++it) {
    int k = it->first;
    int handle = it->second%3;
    int id_e = it->second/3;
    element_edge_handle[std::make_pair(k+k0, id_e)] = handle;
  }

  // for(auto it = element_edge_handle.begin(); it!=element_edge_handle.end();++it) {
  //   std::cout << it->first.first << "\t"
  //             << it->first.second << ") = \t"
  //             << it->second << std::endl;
  //   }
}
void Extension::tag_extension_edges(const MacroElement& macro, const CHyperFace& b) {


  element_edge_handle.clear();
  const FESpace2& Vh(macro.Vh);
  mapMacro[&Vh] = &macro;
  int k0 = k_begin(Vh);

  for(auto it = macro.macro_element.begin(); it!=macro.macro_element.end();++it) {
    for( int ie=0; ie<it->second.inner_edge.size(); ++ie) {

      int k = it->second.inner_edge[ie].first;
      int e = it->second.inner_edge[ie].second;

      element_edge_handle[std::make_pair(k+k0, e)] = extension;
    }
  }

}

void Extension::tag_exhaust_edges(const MacroElement& macro) {

  const FESpace2& Vh(macro.Vh);
  mapMacro[&Vh] = &macro;
  element_edge_handle.clear();
  // get all the dof of small Elements
  std::map<std::pair<int, int>, int>& df2fix(element_edge_handle) ;
  std::set<int> exhaust_element;
  // vector<std::pair<int,int>> small_K_temp = idx_small_K;

  // vector<int> small_K_temp(small_element.size());
  // int ii = 0;
  // for(auto it=small_element.begin(); it!= small_element.end();++it) {
  //   small_K_temp[ii++] = it->second.index;
  // }


  int artificial_good_df = 0;
  // First set the trivial dof (good and trivial exhaust)
  for (auto it=macro.small_element.begin();it!=macro.small_element.end();++it) {
    int k = it->second.index;//.first;

    // need to check this loop using map
    // also need to uncomment erase part a bit further
    const FElement2& FK(Vh[k]);
    int the_domain = FK.whichDomain();
    bool is_exhaust = false;
    int count = 0;
    int count_exhaust = 0;
    int exhaustEdge = 0;
    int notGoodEdge= -1;

    // LOOP OVER EDGES
    for(int e=0;e<3;++e) {
      // int df = FK.loc2glb(e);
      int je = e;
      int kn = Vh.getNeighborElement(k, je, the_domain);

      // CHECK IF EDGE ALREADY SEEN
      auto it = df2fix.find(make_pair(kn, je));
      bool edge_seen = (it != df2fix.end());

      if(edge_seen) continue;
      if(kn == -1) {  // trivial exhaust
        df2fix[std::make_pair(k, e)] = exhaust;
        notGoodEdge = e;
        exhaustEdge = e;
        is_exhaust = true;
        count_exhaust++;
      }
      else if(macro.isRootFat(kn)) {  // triavial good
        df2fix[std::make_pair(k, e)] = good;
        count++;
      }
      else {  // default value
        df2fix[std::make_pair(k, e)] = extension;
        notGoodEdge = e;
      }
    }
    // set element as exhaust
    if(is_exhaust) {
       if(count != 2) exhaust_element.insert(k);
    }
    if(count_exhaust == 2 ) {
      // boundary element
      df2fix[std::make_pair(k, exhaustEdge)] = extension;
    }
    if(count == 2 ) {
      df2fix[std::make_pair(k, notGoodEdge)] = good;
      artificial_good_df += 1;
    }
  }
  // ALL ELEMENTS MUST BECOME EXHAUST
  while (exhaust_element.size() < macro.small_element.size() - artificial_good_df) {
    // LOOP OVER SMALL ELEMENTS
    for(auto it_small = macro.small_element.begin(); it_small!= macro.small_element.end();++it_small){
      int k = it_small->second.index;

      // CHECK IF ALREADY EXHAUST
      auto it_exhaust = exhaust_element.find(k);
      if(it_exhaust != exhaust_element.end()) continue;

      // NOT AN EXHAUST ELEMENT
      const FElement2& FK(Vh[k]);
      int the_domain = FK.whichDomain();
      bool is_exhaust = false;

      //CHECK THE NEIGHBOR TO GET THE EXHAUST EDGE
      for(int e=0;e<3;++e) {
        int je = e;
        int kn = Vh.getNeighborElement(k, je, the_domain);

        it_exhaust = exhaust_element.find(kn);

        bool neighIsExhaust = (it_exhaust != exhaust_element.end() ); // is neighbor exhaust?
        if(neighIsExhaust) {  // if my neighbor is exhaust
          df2fix[std::make_pair(k, e)] = exhaust;
          exhaust_element.insert(k);
          break;  // can stop now, only one exhaust edge
        }
      }
    }
  }
}

void Extension::make_S(){
  typedef typename Mesh2::Partition Partition;
  typedef CutData2 CutData;

  for(int i=0; i<problem.nDoF; ++i){
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
    if(handle == exhaust){
      //
      // do_exhaust
      // do_extension(it);

    }
  }
  std::cout << " removing \t" << dof2rm.size() << "  dof " << std::endl;
}


void Extension::do_extension(const std::map<std::pair<int,int>,int>::const_iterator& it){
  int idxG_Ks = it->first.first;
  int k_begin;
  const FESpace2& Vh(get_space_from_idxK(idxG_Ks, k_begin));
  int idx_Ks = idxG_Ks - k_begin;
  int id_e = it->first.second;
  int handle = it->second;
  const GMacro& macro(get_macro(Vh));
  int idx_Kf = macro.getIndexRootElement(idx_Ks);
  // getchar();
  const FElement2& Ks(Vh[idx_Ks]);
  const FElement2& Kf(Vh[idx_Kf]);


  int ic0 = 0;

  for(int ic=0;ic< Kf.tfe->Sub_ToFE.size();++ic){
    if(Kf.tfe->Sub_ToFE[ic] == &DataFE<Mesh2>::RT0){
      do_extension_RT0(Ks,Kf,id_e,ic);
    }
    else if (Kf.tfe->Sub_ToFE[ic] == &DataFE<Mesh2>::BDM1){
      // do_extension_BDM1(Ks,Kf,id_e,ic);
    }
    else if (Kf.tfe->Sub_ToFE[ic] == &DataFE<Mesh2>::RT1){
      // do_extension_RT1(Ks,Kf,id_e,ic0);
    }
    else if (Kf.tfe->Sub_ToFE[ic] == &DataFE<Mesh2>::P0){
      do_extension_P0(Ks,Kf,ic0);
    }
    else if (Kf.tfe->Sub_ToFE[ic] == &DataFE<Mesh2>::P1dc){
      // do_extension_P1dc(Ks,Kf,ic0);
    }
    else {
      assert(0);
    }
    ic0 += Kf.tfe->Sub_ToFE[ic]->N;
  }
}



void Extension::do_extension_P0   (const FElement2& Ks, const FElement2& Kf, int ic){
  int idx0 = Kf.dfcbegin(ic);   // the pressure index
  int ig0 = problem.mapIdx0[&Ks.Vh];

  if (S[std::make_pair(ig0+Ks.loc2glb(idx0),ig0+Ks.loc2glb(idx0))] == 0) return;
  S[std::make_pair(ig0+Ks.loc2glb(idx0),ig0+Ks.loc2glb(idx0))] = 0;
  St[std::make_pair(ig0+Ks.loc2glb(idx0),ig0+Ks.loc2glb(idx0))] = 0;
  S[std::make_pair(ig0+Ks.loc2glb(idx0),ig0+Kf.loc2glb(idx0))] = 1;
  St[std::make_pair(ig0+Kf.loc2glb(idx0),ig0+Ks.loc2glb(idx0))] = 1;

  dof2rm.insert(ig0+Ks.loc2glb(idx0));

}
void Extension::do_extension_RT0 (const FElement2& Ks, const FElement2& Kf, int id_e, int ic){
  int ndof = Kf.NbDoF();
  KNM<double> val(1, ndof);
  evaluate_dofRT0(Ks, id_e, Kf, val);
  int id_df = id_e;
  int ig0 = problem.mapIdx0[&Ks.Vh];

  S [std::make_pair(ig0+Ks.loc2glb(id_df),ig0+Ks.loc2glb(id_df))] = 0.;        // 1st component small element
  St[std::make_pair(ig0+Ks.loc2glb(id_df),ig0+Ks.loc2glb(id_df))] = 0.;        // 1st component small element
  for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
    S[std::make_pair(ig0+Ks.loc2glb(id_df),ig0+Kf.loc2glb(i))]  = val(0,i);        // 1st component small element
    St[std::make_pair(ig0+Kf.loc2glb(i),ig0+Ks.loc2glb(id_df))] = val(0,i);        // 1st component small element
  }

  dof2rm.insert(ig0+Ks.loc2glb(id_df));
}
void Extension::evaluate_dofRT0(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val) {

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
// void MacroElement::do_extension_P1dc (const FElement2& Ks, const FElement2& Kf, int ic){
//   int idx0 = Kf.dfcbegin(ic);   // the pressure index
//   if (S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] == 0) return;
//
//   int ndof = Kf.NbDoF();
//   KNM<double> val(3, 3);
//   evaluate_dofP1dc(Kf, Ks, val, ic);
//   for(int id_df = Ks.dfcbegin(ic),df=0 ; id_df < Ks.dfcend(ic); ++id_df,++df) {
//     S [std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
//     St[std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
//     for(int i = Kf.dfcbegin(ic),j=0; i < Kf.dfcend(ic); ++i,++j) {
//       S[std::make_pair(Ks.loc2glb(id_df),Kf.loc2glb(i))]  = val(df,j);        // 1st component small element
//       St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_df))] = val(df,j);        // 1st component small element
//     }
//     dof2rm.insert(Ks.loc2glb(id_df));
//   }
// }
// void MacroElement::do_extension_BDM1(const FElement2& Ks, const FElement2& Kf, int id_e, int ic){
//   int ndofOnEdge = Kf.tfe->ndfonEdge;
//   int ndof = Kf.NbDoF();
//   KNM<double> val(ndofOnEdge, ndof);
//   evaluate_dofBDM1(Ks, id_e, Kf, val);
//
//   for(int df = 0; df < ndofOnEdge; ++df) {
//     // if(  df==1) {
//       int id_df = ndofOnEdge*id_e + df;
//       S [std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))]  = 0.;        // 1st component small element
//       St[std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
//       for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
//         // for(int i = 0; i < 6; ++i) {
//         // if( (df==0 && i%2==0) || (df==1 && i%2==1) ) {
//         // if(df==1 && i%2==1 ) {
//         // if(  df==1) {
//         S[std::make_pair(Ks.loc2glb(id_df),Kf.loc2glb(i))]  = val(df,i);        // 1st component small element
//         St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_df))] = val(df,i);        // 1st component small element
//         // }
//
//       }
//       // getchar();
//       dof2rm.insert(Ks.loc2glb(id_df));
//     // }
//   }
//
// }
// void MacroElement::do_extension_RT1(const FElement2& Ks, const FElement2& Kf, int id_e, int ic){
//   int ndofOnEdge = Kf.tfe->ndfonEdge;
//   int ndof = Kf.NbDoF();
//   KNM<double> val(ndofOnEdge+2, ndof); // 2 bubble bf
//   evaluate_dofBDM1(Ks, id_e, Kf, val);
//
//   for(int df = 0; df < ndofOnEdge; ++df) {
//     // if(  df==1) {
//       int id_df = ndofOnEdge*id_e + df;
//       S [std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))]  = 0.;        // 1st component small element
//       St[std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
//       for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
//         S[std::make_pair(Ks.loc2glb(id_df),Kf.loc2glb(i))]  = val(df,i);        // 1st component small element
//         St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_df))] = val(df,i);        // 1st component small element
//       }
//       dof2rm.insert(Ks.loc2glb(id_df));
//     // }
//   }
//
//   for(int df = 0; df < 2; ++df) {
//     int id_df = 3*ndofOnEdge + df;
//     S [std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
//     St[std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
//     for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
//       S[std::make_pair(Ks.loc2glb(id_df),Kf.loc2glb(i))]  = val(df+ndofOnEdge,i);        // 1st component small element
//       St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_df))] = val(df+ndofOnEdge,i);        // 1st component small element
//     }
//     dof2rm.insert(Ks.loc2glb(id_df));
//   }
//
// }



// void MacroElement::evaluate_dofP1dc(const FElement2& FKs, const FElement2& FKf, Rnm& val, int ic) {
//
//   val = 0.;
//   KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
//   What_d Fop = Fwhatd(0);
//
//   for(int v=0;v<3;++v) {
//
//     R2 mip_Ks   = FKs.T.at(v);
//     R2 ip_KfHat = FKf.T.toKref(mip_Ks);
//
//     FKf.BF(Fop_D0, ip_KfHat, bf);
//     for(int i = FKs.dfcbegin(ic), j=0; i < FKs.dfcend(ic); ++i,++j) {
//       val(v,j) += bf(i,ic,0);
//     }
//   }
// }

// void MacroElement::evaluate_dofBDM1(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val) {
//   val = 0.;
//   KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
//   What_d Fop = Fwhatd(0);
//   double meas = FKs.T.lenEdge(e);
//   R2 normal = -FKs.T.Edge(e).perp(); // contain the mesure of edge
//   // R2 normal = FKs.T.EdgeOrientation(e)*FKs.T.N(e);
//
//
//
//   for(int iq=0;iq<QF.getNbrOfQuads();++iq) {
//
//     QuadraturePoint1d ip_1d(QF[iq]);
//     R2 ip_KsHat = FKs.T.toKref(ip_1d, e);
//     R2 mip_Ks   = FKs.map(ip_KsHat);
//     R2 ip_KfHat = FKf.T.toKref(mip_Ks);
//
//     FKf.BF(Fop_D0, ip_KfHat, bf);
//
//     for(int i=0;i<6;++i) {
//       val(0,i) += FKs.T.EdgeOrientation(e)*ip_1d.getWeight()*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
//       val(1,i) += (-6*ip_1d.x+3)          *ip_1d.getWeight()*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
//     }
//   }
//
// }
// void MacroElement::evaluate_dofRT1(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val) {
//   val = 0.;
//   int ndf = FKf.NbDoF();
//   KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
//   What_d Fop = Fwhatd(0);
//   double meas = FKs.T.lenEdge(e);
//   int eOrientation = FKs.T.EdgeOrientation(e);
//   // R2 normal = FKs.T.EdgeOrientation(e)*FKs.T.N(e);
//   R2 normal = -FKs.T.Edge(e).perp();//*eOrientation; // contain the mesure of edge
//
//   for(int iq=0;iq<QF.getNbrOfQuads();++iq) {
//
//     QuadraturePoint1d ip_1d(QF[iq]);
//     R2 ip_KsHat = FKs.T.toKref(ip_1d, e);
//     R2 mip_Ks   = FKs.map(ip_KsHat);
//     R2 ip_KfHat = FKf.T.toKref(mip_Ks);
//
//     FKf.BF(Fop_D0, ip_KfHat, bf);
//
//     // evaluate Fat_bf in small_dof : on edge
//     R l0 = QF[iq].x, l1 = 1 - QF[iq].x;
//     R p0 = (2 * l0 - l1) * 2;     // poly othogonaux to \lambda_1
//     R p1 = (2 * l1 - l0) * 2;     // poly othogonaux to \lambda_0
//     R lambda1 = eOrientation * p0 * QF[iq].a;    // [some quadrature function?]
//     R lambda0 = eOrientation * p1 * QF[iq].a;    //
//     if(eOrientation < 0) {
//       Exchange(lambda1, lambda0);
//     }
//     for(int i=0;i<ndf;++i) {
//       val(0,i) += lambda0*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
//       val(1,i) += lambda1*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
//     }
//
//
//   }
//
//   // evaluate Fat_bf in small_bubble_dof
//   R2 B[2] = {FKs.T.Edge(1), FKs.T.Edge(2)};
//   B[0] = B[0].perp();
//   B[1] = B[1].perp();
//   for(int iq=0;iq<QFK.getNbrOfQuads();++iq) {
//
//     QuadraturePoint2d ip_2d(QFK[iq]);
//     R2 mip_Ks   = FKs.map(ip_2d);
//     R2 ip_KfHat = FKf.T.toKref(mip_Ks);
//
//     FKf.BF(Fop_D0, ip_KfHat, bf);
//     double w = QFK[iq].a * 0.5;
//     for(int i=0;i<ndf;++i) {
//       val(2,i) += w*(bf(i,0,0)*B[0].x + bf(i,1,0)*B[0].y) ;
//       val(3,i) += w*(bf(i,0,0)*B[1].x + bf(i,1,0)*B[1].y) ;
//     }
//   }
// }

void Extension::do_extension() {
  // BUILD THE MATRIX S AND St
  this->make_S();

  // DO MULTIPLY St A S and  St b
  // we dont send
  Rn b(problem.rhs.size());
  this->precond(b);

}

void Extension::solve() {
  if(element_edge_handle.size() == 0) {
    problem.solve(problem.mat, problem.rhs);
    return;
  }

  // BUILD THE MATRIX S AND St
  this->make_S();

  // DO MULTIPLY St A S and  St b
  // we dont send
  Rn b(problem.rhs.size());
  this->precond(b);

  // SOLVE THE LINEAR SYSTEM
  problem.solve(problem.mat, b);

  // RECONSTRUCT THE SOLUTION
  this->reconstruct(b);
}

void Extension::precond(Rn& b) {
  std::map<std::pair<int,int>,double>& A(problem.mat);
  Rn& rhs(problem.rhs);
  assert(b.size() == rhs.size());
  int t0 = MPIcf::Wtime();

  int N = problem.nDoF;
  std::map<std::pair<int,int>,double> R;
  multiply(N, St, A, R);
  A.clear();
  multiply(N, R, S, A);
  multiply(N, N, St, rhs, b);
  // rhs = b;

  std::cout << " initial size \t" << rhs.size() << std::endl;
  // need to remove the bad dof
  removeDF(N, A, b);
  std::cout << " size after \t" << b.size() << std::endl;

  std::cout << " time precond \t" << MPIcf::Wtime() - t0 << std::endl;
}
void Extension::removeDF( int N, std::map<std::pair<int,int>,double>& A, Rn& b){

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
void Extension::reconstruct(Rn& b){
  int N = problem.nDoF;  // original size
  int M = b.size();      // after removed bad dof
  Rn tmp(N);
  multiply(N,M, Pt, b, tmp);  // get back the original size
  multiply(N,N, S, tmp, problem.rhs); // multiply S * b
}
