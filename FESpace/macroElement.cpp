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
    else if (Kf.tfe->Sub_ToFE[ic] == &DataFE<Mesh2>::RT1){
      do_extension_RT1(Ks,Kf,id_e,ic0);
    }
    else if (Kf.tfe->Sub_ToFE[ic] == &DataFE<Mesh2>::P0){
      do_extension_P0(Ks,Kf,ic0);
    }
    else if (Kf.tfe->Sub_ToFE[ic] == &DataFE<Mesh2>::P1dc){
      do_extension_P1dc(Ks,Kf,ic0);
    }
    else {
      assert(0);
    }
    ic0 += Kf.tfe->Sub_ToFE[ic]->N;
  }
}
void MacroElement::do_extension_P0   (const FElement2& Ks, const FElement2& Kf, int ic){
  int idx0 = Kf.dfcbegin(ic);   // the pressure index
  if (S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] == 0) return;
  S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
  St[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] = 0;
  S[std::make_pair(Ks.loc2glb(idx0),Kf.loc2glb(idx0))] = 1;
  St[std::make_pair(Kf.loc2glb(idx0),Ks.loc2glb(idx0))] = 1;

  dof2rm.insert(Ks.loc2glb(idx0));

}
void MacroElement::do_extension_P1dc (const FElement2& Ks, const FElement2& Kf, int ic){
  int idx0 = Kf.dfcbegin(ic);   // the pressure index
  if (S[std::make_pair(Ks.loc2glb(idx0),Ks.loc2glb(idx0))] == 0) return;

  int ndof = Kf.NbDoF();
  KNM<double> val(3, 3);
  evaluate_dofP1dc(Kf, Ks, val, ic);
  for(int id_df = Ks.dfcbegin(ic),df=0 ; id_df < Ks.dfcend(ic); ++id_df,++df) {
    S [std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
    St[std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
    for(int i = Kf.dfcbegin(ic),j=0; i < Kf.dfcend(ic); ++i,++j) {
      S[std::make_pair(Ks.loc2glb(id_df),Kf.loc2glb(i))]  = val(df,j);        // 1st component small element
      St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_df))] = val(df,j);        // 1st component small element
    }
    dof2rm.insert(Ks.loc2glb(id_df));
  }
}
void MacroElement::do_extension_RT0 (const FElement2& Ks, const FElement2& Kf, int id_e, int ic){
  int ndof = Kf.NbDoF();
  KNM<double> val(1, ndof);
  evaluate_dofRT0(Ks, id_e, Kf, val);
  int id_df = id_e;
  S [std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
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
void MacroElement::do_extension_RT1(const FElement2& Ks, const FElement2& Kf, int id_e, int ic){
  int ndofOnEdge = Kf.tfe->ndfonEdge;
  int ndof = Kf.NbDoF();
  KNM<double> val(ndofOnEdge+2, ndof); // 2 bubble bf
  evaluate_dofBDM1(Ks, id_e, Kf, val);

  for(int df = 0; df < ndofOnEdge; ++df) {
    // if(  df==1) {
      int id_df = ndofOnEdge*id_e + df;
      S [std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))]  = 0.;        // 1st component small element
      St[std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
      for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
        S[std::make_pair(Ks.loc2glb(id_df),Kf.loc2glb(i))]  = val(df,i);        // 1st component small element
        St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_df))] = val(df,i);        // 1st component small element
      }
      dof2rm.insert(Ks.loc2glb(id_df));
    // }
  }

  for(int df = 0; df < 2; ++df) {
    int id_df = 3*ndofOnEdge + df;
    S [std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
    St[std::make_pair(Ks.loc2glb(id_df),Ks.loc2glb(id_df))] = 0.;        // 1st component small element
    for(int i = Kf.dfcbegin(ic); i < Kf.dfcend(ic); ++i) {
      S[std::make_pair(Ks.loc2glb(id_df),Kf.loc2glb(i))]  = val(df+ndofOnEdge,i);        // 1st component small element
      St[std::make_pair(Kf.loc2glb(i),Ks.loc2glb(id_df))] = val(df+ndofOnEdge,i);        // 1st component small element
    }
    dof2rm.insert(Ks.loc2glb(id_df));
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
    if(handle == exhaust){
      //
      // do_exhaust

    }
  }
  std::cout << " removing \t" << dof2rm.size() << "  dof " << std::endl;
}

void MacroElement::evaluate_dofP1dc(const FElement2& FKs, const FElement2& FKf, Rnm& val, int ic) {

  val = 0.;
  KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
  What_d Fop = Fwhatd(0);

  for(int v=0;v<3;++v) {

    R2 mip_Ks   = FKs.T.at(v);
    R2 ip_KfHat = FKf.T.toKref(mip_Ks);

    FKf.BF(Fop_D0, ip_KfHat, bf);
    for(int i = FKs.dfcbegin(ic), j=0; i < FKs.dfcend(ic); ++i,++j) {
      val(v,j) += bf(i,ic,0);
    }
  }
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

  // val = 0.;
  // KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
  // What_d Fop = Fwhatd(0);
  //
  // for(int eF=0;eF<3;++eF) {
  //
  //   double meas = FKf.T.lenEdge(eF);
  //   R2 normal = FKf.T.EdgeOrientation(eF)*FKf.T.N(eF);
  //   for(int iq=0;iq<QF.getNbrOfQuads();++iq) {
  //
  //     QuadraturePoint1d ip_1d(QF[iq]);
  //     R2 ip_KsHat = FKf.T.toKref(ip_1d, eF);
  //     R2 mip_Ks   = FKf.map(ip_KsHat);
  //     R2 ip_KfHat = FKs.T.toKref(mip_Ks);
  //
  //     FKs.BF(Fop_D0, ip_KfHat, bf);
  //
  //     val(0,eF) += meas*ip_1d.getWeight()*(bf(e,0,0)*normal.x + bf(e,1,0)*normal.y) ;
  //
  //   }
  // }
}
void MacroElement::evaluate_dofBDM1(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val) {
  val = 0.;
  KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
  What_d Fop = Fwhatd(0);
  double meas = FKs.T.lenEdge(e);
  R2 normal = -FKs.T.Edge(e).perp(); // contain the mesure of edge
  // R2 normal = FKs.T.EdgeOrientation(e)*FKs.T.N(e);



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
void MacroElement::evaluate_dofRT1(const FElement2& FKs, int e, const FElement2& FKf, Rnm& val) {
  val = 0.;
  int ndf = FKf.NbDoF();
  KNMK<double> bf(FKf.NbDoF(),FKf.N,1); //  the value for basic fonction
  What_d Fop = Fwhatd(0);
  double meas = FKs.T.lenEdge(e);
  int eOrientation = FKs.T.EdgeOrientation(e);
  // R2 normal = FKs.T.EdgeOrientation(e)*FKs.T.N(e);
  R2 normal = -FKs.T.Edge(e).perp();//*eOrientation; // contain the mesure of edge

  for(int iq=0;iq<QF.getNbrOfQuads();++iq) {

    QuadraturePoint1d ip_1d(QF[iq]);
    R2 ip_KsHat = FKs.T.toKref(ip_1d, e);
    R2 mip_Ks   = FKs.map(ip_KsHat);
    R2 ip_KfHat = FKf.T.toKref(mip_Ks);

    FKf.BF(Fop_D0, ip_KfHat, bf);

    // evaluate Fat_bf in small_dof : on edge
    R l0 = QF[iq].x, l1 = 1 - QF[iq].x;
    R p0 = (2 * l0 - l1) * 2;     // poly othogonaux to \lambda_1
    R p1 = (2 * l1 - l0) * 2;     // poly othogonaux to \lambda_0
    R lambda1 = eOrientation * p0 * QF[iq].a;    // [some quadrature function?]
    R lambda0 = eOrientation * p1 * QF[iq].a;    //
    if(eOrientation < 0) {
      Exchange(lambda1, lambda0);
    }
    for(int i=0;i<ndf;++i) {
      val(0,i) += lambda0*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
      val(1,i) += lambda1*(bf(i,0,0)*normal.x + bf(i,1,0)*normal.y) ;
    }


  }

  // evaluate Fat_bf in small_bubble_dof
  R2 B[2] = {FKs.T.Edge(1), FKs.T.Edge(2)};
  B[0] = B[0].perp();
  B[1] = B[1].perp();
  for(int iq=0;iq<QFK.getNbrOfQuads();++iq) {

    QuadraturePoint2d ip_2d(QFK[iq]);
    R2 mip_Ks   = FKs.map(ip_2d);
    R2 ip_KfHat = FKf.T.toKref(mip_Ks);

    FKf.BF(Fop_D0, ip_KfHat, bf);
    double w = QFK[iq].a * 0.5;
    for(int i=0;i<ndf;++i) {
      val(2,i) += w*(bf(i,0,0)*B[0].x + bf(i,1,0)*B[0].y) ;
      val(3,i) += w*(bf(i,0,0)*B[1].x + bf(i,1,0)*B[1].y) ;
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
  std::vector<std::pair<int, int>> df2fix;
  element_edge_handle.clear();
  // First set the trivial dof (good and trivial exhaust)
  for(auto it = small_element.begin(); it!=small_element.end();++it) {
    int k = it->second.index;
    const FElement2& FK(Vh[k]);
    int the_domain = FK.whichDomain();
    int pos_k = it->second.chain_position;

    for(int e=0;e<3;++e) {
      int je = e;
      int kn = Vh.getNeighborElement(k, je, the_domain);
      if(isRootFat(kn) ) {  // triavial good
        df2fix.push_back(std::make_pair(k, 3*e+good));
        continue;
      }
      else{
        int pos_kn = (kn == -1) ? pos_k + 1 : small_element[kn].chain_position;
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
    small_K_temp[ii++] = it->second.index;
  }



  // First set the trivial dof (good and trivial exhaust)
  // for(auto it = idx_small_K.begin(); it!=idx_small_K.end();++it) {
  //   int k = it->first;
  // for (int i=idx_small_K.size()-1;i>=0;--i) {
  // for (int i=small_element.size()-1;i>=0;--i) {
    // int k = idx_small_K[i].first;
    for (auto it=small_element.begin();it!=small_element.end();++it) {

    // int k = small_element[i].index;//.first;
    int k = it->second.index;//.first;
    assert(0);
    // need to check this loop using map
    // also need to uncomment erase part a bit further
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
    // if(is_exhaust) {
    //    if(count != 2) exhaust_element.insert(k);
    //    if(count != 2) small_K_temp.erase(small_K_temp.begin()+i);
    // }
    // if(count == 2) {
    //   int df = FK.loc2glb(notGoodEdge);
    //   // df2fix.erase(df);
    //   // df2fix[df] = std::make_pair(k, 3*notGoodEdge+good);
    //   small_K_temp.erase(small_K_temp.begin()+i);
    // }
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

     // Do we care about having unique edge (choose k<kn) or we
    // we just use macro to get element and integrate on faces
  // create all root elements
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

  // for( auto it=macro_element.begin(); it!=macro_element.end(); ++it) {
  //
  //   // for(int i=0;i<it->second.inner_edge.size();++i) {
  //   //   int idxC = it->second.inner_edge[i].first;
  //   //   int ie = it->second.inner_edge[i].second;
  //   //   int k = interface.idxElementOfFace(idxC);
  //   //   int je = ie;
  //   //   int kn = Th.ElementAdj(k, je);
  //   //   assert(k < kn);
  //   // }
  //
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
