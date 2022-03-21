


template<typename M>
void BaseFEM<M>::addToMatrix(const ItemVF<Rd::d>& VFi, const FElement& FKu, const FElement& FKv, const RNMK_& fu, const RNMK_& fv, double Cint) {

  for(int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
    for(int j = FKu.dfcbegin(VFi.cu); j < FKu.dfcend(VFi.cu); ++j) {
      this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j)) +=  Cint * VFi.c * fv(i,VFi.cv,VFi.dv)*fu(j,VFi.cu,VFi.du);
    }
  }
}

template<typename M>
void BaseFEM<M>::addToRHS(const ItemVF<Rd::d>& VFi, const FElement& FKv, const RNMK_& fv, double Cint) {
  for(int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
      (*this)(FKv.loc2glb(i)) =1;//+=   Cint * VFi.c * fv(i,VFi.cv,VFi.dv);
  }
}

//Fun_h are always evaluated on backMesh


// INTEGRATION ON FULL ELEMENT
template<typename M>
void BaseCutFEM<M>::addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh& Th) {
  assert(!VF.isRHS());
  for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {

    if(Th.isCut(k))  BaseCutFEM<M>::addElementContribution(VF, k);
    else             BaseFEM<M>::addElementContribution(VF, k);

    this->addLocalContribution();
  }
}

template<typename M>
void BaseCutFEM<M>::addLinear(const ListItemVF<Rd::d>& VF, const CutMesh& Th) {
  assert(VF.isRHS());
  for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {

    if(Th.isCut(k))  BaseCutFEM<M>::addElementContribution(VF, k);
    else             BaseFEM<M>::addElementContribution(VF, k);

    this->addLocalContribution();
  }
}


template<typename M>
void BaseFEM<M>::addElementContribution(const ListItemVF<Rd::d>& VF, const int k) {

  // CHECK IF IT IS FOR RHS OR MATRIX
  bool to_rhs = VF.isRHS();

  // Compute parameter coonected to the mesh.
  // on can take the one from the first test function
  const FESpace& Vh(*VF[0].fespaceV);
  const FElement& FK(Vh[k]);
  const Element &  K(FK.T);
  double meas = K.get_measure();
  double h    = K.get_h();
  int domain  = FK.get_domain();
  int kb      = Vh.idxElementInBackMesh(k);

  // GET THE QUADRATURE RULE
  const QF& qf(this->get_quadrature_formular_K());


  // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
  for(int l=0; l<VF.size();++l) {
    if(!VF[l].on(domain)) continue;

    // FINTE ELEMENT SPACES && ELEMENTS
    const FESpace& Vhv(VF.get_spaceV(l));
    const FESpace& Vhu(VF.get_spaceU(l));
    const FElement& FKv(Vhv[k]);
    const FElement& FKu(Vhu[k]);
    this->initIndex(FKu, FKv);

    // BF MEMORY MANAGEMENT -
    bool same = (&Vhu == &Vhv);
    int lastop = getLastop(VF[l].du, VF[l].dv);
    RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    RNMK_ fu(this->databf_+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    // COMPUTE COEFFICIENT
    double coef  = VF[l].computeCoefElement(h,meas,meas,meas,domain);

    // LOOP OVER QUADRATURE IN SPACE
    for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
      typename QF::QuadraturePoint ip(qf[ipq]);
      const Rd mip = K.map(ip);
      const double Cint = meas * ip.getWeight();

      // EVALUATE THE BASIS FUNCTIONS
      FKv.BF(Fop,ip, fv);
      if(!same) FKu.BF(Fop, ip, fu);
      //   VF[l].applyFunNL(fu,fv);

      // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
      coef *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip);
      coef *= Cint;

      if(VF.isRHS()) this->addToRHS(VF[l], FKv, fv, coef);
      else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, coef);

    }
  }
}


template<typename M>
void BaseCutFEM<M>::addElementContribution(const ListItemVF<Rd::d>& VF, const int k, const bool extend, const R epsE) {

  // GET CUT AND COMPUTE PARAMETERS
  const FESpace& Vh(VF.get_spaceV(0));
  const CutMesh& Th(Vh.get_mesh());
  const Cut_Part<Element> cutK(Th.get_cut_part(k));
  const FElement& FK(Vh[k]);
  const Element &  K(FK.T);

  if(cutK.multi_interface()){ assert(0);} // not handled yet

  double meas_cut = cutK.get_measure();
  double meas = K.get_measure();
  double h    = K.get_h();
  int domain  = FK.get_domain();
  int kb      = Vh.idxElementInBackMesh(k);

  // ElementSignEnum domPart = cutK.what_part(domain);
  // if(extend) {
  //   assert(domain != -1); // [if for some reason gets called in standard no cut pb]
  //   int otherDom = (domain == 0); // [gets the other domain]
  //   domPart = cutK.what_part(otherDom);
  // }

  // GET THE QUADRATURE RULE
  const QF& qf(this->get_quadrature_formular_K());


  // LOOP OVER ELEMENTS IN THE CUT
  for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){

    // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
    for(int l=0; l<VF.size();++l) {
      if(!VF[l].on(domain)) continue;

      // FINTE ELEMENT SPACES && ELEMENTS
      const FESpace& Vhv(VF.get_spaceV(l));
      const FESpace& Vhu(VF.get_spaceU(l));
      const FElement& FKv(Vhv[k]);
      const FElement& FKu(Vhu[k]);
      this->initIndex(FKu, FKv);

      // BF MEMORY MANAGEMENT -
      bool same = (&Vhu == &Vhv);
      int lastop = getLastop(VF[l].du, VF[l].dv);
      RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
      RNMK_ fu(this->databf_+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
      What_d Fop = Fwhatd(lastop);

      // COMPUTE COEFFICIENT
      double coef  = VF[l].computeCoefElement(h,meas_cut,meas, meas_cut, domain);


      // LOOP OVER QUADRATURE IN SPACE
      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        typename QF::QuadraturePoint ip(qf[ipq]);
        Rd mip = cutK.toPhysicalElement(it, ip);     // to the physical cut part
        Rd cut_ip = K.toReferenceElement(mip); // back to the cut part in reference element
        double Cint = meas_cut * ip.getWeight();

        // EVALUATE THE BASIS FUNCTIONS
        FKv.BF(Fop,cut_ip, fv);
        if(!same) FKu.BF(Fop, cut_ip, fu);
        //   VF[l].applyFunNL(fu,fv);

        // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
        coef *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip);
        coef *= Cint;

        if(VF.isRHS()) this->addToRHS(VF[l], FKv, fv, coef);
        else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, coef);


      }
    }
  }

}




// INTEGRATION ON INNER FACE
template<typename M>
void BaseCutFEM<M>::addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh& Th, const CHyperFace& b) {
  assert(!VF.isRHS());
  for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {
    for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces

      int jfac = ifac;
      int kn = Th.ElementAdj(k, jfac);
      // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
      if(kn == -1 || kn < k) continue;

      std::pair<int,int> e1 = std::make_pair(k,ifac);
      std::pair<int,int> e2 = std::make_pair(kn,jfac);
      // CHECK IF IT IS A CUT EDGE
      if(Th.isCutFace(k, ifac, 0)) BaseCutFEM<M>::addEdgeContribution(VF, e1, e2);
      else BaseFEM<M>::addEdgeContribution(VF, e1, e2);
    }
    this->addLocalContribution();
  }
}



template<typename M>
void BaseFEM<M>::addEdgeContribution(const ListItemVF<Rd::d>& VF, const std::pair<int,int>& e1, const std::pair<int,int>& e2) {

  // typedef typename QFB::QuadraturePoint QuadraturePoint;
  // typedef typename FElement::RdHatBord RdHatBord;
  // typedef typename Mesh::Element Element;
  //
  // assert(VF[0].fespaceU && VF[0].fespaceV);
  // const FESpace& Vhu = *VF[0].fespaceU;
  // const FESpace& Vhv = *VF[0].fespaceV;
  // const FElement& FK(Vhu[k]);
  // const Element & K = FK.T;                  // the triangle
  // int the_domain = FK.whichDomain();
  //
  // int ifacn = ifac;
  // int kn = Vhu.getNeighborElement(k, ifacn, the_domain);
  // if(kn == -1) return;         // border edge
  // if(k > kn) return;           // only compute one time
  //
  // const Rd normal = FK.T.N(ifac);
  // const R meas    = FK.T.mesureBord(ifac);
  // const R measK   = FK.getMeasure();
  // const R h       = FK.T.lenEdge(ifac);//meas;
  //
  // const FElement& FKn(Vhu[kn]);
  // const R meas_macroK  = measK + FKn.getMeasure();
  //
  // for(int l=0; l<VF.size();++l) {
  //
  //   assert(VF[l].fespaceU == VF[l].fespaceV);
  //   int lastop = getLastop(VF[l].du, VF[l].dv);
  //   double coef = BaseProblem<M>::computeCoefEdge(VF[l],h,meas,meas_macroK,measK,the_domain) * VF[l].getCoef(normal) ;
  //
  //   const int ku = (VF[l].domu == 0)? k : kn;
  //   const int kv = (VF[l].domv == 0)? k : kn;
  //   bool same = (ku == kv);
  //
  //   const FElement& FKu(Vhu[ku]);
  //   const FElement& FKv(Vhv[kv]);
  //   this->initIndex(FKu, FKv);
  //   int kuback = Vh->idxElementInBackMesh(ku);
  //   int kvback = Vh->idxElementInBackMesh(kv);
  //
  //   RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
  //   RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
  //   What_d Fop = Fwhatd(lastop);
  //
  //   for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
  //
  //     QuadraturePoint ip(qfb[ipq]); // integration point
  //     const Rd mip = K(K.toKref((RdHatBord)ip, ifac));
  //     const R Cint = meas * ip.getWeight();
  //
  //     FKu.BF(Fop,FKu.T.toKref(mip), fu); // need point in local reference element
  //     if(!same)FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element
  //     double val = VF[l].fx_backMesh_U(kuback,the_domain, mip, normal)
  //                 *VF[l].fx_backMesh_V(kvback,the_domain, mip, normal);
  //     double Cst = Cint * VF[l].c  * val * coef;
  //
  //
  //     for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
  //       for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
  //         this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j)) += Cst* fu(j,VF[l].cu,VF[l].du)*fv(i,VF[l].cv,VF[l].dv);
  //       }
  //     }
  //   }
  // }
  // this->resetIndex();
  // this->addLocalContribution();
}

template<typename M>
void BaseCutFEM<M>::addEdgeContribution(const ListItemVF<Rd::d>& VF, const std::pair<int,int>& e1, const std::pair<int,int>& e2)  {
//   typedef typename Mesh::Partition Partition;
//   typedef typename QFB::QuadraturePoint QuadraturePoint;
//   typedef typename FElement::RdHatBord RdHatBord;
//   typedef typename Mesh::Element Element;
//   assert(VF[0].fespaceU && VF[0].fespaceV);
//   const FESpace& Vhu = *VF[0].fespaceU;
//   const FESpace& Vhv = *VF[0].fespaceV;
//   const FElement& FK(Vhu[k]);
//   const int kb = Vhu.idxElementInBackMesh(k);
//   CutData cutData(Vh->getInterface(0).getCutData(kb));
//
//   if(!Vhu.isCutSpace() || cutData.face_isnot_cut(ifac)){
//   // if(cutData.edge_isnot_cut(ifac)){
//     BaseProblem<M>::addElementMatEdge(VF, k, ifac);
//     return;
//   }
//   double ss = 0.;
//   // for(int the_domain=0;the_domain<2;++the_domain){
//   int the_domain = FK.whichDomain();
//   int ifacn = ifac;
//   int kn = Vhu.getNeighborElement(k, ifacn, the_domain);
//
//
//   if(kn == -1) return;         // border edge
//   if(k > kn) return;           // only compute one time
//
//   const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
//
//   Rd normal = FK.T.N(ifac);
//   const R meas_cut = cutK.mesure(the_domain);
//   // const R measK = cutK.mesure(the_domain);
//
//   const R h = FK.T.lenEdge(ifac);
//   const R meas = cutK.mesureEdge(ifac, the_domain);
//   ss += meas;
//
//   const FElement& FKn(Vhu[kn]);
//   const int kbn = Vhu.idxElementInBackMesh(kn);
//   CutData cutDatan(Vh->getInterface(0).getCutData(kbn));
//   const Partition& cutKn =  Partition(FKn.T, cutDatan);  // build the cut
//
//   double measK = cutK.mesure(the_domain) + cutKn.mesure(the_domain);
//
//
//   for(int l=0; l<VF.size();++l) {
//
//     assert(VF[l].fespaceU == VF[l].fespaceV);
//     int lastop = getLastop(VF[l].du, VF[l].dv);
//
//     // R coef = BaseProblem<M>::computeCoef(VF[l],h,meas, measK,the_domain) * VF[l].getCoef(normal);
//     R coef = BaseCutProblem<M>::computeCoefEdge(VF[l],h,meas, measK,meas_cut, the_domain) * VF[l].getCoef(normal);
//
//     const int ku = (VF[l].domu == 0)? k : kn;
//     const int kv = (VF[l].domv == 0)? k : kn;
//     const FElement& FKu(Vhu[ku]);
//     const FElement& FKv(Vhv[kv]);
//     this->initIndex(FKu, FKv);
//     int kuback = Vh->idxElementInBackMesh(ku);
//     int kvback = Vh->idxElementInBackMesh(kv);
//
//     bool same = (ku == kv );
//     What_d Fop = Fwhatd(lastop);
//     RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
//     RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
//
//
//     for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
//
//       QuadraturePoint ip(qfb[ipq]); // integration point
//       const Rd mip = cutK.toEdge(ifac, (RdHatBord)ip, the_domain); // mip is in the global edge
//       const R Cint = meas * ip.getWeight();
//
//       FKu.BF(Fop,FKu.T.toKref(mip), fu); // need point in local reference element
//       if(!same)FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element
//
//       double val = VF[l].fx_backMesh_U(kuback,the_domain, mip, normal)
//                   *VF[l].fx_backMesh_V(kvback,the_domain, mip, normal);
//       double Cst = Cint * VF[l].c  * val * coef;
//
//       for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
//         for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
//           this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j))  +=  Cst * fv(i,VF[l].cv,VF[l].dv) * fu(j,VF[l].cu,VF[l].du);
//         }
//       }
//     }
//   }
// // }
//   this->resetIndex();
//   this->addLocalContribution();

}
