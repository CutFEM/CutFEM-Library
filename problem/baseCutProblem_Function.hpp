


template<typename M>
void BaseFEM<M>::addToMatrix(const ItemVF<Rd::d>& VFi, const FElement& FKu, const FElement& FKv, const RNMK_& fu, const RNMK_& fv, double Cint) {

  for(int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
    for(int j = FKu.dfcbegin(VFi.cu); j < FKu.dfcend(VFi.cu); ++j) {
      this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j)) +=  Cint * fv(i,VFi.cv,VFi.dv)*fu(j,VFi.cu,VFi.du);
    }
  }
}

template<typename M>
void BaseFEM<M>::addToRHS(const ItemVF<Rd::d>& VFi, const FElement& FKv, const RNMK_& fv, double Cint) {
  for(int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
      (*this)(FKv.loc2glb(i)) =1;//+=   Cint * fv(i,VFi.cv,VFi.dv);
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
  double meas = K.measure();
  double h    = K.get_h();
  int domain  = FK.get_domain();
  int kb      = Vh.idxElementInBackMesh(k);

  // GET THE QUADRATURE RULE
  const QF& qf(this->get_quadrature_formular_K());


  // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
  for(int l=0; l<VF.size();++l) {
    // if(!VF[l].on(domain)) continue;

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
      double Cint = meas * ip.getWeight();

      // EVALUATE THE BASIS FUNCTIONS
      FKv.BF(Fop,ip, fv);
      if(!same) FKu.BF(Fop, ip, fu);
      //   VF[l].applyFunNL(fu,fv);

      // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
      Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip);
      Cint *= coef * VF[l].c;

      if(VF.isRHS()) this->addToRHS(VF[l], FKv, fv, Cint);
      else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);

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

  double meas = K.measure();
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
  const QF& qf(this->get_quadrature_formular_cutK());


  // LOOP OVER ELEMENTS IN THE CUT
  for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){

    double meas_cut = cutK.measure(it);

    // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
    for(int l=0; l<VF.size();++l) {
      // if(!VF[l].on(domain)) continue;

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
        Rd mip = cutK.mapToPhysicalElement(it, ip);     // to the physical cut part
        Rd cut_ip = K.mapToReferenceElement(mip); // back to the cut part in reference element
        double Cint = meas_cut * ip.getWeight();

        // EVALUATE THE BASIS FUNCTIONS
        FKv.BF(Fop,cut_ip, fv);
        if(!same) FKu.BF(Fop, cut_ip, fu);
        //   VF[l].applyFunNL(fu,fv);

        // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
        Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip);
        Cint *= coef * VF[l].c;

        if(VF.isRHS()) this->addToRHS(VF[l], FKv, fv, Cint);
        else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);


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
void BaseCutFEM<M>::addLinear(const ListItemVF<Rd::d>& VF, const CutMesh& Th, const CHyperFace& b) {
  assert(VF.isRHS());
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
  }
}



template<typename M>
void BaseFEM<M>::addEdgeContribution(const ListItemVF<Rd::d>& VF, const std::pair<int,int>& e1, const std::pair<int,int>& e2) {

  typedef typename FElement::RdHatBord RdHatBord;

  // CHECK IF IT IS FOR RHS OR MATRIX
  // CONVENTION ki < kj
  bool to_rhs = VF.isRHS();
  int ki = e1.first, ifac = e1.second;
  int kj = e2.first, jfac = e2.second;
  // Compute parameter coonected to the mesh.
  // on can take the one from the first test function
  const FESpace& Vh(*VF[0].fespaceV);
  const FElement& FKi(Vh[ki]);
  const FElement& FKj(Vh[kj]);
  const Element & Ki(FKi.T);
  const Element & Kj(FKj.T);
  double measK= Ki.measure() + Kj.measure();
  double meas = Ki.mesureBord(ifac);
  double h    = 0.5*(Ki.get_h() + Kj.get_h());
  int domain  = FKi.get_domain();
  Rd normal   = Ki.N(ifac);

  // GET THE QUADRATURE RULE
  const QFB& qfb(this->get_quadrature_formular_dK());

  // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
  for(int l=0; l<VF.size();++l) {
    // if(!VF[l].on(domain)) continue;

    // FINTE ELEMENT SPACES && ELEMENTS
    const FESpace& Vhv(VF.get_spaceV(l));
    const FESpace& Vhu(VF.get_spaceU(l));
    assert(Vhv.get_nb_element() == Vhu.get_nb_element());
    const int kv = VF[l].onWhatElementIsTestFunction(ki,kj);
    const int ku = VF[l].onWhatElementIsTrialFunction(ki,kj);

    int kbv = Vhv.idxElementInBackMesh(kv);
    int kbu = Vhu.idxElementInBackMesh(ku);

    const FElement& FKu(Vhv[ku]);
    const FElement& FKv(Vhu[kv]);
    this->initIndex(FKu, FKv);

    // BF MEMORY MANAGEMENT -
    bool same = (VF.isRHS() || (&Vhu == &Vhv && ku == kv));
    int lastop = getLastop(VF[l].du, VF[l].dv);
    RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop);
    RNMK_ fu(this->databf_+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop);
    What_d Fop = Fwhatd(lastop);

    // COMPUTE COEFFICIENT
    double coef = VF[l].computeCoefElement(h,meas,measK,measK,domain) ;
    coef *= VF[l].computeCoefFromNormal(normal);

    // LOOP OVER QUADRATURE IN SPACE
    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
      typename QFB::QuadraturePoint ip(qfb[ipq]);
      const Rd ip_edge = Ki.mapToReferenceElement((RdHatBord)ip, ifac);
      const Rd mip = Ki.mapToPhysicalElement(ip_edge);
      double Cint = meas * ip.getWeight();

      // EVALUATE THE BASIS FUNCTIONS
      FKv.BF(Fop,ip_edge, fv);
      if(!same) FKu.BF(Fop, ip_edge, fu);

      Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kbu,kbv), std::make_pair(domain,domain), mip, normal);
      Cint *= coef * VF[l].c;

      if(VF.isRHS()) this->addToRHS(VF[l], FKv, fv, Cint);
      else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);

    }
  }
}

template<typename M>
void BaseCutFEM<M>::addEdgeContribution(const ListItemVF<Rd::d>& VF, const std::pair<int,int>& e1, const std::pair<int,int>& e2)  {

  typedef typename FElement::RdHatBord RdHatBord;

  // CHECK IF IT IS FOR RHS OR MATRIX
  // CONVENTION ki < kj
  bool to_rhs = VF.isRHS();
  int ki = e1.first, ifac = e1.second;
  int kj = e2.first, jfac = e2.second;
  // Compute parameter coonected to the mesh.
  // on can take the one from the first test function
  const FESpace& Vh(*VF[0].fespaceV);
  const CutMesh& Th(Vh.get_mesh());
  const FElement& FKi(Vh[ki]);
  const FElement& FKj(Vh[kj]);
  const Element & Ki(FKi.T);
  const Element & Kj(FKj.T);

  const Cut_Part<Element> cutKi(Th.get_cut_part(ki));
  const Cut_Part<Element> cutKj(Th.get_cut_part(kj));
  typename Element::Face face;
  const Cut_Part<typename Element::Face> cutFace(Th.get_cut_face(face, ki, ifac));

  double measK = cutKi.measure() + cutKj.measure();
  double measF = Ki.mesureBord(ifac);
  double h     = 0.5*(Ki.get_h() + Kj.get_h());
  int domain   = FKi.get_domain();
  Rd normal    = Ki.N(ifac);

  // GET THE QUADRATURE RULE
  const QFB& qfb(this->get_quadrature_formular_cutFace());

  // LOOP OVER ELEMENTS IN THE CUT
  for(auto it = cutFace.element_begin();it != cutFace.element_end(); ++it){

    double meas  = cutFace.measure(it);

    for(int l=0; l<VF.size();++l) {
      // if(!VF[l].on(domain)) continue;

      // FINTE ELEMENT SPACES && ELEMENTS
      const FESpace& Vhv(VF.get_spaceV(l));
      const FESpace& Vhu(VF.get_spaceU(l));
      assert(Vhv.get_nb_element() == Vhu.get_nb_element());
      const int kv = VF[l].onWhatElementIsTestFunction (ki,kj);
      const int ku = VF[l].onWhatElementIsTrialFunction(ki,kj);

      // std::cout << VF[l].domv << std::endl;

      int kbv = Vhv.idxElementInBackMesh(kv);
      int kbu = Vhu.idxElementInBackMesh(ku);

      const FElement& FKu(Vhv[ku]);
      const FElement& FKv(Vhu[kv]);
      this->initIndex(FKu, FKv);

      // BF MEMORY MANAGEMENT -
      bool same = (VF.isRHS() || (&Vhu == &Vhv && ku == kv));
      int lastop = getLastop(VF[l].du, VF[l].dv);
      RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop);
      RNMK_ fu(this->databf_+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop);
      What_d Fop = Fwhatd(lastop);

      // COMPUTE COEFFICIENT && NORMAL
      double coef = VF[l].computeCoefElement(h,meas,measK,measK,domain) ;
      coef *= VF[l].computeCoefFromNormal(normal);


      // LOOP OVER QUADRATURE IN SPACE
      for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
        typename QFB::QuadraturePoint ip(qfb[ipq]);
        const Rd mip = cutFace.mapToPhysicalElement(it, (RdHatBord)ip);
        const Rd cut_ip = Ki.mapToReferenceElement(mip);
        double Cint = meas * ip.getWeight();

        // EVALUATE THE BASIS FUNCTIONS
        FKv.BF(Fop,cut_ip, fv);
        if(!same) FKu.BF(Fop, cut_ip, fu);
        //   VF[l].applyFunNL(fu,fv);

        Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kbu,kbv), std::make_pair(domain,domain), mip, normal);
        Cint *= coef * VF[l].c;

        if(VF.isRHS()) this->addToRHS(   VF[l], FKv, fv, Cint);
        else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);

      }
    }
  }
}
