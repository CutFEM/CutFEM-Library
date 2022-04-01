
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
        // Cint = coef;
        if(VF.isRHS()) this->addToRHS(VF[l], FKv, fv, Cint);
        else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);

      }
    }
  }

}




// INTEGRATION ON INNER FACE
template<typename M>
void BaseCutFEM<M>::addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh& Th, const CFacet& b) {
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
      if(Th.isCutFace(k, ifac, 0)) BaseCutFEM<M>::addFaceContribution(VF, e1, e2);
      else BaseFEM<M>::addFaceContribution(VF, e1, e2);
    }
    this->addLocalContribution();
  }
}

template<typename M>
void BaseCutFEM<M>::addLinear(const ListItemVF<Rd::d>& VF, const CutMesh& Th, const CFacet& b) {
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

      // if(Th.isCutFace(k, ifac, 0)) std::cout << " cut edge" << std::endl;
      if(Th.isCutFace(k, ifac, 0)) BaseCutFEM<M>::addFaceContribution(VF, e1, e2);
      else BaseFEM<M>::addFaceContribution(VF, e1, e2);
    }
  }
}

template<typename M>
void BaseCutFEM<M>::addFaceContribution(const ListItemVF<Rd::d>& VF, const std::pair<int,int>& e1, const std::pair<int,int>& e2)  {

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

      // FINITE ELEMENT SPACES && ELEMENTS
      const FESpace& Vhv(VF.get_spaceV(l));
      const FESpace& Vhu(VF.get_spaceU(l));
      assert(Vhv.get_nb_element() == Vhu.get_nb_element());
      const int kv = VF[l].onWhatElementIsTestFunction (ki,kj);
      const int ku = VF[l].onWhatElementIsTrialFunction(ki,kj);

      // std::cout << VF[l].domv << std::endl;

      int kbv = Vhv.idxElementInBackMesh(kv);
      int kbu = Vhu.idxElementInBackMesh(ku);

      const FElement& FKu(Vhu[ku]);
      const FElement& FKv(Vhv[kv]);
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
        // const Rd cut_ip = Ki.mapToReferenceElement(mip);
        double Cint = meas * ip.getWeight();

        // EVALUATE THE BASIS FUNCTIONS
        FKv.BF(Fop,FKv.T.mapToReferenceElement(mip), fv);
        if(!same) FKu.BF(Fop, FKu.T.mapToReferenceElement(mip), fu);
        //   VF[l].applyFunNL(fu,fv);

        Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kbu,kbv), std::make_pair(domain,domain), mip, normal);
        Cint *= coef * VF[l].c;

        if(VF.isRHS()) this->addToRHS(   VF[l], FKv, fv, Cint);
        else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);

      }
    }
  }
}



// INTEGRATION ON BOUNDARY
template<typename M>
void BaseCutFEM<M>::addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh& cutTh, const CBorder& b) {
  assert(!VF.isRHS());
  for(int idx_be=cutTh.first_boundary_element(); idx_be<cutTh.last_boundary_element(); idx_be+= cutTh.next_boundary_element()) {

    int ifac;
    const int kb = cutTh.Th.BoundaryElement(idx_be, ifac);
    std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);

    const Element & K(cutTh.Th[kb]);
    const BorderElement & BE(cutTh.be(idx_be));

    // CHECK IF IT IS A CUT EDGE
    if(cutTh.isCutFace(idxK[0], ifac, 0)) BaseCutFEM<M>::addBorderContribution(VF, K, BE, ifac);
    else BaseFEM<M>::addBorderContribution(VF, K, BE,ifac);
  }
  this->addLocalContribution();

}

template<typename M>
void BaseCutFEM<M>::addLinear(const ListItemVF<Rd::d>& VF, const CutMesh& cutTh, const CBorder& b) {
  assert(VF.isRHS());
  for(int idx_be=cutTh.first_boundary_element(); idx_be<cutTh.last_boundary_element(); idx_be+= cutTh.next_boundary_element()) {

    int ifac;
    const int kb = cutTh.Th.BoundaryElement(idx_be, ifac);
    std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb,-1);

    const Element & K(cutTh.Th[kb]);
    const BorderElement & BE(cutTh.be(idx_be));


    // CHECK IF IT IS A CUT EDGE
    if(cutTh.isCutFace(idxK[0], ifac, 0)) BaseCutFEM<M>::addBorderContribution(VF, K, BE, ifac);
    else BaseFEM<M>::addBorderContribution(VF, K, BE, ifac);
  }
}

template<typename M>
void BaseCutFEM<M>::addBorderContribution(const ListItemVF<Rd::d>& VF, const Element& K,const BorderElement& BE, int ifac)  {

  typedef typename FElement::RdHatBord RdHatBord;

  // Compute parameter connected to the mesh.
  double measK = K.measure();
  double h     = K.get_h();
  Rd normal    = K.N(ifac);

  // U and V HAS TO BE ON THE SAME MESH
  const FESpace& Vh(VF.get_spaceV(0));
  const CutMesh& Th(Vh.get_mesh());
  int kb = Vh.Th(K);
  std::vector<int> idxK = Vh.idxAllElementFromBackMesh(kb,-1);
  assert(idxK.size() == 2);

  // GET THE QUADRATURE RULE
  const QFB& qfb(this->get_quadrature_formular_cutFace());


  for(int e=0;e<2;++e) {

    int k = idxK[e];
    typename Element::Face face;
    const Cut_Part<typename Element::Face> cutFace(Th.get_cut_face(face, k, ifac));
    const Cut_Part<Element> cutK(Th.get_cut_part(k));

    double meas_cut  = cutK.measure();

    // LOOP OVER ELEMENTS IN THE CUT
    for(auto it = cutFace.element_begin();it != cutFace.element_end(); ++it){

      double meas  = cutFace.measure(it);

      for(int l=0; l<VF.size();++l) {
        // if(!VF[l].on(domain)) continue;

        // FINITE ELEMENT SPACES && ELEMENTS
        const FESpace& Vhv(VF.get_spaceV(l));
        const FESpace& Vhu(VF.get_spaceU(l));
        assert(Vhv.get_nb_element() == Vhu.get_nb_element());
        bool same = (VF.isRHS() || (&Vhu == &Vhv));
        const FElement& FKu(Vhu[k]);
        const FElement& FKv(Vhv[k]);
        int domain = FKv.get_domain();
        this->initIndex(FKu, FKv);


        // BF MEMORY MANAGEMENT -
        int lastop = getLastop(VF[l].du, VF[l].dv);
        RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop);
        RNMK_ fu(this->databf_+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop);
        What_d Fop = Fwhatd(lastop);

        // COMPUTE COEFFICIENT && NORMAL
        double coef = VF[l].computeCoefElement(h,meas,measK,meas_cut,domain) ;
        coef *= VF[l].computeCoefFromNormal(normal);


        // LOOP OVER QUADRATURE IN SPACE
        for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
          typename QFB::QuadraturePoint ip(qfb[ipq]);
          const Rd mip = cutFace.mapToPhysicalElement(it, (RdHatBord)ip);
          const Rd cut_ip = K.mapToReferenceElement(mip);
          double Cint = meas * ip.getWeight();

          // EVALUATE THE BASIS FUNCTIONS
          FKv.BF(Fop,cut_ip, fv);
          if(!same) FKu.BF(Fop, cut_ip, fu);
          //   VF[l].applyFunNL(fu,fv);

          Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, normal);
          Cint *= coef * VF[l].c;


          if(VF.isRHS()) this->addToRHS(   VF[l], FKv, fv, Cint);
          else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);

        }
      }
    }
  }
}




// INTEGRATION ON INTERFACE
template<typename M>
void BaseCutFEM<M>::addBilinear(const ListItemVF<Rd::d>& VF, const Interface<M>& gamma,list<int> label) {
  assert(!VF.isRHS());
  bool all_label = (label.size() == 0);

  for(int iface=gamma.first_element(); iface<gamma.last_element(); iface+=gamma.next_element()) {
    const typename Interface<M>::Face& face = gamma[iface];  // the face
    if(util::contain(label, face.lab) || all_label) {

      addInterfaceContribution(VF, gamma, iface);

    }
  }

  this->addLocalContribution();
}

template<typename M>
void BaseCutFEM<M>::addLinear(const ListItemVF<Rd::d>& VF, const Interface<M>& gamma,list<int> label) {
  assert(VF.isRHS());
  bool all_label = (label.size() == 0);

  for(int iface=gamma.first_element(); iface<gamma.last_element(); iface+=gamma.next_element()) {
    const typename Interface<M>::Face& face = gamma[iface];  // the face
    if(util::contain(label, face.lab) || all_label) {

      addInterfaceContribution(VF, gamma, iface);

    }
  }
}

template<typename M>
void BaseCutFEM<M>::addInterfaceContribution(const ListItemVF<Rd::d>& VF, const Interface<M>& interface, int ifac) {
  typedef typename FElement::RdHatBord RdHatBord;

  // GET IDX ELEMENT CONTAINING FACE ON backMes
  const int kb = interface.idxElementOfFace(ifac);
  // const typename Interface<M>::Face& face = interface[ifac];
  const Element & K(interface.get_element(kb));
  double measK = K.measure();
  double h     = K.get_h() ;
  double meas  = interface.measure(ifac);

  const Rd normal(-interface.normal(ifac));

  // GET THE QUADRATURE RULE
  const QFB& qfb(this->get_quadrature_formular_cutFace());

  for(int l=0; l<VF.size();++l) {
    // if(!VF[l].on(domain)) continue;

    // FINITE ELEMENT SPACES && ELEMENTS
    const FESpace& Vhv(VF.get_spaceV(l));
    const FESpace& Vhu(VF.get_spaceU(l));
    bool same = (VF.isRHS() || (&Vhu == &Vhv));

    // NOT OPTIMAL IF DOMAIN IS KNOWN!!!!
    std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
    std::vector<int> idxU = (same)?idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());

    int kv = VF[l].onWhatElementIsTestFunction (idxV);
    int ku = VF[l].onWhatElementIsTrialFunction(idxU);
    int kbv = Vhv.idxElementInBackMesh(kv);
    int kbu = Vhu.idxElementInBackMesh(ku);

    const FElement& FKu(Vhu[ku]);
    const FElement& FKv(Vhv[kv]);
    int domu = FKu.get_domain();
    int domv = FKv.get_domain();
    this->initIndex(FKu, FKv);


    // BF MEMORY MANAGEMENT -
    int lastop = getLastop(VF[l].du, VF[l].dv);
    RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop);
    RNMK_ fu(this->databf_+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop);
    What_d Fop = Fwhatd(lastop);

    // COMPUTE COEFFICIENT && NORMAL
    double coef = VF[l].computeCoefInterface(h,meas,measK,measK,domu, domv) ;
    coef *= VF[l].computeCoefFromNormal(normal);

    // LOOP OVER QUADRATURE IN SPACE
    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = interface.mapToPhysicalFace(ifac,(RdHatBord)ip);
      const Rd face_ip = K.mapToReferenceElement(mip);
      double Cint = meas * ip.getWeight();

      // EVALUATE THE BASIS FUNCTIONS
      FKv.BF(Fop,face_ip, fv);
      if(!same) FKu.BF(Fop, face_ip, fu);


      Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kbu,kbv), std::make_pair(domu, domv), mip, normal);
      Cint *= coef * VF[l].c;

      if(VF.isRHS()) this->addToRHS(   VF[l], FKv, fv, Cint);
      else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);

    }

  }
  // getchar();
}


// FACE STABILIZATION
template<typename M>
void BaseCutFEM<M>::addFaceStabilization(const ListItemVF<Rd::d>& VF, const CutMesh& Th) {
  assert(!VF.isRHS());
  for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {
    if(!Th.isCut(k)) continue;
    for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces

      int jfac = ifac;
      int kn = Th.ElementAdj(k, jfac);
      // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
      if(kn == -1 || kn < k) continue;

      std::pair<int,int> e1 = std::make_pair(k,ifac);
      std::pair<int,int> e2 = std::make_pair(kn,jfac);
      BaseFEM<M>::addFaceContribution(VF, e1, e2);

    }
    this->addLocalContribution();
  }
}

template<typename M>
void BaseCutFEM<M>::addFaceStabilization(const ListItemVF<Rd::d>& VF, const CutMesh& Th, const MacroElement<M>& macro) {

  for(auto me=macro.macro_element.begin(); me!=macro.macro_element.end();++me) {
    for(auto it=me->second.inner_edge.begin(); it!=me->second.inner_edge.end();++it){

      int k = it->first;
      int ifac = it->second;
      int jfac = ifac;
      int kn = Th.ElementAdj(k, jfac);

      std::pair<int,int> e1 = std::make_pair(k,ifac);
      std::pair<int,int> e2 = std::make_pair(kn,jfac);

      BaseFEM<M>::addFaceContribution(VF, e1, e2);
    }
    this->addLocalContribution();
  }
}








//dd
