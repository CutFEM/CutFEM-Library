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
      (*this)(FKv.loc2glb(i)) += Cint * fv(i,VFi.cv,VFi.dv);
  }
}

//Fun_h are always evaluated on backMesh
template<typename M>
void BaseFEM<M>::addDiagonal(const FESpace& Qh, double val){
  int idxBegin = this->mapIdx0_[&Qh];
  int idxEnd   = Qh.nbDoF;
  for(int i=idxBegin;i<idxEnd ;++i){
    (*this->pmat_)[std::make_pair(i,i)] += val;
  }
}
template<typename M>
void BaseFEM<M>::setDiagonal(const FESpace& Qh, double val){
  int idxBegin = this->mapIdx0_[&Qh];
  int idxEnd   = Qh.nbDoF;
  for(int i=idxBegin;i<idxEnd ;++i){
    (*this->pmat_)[std::make_pair(i,i)] = val;
  }
}



// INTEGRATION ON FULL ELEMENT
template<typename Mesh>
void BaseFEM<Mesh>::addBilinear(const ListItemVF<Rd::d>& VF, const Mesh& Th) {
  assert(!VF.isRHS());
  for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {

    BaseFEM<Mesh>::addElementContribution(VF, k);

    this->addLocalContribution();
  }
}
template<typename Mesh>
void BaseFEM<Mesh>::addLinear(const ListItemVF<Rd::d>& VF, const Mesh& Th) {
  assert(VF.isRHS());
  for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {

    BaseFEM<Mesh>::addElementContribution(VF, k);
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

      // Cint = coef;
      if(VF.isRHS()) this->addToRHS(VF[l], FKv, fv, Cint);
      else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);

    }
  }
}



// INTEGRATION ON INNER FACE
template<typename Mesh>
void BaseFEM<Mesh>::addBilinear(const ListItemVF<Rd::d>& VF, const Mesh& Th, const CFacet& b) {
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
      BaseFEM<Mesh>::addFaceContribution(VF, e1, e2);
    }
    this->addLocalContribution();
  }
}
template<typename Mesh>
void BaseFEM<Mesh>::addLinear(const ListItemVF<Rd::d>& VF, const Mesh& Th, const CFacet& b) {
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
      BaseFEM<Mesh>::addFaceContribution(VF, e1, e2);
    }
  }
}
template<typename M>
void BaseFEM<M>::addFaceContribution(const ListItemVF<Rd::d>& VF, const std::pair<int,int>& e1, const std::pair<int,int>& e2) {

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

    const FElement& FKu(Vhu[ku]);
    const FElement& FKv(Vhv[kv]);
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
      FKv.BF(Fop,FKv.T.mapToReferenceElement(mip), fv);
      if(!same) FKu.BF(Fop,FKu.T.mapToReferenceElement(mip), fu);

      Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kbu,kbv), std::make_pair(domain,domain), mip, normal);
      Cint *= coef * VF[l].c;
      if(VF.isRHS()) this->addToRHS(VF[l], FKv, fv, Cint);
      else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);

    }
  }
}


// INTEGRATION ON BOUNDARY
template<typename Mesh>
void BaseFEM<Mesh>::addBilinear(const ListItemVF<Rd::d>& VF, const Mesh& Th, const CBorder& b) {
  assert(!VF.isRHS());
  for(int idx_be=Th.first_boundary_element(); idx_be<Th.last_boundary_element(); idx_be+= Th.next_boundary_element()) {

    int ifac;
    const int kb = Th.BoundaryElement(idx_be, ifac);
    const Element & K(Th[kb]);
    const BorderElement & BE(Th.be(idx_be));

    BaseFEM<Mesh>::addBorderContribution(VF, K, BE, ifac);
  }
  this->addLocalContribution();

}
template<typename Mesh>
void BaseFEM<Mesh>::addLinear(const ListItemVF<Rd::d>& VF, const Mesh& Th, const CBorder& b) {
  assert(VF.isRHS());
  for(int idx_be=Th.first_boundary_element(); idx_be<Th.last_boundary_element(); idx_be+= Th.next_boundary_element()) {

    int ifac;
    const int kb = Th.BoundaryElement(idx_be, ifac);

    const Element & K(Th[kb]);
    const BorderElement & BE(Th.be(idx_be));


    // CHECK IF IT IS A CUT EDGE
    BaseFEM<Mesh>::addBorderContribution(VF, K, BE, ifac);
  }
}
template<typename M>
void BaseFEM<M>::addBorderContribution(const ListItemVF<Rd::d>& VF, const Element& K,const BorderElement& BE, int ifac) {

  typedef typename FElement::RdHatBord RdHatBord;

  // Compute parameter connected to the mesh.
  double measK= K.measure();
  double meas = K.mesureBord(ifac);
  double h    = K.get_h();
  Rd normal   = K.N(ifac);


  // U and V HAS TO BE ON THE SAME MESH
  const FESpace& Vh(VF.get_spaceV(0));
  int kb = Vh.Th(K);
  std::vector<int> idxK = Vh.idxAllElementFromBackMesh(kb,-1);
  assert(idxK.size() == 1);
  // if(idxK.size() != 1) return;
  int k = VF[0].onWhatElementIsTestFunction (idxK);

  // GET THE QUADRATURE RULE
  const QFB& qfb(this->get_quadrature_formular_dK());

  // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
  for(int l=0; l<VF.size();++l) {
    // if(!VF[l].on(domain)) continue;

    // FINTE ELEMENT SPACES && ELEMENTS
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

    // COMPUTE COEFFICIENT
    double coef = VF[l].computeCoefElement(h,meas,measK,measK,domain) ;
    coef *= VF[l].computeCoefFromNormal(normal);

    // LOOP OVER QUADRATURE IN SPACE
    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
      typename QFB::QuadraturePoint ip(qfb[ipq]);
      const Rd mip = BE.mapToPhysicalElement((RdHatBord)ip);
      const Rd ip_edge = K.mapToReferenceElement(mip);
      double Cint = meas * ip.getWeight();

      // EVALUATE THE BASIS FUNCTIONS
      FKv.BF(Fop,ip_edge, fv);
      if(!same) FKu.BF(Fop, ip_edge, fu);

      Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, normal);
      Cint *= coef * VF[l].c;

      if(VF.isRHS()) this->addToRHS(VF[l], FKv, fv, Cint);
      else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);

    }
  }
}
