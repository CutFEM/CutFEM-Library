


template<typename M>
void BaseProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF, const TimeSlab& In) {

  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

    for(int itq=0;itq<qTime.n;++itq) {
      GTime::iter_step_quadrature = itq;               // define which time quad point
      addElementMat(VF, k, In);
    }
    this->addLocalContribution();
  }
}

template<typename M>
void BaseProblem<M>::addElementMat(const ListItemVF<Rd::d>& VF, const int k, const TimeSlab& In) {

  typedef typename QF::QuadraturePoint QuadraturePoint;

  const FElement& FK((*Vh)[k]);
  const R meas = FK.getMeasure();
  const R h = FK.T.lenEdge(0);
  GQuadraturePoint<R1> tq((qTime)[GTime::iter_step_quadrature]);

  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq.x, basisFunTime);                  // compute time basic funtions

  for(int l=0; l<VF.size();++l) {
    bool same = (VF[l].fespaceU == VF[l].fespaceV);
    int lastop = getLastop(VF[l].du, VF[l].dv);

    const int ku = VF[l].fespaceU->idxElementFromBackMesh(k, VF[l].domu);
    const int kv = VF[l].fespaceV->idxElementFromBackMesh(k, VF[l].domv);
    const FElement& FKu((*VF[l].fespaceU)[ku]);
    const FElement& FKv((*VF[l].fespaceV)[kv]);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    R coef = computeCoef(VF[l],h,meas, meas);
    this->initIndex(FKu, FKv);

    for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qf[ipq]); // integration point
      const Rd mip = FK.map(ip);
      const R Cint = meas * ip.getWeight() * In.T.mesure()* tq.a;

      FKu.BF(Fop, ip, fu); // need point in local reference element
      if(!same) FKv.BF(Fop,ip, fv);

      double val = VF[l].fxu(ku, mip, In.map(tq)) ;

      for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
        for(int jt=In.dfcbegin(0); jt<In.dfcend(0); ++jt) {
          const R tval = basisFunTime(it,0,VF[l].dtv)*basisFunTime(jt,0,VF[l].dtu);
          for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
            for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
              // (*this)(FKv.loc2glb(i,it),FKu.loc2glb(j,jt)) +=  Cint * coef * tval * val * VF[l].c * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du) ;
              this->addToLocalContribution(FKv.loc2glb(i,it),FKu.loc2glb(j,jt))+=  Cint * coef * tval * val * VF[l].c * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du) ;
            }
          }
        }
      }
    }
  }
}



template<typename M>
void BaseProblem<M>::addLinear(const ListItemVF<Rd::d>& VF, const TimeSlab& In) {

  for(int itq=0;itq<qTime.n;++itq) {
    GTime::iter_step_quadrature = itq;               // define which time quad point

    for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

      addElementRHS(VF, k, In);

    }
  }
}

template<typename M>
void BaseProblem<M>::addElementRHS(const ListItemVF<Rd::d>& VF, const int k, const TimeSlab& In) {
assert(0);
}

template<typename M>
void BaseProblem<M>::addStrongBC(std::list<ExpressionFunFEM<typename typeMesh<Rd::d>::Mesh>> gh, const TimeSlab& In, list<int> label) {

  typedef typename Mesh::BorderElement BorderElement;
  bool all_label = (label.size() == 0);

  for( int ifac = Vh->first_boundary_element(); ifac < Vh->last_boundary_element(); ifac+=Vh->next_boundary_element()) {
    const BorderElement & face(Vh->Th.be(ifac)); // The face
    if(contain(label, face.lab) || all_label) {
      for(auto it = gh.begin(); it != gh.end(); ++it) {
        setElementStrongBC(ifac, In, *it);
      }
    }
  }
}

template<typename M>
void BaseProblem<M>::addStrongBC(const ExpressionVirtual& gh, const TimeSlab& In, list<int> label) {

  typedef typename Mesh::BorderElement BorderElement;
  bool all_label = (label.size() == 0);

  for( int ifac = Vh->first_boundary_element(); ifac < Vh->last_boundary_element(); ifac+=Vh->next_boundary_element()) {
    const BorderElement & face(Vh->Th.be(ifac)); // The face
    if(contain(label, face.lab) || all_label) {
      setElementStrongBC(ifac, In, gh);
    }
  }
}


template<typename M>
void BaseProblem<M>::setElementStrongBC(int ifac, const TimeSlab& In, const ExpressionVirtual& gh) {

  typedef typename Mesh::BorderElement BorderElement;
  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;


  const BorderElement & BE(Vh->Th.be(ifac)); // The face
  int ifaceK; // index of face of triangle corresp to edge (0,1,2)
  const int kb = Vh->Th.BoundaryElement(ifac, ifaceK); // index of element (triangle), ifaceK gets modified inside
  const int k = Vh->idxElementFromBackMesh(kb, 0);
  const FElement& FK((*Vh)[k]);


  for(int i = FK.dfcbegin(gh.cu); i < FK.dfcend(gh.cu); ++i) {

    // std::cout << " dof id local " << i << std::endl;
    // std::cout << " df on What " << FK.DFOnWhat(i) << std::endl;
    // std::cout << " id face loc" << ifaceK << std::endl;

    if(Element::onWhatBorder[ifaceK][FK.DFOnWhat(i)]) {
      for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
        long df = FK.loc2glb(i,it);

        (*this)(df,df) = tgv         ;
        (*this)(df)    = tgv * gh(df);
      }
    }
  }
}




template<typename M>
void BaseProblem<M>::addBilinearFormBorder(const ListItemVF<Rd::d>& VF, const TimeSlab& In, list<int> label) {
  typedef typename Mesh::BorderElement BorderElement;
  bool all_label = (label.size() == 0);

  for( int ifac = Vh->first_boundary_element(); ifac < Vh->last_boundary_element(); ifac+=Vh->next_boundary_element()) {
    const BorderElement & face(Vh->Th.be(ifac)); // The face

    if(contain(label, face.lab) || all_label) {

      for(int itq=0;itq<qTime.n;++itq) {
        GTime::iter_step_quadrature = itq;               // define which time quad point

        addElementMatBorder(VF, ifac, In);
      }
      this->addLocalContribution();
    }
  }
}

template<typename M>
void BaseProblem<M>::addElementMatBorder(const ListItemVF<Rd::d>& VF, const int ifac, const TimeSlab& In) {

  typedef typename Mesh::BorderElement BorderElement;
  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;


  const BorderElement & BE(Vh->Th.be(ifac)); // The face
  int ifaceK; // index of face of triangle corresp to edge (0,1,2)
  const int kb = Vh->Th.BoundaryElement(ifac, ifaceK); // index of element (triangle), ifaceK gets modified inside
  const int k = Vh->idxElementFromBackMesh(kb, 0);

  const FElement& FK((*Vh)[k]);
  Rd normal = FK.T.N(ifaceK);
  const R meas = BE.mesure();
  const R measK = FK.getMeasure();
  const R h = FK.T.lenEdge(0);

  int nb_face_onB = 0;
  for(int i=0;i<Element::nea;++i){
    int ib = i;
    if(Vh->Th.ElementAdj(kb,ib) == -1) nb_face_onB += 1;
  }
  assert(nb_face_onB > 0);
  double measOnB = nb_face_onB*meas;

  const int itq = GTime::iter_step_quadrature;
  GQuadraturePoint<R1> tq((qTime)[itq]);
  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq.x, basisFunTime);

  for(int l=0; l<VF.size();++l) {
    bool same = (VF[l].fespaceU == VF[l].fespaceV);
    int lastop = getLastop(VF[l].du, VF[l].dv);
    What_d Fop = Fwhatd(lastop);

    const int ku = VF[l].fespaceU->idxElementFromBackMesh(kb,0);
    const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb,0);
    const FElement& FKu((*VF[l].fespaceU)[ku]);
    const FElement& FKv((*VF[l].fespaceV)[kv]);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction

    R coef = computeCoef(VF[l],h,measOnB,measK);
    this->initIndex(VF[l].fespaceU, VF[l].fespaceV);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = BE((RdHatBord)ip); // mip is in the global edge
      const R Cint = meas * ip.getWeight() * In.T.mesure()* tq.a;

      FKu.BF(Fop, FKu.T.toKref(mip), fu); // need point in local reference element
      if(!same) FKv.BF(Fop,FKv.T.toKref(mip), fv);

      double cst_normal = VF[l].getCoef(normal);
      double Cst = Cint * coef * VF[l].c * cst_normal;

      for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
        for(int jt=In.dfcbegin(0); jt<In.dfcend(0); ++jt) {
          const R tval = basisFunTime(it,0,VF[l].dtv)*basisFunTime(jt,0,VF[l].dtu);
          for(int i = FK.dfcbegin(VF[l].cv); i < FK.dfcend(VF[l].cv); ++i) {
            for(int j = FK.dfcbegin(VF[l].cu); j < FK.dfcend(VF[l].cu); ++j) {
              // (*this)(FK.loc2glb(i,it),FK.loc2glb(j,jt)) +=  Cst* tval * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du);
              this->addToLocalContribution(FKv.loc2glb(i,it),FKu.loc2glb(j,jt)) += Cst* tval * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du);
            }
          }
        }
      }
    }
  }
}


template<typename M>
void BaseProblem<M>::addLinearFormBorder(const ListItemVF<Rd::d>& VF, const TimeSlab& In, list<int> label) {
  typedef typename Mesh::BorderElement BorderElement;
  bool all_label = (label.size() == 0);


  for(int itq=0;itq<qTime.n;++itq) {
    GTime::iter_step_quadrature = itq;               // define which time quad point

    for( int ifac = Vh->first_boundary_element(); ifac < Vh->last_boundary_element(); ifac+=Vh->next_boundary_element()) {
      const BorderElement & face(Vh->Th.be(ifac)); // The face
      if(contain(label, face.lab) || all_label) {
        addElementRHSBorder(VF, ifac, In);
      }
    }
  }
}


template<typename M>
void BaseProblem<M>::addElementRHSBorder(const ListItemVF<Rd::d>& VF, const int ifac, const TimeSlab& In) {

  typedef typename Mesh::BorderElement BorderElement;
  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;

  int domain = 0; // ok for normal fem
  const BorderElement & BE(Vh->Th.be(ifac)); // The face
  int ifaceK; // index of face of triangle corresp to edge (0,1,2)
  const int kb = Vh->Th.BoundaryElement(ifac, ifaceK); // index of element (triangle), ifaceK gets modified inside
  const int k = Vh->idxElementFromBackMesh(kb, domain); // index in Omega0

  GQuadraturePoint<R1> tq((qTime)[GTime::iter_step_quadrature]);
  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq.x, basisFunTime);                  // compute time basic funtions

  const FElement& FK((*Vh)[k]);
  Rd normal = FK.T.N(ifaceK);
  const R meas = BE.mesure();
  const R measK = FK.getMeasure();
  const R h = FK.T.lenEdge(0);

  int nb_face_onB = 0;
  for(int i=0;i<Element::nea;++i){
    int ib = i;
    if(Vh->Th.ElementAdj(kb,ib) == -1) nb_face_onB += 1;
  }
  // std::cout << nb_face_onB << std::endl;
  assert(nb_face_onB > 0);
  double measOnB = nb_face_onB*meas;


  for(int l=0; l<VF.size();++l) {
    int lastop = getLastop(0, VF[l].dv);
    const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb, 0);
    const FElement& FKv((*VF[l].fespaceV)[kv]);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    R coef = computeCoef(VF[l],h,measOnB,measK);
    this->initIndex(VF[l].fespaceV);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = BE((RdHatBord)ip); // mip is in the global edge
      const R Cint = meas * ip.getWeight() * In.T.mesure()* tq.a;

      FKv.BF(Fop,FK.T.toKref(mip), fv);

      double cst_normal = VF[l].getCoef(normal) ;
      double val_fh = VF[l].fxu_backMesh(kb, domain, mip, In.map(tq));
      double Cst = Cint * coef * cst_normal * val_fh * VF[l].c ;

      for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
        const R tval = basisFunTime(it,0,VF[l].dtv);

        for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
          (*this)(FKv.loc2glb(i,it)) +=   Cst * tval * fv(i,VF[l].cv,VF[l].dv);
        }
      }
    }
  }
}



template<typename M>
void BaseProblem<M>::addEdgeIntegral(const ListItemVF<Rd::d>& VF, const TimeSlab& In) {
  typedef typename Mesh::Element Element;
  const FESpace& Sh =(VF[0].fespaceU)? *VF[0].fespaceU : *Vh;

  for(int k=Sh.first_element(); k<Sh.last_element(); k+= Sh.next_element()) {
    for(int itq=0;itq<qTime.n;++itq) {
      GTime::iter_step_quadrature = itq;               // define which time quad point

      for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces

        addElementMatEdge(VF, k, ifac, In);

      }
    }
    this->addLocalContribution();
  }
}

template<typename M>
void BaseProblem<M>::addElementMatEdge(const ListItemVF<Rd::d>& VF, const int k, const int ifac, const TimeSlab& In) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;

  assert(VF[0].fespaceU && VF[0].fespaceV);
  const FESpace& Vhu = *VF[0].fespaceU;
  const FESpace& Vhv = *VF[0].fespaceV;

  GQuadraturePoint<R1> tq((qTime)[GTime::iter_step_quadrature]);
  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq.x, basisFunTime);

  const FElement& FK(Vhu[k]);
  const Element & K = FK.T;                  // the triangle
  int the_domain = FK.whichDomain();

  int ifacn = ifac;
  int kn = Vhu.getNeighborElement(k, ifacn, the_domain);
  if(kn == -1) return;
  if(k > kn) return;      // only compute one time

  const Rd normal = FK.T.N(ifac);
  const R meas    = FK.T.mesureBord(ifac);
  const R measK   = FK.getMeasure();
  const R h       = FK.T.lenEdge(ifac);

  for(int l=0; l<VF.size();++l) {

    assert(VF[l].fespaceU == VF[l].fespaceV);
    int lastop = getLastop(VF[l].du, VF[l].dv);
    double coef = BaseProblem<M>::computeCoef(VF[l],h,meas,measK,the_domain) * VF[l].getCoef(normal) ;

    const int ku = (VF[l].domu == 0)? k : kn;
    const int kv = (VF[l].domv == 0)? k : kn;
    bool same = (ku == kv);

    const FElement& FKu(Vhu[ku]);
    const FElement& FKv(Vhv[kv]);
    this->initIndex(FKu, FKv);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = K(K.toKref((RdHatBord)ip, ifac));
      const R Cint = meas * ip.getWeight() * In.T.mesure()* tq.a;

      FKu.BF(Fop,FKu.T.toKref(mip), fu);//basisFunMat); // need point in local reference element
      if(!same) FKv.BF(Fop,FKv.T.toKref(mip), fv);

      for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
        for(int jt=In.dfcbegin(0); jt<In.dfcend(0); ++jt) {
          const R tval = basisFunTime(it,0,VF[l].dtv)*basisFunTime(jt,0,VF[l].dtu);
          for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
            for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
              // (*this)(FKv.loc2glb(i,it),FKu.loc2glb(j,jt))  += Cint * coef * tval * VF[l].c * fu(j,VF[l].cu,VF[l].du)*fv(i,VF[l].cv,VF[l].dv);
              this->addToLocalContribution(FKv.loc2glb(i,it),FKu.loc2glb(j,jt))  += Cint * coef * tval * VF[l].c * fu(j,VF[l].cu,VF[l].du)*fv(i,VF[l].cv,VF[l].dv);
            }
          }
        }
      }
    }
  }
}

// template<typename M>
// void BaseProblem<M>::addElementMatEdge(const ListItemVF<Rd::d>& VF, const int k, const TimeSlab& In) {
//
//   typedef typename QFB::QuadraturePoint QuadraturePoint;
//   typedef typename FElement::RdHatBord RdHatBord;
//   typedef typename Mesh::Element Element;
//
//   const FESpace& Sh =(VF[0].fespaceU)? *VF[0].fespaceU : *Vh;
//   this->initIndex(VF[0].fespaceU, VF[0].fespaceV);
//   GQuadraturePoint<R1> tq((qTime)[GTime::iter_step_quadrature]);
//
//   KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
//   In.BF(tq.x, basisFunTime);
//
//   const FElement& FK(Sh[k]);
//   const Element & K = FK.T;                  // the triangle
//   const R measK = FK.getMeasure();
//
//   for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces
//     int ifacn = ifac;
//     int kn = Sh.Th.ElementAdj(k,ifacn);
//     if(kn == -1) continue;
//     if(k > kn) continue;      // only compute one time
//
//     const R h = FK.T.lenEdge(ifac);
//
//     const FElement & FKn(Sh[kn]);                     // the neighboor finite element
//     Rd normal = (k < kn)? FK.T.N(ifac) : FKn.T.N(ifacn);
//
//     const R meas = FK.T.mesureBord( ifac);
//
//     for(int l=0; l<VF.size();++l) {
//
//       assert(VF[l].fespaceU == VF[l].fespaceV);
//       int lastop = getLastop(0, VF[l].dv);
//
//       double coef = computeCoef(VF[l],h,meas, measK) * VF[l].getCoef(normal);
//
//       const int ku = (VF[l].domu == 0)? k : kn;
//       const int kv = (VF[l].domv == 0)? k : kn;
//       bool same = (ku == kv);
//       const FElement& FKu(Sh[ku]);
//       const FElement& FKv(Sh[kv]);
//       this->initIndex(FKu, FKv);
//
//       RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
//       RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
//       What_d Fop = Fwhatd(lastop);
//
//       for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
//
//         QuadraturePoint ip(qfb[ipq]); // integration point
//         const Rd mip = K(K.toKref((RdHatBord)ip, ifac));
//         const R Cint = meas * ip.getWeight() * In.T.mesure()* tq.a;
//
//         FKu.BF(Fop,FKu.T.toKref(mip), fu);//basisFunMat); // need point in local reference element
//         if(!same) FKv.BF(Fop,FKv.T.toKref(mip), fv);
//
//         for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
//           for(int jt=In.dfcbegin(0); jt<In.dfcend(0); ++jt) {
//             const R tval = basisFunTime(it,0,VF[l].dtv)*basisFunTime(jt,0,VF[l].dtu);
//             for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
//               for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
//                 // (*this)(FKv.loc2glb(i,it),FKu.loc2glb(j,jt))  += Cint * coef * tval * VF[l].c * fu(j,VF[l].cu,VF[l].du)*fv(i,VF[l].cv,VF[l].dv);
//                 this->addToLocalContribution(FKv.loc2glb(i,it),FKu.loc2glb(j,jt))  += Cint * coef * tval * VF[l].c * fu(j,VF[l].cu,VF[l].du)*fv(i,VF[l].cv,VF[l].dv);
//
//               }
//             }
//           }
//         }
//       }
//     }
//   }
// }
