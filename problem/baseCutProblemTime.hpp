template<typename M>
void BaseCutProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF, const TimeSlab& In) {

  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

    for(int itq=0;itq<qTime.n;++itq) {
      GTime::iter_step_quadrature = itq;               // define which time quad point

      addElementMat(VF, k, In);

    }
    this->addLocalContribution();
  }
}

template<typename M>
void BaseCutProblem<M>::addElementMat(const ListItemVF<Rd::d>& VF, const int k, const TimeSlab& In) {
  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;

  const int itq = GTime::iter_step_quadrature;

  const FElement& FK((*Vh)[k]);
  const int domain = FK.whichDomain();

  const int kb = Vh->Th(FK.T);
  CutData cutData(Vh->getInterface(itq).getCutData(kb));     // get the cut data
  const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
  ElementSignEnum the_part = cutK.what_part(domain);
  const R h = FK.T.lenEdge(0);

  GQuadraturePoint<R1> tq((qTime)[itq]);
  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq.x, basisFunTime);                  // compute time basic funtions

  for(typename Partition::const_element_iterator itpart = cutK.element_begin(the_part);
  itpart != cutK.element_end(the_part); ++itpart){

    const R meas = cutK.mesure(itpart);

    for(int l=0; l<VF.size();++l) {
      if(!VF[l].on(domain)) continue;

      int lastop = getLastop(VF[l].du, VF[l].dv);

      const int ku = (VF[l].domu != -1)? VF[l].fespaceU->idxElementFromBackMesh(kb, VF[l].domu) : k;
      const int kv = (VF[l].domv != -1)? VF[l].fespaceV->idxElementFromBackMesh(kb, VF[l].domv) : k;
      const FElement& FKu((*VF[l].fespaceU)[ku]);
      const FElement& FKv((*VF[l].fespaceV)[kv]);
      bool same = (VF[l].fespaceU==VF[l].fespaceV);
      RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
      RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
      What_d Fop = Fwhatd(lastop);

      R coef = computeCoef(VF[l],domain,cutK);
      this->initIndex(FKu,FKv);

      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = cutK.toK(itpart, ip);
        Rd cut_ip = FKv.T.toKref(mip); // in Kref
        const R Cint = meas * ip.getWeight() * In.T.mesure()* tq.a;

        FKu.BF(Fop, cut_ip, fu); // need point in local reference element
        if(!same) FKv.BF(Fop,cut_ip, fv);

        double val_fh = VF[l].fxu_backMesh(kb, domain, mip, In.map(tq));

        for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
          for(int jt=In.dfcbegin(0); jt<In.dfcend(0); ++jt) {
            const R tval = basisFunTime(it,0,VF[l].dtv)*basisFunTime(jt,0,VF[l].dtu);
            for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
              for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
                // (*this)(FKv.loc2glb(i,it),FKu.loc2glb(j,jt)) +=  Cint * coef * val_fh * tval * VF[l].c * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du) ;
                this->addToLocalContribution(FKv.loc2glb(i,it),FKu.loc2glb(j,jt)) += Cint * coef * val_fh * tval * VF[l].c * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du) ;
              }
            }
          }
        }
      }
    }
  }
}

template<typename M>
void BaseCutProblem<M>::addLinear(const ListItemVF<Rd::d>& VF, const TimeSlab& In) {

  for(int itq=0;itq<qTime.n;++itq) {
    GTime::iter_step_quadrature = itq;               // define which time quad point

    for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

      addElementRHS(VF, k, In);

    }
  }
}

template<typename M>
void BaseCutProblem<M>::addElementRHS(const ListItemVF<Rd::d>& VF, const int k, const TimeSlab& In) {

  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;

  const FElement& FK((*Vh)[k]);
  const int domain = FK.whichDomain();
  const R h = FK.T.lenEdge(0);
  const int itq = GTime::iter_step_quadrature;

  const int kb = Vh->Th(FK.T);
  CutData cutData(Vh->getInterface(itq).getCutData(kb));     // get the cut data
  const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
  ElementSignEnum the_part = cutK.what_part(domain);

  GQuadraturePoint<R1> tq((qTime)[itq]);
  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq.x, basisFunTime);                  // compute time basic funtions

  for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
  it != cutK.element_end(the_part); ++it){

    const R meas = cutK.mesure(it);

    for(int l=0; l<VF.size();++l) {
      if(!VF[l].on(domain)) continue;

      int lastop = getLastop(0, VF[l].dv);

      const int kv = (VF[l].domv != -1)? VF[l].fespaceV->idxElementFromBackMesh(kb, VF[l].domv) : k;
      const FElement& FKv((*VF[l].fespaceV)[kv]);
      RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
      What_d Fop = Fwhatd(lastop);

      R coef = computeCoef(VF[l],domain,cutK);
      this->initIndex( VF[l].fespaceV);

      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = cutK.toK(it, ip);
        Rd cut_ip = FKv.T.toKref(mip); // in Kref

        const R Cint = meas * ip.getWeight() * In.T.mesure()* tq.a;

        FKv.BF(Fop, cut_ip, fv); // need point in local reference element

        double val_fh = VF[l].fxu_backMesh(kb, domain, mip, In.map(tq));
        double Cst = Cint * coef * VF[l].c * val_fh;


        for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
          const R tval = basisFunTime(it,0,VF[l].dtv);
          for(int i = FKv.dfcbegin(VF[l].cv); i < FK.dfcend(VF[l].cv); ++i) {
            (*this)(FKv.loc2glb(i,it)) +=  Cst * tval * fv(i,VF[l].cv,VF[l].dv);

          }
        }
      }
    }
  }
}


template<typename M>
void BaseCutProblem<M>::addFaceStabilization(const ListItemVF<Rd::d>& VF, const TimeSlab& In) {
  typedef typename Mesh::Element Element;
  const FESpace& Sh =(VF[0].fespaceU)? *VF[0].fespaceU : *Vh;

  for(int k=Sh.first_element(); k<Sh.last_element(); k+= Sh.next_element()) {

    if(!Sh.isCut(k)) continue;

    for(int itq=0;itq<qTime.n;++itq) {
      GTime::iter_step_quadrature = itq;

      for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces

        BaseProblem<M>::addElementMatEdge(VF, k, ifac, In);
        // addElementMatEdgeCut(VF, k, In);
      }
    }
    this->addLocalContribution();
  }
}



template<typename M>
void BaseCutProblem<M>::addElementMatEdge(const ListItemVF<Rd::d>& VF, const int k, const int ifac, const TimeSlab& In) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::Partition Partition;


  const int itq = GTime::iter_step_quadrature;
  assert(VF[0].fespaceU && VF[0].fespaceV);
  const FESpace& Vhu = *VF[0].fespaceU;
  const FESpace& Vhv = *VF[0].fespaceV;
  const FElement& FK(Vhu[k]);
  const int kb = Vh->Th(FK.T);
  CutData cutData(Vh->getInterface(itq).getCutData(kb));

  if(!Vh->isCutSpace() || cutData.edge_isnot_cut(ifac)){
    BaseProblem<M>::addElementMatEdge(VF, k, ifac, In);
    return;
  }

  int the_domain = FK.whichDomain();
  int ifacn = ifac;
  int kn = Vhu.getNeighborElement(k, ifacn, the_domain);

  if(kn == -1) return;         // border edge
  if(k > kn) return;           // only compute one time

  const Partition& cutK =  Partition(FK.T, cutData);  // build the cut

  const Rd normal = FK.T.N(ifac);
  const R meas    = cutK.mesureEdge(ifac, the_domain);
  const R measK   = cutK.mesure(the_domain);
  const R h       = FK.T.lenEdge(ifac);

  GQuadraturePoint<R1> tq((qTime)[itq]);
  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq.x, basisFunTime);

  for(int l=0; l<VF.size();++l) {
    assert(VF[l].fespaceU == VF[l].fespaceV);
    int lastop = getLastop(VF[l].du, VF[l].dv);

    R coef = BaseProblem<M>::computeCoef(VF[l],h,meas, measK,the_domain)*VF[l].getCoef(normal);

    const int ku = (VF[l].domu == 0)? k : kn;
    const int kv = (VF[l].domv == 0)? k : kn;
    const FElement& FKu(Vhu[ku]);
    const FElement& FKv(Vhv[kv]);
    this->initIndex(FKu, FKv);

    bool same = (k == kn );
    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      // const Rd mip = K(K.toKref((RdHatBord)ip, ifac));
      const Rd mip = cutK.toEdge(ifac, (RdHatBord)ip, the_domain); // mip is in the global edge
      const R Cint = meas * ip.getWeight() * In.T.mesure()* tq.a * coef;

      FKu.BF(Fop,FKu.T.toKref(mip), fu);//basisFunMat); // need point in local reference element
      if(!same) FKv.BF(Fop,FKv.T.toKref(mip), fv);

      for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
        for(int jt=In.dfcbegin(0); jt<In.dfcend(0); ++jt) {
          const R tval = basisFunTime(it,0,VF[l].dtv)*basisFunTime(jt,0,VF[l].dtu);
          for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
            for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
              // (*this)(FKv.loc2glb(i,it),FKu.loc2glb(j,jt))  +=  Cint * tval * VF[l].c * fv(i,VF[l].cv,VF[l].dv) * fu(j,VF[l].cu,VF[l].du);
              this->addToLocalContribution(FKv.loc2glb(i,it),FKu.loc2glb(j,jt)) +=  Cint * tval * VF[l].c * fv(i,VF[l].cv,VF[l].dv) * fu(j,VF[l].cu,VF[l].du);
            }
          }
        }
      }
    }
  }
}



template<typename M>
void BaseCutProblem<M>::addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, double val,  const TimeSlab& In) {

  int ndf = this->rhs.size();
  this->rhs.resize(ndf+1);
  this->rhs(ndf) = MPIcf::IamMaster()*val;
  // this->nDoF += 1;
  const int nbDof = (*Vh)[0].NbDoF(); // same local dof for every element
  int lastop = 0;
  for(int i=0;i<VF.size();++i) lastop = Max(lastop, VF[i].du, VF[i].dv);
  lastop += 1;

  KNMK<double> basisFunMat(nbDof,Vh->N,lastop);
  for(int itq=0;itq<qTime.n;++itq) {
    GTime::iter_step_quadrature = itq;               // define which time quad point
    for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {
      addElementLagrange(basisFunMat,VF, k, In);
    }
  }
}

template<typename M>
void BaseCutProblem<M>::addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, double val, int itq, const TimeSlab& In) {

  int ndf = this->rhs.size();
  this->rhs.resize(ndf+1);
  this->rhs(ndf) = MPIcf::IamMaster()*val;
  GTime::iter_step_quadrature = itq;               // define which time quad point
  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

    addElementLagrange(VF, k, In);

  }
}

template<typename M>
void BaseCutProblem<M>::addElementLagrange(const ListItemVF<Rd::d>& VF , const int k, const TimeSlab& In){

  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;
  int nend = this->rhs.size()-1;
  const int itq = GTime::iter_step_quadrature;

  const FElement& FK((*Vh)[k]);
  const int domain = FK.whichDomain();

  const int kb = Vh->Th(FK.T);
  CutData cutData(Vh->getInterface(itq).getCutData(kb));     // get the cut data
  const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
  ElementSignEnum the_part = cutK.what_part(domain);

  GQuadraturePoint<R1> tq((qTime)[itq]);
  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq.x, basisFunTime);

  const R h = FK.T.lenEdge(0);

  for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
  it != cutK.element_end(the_part); ++it){

    const R meas = cutK.mesure(it);

    for(int l=0; l<VF.size();++l) {
      if(VF[l].on(domain)){
        const FESpace& Wh(*VF[l].fespaceV);

        int lastop = getLastop(VF[l].du, VF[l].dv);
        What_d Fop = Fwhatd(lastop);
        RNMK_ fv(this->databf,FK.NbDoF(),FK.N,lastop); //  the value for basic fonction
        R coef = computeCoef(VF[l],domain,cutK);

        for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

          QuadraturePoint ip(qf[ipq]); // integration point
          Rd mip = cutK.toK(it, ip);
          const R Cint = meas * ip.getWeight();

          FK.BF(Fop,FK.T.toKref(mip), fv);

          for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
            for(int j = FK.dfcbegin(VF[l].cv); j < FK.dfcend(VF[l].cv); ++j) {
              (*this)(nend, FK.loc2glb(j,it)+this->mapIdx0[&Wh]) +=  Cint * coef * VF[l].c * fv(j,VF[l].cv,VF[l].dv)*basisFunTime(it,0,VF[l].dtv);
              (*this)(FK.loc2glb(j,it)+this->mapIdx0[&Wh], nend) +=  Cint * coef * VF[l].c * fv(j,VF[l].cv,VF[l].dv)*basisFunTime(it,0,VF[l].dtv);
            }
          }
        }
      }
    }
  }
}
