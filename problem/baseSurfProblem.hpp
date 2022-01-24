#ifndef _BASE_SURFACE_PROBLEM_HPP
#define _BASE_SURFACE_PROBLEM_HPP



static R2 operator*(const KNM_<double>& mat, const R2 v){
  R2 val;
  val.x = mat(0,0)*v.x + mat(0,1)*v.y;
  val.y = mat(1,0)*v.x + mat(1,1)*v.y;
  return val;
}
static R3 operator*(const KNM_<double>& mat, const R3 v){
  R3 val;
  val.x = mat(0,0)*v.x + mat(0,1)*v.y + mat(0,2)*v.z;
  val.y = mat(1,0)*v.x + mat(1,1)*v.y + mat(1,2)*v.z;
  val.z = mat(2,0)*v.x + mat(2,1)*v.y + mat(2,2)*v.z;
  return val;
}
static R determinant(const KNM_<double>& a) {
  R sum = 0.;
  if(a.N() == 2 && a.M() == 2) {
    sum = a(0,0)*a(1,1) - a(1,0)*a(0,1);
  }
  else if (a.N() == 3 && a.M() == 3){
    sum = a(0,0)*( a(1,1)*a(2,2) - a(1,2)*a(2,1) )
        - a(0,1)*( a(1,0)*a(2,2) - a(2,0)*a(1,2) )
        + a(0,2)*( a(1,0)*a(2,1) - a(2,0)*a(1,1));
  }
  else  assert(0);

  assert(std::fabs(sum) > 1e-16);
  return std::fabs(sum);
}


/*
        INTEGRATION SPACE  \int_\Omega f(x) dx
*/
template<typename M>
void BaseProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF, const Interface& gamma,list<int> label,const Mapping& mapping) {

  bool all_label = (label.size() == 0);
  for(int iface=gamma.first_element(); iface<gamma.last_element(); iface+=gamma.next_element()) {
    const typename Interface::FaceIdx& face = gamma[iface];  // the face
    if(contain(label, face.lab) || all_label) {
      addElementMat(VF, gamma, iface,mapping);
    }
  }
}

template<typename M>
void BaseProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF, const Interface& gamma,const GMacro& macro, list<int> label,const Mapping& mapping) {

  bool all_label = (label.size() == 0);
  for(auto it=macro.macro_element.begin(); it!=macro.macro_element.end();++it) {
    for(int i=0;i<it->second.idx_element.size();++i) {

      int iface = it->second.idx_element[i];
      const typename Interface::FaceIdx& face = gamma[iface];  // the face
      if(contain(label, face.lab) || all_label) {
        addElementMat(VF, gamma, iface,mapping);
      }
    }
  }
}



template<typename M>
void BaseProblem<M>::addElementMat(const ListItemVF<Rd::d>& VF,const Interface& interface, const int iface,const Mapping& mapping) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Partition Partition;
  typedef typename TypeCutData<Rd::d>::CutData CutData;

  KNM<double> invJ(Rd::d, Rd::d);

  const int kb = interface.idxElementOfFace(iface);   // idx on
  const int km = mapping.idxElementFromBackMesh(kb);
  const typename Interface::FaceIdx& face = interface[iface];  // the face
  const R meas = interface.computeDx(face).norm();
  const double h = meas;
  const Rd linear_normal(-interface.normal(iface));
  assert(fabs(linear_normal.norm()-1)<1e-14);

  CutData cutData(interface.getCutData(kb));     // get the cut data
  const Partition& cutK =  Partition(Vh->backSpace->Th[kb], cutData);

  for(int l=0; l<VF.size();++l) {
    bool same = (VF[l].fespaceU==VF[l].fespaceV);

    int lastop = getLastop(VF[l].du, VF[l].dv);

    const int ku = VF[l].fespaceU->idxElementFromBackMesh(kb,VF[l].domu);
    const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb,VF[l].domv);

    const FElement& FKu((*VF[l].fespaceU)[ku]);
    const FElement& FKv((*VF[l].fespaceV)[kv]);
    double measK = FKu.getMeasure();

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    R coef = computeCoefInterface(VF[l],h,meas,measK );

    this->initIndex(VF[l].fespaceU, VF[l].fespaceV);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = interface.mapToFace(face,(RdHatBord)ip);
      mapping.computeInverseJacobian(km, mip, invJ);


      Rd normal = invJ*linear_normal;
      double DetJ = 1./determinant(invJ);

      const R Cint = meas * ip.getWeight()*DetJ*normal.norm();
      normal = normal / normal.norm();

      double Cnormal = VF[l].getCoef(normal);

      FKu.BF(Fop,FKu.T.toKref(mip), fu);//basisFunMat); // need point in local reference element
      if(!same) FKv.BF(Fop,FKv.T.toKref(mip), fv);

      mapping.transform(FKu, fu, invJ);
      if(!same) mapping.transform(FKv, fv, invJ);
      R val_fh = VF[l].fxu_backMesh(kb, 0, mip, normal);

      for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
        for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
          // (*this)(FKv.loc2glb(i),FKu.loc2glb(j)) += Cint * coef *  VF[l].c * Cnormal * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du);
          this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j)) += Cint * coef * val_fh * VF[l].c * Cnormal * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du);
        }
      }
    }

    this->resetIndex();
    this->addLocalContribution();

  }
}

template<typename M>
void BaseProblem<M>::addLinear(const ListItemVF<Rd::d>& VF, const Interface& gamma,list<int> label,const Mapping& mapping) {
  bool all_label = (label.size() == 0);
  for(int iface=gamma.first_element(); iface<gamma.last_element(); iface+=gamma.next_element()) {
    const typename Interface::FaceIdx& face = gamma[iface];  // the face
    if(contain(label, face.lab) || all_label) {
      addElementRHS(VF, gamma, iface,mapping);
    }
  }
}

template<typename M>
void BaseProblem<M>::addElementRHS(const ListItemVF<Rd::d>& VF, const Interface& interface, const int iface,const Mapping& mapping) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Partition Partition;
  typedef typename TypeCutData<Rd::d>::CutData CutData;

  KNM<double> invJ(Rd::d, Rd::d);
  const int kb = interface.idxElementOfFace(iface);   // idx on
  const int km = mapping.idxElementFromBackMesh(kb);
  const typename Interface::FaceIdx& face = interface[iface];  // the face
  const R meas = interface.computeDx(face).norm();
  const double h = meas;
  const Rd linear_normal(-interface.normal(iface));
  assert(fabs(linear_normal.norm()-1)<1e-14);

  CutData cutData(interface.getCutData(kb));     // get the cut data
  const Partition& cutK =  Partition(Vh->backSpace->Th[kb], cutData);

  for(int l=0; l<VF.size();++l) {
    int lastop = getLastop(VF[l].dv, VF[l].dv);
    // assert(VF[l].fespaceU && VF[l].fespaceV);
    assert(VF[l].fespaceV);
    // const int kf = VF[l].fespaceU->idxElementFromBackMesh(kb,VF[l].domu);
    const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb,VF[l].domv);
    const FElement& FKv((*VF[l].fespaceV)[kv]);
    double measK = FKv.getMeasure();
    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    R coef = computeCoefInterface(VF[l],h,meas, measK);

    this->initIndex(VF[l].fespaceV);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = interface.mapToFace(face,(RdHatBord)ip);

      mapping.computeInverseJacobian(km, mip, invJ);
      Rd normal = invJ*linear_normal;
      double DetJ = 1./determinant(invJ);
      const R Cint = meas * ip.getWeight()*DetJ*normal.norm();
      normal = normal / normal.norm();
      double Cnormal = VF[l].getCoef(normal);

      FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element
      mapping.transform(FKv, fv, invJ);
      R val_fh = VF[l].fx_backMesh_V(kb, VF[l].domu, mip, normal)
      //VF[l].fxU(kf, mip, normal)
      *VF[l].fx_backMesh_V(kb, VF[l].domv, mip, normal);
      double Cst = Cint * coef * VF[l].c * Cnormal * val_fh;
      for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
        (*this)(FKv.loc2glb(i)) +=  Cst * fv(i,VF[l].cv,VF[l].dv) ;
      }
    }
    this->resetIndex();
  }

}

template<typename M>
void BaseProblem<M>::addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, const Interface& gamma, double val,const Mapping& mapping) {

  int ndf = rhs.size();
  rhs.resize(ndf+1);
  rhs(ndf) = val;
  nDoF += 1;

  for(int iface=gamma.first_element(); iface<gamma.last_element(); iface+=gamma.next_element()) {

    addElementLagrange(VF, gamma, iface,mapping);

  }
}

template<typename M>
void BaseProblem<M>::addElementLagrange(const ListItemVF<Rd::d>& VF, const Interface& interface, const int iface,const Mapping& mapping) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;

  KNM<double> invJ(Rd::d, Rd::d);
  const int kb = interface.idxElementOfFace(iface);   // idx on
  const int km = mapping.idxElementFromBackMesh(kb);
  const typename Interface::FaceIdx& face = interface[iface];  // the face
  const R meas = interface.computeDx(face).norm();
  const double h = meas;
  const Rd linear_normal(-interface.normal(iface));
  assert(fabs(linear_normal.norm()-1)<1e-14);
  int nend = rhs.size()-1;


  for(int l=0; l<VF.size();++l) {
    int lastop = getLastop(VF[l].dv, 0);

    const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb);
    const FElement& FKv((*VF[l].fespaceV)[kv]);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    this->initIndex(VF[l].fespaceV);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = interface.mapToFace(face,(RdHatBord)ip);
      mapping.computeInverseJacobian(km, mip, invJ);
      Rd normal = invJ*linear_normal;
      double DetJ = 1./determinant(invJ);
      const R Cint = meas * ip.getWeight()*DetJ*normal.norm();

      FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element

      mapping.transform(FKv, fv, invJ);


      for(int j = FKv.dfcbegin(VF[l].cv); j < FKv.dfcend(VF[l].cv); ++j) {
        (*this)(nend, FKv.loc2glb(j)) +=  Cint *  VF[l].c *fv(j,VF[l].cv,VF[l].dv);
        (*this)(FKv.loc2glb(j), nend) +=  Cint *  VF[l].c *fv(j,VF[l].cv,VF[l].dv);

      }
    }
    this->resetIndex();
  }
}

template<typename M>
void BaseProblem<M>::addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, const Interface& gamma, const TimeSlab& In, double tq, double val,const Mapping& mapping) {

  int ndf = rhs.size();
  rhs.resize(ndf+1);
  rhs(ndf) = val;
  nDoF += 1;

  for(int iface=gamma.first_element(); iface<gamma.last_element(); iface+=gamma.next_element()) {

    addElementLagrange(VF, gamma,iface, In, tq, mapping);

  }
}

template<typename M>
void BaseProblem<M>::addElementLagrange(const ListItemVF<Rd::d>& VF, const Interface&  interface, const int iface, const TimeSlab& In, double tq, const Mapping& mapping) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;

  KNM<double> invJ(Rd::d, Rd::d);
  const int kb = interface.idxElementOfFace(iface);   // idx on
  const int km = mapping.idxElementFromBackMesh(kb);
  const typename Interface::FaceIdx& face = interface.getFace(iface); ;  // the face
  const R meas = interface.computeDx(face).norm();
  const double h = meas;
  const Rd linear_normal(-interface.normal(iface));
  assert(fabs(linear_normal.norm()-1)<1e-14);
  int nend = rhs.size()-1;

  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq, basisFunTime);                  // compute time basic funtions


  for(int l=0; l<VF.size();++l) {
    int lastop = getLastop(VF[l].dv, 0);

    const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb);
    const FElement& FKv((*VF[l].fespaceV)[kv]);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    this->initIndex(VF[l].fespaceV);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = interface.mapToFace(face,(RdHatBord)ip);
      mapping.computeInverseJacobian(km, mip, invJ);
      Rd normal = invJ*linear_normal;
      double DetJ = 1./determinant(invJ);
      const R Cint = meas * ip.getWeight()*DetJ*normal.norm();

      FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element

      mapping.transform(FKv, fv, invJ);

      for(int jt=In.dfcbegin(0); jt<In.dfcend(0); ++jt) {
        const R tval = basisFunTime(jt,0,VF[l].dtv);
        for(int j = FKv.dfcbegin(VF[l].cv); j < FKv.dfcend(VF[l].cv); ++j) {
          (*this)(nend, FKv.loc2glb(j,jt)) +=  Cint *  tval * VF[l].c *fv(j,VF[l].cv,VF[l].dv);
          (*this)(FKv.loc2glb(j,jt), nend) +=  Cint *  tval * VF[l].c *fv(j,VF[l].cv,VF[l].dv);
        }
      }
    }
    this->resetIndex();
  }

}

template<typename M>
void BaseProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF, const Interface& gamma,const CNode& nodeEval) {

  for(int iface=gamma.first_element(); iface<gamma.last_element(); iface+=gamma.next_element()) {
    addElementMat(VF, gamma, iface,nodeEval);
  }
}

template<typename M>
void BaseProblem<M>::addElementMat(const ListItemVF<Rd::d>& VF,const Interface& interface, const int iface0,const CNode& nodeEval) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Partition Partition;
  typedef typename TypeCutData<Rd::d>::CutData CutData;
  const Mesh& Th(*interface.backMesh);

  const int kb = interface.idxElementOfFace(iface0);   // idx on
  const typename Interface::FaceIdx& face0 = interface[iface0];  // the face

  CutData cutData(interface.getCutData(kb));     // get the cut data
  const Partition& cutK =  Partition(Vh->backSpace->Th[kb], cutData);


  for(int inode=0;inode<2;++inode) {

    int idx_node = face0[inode];
    int e0 = interface.edge_of_node_[idx_node];
    int e1 = e0;
    int knb = Th.ElementAdj(kb,e1);
    const Rd mip = interface(idx_node);

    if(kb > knb) continue;  //only does evaluation once;

    int iface1 = interface.face_of_element_.find(knb)->second;
    const typename Interface::FaceIdx& face1 = interface[iface1];  // the face
    int inode1 = (interface.edge_of_node_[face1[0]]==e1)? 0 : 1;

    // here the normal can be inward or outward since it is only squared in the gradient
    // then we can define the normal that will give the tangent to be the outward conormal
    Rd tangent0 = interface(face0[inode] ) - interface(face0[(inode +1)%2]);
    Rd tangent1 = interface(face1[inode1]) - interface(face1[(inode1+1)%2]);
    tangent0 /= tangent0.norm();
    tangent1 /= tangent1.norm();
    const Rd normal0(tangent0.y, -tangent0.x);
    const Rd normal1(tangent1.y, -tangent1.x);

    // assert( (n1, tangent0) >0);
    // assert( (n2, tangent1) >0);

    for(int l=0; l<VF.size();++l) {
      assert(!VF[l].fespaceU->isCutSpace());
      assert(VF[l].fespaceU==VF[l].fespaceV);
      int lastop = getLastop(VF[l].du, VF[l].dv);
      this->initIndex(VF[l].fespaceU, VF[l].fespaceV);

      const int kbu = (VF[l].domu == 0)? kb : knb;
      const int kbv = (VF[l].domv == 0)? kb : knb;
      const Rd normalu = (VF[l].domu == 0)? normal0 : normal1;
      const Rd normalv = (VF[l].domv == 0)? normal0 : normal1;

      const int ku = VF[l].fespaceU->idxElementFromBackMesh(kbu);
      const int kv = VF[l].fespaceV->idxElementFromBackMesh(kbv);
      bool same = (ku == kv);

      const FElement& FKu((*VF[l].fespaceU)[ku]);
      const FElement& FKv((*VF[l].fespaceV)[kv]);

      RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
      RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
      What_d Fop = Fwhatd(lastop);

      FKu.BF(Fop,FKu.T.toKref(mip), fu);//basisFunMat); // need point in local reference element
      if(!same) FKv.BF(Fop,FKv.T.toKref(mip), fv);

      R coef = computeCoefInterface(VF[l],0,0,0 );

      R val_fh = VF[l].fx_backMesh_U(kbu, 0, mip, normalu)
                *VF[l].fx_backMesh_V(kbv, 0, mip, normalv);
      double Cnormal = VF[l].getCoefU(normalu)*VF[l].getCoefV(normalv);

      double Cst = Cnormal*coef*val_fh*VF[l].c;

      for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
        for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
          (*this)(FKv.loc2glb(i),FKu.loc2glb(j)) += Cst*fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du);
        }
      }
      this->resetIndex();
    }
  }
}


template<typename M>
void BaseProblem<M>::addLinear(const ListItemVF<Rd::d>& VF, const Interface& gamma, const CBorder& b, list<int> label ) {
  typedef typename Mesh::BorderElement BorderElement;
  bool all_label = (label.size() == 0);

  for( int ifac = Vh->first_boundary_element(); ifac < Vh->last_boundary_element(); ifac+=Vh->next_boundary_element()) {
    const BorderElement & face(Vh->Th.be(ifac));
    if(contain(label, face.lab) || all_label) {
      int ifaceK;
      const int kb = Vh->Th.BoundaryElement(ifac, ifaceK);
      if(gamma.face_of_element_.find(kb) != gamma.face_of_element_.end()){
        // std::cout << "hey " << face[0] << "\t" << face[1] << std::endl;

        addElementRHSBorder(VF, gamma, ifac);


      }
    }
  }
}

template<typename M>
void BaseProblem<M>::addElementRHSBorder(const ListItemVF<Rd::d>& VF,const Interface& interface, const int iface_B) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Partition Partition;
  typedef typename TypeCutData<Rd::d>::CutData CutData;
  typedef typename Mesh::BorderElement BorderElement;

  const Mesh& Th(*interface.backMesh);


  // GET THE BOUNDARY FACE
  int ifaceK;
  const int kb = Th.BoundaryElement(iface_B, ifaceK);
  const BorderElement & BE(Th.be(iface_B));
  // Segment<Rd> segment(BE[0], BE[1]);

  // GET THE INTERFACE FACE
  int iface_G = interface.idxFaceOfElement(kb);
  const typename Interface::FaceIdx& face = interface[iface_G];

  // FIND THE POINT OF THE INTERFACE ON THE BOUNDARY
  int inode = (interface.edge_of_node_[face[0]]==ifaceK)?0:1;
  int idx_node = face[inode];
  const Rd mip = interface(idx_node);


  for(int l=0; l<VF.size();++l) {
    assert(!VF[l].fespaceV->isCutSpace());
    assert(VF[l].fespaceV);
    int lastop = getLastop(0,VF[l].dv);
    this->initIndex( VF[l].fespaceV);

    // THE ELEMENT IN THE SURFACE MESH
    const FESpace& Sh(*VF[l].fespaceU);
    const int k = Sh.idxElementFromBackMesh(kb);
    const FElement& FK(Sh[k]);

    // THE NORMAL
    // Rd normal = FK.T.N(ifaceK);
    //
    // // EVALUATE THE BASIS FUNCTION
    // RNMK_ fv(this->databf,FK.NbDoF(),FK.N,lastop); //  the value for basic fonction
    // What_d Fop = Fwhatd(lastop);
    // FK.BF(Fop,FK.T.toKref(mip), fv);
    //
    // // COMPUTE COEFFICIENTS
    // R coef = computeCoefInterface(VF[l],0,0,0 );
    // R val_fh = VF[l].fx_backMesh_U(kb, 0, mip, normal)
    // *VF[l].fx_backMesh_V(kb, 0, mip, normal);
    // double Cnormal = VF[l].getCoefU(normal)*VF[l].getCoefV(normal);
    // double Cst = Cnormal*coef*val_fh*VF[l].c;
    //
    // // FILL THE MATRIX
    // for(int i = FK.dfcbegin(VF[l].cv); i < FK.dfcend(VF[l].cv); ++i) {
    //   for(int j = FK.dfcbegin(VF[l].cu); j < FK.dfcend(VF[l].cu); ++j) {
    //     (*this)(FK.loc2glb(i),FK.loc2glb(j)) += Cst*fv(i,VF[l].cv,VF[l].dv)*fv(j,VF[l].cu,VF[l].du);
    //   }
    // }
    this->resetIndex();
  }
}


/*
        INTEGRATION SPACE AND TIME  \int_\Omega(t) f(x,t) dx dt
*/
template<typename M>
void BaseProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, int itq, const TimeSlab& In, KN<const Mapping*>  vector_mapping) {

  assert(gamma.size() >= qTime.n);

  GTime::iter_step_quadrature = itq;               // define which time quad point
  const Mapping& mapping = (vector_mapping.size()>0)? *(vector_mapping[itq]) : DataMapping<Mesh>::Id;
  for(int iface=gamma[itq]->first_element(); iface<gamma[itq]->last_element(); iface+=gamma[itq]->next_element()) {
    addElementMat(VF, gamma, iface, In, mapping);

    this->addLocalContribution();

  }
}

template<typename M>
void BaseProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, const TimeSlab& In, KN<const Mapping*>  mapping) {
  assert(gamma.size() >= qTime.n);

  for(int itq=0;itq<qTime.n;++itq) {
    addBilinear(VF, gamma, itq, In, mapping);
  }
}

template<typename M>
void BaseProblem<M>::addElementMat(const ListItemVF<Rd::d>& VF, const TimeInterface<M>&  interface, const int iface, const TimeSlab& In, const Mapping& mapping) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Partition Partition;
  typedef typename TypeCutData<Rd::d>::CutData CutData;

  const int itq = GTime::iter_step_quadrature;

  KNM<double> invJ(Rd::d, Rd::d);

  const int kb = interface(itq)->idxElementOfFace(iface);   // idx on
  const int km = mapping.idxElementFromBackMesh(kb);
  const typename Interface::FaceIdx& face = interface(itq)->getFace(iface); ;  // the face
  const R meas = interface(itq)->computeDx(face).norm();
  const double h = meas;
  const Rd linear_normal(-interface(itq)->normal(iface));
  assert(fabs(linear_normal.norm()-1)<1e-14);

  GQuadraturePoint<R1> tq((qTime)[itq]);
  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq.x, basisFunTime);                  // compute time basic funtions

  CutData cutData(interface(itq)->getCutData(kb));     // get the cut data
  const Partition& cutK =  Partition(Vh->backSpace->Th[kb], cutData);

  for(int l=0; l<VF.size();++l) {
    bool same = (VF[l].fespaceU==VF[l].fespaceV);

    int lastop = getLastop(VF[l].du, VF[l].dv);

    const int ku = VF[l].fespaceU->idxElementFromBackMesh(kb, VF[l].domu);
    const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb, VF[l].domv);
    const FElement& FKu((*VF[l].fespaceU)[ku]);
    const FElement& FKv((*VF[l].fespaceV)[kv]);
    double measK = FKu.getMeasure();


    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    R coef = computeCoefInterface(VF[l],h,meas,measK );

    this->initIndex(VF[l].fespaceU, VF[l].fespaceV);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = interface(itq)->mapToFace(face,(RdHatBord)ip);
      mapping.computeInverseJacobian(km, mip, invJ);

      Rd normal = invJ*linear_normal;
      double DetJ = 1./determinant(invJ);

      const R Cint = meas * ip.getWeight()*DetJ*normal.norm() * In.T.mesure()* tq.a;
      normal = normal / normal.norm();

      double Cnormal = VF[l].getCoef(normal);

      double Cst = coef * Cnormal * VF[l].c * Cint;

      FKu.BF(Fop,FKu.T.toKref(mip), fu);//basisFunMat); // need point in local reference element
      if(!same) FKv.BF(Fop,FKv.T.toKref(mip), fv);

      double val_fh = VF[l].fxu_backMesh(kb, 0, mip, In.map(tq),normal);

      mapping.transform(FKu, fu, invJ);
      if(!same) mapping.transform(FKv, fv, invJ);
      for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
        for(int jt=In.dfcbegin(0); jt<In.dfcend(0); ++jt) {
          const R tval = basisFunTime(it,0,VF[l].dtv)*basisFunTime(jt,0,VF[l].dtu);
          for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
            for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
              // (*this)(FKv.loc2glb(i,it),FKu.loc2glb(j,jt)) += Cst * tval * val_fh* fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du);
              this->addToLocalContribution(FKv.loc2glb(i,it),FKu.loc2glb(j,jt))+=  Cst * tval * val_fh* fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du);
            }
          }
        }
      }
    }
  }
}


template<typename M>
void BaseProblem<M>::addLinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, int itq, const TimeSlab& In, KN<const Mapping*> vector_mapping) {
  assert(gamma.size() >= qTime.n);
  GTime::iter_step_quadrature = itq;               // define which time quad point
  const Mapping& mapping = (vector_mapping.size()>0)? *(vector_mapping[itq]) : DataMapping<Mesh>::Id;

  for(int iface=gamma[itq]->first_element(); iface<gamma[itq]->last_element(); iface+=gamma[itq]->next_element()) {
    addElementRHS(VF, gamma, iface, In,mapping);
  }
}

template<typename M>
void BaseProblem<M>::addLinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, const TimeSlab& In, KN<const Mapping*> mapping) {
  assert(gamma.size() >= qTime.n);
  for(int itq=0;itq<qTime.n;++itq) {
    addLinear(VF, gamma, itq, In, mapping);
  }
}


template<typename M>
void BaseProblem<M>::addElementRHS(const ListItemVF<Rd::d>& VF, const TimeInterface<M>&  interface, const int iface, const TimeSlab& In, const Mapping& mapping) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Partition Partition;
  typedef typename TypeCutData<Rd::d>::CutData CutData;

  const int itq = GTime::iter_step_quadrature;

  KNM<double> invJ(Rd::d, Rd::d);
  const int kb = interface(itq)->idxElementOfFace(iface);   // idx on
  const int km = mapping.idxElementFromBackMesh(kb);
  const typename Interface::FaceIdx& face = interface(itq)->getFace(iface); ;  // the face
  const R meas = interface(itq)->computeDx(face).norm();
  const double h = meas;
  const Rd linear_normal(-interface(itq)->normal(iface));
  assert(fabs(linear_normal.norm()-1)<1e-14);

  GQuadraturePoint<R1> tq((qTime)[itq]);
  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq.x, basisFunTime);                  // compute time basic funtions

  CutData cutData(interface(itq)->getCutData(kb));     // get the cut data
  const Partition& cutK =  Partition(Vh->backSpace->Th[kb], cutData);


  for(int l=0; l<VF.size();++l) {
    int lastop = getLastop(VF[l].dv, VF[l].dv);

    const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb, VF[l].domv);
    const FElement& FKv((*VF[l].fespaceV)[kv]);
    double measK = FKv.getMeasure();
    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    R coef = computeCoefInterface(VF[l],h,meas,measK );

    this->initIndex(VF[l].fespaceV);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = interface(itq)->mapToFace(face,(RdHatBord)ip);

      mapping.computeInverseJacobian(km, mip, invJ);
      Rd normal = invJ*linear_normal;
      double DetJ = 1./determinant(invJ);
      const R Cint = meas * ip.getWeight()*DetJ*normal.norm() * In.T.mesure()* tq.a;
      normal = normal / normal.norm();
      double Cnormal = VF[l].getCoef(normal);

      FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element
      mapping.transform(FKv, fv, invJ);
      R val_fh = VF[l].fxu_backMesh(kb, 0, mip, In.map(tq), normal);

      double Cst = coef * Cnormal * VF[l].c * Cint * val_fh;

      for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
        const R tval = basisFunTime(it,0,VF[l].dtv);
        for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
          (*this)(FKv.loc2glb(i,it)) +=  Cst * tval * fv(i,VF[l].cv,VF[l].dv);
        }
      }
    }
    this->resetIndex();
  }
}


template<typename M>
void BaseProblem<M>::addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, const TimeSlab& In, double val, KN<const Mapping*> vector_mapping) {
  assert(gamma.size() >= qTime.n);

  int ndf = rhs.size();
  rhs.resize(ndf+1);
  rhs(ndf) = val;
  nDoF += 1;

  for(int itq=0;itq<qTime.n;++itq) {
    GTime::iter_step_quadrature = itq;               // define which time quad point
    const Mapping& mapping = (vector_mapping.size()>0)? *(vector_mapping[itq]) : DataMapping<Mesh>::Id;

    for(int iface=gamma[itq]->first_element(); iface<gamma[itq]->last_element(); iface+=gamma[itq]->next_element()) {

      addElementLagrange(VF, gamma, iface, In,mapping);
    }
  }
}

template<typename M>
void BaseProblem<M>::addElementLagrange(const ListItemVF<Rd::d>& VF, const TimeInterface<M>&  interface, const int iface, const TimeSlab& In, const Mapping& mapping) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  const int itq = GTime::iter_step_quadrature;

  KNM<double> invJ(Rd::d, Rd::d);
  const int kb = interface(itq)->idxElementOfFace(iface);   // idx on
  const int km = mapping.idxElementFromBackMesh(kb);
  const typename Interface::FaceIdx& face = interface(itq)->getFace(iface); ;  // the face
  const R meas = interface(itq)->computeDx(face).norm();
  const double h = meas;
  const Rd linear_normal(-interface(itq)->normal(iface));
  assert(fabs(linear_normal.norm()-1)<1e-14);
  int nend = rhs.size()-1;

  GQuadraturePoint<R1> tq((qTime)[itq]);
  KNMK<double> basisFunTime(In.NbDoF(),1,op_dz+1);
  In.BF(tq.x, basisFunTime);                  // compute time basic funtions


  for(int l=0; l<VF.size();++l) {
    int lastop = getLastop(VF[l].dv, 0);

    const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb);
    const FElement& FKv((*VF[l].fespaceV)[kv]);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    this->initIndex(VF[l].fespaceV);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = interface(itq)->mapToFace(face,(RdHatBord)ip);
      mapping.computeInverseJacobian(km, mip, invJ);
      Rd normal = invJ*linear_normal;
      double DetJ = 1./determinant(invJ);
      const R Cint = meas * ip.getWeight()*DetJ*normal.norm() * In.T.mesure()* tq.a;

      FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element

      mapping.transform(FKv, fv, invJ);

      for(int jt=In.dfcbegin(0); jt<In.dfcend(0); ++jt) {
        const R tval = basisFunTime(jt,0,VF[l].dtv);
        for(int j = FKv.dfcbegin(VF[l].cv); j < FKv.dfcend(VF[l].cv); ++j) {
          (*this)(nend, FKv.loc2glb(j,jt)) +=  Cint *  tval * VF[l].c *fv(j,VF[l].cv,VF[l].dv);
          (*this)(FKv.loc2glb(j,jt), nend) +=  Cint *  tval * VF[l].c *fv(j,VF[l].cv,VF[l].dv);
        }
      }
    }
    this->resetIndex();
  }

}






#endif
