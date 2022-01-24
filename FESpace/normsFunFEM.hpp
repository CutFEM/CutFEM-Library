template<typename M>
double L2normSurf_2(const ExpressionVirtual& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i, double t),double tt,const GFESpace<M>& Vh) {
  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QFB QFB;
  typedef typename FElement::Rd Rd;
  typedef typename Mesh::Partition Partition;
  typedef GenericInterface<Mesh> Interface;
  typedef typename QFB::QuadraturePoint QuadraturePoint;

  const Interface& interface(Vh.getInterface(0));
  const QFB& qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
  What_d Fop = Fwhatd(op_id);

  double val = 0.;

  for(int iface=interface.first_element(); iface<interface.last_element(); iface+=interface.next_element()) {

    const int kb = interface.idxElementOfFace(iface);   // idx on backMesh
    const typename Interface::Face& face = interface.getFace(iface);  // the face
    const R meas = interface.computeDx(face).norm();

    int k = Vh.idxElementFromBackMesh(kb,0);
    const FElement& FK(Vh[k]);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = interface.mapToFace(face,(typename FElement::RdHatBord)ip);
      const R Cint = meas * ip.getWeight();

      double a = fh.eval(k, mip, tt) - fex(mip, fh.cu, tt);
      val += Cint * a * a;

    }
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;

}

template<typename M>
double L2normSurf( const FunFEM<M>& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i, double t), double tt,int c0=0, int num_comp = GFESpace<M>::FElement::Rd::d) {

    const GFESpace<M>& Vh(*fh.Vh);

    double val = 0;
    for(int i=c0;i<num_comp+c0;++i) {
      ExpressionFunFEM<M> ui(fh, i, op_id);
      val += L2normSurf_2(ui,fex, tt,Vh);
    }
    return sqrt(val);
  }



template<typename M>
double L2norm_2(const ExpressionVirtual& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i), const GFESpace<M>& Vh) {
  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QF QF;
  typedef typename FElement::Rd Rd;
  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;
  typedef typename TypeCutData<Rd::d>::CutData CutData;

  const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
  What_d Fop = Fwhatd(op_id);
  double val = 0.;

  for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {

    const FElement& FK(Vh[k]);

    const R meas = FK.getMeasure();
    for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qf[ipq]); // integration point
      Rd mip = FK.map(ip);
      const R Cint = meas * ip.getWeight();
      double a = fh.eval(k, mip) - fex(mip, fh.cu);

      val += Cint * a * a;
    }
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;
}

template<typename M>
double L2norm_2(const ExpressionVirtual& fh, const GFESpace<M>& Vh) {
  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QF QF;
  typedef typename FElement::Rd Rd;
  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;
  typedef typename TypeCutData<Rd::d>::CutData CutData;

  const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
  What_d Fop = Fwhatd(op_id);
  double val = 0.;

  for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {

    const FElement& FK(Vh[k]);

    const R meas = FK.getMeasure();
    for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qf[ipq]); // integration point
      Rd mip = FK.map(ip);
      const R Cint = meas * ip.getWeight();
      double a = fh.eval(k, mip);
      val += Cint * a * a;
    }
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;
}


template<typename M>
double L2normCut_2(const ExpressionVirtual& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i, int dom),int domain, const GFESpace<M>& Vh, const MacroElement* macro) {
  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QF QF;
  typedef typename FElement::Rd Rd;
  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;
  typedef typename TypeCutData<Rd::d>::CutData CutData;

  const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));

  What_d Fop = Fwhatd(op_id);

  double val = 0.;

  for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {

    const FElement& FK(Vh[k]);
    if(domain != FK.whichDomain()) continue;

    const int kb = Vh.Th(FK.T);
    CutData cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
    const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
    ElementSignEnum the_part = cutK.what_part(domain);

    int kk = k;
    // if(macro){
    //   if(!macro->isRootFat(k)) {
    //     kk = macro->getIndexRootElement(k);
    //   }
    // }

    for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
    it != cutK.element_end(the_part); ++it){

      const R meas = cutK.mesure(it);

      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = cutK.toK(it, ip);
        const R Cint = meas * ip.getWeight();

        double a = fh.eval(kk, mip) - fex(mip, fh.cu, domain);

        val += Cint * a * a;

      }
    }

  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;

}

template<typename M>
double L2normCut_2(const ExpressionVirtual& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i, int dom), const GFESpace<M>& Vh, const MacroElement* macro) {

  return L2normCut_2(fh,fex, 0, Vh, macro) + L2normCut_2(fh,fex,1, Vh, macro);
}



template<typename M>
double L2norm( const FunFEM<M>& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i),int c0=0, int num_comp = GFESpace<M>::FElement::Rd::d) {

  const GFESpace<M>& Vh(*fh.Vh);

  double val = 0;
  for(int i=c0;i<num_comp+c0;++i) {
    ExpressionFunFEM<M> ui(fh, i, op_id);
    val += L2norm_2(ui,fex, Vh);
  }
  return sqrt(val);
}

template<typename M>
double L2norm( const ExpressionVirtual& fh, R (fex)(const typename GFESpace<M>::FElement::Rd, int i), const GFESpace<M>& Vh) {

    double val = L2norm_2(fh,fex, Vh);

  return sqrt(val);
}

template<typename M>
double L2norm( const FunFEM<M>& fh,int c0=0, int num_comp = GFESpace<M>::FElement::Rd::d) {

  const GFESpace<M>& Vh(*fh.Vh);

  double val = 0;
  for(int i=c0;i<num_comp+c0;++i) {
    ExpressionFunFEM<M> ui(fh, i, op_id);
    val += L2norm_2(ui, Vh);
  }
  return sqrt(val);
}

template<typename M>
double L2norm( const ExpressionVirtual& fh, const GFESpace<M>& Vh) {

    double val = L2norm_2(fh, Vh);

  return sqrt(val);
}


template<typename M>
double L2normCut( const FunFEM<M>& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i, int dom),int c0=0, int num_comp = GFESpace<M>::FElement::Rd::d, const MacroElement* macro=nullptr) {

  const GFESpace<M>& Vh(*fh.Vh);

  double val = 0;
  for(int i=c0;i<num_comp+c0;++i) {
    ExpressionFunFEM<M> ui(fh, i, op_id);
    val += L2normCut_2(ui,fex, Vh, macro);
  }
  return sqrt(val);
}

template<typename M>
double L2normCut( const ExpressionVirtual& fh, R (fex)(const typename GFESpace<M>::FElement::Rd, int i, int dom), const GFESpace<M>& Vh, const MacroElement* macro=nullptr) {

    double val = L2normCut_2(fh,fex, Vh, macro);

  return sqrt(val);
}


template<typename M>
double maxNorm(const ExpressionVirtual& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i),const GFESpace<M>& Vh) {

  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QF QF;
  typedef typename FElement::Rd Rd;
  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;
  typedef typename TypeCutData<Rd::d>::CutData CutData;

  const QF& qf(*QF_Simplex<typename FElement::RdHat>(0));
  What_d Fop = Fwhatd(op_id);
  double val = 0.;

  for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {

    const FElement& FK(Vh[k]);
    const R meas = FK.getMeasure();
    for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qf[ipq]); // integration point
      Rd mip = FK.map(ip);
      const R Cint = meas * ip.getWeight();

      val = max(val, fabs(fh.eval(k, mip)-fex(mip, fh.cu)));
    }
  }

  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_MAX);

  return val_receive;
}

template<typename M>
double maxNormCut(const ExpressionVirtual& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i, int dom),int domain,const GFESpace<M>& Vh) {

    typedef M Mesh;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename FElement::Rd Rd;
    typedef typename Mesh::Partition Partition;
    typedef typename QF::QuadraturePoint QuadraturePoint;
    typedef typename TypeCutData<Rd::d>::CutData CutData;

    const QF& qf(*QF_Simplex<typename FElement::RdHat>(0));

    What_d Fop = Fwhatd(op_id);

    double val = 0.;

    for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {

      const FElement& FK(Vh[k]);
      if(domain != FK.whichDomain()) continue;
      // if(!Vh.isCut(k)) continue; // [if included; only checks cut elements]

      const int kb = Vh.Th(FK.T);
      CutData cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
      const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
      ElementSignEnum the_part = cutK.what_part(domain);

      for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
      it != cutK.element_end(the_part); ++it){

        for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

          QuadraturePoint ip(qf[ipq]); // integration point
          Rd mip = cutK.toK(it, ip);

          val = max(val, fabs(fh.eval(k, mip)-fex(mip, fh.cu, domain)));

        }
      }
    }
    double val_receive = 0;
    MPIcf::AllReduce(val, val_receive, MPI_MAX);

    return val_receive;
}

template<typename M>
double maxNormCut(const ExpressionVirtual& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i, int dom),const GFESpace<M>& Vh) {
  return max(maxNormCut(fh,fex,0,Vh), maxNormCut(fh,fex,1,Vh));
}




/*---------------------------- cut local L2 norm ----------------------------
Integrates only over elements whose edges are not part of the stabilised edges,
thereby neglecting elements messed up by the mcdonald stab
*/
template<typename M>
double L2normCutLoc_2(const ExpressionVirtual& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i, int dom),int domain, const GFESpace<M>& Vh, const MacroElement& bigMac ) {
  typedef Mesh2 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QF QF;
  typedef typename FElement::Rd Rd;
  typedef typename Mesh::Partition Partition;
  typedef typename Mesh::Element Element; // [Needed for neighbor element]
  typedef typename QF::QuadraturePoint QuadraturePoint;
  typedef typename TypeCutData<Rd::d>::CutData CutData;


  const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
  What_d Fop = Fwhatd(op_id);
  double val = 0.;

  vector<int> elements_to_integrate(Vh.NbElement(),1);
  // std::vector<std::pair<int,int>>::iterator it;
  // 0. Mark all elements as to integrate
  // 1. Loop over all elements
  // 2. Remove from elements_to_integrate when isSmall or neighbor
  // 3. Standard L2 error, but only over elements_to_integrate


assert(0);
  // for (auto it=bigMac.edges_to_stabilize.begin(); it<bigMac.edges_to_stabilize.end(); it++) {
  //
  //   int k = it->first;
  //   int ifac = it->second;
  //   //
  //   int ifacn = ifac;
  //   // for(int ifac=0;ifac<3;++ifac){
  //     // int ifacn = ifac;
  //
  //     int k_back  = Vh.Th(Vh[k].T);
  //     int kn_back = Vh.Th.ElementAdj(k_back,ifacn);
  //     if(kn_back == -1) kn_back = k_back;
  //     int kn = Vh.idxElementFromBackMesh(kn_back, domain);
  //     if(kn == -1) kn = k;
  //     elements_to_integrate[k] = 0;
  // }



  for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
    const FElement& FK(Vh[k]);
    if(domain != FK.whichDomain()) continue;

    if(elements_to_integrate[k]==0) continue;

    const int kb = Vh.Th(FK.T);
    CutData cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
    const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
    ElementSignEnum the_part = cutK.what_part(domain);
    double locV = 0;

    for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
    it != cutK.element_end(the_part); ++it){

      const R meas = cutK.mesure(it);
      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = cutK.toK(it, ip);
        const R Cint = meas * ip.getWeight();

        double a = fh.eval(k, mip) - fex(mip, fh.cu, domain);

        val += Cint * a * a;
        locV += Cint * a * a;
      }
    }
    // std::cout << k << "\t" << locV << "\t"<< val << std::endl;
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;
}

template<typename M>
double L2normCutLoc_2(const ExpressionVirtual& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i, int dom), const GFESpace<M>& Vh, const MacroElement& bigMac) {
  return L2normCutLoc_2(fh,fex,0,Vh,bigMac) + L2normCutLoc_2(fh,fex,1,Vh,bigMac);
}

// Several components, ie vector valued L2norm
template<typename M>
double L2normCutLoc( const FunFEM<M>& fh,R (fex)(const typename GFESpace<M>::FElement::Rd, int i, int dom), int c0, int num_comp , const MacroElement& bigMac) {

  const GFESpace<M>& Vh(*fh.Vh);

  double val = 0;
  for(int i=c0;i<num_comp+c0;++i) {
    ExpressionFunFEM<M> ui(fh, i, op_id);
    val += L2normCutLoc_2(ui,fex,Vh,bigMac);
  }
  return sqrt(val);
}

// Only one component
template<typename M>
double L2normCutLoc( const ExpressionVirtual& fh, R (fex)(const typename GFESpace<M>::FElement::Rd, int i, int dom), const GFESpace<M>& Vh, const MacroElement& bigMac) {
  double val = L2normCutLoc_2(fh,fex,Vh,bigMac);
  return sqrt(val);
}
