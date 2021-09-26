// #ifndef REINIT_TPP_
// #define REINIT_TPP_

template<class Mesh>
void CReinitialization<Mesh>::perform(Fun& levelSet_P1){

  if(ON_OFF == "ON"){
    const FESpace& Vh(*levelSet_P1.Vh);
    Interface interface(Vh.Th, levelSet_P1.v);
    double areaInitial = computeArea(*levelSet_P1.Vh, interface);
    // std::cout << " area before reinitialization ->\t" << areaInitial  << std::endl;

    Reinitialization<Mesh> reinitialization(levelSet_P1, interface, number_iteration, epsilon_diffusion, ODE_method, dt);

    Interface interface0(Vh.Th, levelSet_P1.v);
    double areaBefore = computeArea(*levelSet_P1.Vh, interface0);
    // std::cout << " area before correction      ->\t" << areaBefore  << std::endl;

    if(mass_correction == "ON") {

      R delta = correction(levelSet_P1, areaInitial);

      Interface interface1(Vh.Th, levelSet_P1.v);
      double areaFinal = computeArea(*levelSet_P1.Vh, interface1);
      // std::cout << " area after correction      ->\t" << areaFinal  << std::endl;
    }
  }
}

template<class Mesh>
void CReinitialization<Mesh>::perform(Fun& levelSet_k, Fun& levelSet_P1){

  if(ON_OFF == "ON"){
    const FESpace& Vh(*levelSet_P1.Vh);
    Interface interface(Vh.Th, levelSet_P1.v);
    double areaInitial = computeArea(*levelSet_P1.Vh, interface);
    // std::cout << " area before reinitialization ->\t" << areaInitial  << std::endl;

    Reinitialization<Mesh> reinitialization(levelSet_k, interface, number_iteration, epsilon_diffusion, ODE_method, dt);
    projection(levelSet_k, levelSet_P1);

    Interface interface0(Vh.Th, levelSet_P1.v);
    double areaBefore = computeArea(*levelSet_P1.Vh, interface0);
    // std::cout << " area before correction      ->\t" << areaBefore  << std::endl;

    if(mass_correction == "ON") {

      R delta = correction(levelSet_P1, areaInitial);

      Interface interface1(Vh.Th, levelSet_P1.v);
      double areaFinal = computeArea(*levelSet_P1.Vh, interface1);
      // std::cout << " area after correction      ->\t" << areaFinal  << std::endl;

      levelSet_k.v += delta;
    }
  }
}

template<class Mesh>
double CReinitialization<Mesh>::correction(Fun& levelSet, const double areaInitial) {

  const FESpace& Vh(*levelSet.Vh);
  double dA = precision_correction;

  // find first Step
  Interface newInterface0(Vh.Th, levelSet.v);
  double area0 = computeArea(Vh, newInterface0);

  const double diffInitial = area0 - areaInitial;
  double diff0 = area0 - areaInitial;

  if (fabs(diff0) < dA) { return 0; }
  double dz0 = 0;
  double dz1 = diff0; //approximation
  double diff1 = diff0;
  int nn = 0;
  while(diff1*diff0 > 0){
    levelSet.v += dz1;
    Interface newInterface1(Vh.Th, levelSet.v);
    double area1 = computeArea(Vh, newInterface1);
    diff1 = area1 - areaInitial;
    nn += 1;
    if (fabs(diff1) < dA) { return dz1; }

  }
  dz1 = nn*dz1;
  levelSet.v -= dz1;



  for(int i=0;i<max_iteration_correction;++i){
    double dz2 = 0.5*(dz0+dz1);
    levelSet.v += dz2;
    Interface newInterface1(Vh.Th, levelSet.v);
    double area2 = computeArea(Vh, newInterface1);
    double diff2 = area2 - areaInitial;
    if(fabs(diff2) < dA) {
      std::cout << " Mass correction error \t" << diff2 << std::endl;
      return dz2;
    }
    if(i == max_iteration_correction-1 && fabs(diffInitial) > fabs(diff2)) {
      std::cout << " Correction of levelSet not reached the required precision ("<< dA << ")" << std::endl;
      std::cout << " Mass correction error \t" << diff2 << std::endl;
      return dz2;
    }
    levelSet.v -= dz2;
    if(diff0*diff2 < 0) {dz1 = dz2; diff1 = diff2;}
    else{ dz0 = dz2; diff0 = diff2;}

  }
  std::cout << " Correction of levelSet failed!!" << std::endl;

}

template<class Mesh>
double CReinitialization<Mesh>::computeArea(const FESpace& Vh, const Interface& gamma){
  typedef typename FESpace::FElement FElement;
  typedef typename Mesh::Partition Partition;
  typedef typename FElement::Rd Rd;
  typedef typename TypeCutData<Rd::d>::CutData CutData;
  double area = 0.;
  int domain = 1;
  for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
  // for(int k=0; k<Vh.last_element(); k+=1) {

    const FElement& FK(Vh[k]);
    const int kb = Vh.Th(FK.T);
    CutData cutData(gamma.getCutData(kb));     // get the cut data
    const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
    ElementSignEnum the_part = cutK.what_part(domain);

    for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
    it != cutK.element_end(the_part); ++it){

      int signElement = cutK.whatSign(it);
      area += (signElement < 0)? cutK.mesure(it) : 0.;

    }
  }
  double sum_area = 0;
  MPIcf::AllReduce(area, sum_area, MPI_SUM);
  return sum_area;
}




// #endif
