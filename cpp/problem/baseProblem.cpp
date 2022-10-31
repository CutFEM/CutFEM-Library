#include "baseProblem.hpp"


template<>
void BaseCutFEM<Mesh2>::addInterfaceRidgeContribution(const ListItemVF<Rd::d>& VF, const Interface<Mesh2>& interface, int ifac, const TimeSlab* In, int itq, double cst_time) {

  // GET IDX ELEMENT CONTAINING FACE ON backMes
  const Mesh2& Kh(interface.get_mesh());
  const int kb = interface.idxElementOfFace(ifac);
  const typename Interface<Mesh2>::Face& face0 = interface[ifac];
  // const Element & K(interface.get_element(kb));
  // double measK = K.measure();
  // double h     = K.get_h() ;
  // double meas  = interface.measure(ifac);
  const Rd normal0(-interface.normal(ifac));

  for(int inode=0;inode<2;++inode) {

    int idx_node = face0[inode];
    int e0 = interface.edge_of_node_[idx_node];
    int e1 = e0;
    int knb = Kh.ElementAdj(kb,e1);
    const Rd mip = interface(idx_node);

    if(kb > knb) continue;

    // GET NEIGHBOR FACE AND CORRESPONDING INDEX OF THE NODE
    int iface1 = interface.face_of_element_.find(knb)->second;
    const typename Interface<Mesh2>::Face& face1 = interface[iface1];
    int inode1 = (interface.edge_of_node_[face1[0]]==e1)? 0 : 1;
    const Rd normal1(-interface.normal(iface1));


    Rd conormal0 = interface(face0[inode] ) - interface(face0[(inode +1)%2]);
    Rd conormal1 = interface(face1[inode1]) - interface(face1[(inode1+1)%2]);
    conormal0 /= conormal0.norm();
    conormal1 /= conormal1.norm();

    for(int l=0; l<VF.size();++l) {
      // if(!VF[l].on(domain)) continue;

      const int kbv = VF[l].onWhatElementIsTestFunction(kb, knb);
      const int kbu = VF[l].onWhatElementIsTrialFunction(kb, knb);

      // FINITE ELEMENT SPACES && ELEMENTS
      const FESpace& Vhv(VF.get_spaceV(l));
      const FESpace& Vhu(VF.get_spaceU(l));

      bool same = (VF.isRHS() || (&Vhu == &Vhv && kbu == kbv) );

      std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kbv, VF[l].get_domain_test_function());
      std::vector<int> idxU = (same)?idxV : Vhu.idxAllElementFromBackMesh(kbu, VF[l].get_domain_trial_function());
      int kv = VF[l].onWhatElementIsTestFunction (idxV);
      int ku = (same)? kv : VF[l].onWhatElementIsTrialFunction(idxU);

      const FElement& FKu(Vhu[ku]);
      const FElement& FKv(Vhv[kv]);
      this->initIndex(FKu, FKv);
      int domu = FKu.get_domain();
      int domv = FKv.get_domain();

      // DEFINE THE GOOD NORMALS AND CONORMALS
      const Rd normalv = VF[l].onWhatElementIsTestFunction( normal0,normal1);
      const Rd normalu = (same)? normalv : VF[l].onWhatElementIsTrialFunction(normal0,normal1);
      const Rd conormalv = VF[l].onWhatElementIsTestFunction(conormal0,conormal1);
      const Rd conormalu = (same)? conormalv : VF[l].onWhatElementIsTrialFunction(conormal0,conormal1);

      // BF MEMORY MANAGEMENT -
      int lastop = getLastop(VF[l].du, VF[l].dv);
      RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop);
      RNMK_ fu(this->databf_+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop);
      What_d Fop = Fwhatd(lastop);

      // COMPUTE COEFFICIENT && NORMAL
      double coef = cst_time;//VF[l].computeCoefInterface(h,meas,measK,measK,domu, domv) ;
      coef *= VF[l].computeCoefFromNormal(normalu, normalv);
      coef *= VF[l].computeCoefFromConormal(conormalu, conormalv);

      // EVALUATE THE BASIS FUNCTIONS
      FKv.BF(Fop,FKv.T.mapToReferenceElement(mip), fv);
      if(!same) FKu.BF(Fop, FKu.T.mapToReferenceElement(mip), fu);

      coef *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kbu,kbv), std::make_pair(domu, domv), mip, 0., std::make_pair(conormalu, conormalv));
      coef *= VF[l].c;

      if( In ){
        if(VF.isRHS()) this->addToRHS(   VF[l], *In, FKv, fv, coef);
        else           this->addToMatrix(VF[l], *In, FKu, FKv, fu, fv, coef);
      }
      else {
        if(VF.isRHS()) this->addToRHS(   VF[l], FKv, fv, coef);
        else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, coef);
      }
    }

  }
}




template class CutFEM<Mesh2> ;
template class FEM<Mesh2> ;
template class CutFEM<MeshQuad2> ;
template class FEM<MeshQuad2> ;
