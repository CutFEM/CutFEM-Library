#include "iceSheet.hpp"
#include "integral_iceSheet.hpp"


IceSheet::IceSheet(const CutFESpace& vh) :
  ShapeOfLinProblem(vh.NbDoF()),
  Solver(),
  Vh(&vh)
{
  rhs.resize(nDoF+1); rhs = 0.0;


  std::cout <<  "rhow \t " << rhow << "\n"
	    <<  "rhoi \t " << rhoi << "\n"
	    <<  "mu \t " << mu_ << "\n"
	    <<  "A \t " << A << "\n"
	    <<  "C \t " << C << "\n"
    	    <<  "gravity \t " << gravity << "\n"
	    <<  "frhs \t " << frhs(R2(100,100))*rhoi << "\n"
    	    <<  "Pw \t " << Pwater(R2(100,-100)) << std::endl;

}



void IceSheet::init(const CutFESpace& vh) {
  Vh = &vh;
  nDoF = Vh->NbDoF();
  rhs.resize(nDoF+1); rhs = 0.0;

}



void IceSheet::solve(const Marker& marker){

  R t0 = CPUtime();
  this->assembly(marker);
  std::cout << "Time Assembly \t \t " << CPUtime() - t0 << std::endl;

  matlab::Export(this->mat, "mat0.dat");
  matlab::Export(this->rhs, "rhs0.dat");


  // matlab::Export(mat, "../../matlabFiles/mat.dat");
  t0 = CPUtime();
  Solver::solve(mat, rhs);
  std::cout << "Time Stokes Solver \t " << CPUtime() - t0 << std::endl;
}



void IceSheet::assembly(const Marker& marker) {

  double t0 = CPUtime();
  assembly_full(marker);
  // std::cout << "Time full assembly \t" << CPUtime() - t0 << std::endl;
  t0 = CPUtime();

  assembly_surface(marker);
  // // std::cout << "Time surface assembly \t" << CPUtime() - t0 << std::endl;
  // // t0 = CPUtime();
  //
  // stabilization(marker);
  // // std::cout << "Time stabilization \t" << CPUtime() - t0 << std::endl;
  // // t0 = CPUtime();
}


void IceSheet::assembly_full(const Marker& marker) {

  IntegralIceSheet integration(*this); //class that compute integrals

  const int nbDoF = (*Vh)[0].NbDoF();
  KNM<R> ML(nbDoF,nbDoF);             // matrix local contribution
  KN<R> VL(nbDoF);                    // vector local contribution
  KN<R> VP(nbDoF);                    // vector for lagrange multiplier

// loop over elements
  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

    const FElement& FK((*Vh)[k]);             // the element
    const int k_backMesh = Vh->Th(FK.T);      // index in backMesh
    CutData cutData(marker.getCutData(k_backMesh));     // get the cut data
    const Partition& cutK =  Partition(FK.T, cutData);  // build the cut

    integration.locFullAssembly(FK, ML, VL, VP, cutK);  // compute integrals

    // fill the global matrix/vector
    // FK(i) gives the loc2glob index mapping
    for (int i=0;i<nbDoF;++i){
      for (int j=0;j<nbDoF;++j) {
    	(*this)(FK(i),FK(j)) += ML(i,j);
      }
      (*this)(FK(i)) += VL(i);
    }

    // set int_Omega p dx = 0
    for(int jp=FK.dfcbegin(Rd::d); jp<FK.dfcend(Rd::d); ++jp){
      (*this)(nDoF,FK(jp)) += VP(jp);
      (*this)(FK(jp),nDoF) += VP(jp);
    }


  } //end loop over elements
  (*this)(nDoF ,nDoF) = 0;
  (*this)(nDoF) = 0;

  // plot.close();
}





void IceSheet::assembly_surface(const Marker& marker) {

  IntegralIceSheet integration(*this);
  // kappa1 = mu2 / (mu1+mu2);
  // kappa2 = mu1 / (mu1+mu2);

  const int nbDoF = (*Vh)[0].NbDoF();
  KNM<R> ML(nbDoF,nbDoF);
  KN<R> VL(nbDoF);

  const Marker& interface(marker);
  for(int iface=interface.first_element(); iface<interface.last_element();
      iface+= interface.next_element()) {

    const FaceMarker& face = marker.getFace(iface);  // the face
    const int kb = face.k; // idx on backMesh
    const int k = Vh->idxElementFromBackMesh(kb, 0);
    const FElement& FK((*Vh)[k]);

    integration.locSurfAssembly(FK, iface, ML, VL);

    for (int i=0;i<nbDoF;++i){
      for (int j=0;j<nbDoF;++j) {
    	(*this)(FK(i),FK(j)) += ML(i,j);
      }
      (*this)(FK(i)) += VL(i);
    }

  }

}







void IceSheet::stabilization(const Marker & marker) {

  IntegralIceSheet integration(*this);
  integration.initStab(Edu, Edp);
  const int nbDoF = (*Vh)[0].NbDoF();

  KNM<R> ML(nbDoF,nbDoF);
  KNM<R> ML_M(nbDoF,nbDoF);
  KNM<R> MLN(nbDoF,nbDoF);

  const Marker& interface(marker);
  for(int iface=interface.first_element(); iface<interface.last_element();
      iface+= interface.next_element()) {

    const FaceMarker& face = marker.getFace(iface);  // the face
    const int kb = face.k; // idx on backMesh

//     for(int the_domain=0;the_domain<=1;++the_domain) { //loop over sub_domains
    const int k = Vh->idxElementFromBackMesh(kb, 0);
    const FElement& FK((*Vh)[k]);

    for(int ifac = 0; ifac < Element::nea; ++ifac) {         //loop over the faces

      int kn = integration.locStabilization(ML, ML_M, MLN, FK, ifac);
      if(kn==-1) continue;
      const FElement& FKn((*Vh)[kn]);

      for (int i=0;i<nbDoF;++i){
	for (int j=0;j<nbDoF;++j) {
	  (*this)(FK(i),FK(j)) += ML(i,j);
	  (*this)(FK(i),FKn(j)) -= ML_M(i,j);
	}
      }

// 	levelSet.Fvec(Vh->Th(FKn.T), loc_ls);
// 	if(!changeSign(loc_ls)) {            // neighbor not cut
// 	  for (int i=0;i<nbDoF;++i){
// 	    for (int j=0;j<nbDoF;++j) {
// 	      (*this)(FKn(i),FKn(j)) += MLN(i,j);
// 	      (*this)(FKn(i),FK(j))  -= ML_M(j,i);
// 	    }
// 	  }
// 	}
//       }
    }
  }

}
