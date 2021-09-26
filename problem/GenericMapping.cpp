#include "GenericMapping.hpp"


R Mapping2::errorInfty(const Interface2& interface, R (*fun)(const R2, int i)) {

  std::ofstream plot1;
  std::ofstream plot2;
  // plot1.open("initPT.dat", std::ofstream::out);
  plot2.open("defPTP3.dat", std::ofstream::out);

  typedef typename FElement::RdHatBord RdHatBord;
  Fun_h deform(*Vh, deformation_);
  R errMax = 0;

  const int numberOfPoint = 10;
  R deltaX = 1.0/ numberOfPoint;

  for(int iface=0;iface<interface.nbElement();++iface) {
    const typename Interface::FaceIdx& face = interface[iface];
    int k = interface.idxElementOfFace(iface);
    if(interface.backMesh != &(Vh->Th)) k = Vh->idxElementFromBackMesh(k);

    for(int i = 0; i < numberOfPoint; ++i) {
      RdHatBord ip = i*deltaX;
      Rd mip = interface.mapToFace(face, ip);
      Rd defPt;
      for(int c=0;c<Rd::d;++c) defPt[c] = mip[c] + deform.eval(k,mip,c, op_id);
      const R val = fabs(fun(defPt,0));
      errMax = std::max(errMax, val );
      plot1 << mip << std::endl;
      plot2 << defPt<< std::endl;
     }
     // std::cout << interface.mapToFace(face, ip) << "\t" << dd << "\t" << fabs(fun(interface.mapToFace(face, ip),0)) << "\n";
     // std::cout << mip << "\t" << fabs(fun(mip,0))  << std::endl ;
     // std::cout << defPt << "\t" << val << std::endl << std::endl << std::endl;

  }

  R err = errMax;
  std::cout << " Error L_infty = " << err << std::endl;
  return err;
}




class IdMapping2 : public GenericMapping<Mesh2> {
  typedef GenericMapping<Mesh2> Mapping;
 public :
 IdMapping2() : Mapping() { }

 virtual int idxElementFromBackMesh(int kb) const {
   return kb;
 }
 virtual void computeInverseJacobian(int k, Rd mip, KNM_<double>& invJ) const {
   for(int i=0;i<Rd::d;++i){
     for(int j=0;j<Rd::d;++j){
       invJ(i,j) = (i==j);
     }
   }
 }

 virtual void transform(KNMK_<double>& bf, const KNM_<double>& invJ) const {}

 virtual Rd transform(const int k, const Rd x) const { return x;}

};

class IdMapping3 : public GenericMapping<Mesh3> {
  typedef GenericMapping<Mesh3> Mapping;
 public :
 IdMapping3() : Mapping() { }

 virtual int idxElementFromBackMesh(int kb) const {
   return kb;
 }
 virtual void computeInverseJacobian(int k, Rd mip, KNM_<double>& invJ) const {
   for(int i=0;i<Rd::d;++i){
     for(int j=0;j<Rd::d;++j){
       invJ(i,j) = (i==j);
     }
   }
 }

 virtual void transform(KNMK_<double>& bf, const KNM_<double>& invJ) const {}

 virtual Rd transform(const int k, const Rd x) const { return x;}

};



static IdMapping2 IdMapping_2d;
static IdMapping3 IdMapping_3d;
GenericMapping<Mesh2> & IdMapping2d(IdMapping_2d);
GenericMapping<Mesh3> & IdMapping3d(IdMapping_3d);
template<> GenericMapping<Mesh2> & DataMapping<Mesh2>::Id = IdMapping_2d;
template<> GenericMapping<Mesh3> & DataMapping<Mesh3>::Id = IdMapping_3d;




// R Mapping3::errorInfty(const Interface3& interface, R (*fun)(const R3)) const { return 0.0;}

//   typedef typename FElement::QFB  QFB;

//   const QFB & qfb(QuadratureFormular_T_9);

//   // if(interface.nbElement() != Vh.Th.nt) {
//   //   std::cout << interface.nbElement() << "\t" << Vh.Th.nt << std::endl;
//   //   std::cout << " Error need to be computed on the cut mesh created with the same interface " << std::endl;
//   //   std::cout << " The error \"errorInfty \" has not been computed " << std::endl;
//   //   return 0.0;
//   // }

//   R errMax = 0;
//   // const int nve = Vh[0].NbOfNodes();
//   // R loc_lset[nve];

//   // const int numberOfPoint = 10;
//   // R deltaX = 1.0/ numberOfPoint;

//   // const Mesh & Th = Vh.Th;                               // the mesh

//   for(Uint iface=0;iface<interface.nbElement();++iface) {

//     const int k = Vh.idxElementInLocMesh(interface.idxElementOfFace(iface) );
//     const FElement FK(Vh[k]);
//     const typename Interface::FaceIdx& face = interface[iface];  // the face

//   //   for(int i = 0; i < numberOfPoint; ++i) {
//     for(int l=0;l<qfb.n;++l) {
//       const typename QFB::QuadraturePoint ip(qfb[l]);                     // quadrature rule
//       Rd mip = interface.mapToFace(face, ip);

//       mip= deformNode(FK,mip);

//       const R val = fabs(fun(mip));
//       errMax = std::max(errMax, val );
//     }
//   }

//   R err = errMax;
//   return err;
// }


// R Mapping3::errorL2(const Interface3& interface, R (*fun)(const R3)) const {

//   // typedef typename FElement::RdHatBord RdHatBord;
//   typedef typename FElement::QFB  QFB;

//   const QFB & qfb(QuadratureFormular_T_9);
//   // const int nvx = Vh.global.NbOfNodes;
//   const int nbDoF = Vh[0].nbDoF();                    // nb of local dof
//   const What_d whatd = Fop_D0;                 // only phi
//   KNMK<R> w(nbDoF,1,op_id+1);                  // tensor that save bf

//   R errL2 = 0;

//   for(Uint iface=0;iface<interface.nbElement();++iface) {

//     const int k = Vh.idxElementInLocMesh(interface.idxElementOfFace(iface) );
//     const FElement FK(Vh[k]);
//     const typename Interface::FaceIdx& face = interface[iface];  // the face
//     const Rd dx = interface.computeDx(face);            // integration constant



//   //   for(int i = 0; i < numberOfPoint; ++i) {
//     for(int l=0;l<qfb.n;++l) {
//       const typename QFB::QuadraturePoint ip(qfb[l]);                     // quadrature rule
//       Rd mip = interface.mapToFace(face, ip);
//       const R Cint = dx.norm() * ip.a;


//       FK.BF(whatd, FK.T.toKref(mip), w);                // compute basic funtions
//       mip = deformNode(FK,mip);

//       R uk=0;
//       // for(int i=0;i<nbDoF;++i){
//       // 	uk  += w(i,0,op_id) * fun(mip);
//       // 	//	uk  += fun(mip);
//       // }
//       uk  = std::fabs(fun(mip));
//       errL2 += Cint * uk*uk;
//     }
//   }

//   return sqrt(errL2);
// }
