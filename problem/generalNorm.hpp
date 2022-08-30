#ifndef _GENERAL_NORM_
#define _GENERAL_NORM_

#include "../FESpace/FESpace.hpp"
#include "itemVF.hpp"


template<typename Mesh>
double integralCut(const ListItemVF<Mesh::Rd::d>& VF, const FunFEM<Mesh>& f,const FunFEM<Mesh>& g, const ActiveMesh<Mesh>& Th) {

  // TYPEDEF NAMES OF TEMPLATE CLASSES
  typedef GFESpace<Mesh> FESpace;
  typedef typename Mesh::Element  Element;
  typedef typename FESpace::FElement  FElement;
  typedef typename FElement::QF QF;

  const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));

  double val = 0.;

  for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {

    const Cut_Part<Element> cutK(Th.get_cut_part(k,0));
    const Element & K(Th[k]);
    int domain  = Th.get_domain_element(k);
    double measK = K.measure();
    double h = K.hElement();
    int kb = Th.idxElementInBackMesh(k);

    // LOOP OVER PART IN THE DOMAIN
    for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){

      const R meas = cutK.measure(it);

      // LOOP OVER TERMS IN THE WEAK FORMULATION
      for(int l=0; l<VF.size();++l) {

        double coef  = 1.;//VF[l].computeCoefElement(h,meas,measK, meas, domain);

        for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

          typename QF::QuadraturePoint ip(qf[ipq]);
          typename Mesh::Rd mip = cutK.mapToPhysicalElement(it, ip);

          double Cint = meas * ip.getWeight() * coef;
          double val_fh = f.evalOnBackMesh(kb, domain, mip, VF[l].cu, VF[l].du);
          double val_gh = g.evalOnBackMesh(kb, domain, mip, VF[l].cv, VF[l].dv);

          val += Cint * val_fh * val_gh * VF[l].c;
        }
      }
    }
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);
  return sqrt(val_receive);
}

// template<typename Mesh>
// double integralCut(const ExpressionVirtual& f,const ExpressionVirtual& g, const ActiveMesh<Mesh>& Th) {
//
//   // TYPEDEF NAMES OF TEMPLATE CLASSES
//   typedef CutFESpace<Mesh> FESpace;
//   typedef typename FESpace::FElement  FElement;
//   typedef typename FElement::QF QF;
//
//   const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
//
//   double val = 0.;
//
//   for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {
//
//     const Cut_Part<Element> cutK(Th.get_cut_part(k));
//     const Element & K(Th[k]);
//     int domain  = Th.get_domain(k);
//     double measK = K.measure();
//     int kb = Th.idxElementInBackMesh(k);
//
//     // LOOP OVER PART IN THE DOMAIN
//     for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
//
//       const R meas = cutK.mesure(it);
//
//       // LOOP OVER TERMS IN THE WEAK FORMULATION
//       for(int l=0; l<VF.size();++l) {
//
//         double coef  = VF[l].computeCoefElement(h,meas,measK, meas, domain);
//
//         for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//
//           typename QF::QuadraturePoint ip(qf[ipq]);
//           typename Mesh::Rd mip = cutK.mapToPhysicalElement(it, ip);
//
//           double Cint = meas * ip.getWeight() * coef;
//           double val_fh = f.evalOnBackMesh(kb, mip, VF[l].cu, VF[l].du, domain);
//           double val_gh = g.evalOnBackMesh(kb, mip, VF[l].cv, VF[l].dv, domain);
//
//           val += Cint * val_fh * val_gh * VF[l].c;
//         }
//       }
//     }
//   }
//   double val_receive = 0;
//   MPIcf::AllReduce(val, val_receive, MPI_SUM);
//   return sqrt(val_receive);
// }


// template<typename Mesh>
// double integral(const ListItemVF<Mesh::Rd::d>& VF, const FunFEM<Mesh>& f,const FunFEM<Mesh>& g, const GenericInterface<Mesh>& interface) {
//
//   // TYPEDEF NAMES OF TEMPLATE CLASSES
//   typedef CutGFESpace<Mesh> FESpace;
//   typedef typename FESpace::FElement  FElement;
//   typedef typename TypeCutData<Mesh::Rd::d>::CutData CutData;
//   typedef typename Mesh::Partition Partition;
//   typedef typename FElement::QFB QFB;
//   typedef typename FElement::RdHatBord RdHatBord;
//   typedef GenericInterface<Mesh> Interface;
//   typedef typename Mesh::Rd Rd;
//
//
//   const QFB& qf(*QF_Simplex<typename FElement::RdHatBord>(5));
//   double val = 0.;
//
//   for(int iface=interface.first_element(); iface<interface.last_element(); iface+=interface.next_element()) {
//
//     // COMPUTE THE PARTITION OF THE ELEMENT
//     const int kb = interface.idxElementOfFace(iface);   // idx on
//     const typename Interface::FaceIdx& face = interface[iface];  // the face
//     const R meas = interface.computeDx(face).norm();
//     const Rd normal(-interface.normal(iface));
//     const double h = meas;
//       // LOOP OVER TERMS IN THE WEAK FORMULATION
//       for(int l=0; l<VF.size();++l) {
//
//         double measK = (*interface.backMesh)[kb].mesure();
//
//         R coef = VF[l].computeCoefInterface(h,meas,measK );
//         double Cnormal = VF[l].getCoef(normal);
//
//         for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//
//           typename QFB::QuadraturePoint ip(qf[ipq]);
//           Rd mip = interface.mapToFace(face,(RdHatBord)ip);
//
//           double Cint = meas * ip.getWeight() * coef* VF[l].c * Cnormal;
//
//           double val_fh = f.evalOnBackMesh(kb, mip, VF[l].cu, VF[l].du, VF[l].domu);
//           double val_gh = g.evalOnBackMesh(kb, mip, VF[l].cv, VF[l].dv, VF[l].domv);
//           val += Cint * val_fh * val_gh ;
//         }
//       }
//     }
//
//   double val_receive = 0;
//   MPIcf::AllReduce(val, val_receive, MPI_SUM);
//   return sqrt(val_receive);
// }
//
//




















#endif
