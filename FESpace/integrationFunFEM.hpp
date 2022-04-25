#ifndef INTEGRATION_FUNFEM_HPP_
#define INTEGRATION_FUNFEM_HPP_

#include "expression.hpp"
// #include "../common/Interface2dn.hpp"
#include "macroElement.hpp"


//
// /*
//           Space integration of f(x)
// */
// template<typename M>
// double integral(FunFEM<M>& fh, int cu, int domain = -1){
//
//   typedef M Mesh;
//   typedef GFESpace<Mesh> FESpace;
//   typedef typename FESpace::FElement FElement;
//   typedef typename FElement::QF QF;
//   typedef typename FElement::Rd Rd;
//   typedef typename Mesh::Partition Partition;
//   typedef typename QF::QuadraturePoint QuadraturePoint;
//   // typedef typename TypeCutData<Rd::d>::CutData CutData;
//   const FESpace& Vh(*fh.Vh);
//   const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
//
//   What_d Fop = Fwhatd(op_id);
//   double val = 0.;
//
//   for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
//
//     const FElement& FK(Vh[k]);
//
//     if(domain == -1){
//       const R meas = FK.getMeasure();
//       for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//
//         QuadraturePoint ip(qf[ipq]); // integration point
//         Rd mip = FK.map(ip);
//         const R Cint = meas * ip.getWeight();
//         val += Cint * fh.eval(k, mip, cu, op_id);
//       }
//     }
//     else if(domain == FK.whichDomain()) {
//       // if(domain != FK.whichDomain()) continue;
//
//       const int kb = Vh.Th(FK.T);
//       CutData cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
//       const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
//       ElementSignEnum the_part = cutK.what_part(domain);
//
//       for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
//       it != cutK.element_end(the_part); ++it){
//
//         const R meas = cutK.mesure(it);
//
//         for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//
//           QuadraturePoint ip(qf[ipq]); // integration point
//           Rd mip = cutK.toK(it, ip);
//           const R Cint = meas * ip.getWeight();
//
//           val += Cint * fh.eval(k, mip, cu, op_id);
//
//         }
//       }
//     }
//   }
//   double val_receive = 0;
//   MPIcf::AllReduce(val, val_receive, MPI_SUM);
//
//   return val_receive;
// }
//
// template<typename M>
// double integral(const ExpressionFunFEM<M>& fh, int domain) {
//
//   typedef M Mesh;
//   typedef GFESpace<Mesh> FESpace;
//   typedef typename FESpace::FElement FElement;
//   typedef typename FElement::QF QF;
//   typedef typename FElement::Rd Rd;
//   typedef typename Mesh::Partition Partition;
//   typedef typename QF::QuadraturePoint QuadraturePoint;
//   typedef typename TypeCutData<Rd::d>::CutData CutData;
//
//   const FESpace& Vh(fh.getSpace());
//   const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
//
//   What_d Fop = Fwhatd(op_id);
//
//   double val = 0.;
//
//   for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
//
//     const FElement& FK(Vh[k]);
//     if(domain == -1){
//       const R meas = FK.getMeasure();
//       for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//
//         QuadraturePoint ip(qf[ipq]); // integration point
//         Rd mip = FK.map(ip);
//         const R Cint = meas * ip.getWeight();
//         val += Cint * fh.eval(k, mip);
//       }
//     }
//     else if(domain == FK.whichDomain()) {
//       // if(domain != FK.whichDomain()) continue;
//
//       const int kb = Vh.Th(FK.T);
//       CutData cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
//       const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
//       ElementSignEnum the_part = cutK.what_part(domain);
//
//       for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
//       it != cutK.element_end(the_part); ++it){
//
//         const R meas = cutK.mesure(it);
//
//         for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//
//           QuadraturePoint ip(qf[ipq]); // integration point
//           Rd mip = cutK.toK(it, ip);
//           const R Cint = meas * ip.getWeight();
//
//           val += Cint * fh.eval(k, mip);
//
//         }
//       }
//     }
//   }
//   double val_receive = 0;
//   MPIcf::AllReduce(val, val_receive, MPI_SUM);
//
//   return val_receive;
// }
//
// template<typename M>
// double integral(const ExpressionFunFEM<M>& fh) {
//   if(fh.getSpace().isCut()) return integral(fh,0) + integral(fh,1);
//   else return integral(fh,-1) ;
// }

template<typename M>
double integral( const ActiveMesh<M>& Th, const FunFEM<M>& fh, int c0) {
  int nb_dom = Th.get_nb_domain();
  double val = 0.;
  ExpressionFunFEM<M> ui(fh, c0, op_id);
  for(int i=0;i<nb_dom;++i) {
    val += integral(Th, ui, i);
  }
  return val;
}
template<typename M>
double integral( const ActiveMesh<M>& Th, const ExpressionVirtual& fh) {
  int nb_dom = Th.get_nb_domain();
  double val = 0.;
  for(int i=0;i<nb_dom;++i) {
    val += integral(Th, fh, i);
  }
  return val;
}
template<typename M>
double integral(const ActiveMesh<M>& Th, const ExpressionVirtual& fh, int domain) {
  typedef M Mesh;
  typedef typename Mesh::Element Element;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QF QF;
  typedef typename FElement::Rd Rd;
  typedef typename QF::QuadraturePoint QuadraturePoint;

  const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));

  double val = 0.;

  for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {

    if(domain != Th.get_domain_element(k)) continue;

    const Cut_Part<Element> cutK(Th.get_cut_part(k));
    int kb = Th.idxElementInBackMesh(k);
    for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){

      const R meas = cutK.measure(it);

      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = cutK.mapToPhysicalElement(it, ip);
        const R Cint = meas * ip.getWeight();

        val += Cint * fh.evalOnBackMesh(kb, domain, mip);
      }
    }
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;
}



template<typename M>
double integral(FunFEM<M>& fh, const Interface<M>& interface, int cu) {
  return integral(fh, interface, cu, -1);
}
template<typename M>
double integral(FunFEM<M>& fh, const Interface<M>* interface, int cu) {
  return integral(fh, *interface, cu, -1);
}
template<typename M>
double integral(FunFEM<M>& fh, const Interface<M>& interface, int cu, double t) {

  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QFB QFB;
  typedef typename FElement::Rd Rd;
  typedef typename QFB::QuadraturePoint QuadraturePoint;

  if(t > -Epsilon && fh.In){
    assert(fh.In->Pt(0) <=t && t <= fh.In->Pt(1));
  }
  if(t < -Epsilon && fh.In) {
    t = fh.In->Pt(0);
    std::cout << " Use default value In(0) \t -> " << t << std::endl;
  }

  const QFB& qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
  double val = 0.;

  for(int iface=interface.first_element(); iface<interface.last_element(); iface+=interface.next_element()) {
    const int kb = interface.idxElementOfFace(iface);   // idx on backMesh
    const R meas = interface.measure(iface);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = interface.mapToPhysicalFace(iface,(typename FElement::RdHatBord)ip);
      const R Cint = meas * ip.getWeight();

      val += Cint * fh.evalOnBackMesh(kb, 0, mip, t, cu, 0, 0);

    }
  }

  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;
}

// template<typename M>
// double integralSurf(const ExpressionVirtual& fh, const GFESpace<M>& Vh, int it = 0, double t = 0.) {
//
//   typedef M Mesh;
//   typedef GenericInterface<Mesh> Interface;
//   typedef GFESpace<Mesh> FESpace;
//   typedef typename FESpace::FElement FElement;
//   typedef typename FElement::QFB QFB;
//   typedef typename FElement::Rd Rd;
//   typedef typename Mesh::Partition Partition;
//   typedef typename QFB::QuadraturePoint QuadraturePoint;
//
//   const Interface& interface(Vh.getInterface(it));
//   const QFB& qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
//
//   What_d Fop = Fwhatd(op_id);
//   double val = 0.;
//
//   for(int iface=interface.first_element(); iface<interface.last_element(); iface+=interface.next_element()) {
//
//     const int kb = interface.idxElementOfFace(iface);   // idx on backMesh
//     const typename Interface::Face& face = interface.getFace(iface);  // the face
//     const R meas = interface.computeDx(face).norm();
//
//     int k = Vh.idxElementFromBackMesh(kb,0);
//     const FElement& FK(Vh[k]);
//
//     for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
//
//       QuadraturePoint ip(qfb[ipq]); // integration point
//       const Rd mip = interface.mapToFace(face,(typename FElement::RdHatBord)ip);
//       const R Cint = meas * ip.getWeight();
//
//       val += Cint * fh.eval(k, mip, t);
//
//     }
//   }
//   double val_receive = 0;
//   MPIcf::AllReduce(val, val_receive, MPI_SUM);
//
//   return val_receive;
// }


template<typename M>
double integral(FunFEM<M>& fh, const TimeSlab& In, const TimeInterface<M>& gamma, int cu) {

  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QFB QFB;
  typedef typename FElement::Rd Rd;
  typedef typename QFB::QuadraturePoint QuadraturePoint;


  const QFB& qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
  double val = 0.;

  for(int it = 0; it<gamma.size(); ++it){
    const Interface<M>& interface(*gamma(it));
    const QuadratureFormular1d& qTime( gamma.get_quadrature_time());
    GQuadraturePoint<R1> tq((qTime)[it]);
    const double t = In.mapToPhysicalElement(tq);


    for(int iface=interface.first_element(); iface<interface.last_element(); iface+=interface.next_element()) {

      const int kb = interface.idxElementOfFace(iface);   // idx on backMesh
      const R meas = interface.measure(iface);

      for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qfb[ipq]); // integration point
        const Rd mip = interface.mapToPhysicalFace(iface,(typename FElement::RdHatBord)ip);
        const R Cint = meas * ip.getWeight() * In.T.mesure()* tq.a ;
        val += Cint * fh.evalOnBackMesh(kb, 0, mip, t, cu, 0, 0);

      }
    }
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;
}




/*
          Space integration of f(x,t)
*/
// template<typename M>
// double integralSurf(FunFEM<M>& fh, const TimeInterface<M>& interface, int cu, int it, double t) {
//
//   typedef M Mesh;
//   typedef GenericInterface<Mesh> Interface;
//   typedef GFESpace<Mesh> FESpace;
//   typedef typename FESpace::FElement FElement;
//   typedef typename FElement::QFB QFB;
//   typedef typename FElement::Rd Rd;
//   typedef typename Mesh::Partition Partition;
//   typedef typename QFB::QuadraturePoint QuadraturePoint;
//
//   assert(fh.In);
//
//   const FESpace& Vh(*fh.Vh);
//   const Interface& interface(Vh.getInterface(it));
//   const QFB& qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
//
//
//   What_d Fop = Fwhatd(op_id);
//   double val = 0.;
//
//   for(int iface=interface.first_element(); iface<interface.last_element(); iface+=interface.next_element()) {
//
//     const int kb = interface.idxElementOfFace(iface);   // idx on backMesh
//     const typename Interface::Face& face = interface.getFace(iface);  // the face
//     const R meas = interface.computeDx(face).norm();
//
//     int k = Vh.idxElementFromBackMesh(kb,0);
//     const FElement& FK(Vh[k]);
//
//     for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
//
//       QuadraturePoint ip(qfb[ipq]); // integration point
//       const Rd mip = interface.mapToFace(face,(typename FElement::RdHatBord)ip);
//       const R Cint = meas * ip.getWeight();
//
//       val += Cint * fh.eval(k, mip, t, cu, op_id);
//
//     }
//   }
//
//   double val_receive = 0;
//   MPIcf::AllReduce(val, val_receive, MPI_SUM);
//
//   return val_receive;
// }

// template<typename M>
// double integralSurf(const ExpressionVirtual& fh, const GFESpace<M>& Vh, int it, double t) {
//
//   typedef M Mesh;
//   typedef GenericInterface<Mesh> Interface;
//   typedef GFESpace<Mesh> FESpace;
//   typedef typename FESpace::FElement FElement;
//   typedef typename FElement::QFB QFB;
//   typedef typename FElement::Rd Rd;
//   typedef typename Mesh::Partition Partition;
//   typedef typename QFB::QuadraturePoint QuadraturePoint;
//
//   const Interface& interface(Vh.getInterface(it));
//   const QFB& qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
//
//   What_d Fop = Fwhatd(op_id);
//   double val = 0.;
//
//   for(int iface=interface.first_element(); iface<interface.last_element(); iface+=interface.next_element()) {
//
//     const int kb = interface.idxElementOfFace(iface);   // idx on backMesh
//     const typename Interface::Face& face = interface.getFace(iface);  // the face
//     const R meas = interface.computeDx(face).norm();
//
//     int k = Vh.idxElementFromBackMesh(kb,0);
//     const FElement& FK(Vh[k]);
//
//     for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
//
//       QuadraturePoint ip(qfb[ipq]); // integration point
//       const Rd mip = interface.mapToFace(face,(typename FElement::RdHatBord)ip);
//       const R Cint = meas * ip.getWeight();
//
//       val += Cint * fh.eval(k, mip, t);
//
//     }
//   }
//   double val_receive = 0;
//   MPIcf::AllReduce(val, val_receive, MPI_SUM);
//
//   return val_receive;
// }

/*
          Space & Time integration of f(x,t)
*/

#include "normsFunFEM.hpp"





#endif
