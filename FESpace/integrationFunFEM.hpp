#ifndef INTEGRATION_FUNFEM_HPP_
#define INTEGRATION_FUNFEM_HPP_

#include "expression.hpp"
#include "../common/Interface2dn.hpp"
#include "macroElement.hpp"



/*
          Space integration of f(x)
*/
template<typename M>
double integral(FunFEM<M>& fh, int cu, int domain = -1){

  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QF QF;
  typedef typename FElement::Rd Rd;
  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;
  typedef typename TypeCutData<Rd::d>::CutData CutData;
  const FESpace& Vh(*fh.Vh);
  const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));

  What_d Fop = Fwhatd(op_id);
  double val = 0.;

  for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {

    const FElement& FK(Vh[k]);

    if(domain == -1){
      const R meas = FK.getMeasure();
      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = FK.map(ip);
        const R Cint = meas * ip.getWeight();
        val += Cint * fh.eval(k, mip, cu, op_id);
      }
    }
    else if(domain == FK.whichDomain()) {
      // if(domain != FK.whichDomain()) continue;

      const int kb = Vh.Th(FK.T);
      CutData cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
      const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
      ElementSignEnum the_part = cutK.what_part(domain);

      for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
      it != cutK.element_end(the_part); ++it){

        const R meas = cutK.mesure(it);

        for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

          QuadraturePoint ip(qf[ipq]); // integration point
          Rd mip = cutK.toK(it, ip);
          const R Cint = meas * ip.getWeight();

          val += Cint * fh.eval(k, mip, cu, op_id);

        }
      }
    }
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;
}

template<typename M>
double integral(const ExpressionFunFEM<M>& fh, int domain) {

  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QF QF;
  typedef typename FElement::Rd Rd;
  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;
  typedef typename TypeCutData<Rd::d>::CutData CutData;

  const FESpace& Vh(fh.getSpace());
  const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));

  What_d Fop = Fwhatd(op_id);

  double val = 0.;

  for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {

    const FElement& FK(Vh[k]);
    if(domain == -1){
      const R meas = FK.getMeasure();
      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = FK.map(ip);
        const R Cint = meas * ip.getWeight();
        val += Cint * fh.eval(k, mip);
      }
    }
    else if(domain == FK.whichDomain()) {
      // if(domain != FK.whichDomain()) continue;

      const int kb = Vh.Th(FK.T);
      CutData cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
      const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
      ElementSignEnum the_part = cutK.what_part(domain);

      for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
      it != cutK.element_end(the_part); ++it){

        const R meas = cutK.mesure(it);

        for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

          QuadraturePoint ip(qf[ipq]); // integration point
          Rd mip = cutK.toK(it, ip);
          const R Cint = meas * ip.getWeight();

          val += Cint * fh.eval(k, mip);

        }
      }
    }
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;
}

template<typename M>
double integral(const ExpressionFunFEM<M>& fh) {
  if(fh.getSpace().isCut()) return integral(fh,0) + integral(fh,1);
  else return integral(fh,-1) ;
}


template<typename M>
double integral(const ExpressionVirtual& fh, const GFESpace<M>& Vh, int domain) {
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

    if(domain == -1){
      const R meas = FK.getMeasure();
      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = FK.map(ip);
        const R Cint = meas * ip.getWeight();
        val += Cint * fh.eval(k, mip);
      }
    }
    else if(domain == FK.whichDomain()) {
      // if(domain != FK.whichDomain()) continue;

      const int kb = Vh.Th(FK.T);
      CutData cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
      const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
      ElementSignEnum the_part = cutK.what_part(domain);

      for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
      it != cutK.element_end(the_part); ++it){

        const R meas = cutK.mesure(it);

        for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

          QuadraturePoint ip(qf[ipq]); // integration point
          Rd mip = cutK.toK(it, ip);
          const R Cint = meas * ip.getWeight();

          val += Cint * fh.eval(k, mip);

        }
      }
    }
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;

}

template<typename M>
double integral(const ExpressionVirtual& fh, const GFESpace<M>& Vh) {
  if(Vh.isCut()) return integral(fh,Vh,0) + integral(fh,Vh,1);
  else return integral(fh,Vh,-1);
}


template<typename M>
double integralSurf(FunFEM<M>& fh, int cu=0, int it = 0, double t = -1.) {

  typedef M Mesh;
  typedef GenericInterface<Mesh> Interface;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QFB QFB;
  typedef typename FElement::Rd Rd;
  typedef typename Mesh::Partition Partition;
  typedef typename QFB::QuadraturePoint QuadraturePoint;

  if(t > -Epsilon && fh.In){
    assert(fh.In->Pt(0) <=t && t <= fh.In->Pt(1));
  }
  if(t < -Epsilon && fh.In) {
    t = fh.In->Pt(0);
    std::cout << " Use default value In(0) \t -> " << t << std::endl;
  }

  const FESpace& Vh(*fh.Vh);
  const Interface& interface(Vh.getInterface(it));
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

      val += Cint * fh.eval(k, mip, t, cu);

    }
  }

  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;
}

template<typename M>
double integralSurf(const ExpressionVirtual& fh, const GFESpace<M>& Vh, int it = 0, double t = 0.) {

  typedef M Mesh;
  typedef GenericInterface<Mesh> Interface;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QFB QFB;
  typedef typename FElement::Rd Rd;
  typedef typename Mesh::Partition Partition;
  typedef typename QFB::QuadraturePoint QuadraturePoint;

  const Interface& interface(Vh.getInterface(it));
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

      val += Cint * fh.eval(k, mip, t);

    }
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;
}

// /*
//           Space integration of f(x,t)
// */
// template<typename M>
// double integralSurf(FunFEM<M>& fh, int cu, int it, double t) {
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
template<typename M>
double integralSurf(FunFEM<M>& fh, const TimeSlab& In, const QuadratureFormular1d& qTime, int cu=0) {

  typedef M Mesh;
  typedef GenericInterface<Mesh> Interface;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QFB QFB;
  typedef typename FElement::Rd Rd;
  typedef typename Mesh::Partition Partition;
  typedef typename QFB::QuadraturePoint QuadraturePoint;


  const FESpace& Vh(*fh.Vh);
  const QFB& qfb(*QF_Simplex<typename FElement::RdHatBord>(5));

  What_d Fop = Fwhatd(op_id);
  double val = 0.;

  for(int it = 0; it<qTime.n; ++it){
    const Interface& interface(Vh.getInterface(it));
    GQuadraturePoint<R1> tq((qTime)[it]);
    const double t = In.map(tq);
    for(int iface=interface.first_element(); iface<interface.last_element(); iface+=interface.next_element()) {

      const int kb = interface.idxElementOfFace(iface);   // idx on backMesh
      const typename Interface::Face& face = interface.getFace(iface);  // the face
      const R meas = interface.computeDx(face).norm();

      int k = Vh.idxElementFromBackMesh(kb,0);
      const FElement& FK(Vh[k]);

      for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qfb[ipq]); // integration point
        const Rd mip = interface.mapToFace(face,(typename FElement::RdHatBord)ip);
        const R Cint = meas * ip.getWeight() * In.T.mesure()* tq.a ;

        val += Cint * fh.eval(k, mip, t, cu);

      }
    }
  }
  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);

  return val_receive;
}

#include "normsFunFEM.hpp"





#endif
