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
      val += integral(Th, ui, i, 0);
    }
    return val;
}

template<typename M>
double integral( const ActiveMesh<M>& Th, const FunFEM<M>& fh, int c0, int itq) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    ExpressionFunFEM<M> ui(fh, c0, op_id);
    for(int i=0;i<nb_dom;++i) {
      val += integral(Th, ui, i, itq);
    }
    return val;
}

template<typename M>
double integral( const ActiveMesh<M>& Th, const FunFEM<M>& fh, int c0, int domain, int itq) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    ExpressionFunFEM<M> ui(fh, c0, op_id);
      val += integral(Th, ui, domain, itq);
    return val;
}
template<typename M>
double integral( const ActiveMesh<M>& Th, const ExpressionVirtual& fh) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for(int i=0;i<nb_dom;++i) {
      val += integral(Th, fh, i, 0);
    }
    return val;
}
template<typename M>
double integral( const ActiveMesh<M>& Th, const ExpressionVirtual& fh, int itq) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for(int i=0;i<nb_dom;++i) {
      val += integral(Th, fh, i, itq);
    }
    return val;
}
template<typename M>
double integral(const ActiveMesh<M>& Th, const ExpressionVirtual& fh, int domain, int itq) {
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

      if(Th.isInactive(k, itq)) continue;

      const Cut_Part<Element> cutK(Th.get_cut_part(k,itq));
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
double integral( const ActiveMesh<M>& Th, const TimeSlab& In, const FunFEM<M>& fh, int c0, const QuadratureFormular1d& qTime) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    ExpressionFunFEM<M> ui(fh, c0, op_id);
    for(int i=0;i<nb_dom;++i) {
      val += integral(Th, In, ui, qTime, i);
    }
    return val;
}


template<typename M>
double integral(const ActiveMesh<M>& Th, const TimeSlab& In, const ExpressionVirtual& fh, const QuadratureFormular1d& qTime, int domain){
    typedef M Mesh;
    typedef typename Mesh::Element Element;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename FElement::Rd Rd;
    typedef typename QF::QuadraturePoint QuadraturePoint;

    const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
    double val = 0.;

    for(int itq = 0; itq<qTime.n; ++itq){
      GQuadraturePoint<R1> tq((qTime)[itq]);
      const double t = In.mapToPhysicalElement(tq);


      for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {

        if(domain != Th.get_domain_element(k)) continue;
        if(Th.isInactive(k, itq)) continue;


        const Cut_Part<Element> cutK(Th.get_cut_part(k, itq));
        int kb = Th.idxElementInBackMesh(k);
        for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){

          const R meas = cutK.measure(it);

          for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

            QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip = cutK.mapToPhysicalElement(it, ip);
            const R Cint = meas * ip.getWeight()* In.T.mesure()* tq.a ;

            val += Cint * fh.evalOnBackMesh(kb, domain, mip, t);

          }
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






template<typename Mesh>
double integral( const Mesh& Th, const ExpressionVirtual& fh) {
    double val = integral(Th, fh, 0);
    return val;
}

template<typename Mesh>
double integral(const Mesh& Th, const ExpressionVirtual& fh, int itq) {
    typedef typename Mesh::Element Element;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename FElement::Rd Rd;
    typedef typename QF::QuadraturePoint QuadraturePoint;

    const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));

    double val = 0.;
    for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {

      const Element& K(Th[k]);
      const R meas = K.measure();

      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = K.mapToPhysicalElement(ip);
        const R Cint = meas * ip.getWeight();

        val += Cint * fh.evalOnBackMesh(k, 0, mip);
      }
    }
    double val_receive = 0;
    MPIcf::AllReduce(val, val_receive, MPI_SUM);

    return val_receive;
}




// new
template<typename M>
double integral( const ActiveMesh<M>& Th, const TimeSlab& In, const ExpressionVirtual& fh,  const CBorder& b, const QuadratureFormular1d& qTime, list<int> label = {}) {
  
    int nb_dom = Th.get_nb_domain();
    double val = 0.;

    for(int i=0;i<nb_dom;++i) {
      val += integral(Th, In, fh, b, qTime, i, label);
    }
    return val;
}

template<typename M>
double integral(const ActiveMesh<M>& Th, const TimeSlab& In, const ExpressionVirtual& fh, const CBorder& b, const QuadratureFormular1d& qTime, int domain, list<int> label){
    
    typedef M Mesh;
    typedef typename Mesh::Element Element;
    typedef typename Mesh::BorderElement BorderElement;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QFB QFB;
    typedef typename FElement::Rd Rd;
    typedef typename QFB::QuadraturePoint QuadraturePoint;
    

    bool all_label = (label.size() == 0);

    double val = 0.;

    for(int itq = 0; itq<qTime.n; ++itq){
    
        GQuadraturePoint<R1> tq((qTime)[itq]);
        const double t = In.mapToPhysicalElement(tq);


        for(int idx_be=Th.first_boundary_element(); idx_be<Th.last_boundary_element(); idx_be+= Th.next_boundary_element()) {

            int ifac;
            const int kb = Th.Th.BoundaryElement(idx_be, ifac);
            std::vector<int> idxK = Th.idxAllElementFromBackMesh(kb, -1);

            const BorderElement & BE(Th.be(idx_be));

            if (util::contain(label, BE.lab) || all_label) {
                
                for (int i=0; i < idxK.size(); i++) {
            
                int k = idxK[i];
                if (Th.get_domain_element(k) != domain) continue;

                if (Th.isInactive(k, itq)) continue;

                if (Th.isCutFace(k, ifac, itq)) val += integral_dK_cut(fh, Th, In, k, ifac, qTime, itq);
                else val += integral_dK(fh, Th, In, k, BE, ifac, qTime, itq);

                }
            }
          
        }

    }

    double val_receive = 0;
    MPIcf::AllReduce(val, val_receive, MPI_SUM);

    return val_receive;
}


template<typename Mesh>
double integral_dK_cut(const ExpressionVirtual& fh, const ActiveMesh<Mesh>& Th, const TimeSlab& In, int k, int ifac, const QuadratureFormular1d& qTime, int itq) {
    
    typedef typename Mesh::Element Element;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QFB QFB;
    typedef typename FElement::Rd Rd;
    typedef typename FElement::RdHatBord RdHatBord;
    typedef typename QFB::QuadraturePoint QuadraturePoint;
    
    const QFB& qfb(*QF_Simplex<typename FElement::RdHatBord>(5));

    double val = 0.;
    auto tq = qTime[itq];
    double tid = In.map(tq);
    double cst_time = tq.a*In.get_measure();

    int domain = Th.get_domain_element(k);

    const Element & K(Th[k]);
    Rd normal    = K.N(ifac);

    const Cut_Part<Element> cutK(Th.get_cut_part(k, itq));

    typename Element::Face face;
    const Cut_Part<typename Element::Face> cutFace(Th.get_cut_face(face, k, ifac));

    int kb = Th.idxElementInBackMesh(k);

    for(auto it = cutFace.element_begin();it != cutFace.element_end(); ++it){

      const R meas = cutFace.measure(it);

      for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

        typename QFB::QuadraturePoint ip(qfb[ipq]);
        const Rd mip = cutFace.mapToPhysicalElement(it, (RdHatBord)ip);
        double Cint = meas * ip.getWeight();

        val += Cint * fh.evalOnBackMesh(kb, domain, mip, tid, normal) * cst_time;

      }
    }

    return val;
}

template<typename Mesh>
double integral_dK(const ExpressionVirtual& fh, const ActiveMesh<Mesh>& Th, const TimeSlab& In, int k, const typename Mesh::BorderElement& BE, int ifac, const QuadratureFormular1d& qTime, int itq) {
  
    typedef typename Mesh::Element Element;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QFB QFB;
    typedef typename FElement::Rd Rd;
    typedef typename FElement::RdHatBord RdHatBord;
    typedef typename QFB::QuadraturePoint QuadraturePoint;
    typedef typename Mesh::BorderElement BorderElement;

    const QFB& qfb(*QF_Simplex<typename FElement::RdHatBord>(5));

    double val = 0.;
    auto tq = qTime[itq];
    double tid = In.map(tq);
    double cst_time = tq.a*In.get_measure();

    int domain = Th.get_domain_element(k);

    const Element & K(Th[k]);

    Rd normal    = K.N(ifac);

    int kb = Th.idxElementInBackMesh(k);

    double meas = K.mesureBord(ifac);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      typename QFB::QuadraturePoint ip(qfb[ipq]);
      const Rd mip = BE.mapToPhysicalElement((RdHatBord)ip);
      double Cint = meas * ip.getWeight();

      val += Cint * fh.evalOnBackMesh(kb, domain, mip, tid, normal) * cst_time;
    }

    return val;
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


#include "normsFunFEM.hpp"





#endif
