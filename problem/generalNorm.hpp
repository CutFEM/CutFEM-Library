#ifndef _GENERAL_NORM_
#define _GENERAL_NORM_

#include "../FESpace/FESpace.hpp"
#include "itemVF.hpp"

template <typename Mesh>
double integralCut(const ListItemVF<Mesh::Rd::d> &VF, const FunFEM<Mesh> &f, const FunFEM<Mesh> &g,
                   const ActiveMesh<Mesh> &Th) {

    // TYPEDEF NAMES OF TEMPLATE CLASSES
    typedef GFESpace<Mesh> FESpace;
    typedef typename Mesh::Element Element;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QF QF;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(5));

    double val = 0.;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        const Cut_Part<Element> cutK(Th.get_cut_part(k, 0));
        const Element &K(Th[k]);
        int domain   = Th.get_domain_element(k);
        double measK = K.measure();
        double h     = K.hElement();
        int kb       = Th.idxElementInBackMesh(k);

        // LOOP OVER PART IN THE DOMAIN
        for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

            const R meas = cutK.measure(it);

            // LOOP OVER TERMS IN THE WEAK FORMULATION
            for (int l = 0; l < VF.size(); ++l) {

                double coef = 1.; // VF[l].computeCoefElement(h,meas,measK, meas, domain);

                for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

                    typename QF::QuadraturePoint ip(qf[ipq]);
                    typename Mesh::Rd mip = cutK.mapToPhysicalElement(it, ip);

                    double Cint   = meas * ip.getWeight() * coef;
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
// double integralCut(const ExpressionVirtual& f,const ExpressionVirtual& g,
// const ActiveMesh<Mesh>& Th) {
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
//           double val_fh = f.evalOnBackMesh(kb, mip, VF[l].cu, VF[l].du,
//           domain); double val_gh = g.evalOnBackMesh(kb, mip, VF[l].cv,
//           VF[l].dv, domain);
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
// double integral(const ListItemVF<Mesh::Rd::d>& VF, const FunFEM<Mesh>&
// f,const FunFEM<Mesh>& g, const GenericInterface<Mesh>& interface) {

//   // TYPEDEF NAMES OF TEMPLATE CLASSES
//   //typedef CutGFESpace<Mesh> FESpace;
//   typedef typename FESpace::FElement  FElement;
//   typedef typename TypeCutData<Mesh::Rd::d>::CutData CutData;
//   typedef typename Mesh::Partition Partition;
//   typedef typename FElement::QFB QFB;
//   typedef typename FElement::RdHatBord RdHatBord;
//   typedef GenericInterface<Mesh> Interface;
//   typedef typename Mesh::Rd Rd;

//   const QFB& qf(*QF_Simplex<typename FElement::RdHatBord>(5));
//   double val = 0.;

//   for(int iface=interface.first_element(); iface<interface.last_element();
//   iface+=interface.next_element()) {

//     // COMPUTE THE PARTITION OF THE ELEMENT
//     const int kb = interface.idxElementOfFace(iface);   // idx on
//     const typename Interface::FaceIdx& face = interface[iface];  // the face
//     const R meas = interface.computeDx(face).norm();
//     const Rd normal(-interface.normal(iface));
//     const double h = meas;
//       // LOOP OVER TERMS IN THE WEAK FORMULATION
//       for(int l=0; l<VF.size();++l) {

//         double measK = (*interface.backMesh)[kb].mesure();

//         R coef = VF[l].computeCoefInterface(h,meas,measK );
//         double Cnormal = VF[l].getCoef(normal);

//         for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

//           typename QFB::QuadraturePoint ip(qf[ipq]);
//           Rd mip = interface.mapToFace(face,(RdHatBord)ip);

//           double Cint = meas * ip.getWeight() * coef* VF[l].c * Cnormal;

//           double val_fh = f.evalOnBackMesh(kb, mip, VF[l].cu, VF[l].du,
//           VF[l].domu); double val_gh = g.evalOnBackMesh(kb, mip, VF[l].cv,
//           VF[l].dv, VF[l].domv); val += Cint * val_fh * val_gh ;
//         }
//       }
//     }

//   double val_receive = 0;
//   MPIcf::AllReduce(val, val_receive, MPI_SUM);
//   return sqrt(val_receive);
// }

template <typename Mesh>
double integralCut(const ListItemVF<Mesh::Rd::d> &VF, const CutFEM<Mesh> &cutfem, const FunFEM<Mesh> &f,
                   const FunFEM<Mesh> &g, const TimeInterface<Mesh> &gamma, const TimeSlab &In) {

    double val = 0.;

    typedef typename Mesh::Element Element;
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::Rd Rd;
    typedef typename FElement::QF QF;
    typedef typename FElement::QFB QFB;
    typedef typename Mesh::BorderElement BorderElement;

    const QuadratureFormular1d &qTime = gamma.get_quadrature_time();

    for (int itq = 0; itq < qTime.n; ++itq) {
        assert(!VF.isRHS());
        auto tq    = cutfem.get_quadrature_time(itq);
        double tid = In.map(tq);

        KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
        // RNMK_ bf_time(cutfem.databf_time_, In.NbDoF(), 1, op_dz);
        // In.BF(tq.x, bf_time); // compute time basic funtions
        double cst_time = tq.a * In.get_measure();

        for (int iface = gamma[itq]->first_element(); iface < gamma[itq]->last_element();
             iface += gamma[itq]->next_element()) {
            const typename Interface<Mesh>::Face &face = (*gamma[itq])[iface]; // the face

            typedef typename FElement::RdHatBord RdHatBord;

            // GET IDX ELEMENT CONTAINING FACE ON backMes
            const int kb = gamma[itq]->idxElementOfFace(iface);
            const Element &K(gamma[itq]->get_element(kb));
            double measK = K.measure();
            double h     = K.get_h();
            double meas  = gamma[itq]->measure(iface);

            const Rd normal(-gamma[itq]->normal(iface));

            // GET THE QUADRATURE RULE
            // const QFB &qfb(cutfem.get_quadrature_formular_cutFace());
            const QFB &qfb(*QF_Simplex<typename FElement::RdHatBord>(5));

            for (int l = 0; l < VF.size(); ++l) {
                // if(!VF[l].on(domain)) continue;

                // FINITE ELEMENT SPACES && ELEMENTS
                const FESpace &Vhv(VF.get_spaceV(l));
                const FESpace &Vhu(VF.get_spaceU(l));
                bool same = (VF.isRHS() || (&Vhu == &Vhv));

                std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
                std::vector<int> idxU =
                    (same) ? idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());
                int kv = VF[l].onWhatElementIsTestFunction(idxV);
                int ku = VF[l].onWhatElementIsTrialFunction(idxU);

                const FElement &FKu(Vhu[ku]);
                const FElement &FKv(Vhv[kv]);
                int domu = FKu.get_domain();
                int domv = FKv.get_domain();

                // COMPUTE COEFFICIENT && NORMAL
                double coef = VF[l].computeCoefInterface(h, meas, measK, measK, domu, domv);
                coef *= VF[l].computeCoefFromNormal(normal);

                // LOOP OVER QUADRATURE IN SPACE
                for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {

                    typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
                    const Rd mip = gamma[itq]->mapToPhysicalFace(iface, (RdHatBord)ip);
                    double Cint  = meas * ip.getWeight() * cst_time * VF[l].c * coef;

                    double val_fh = f.evalOnBackMesh(kb, 0, mip, tid, VF[l].cu, VF[l].du, VF[l].dtu);
                    double val_gh = g.evalOnBackMesh(kb, 0, mip, tid, VF[l].cv, VF[l].dv, VF[l].dtu);
                    // std::cout << "Cint = " << Cint << ", val_fh = " << val_fh << ", val_gh = " << val_gh <<
                    // std::endl;
                    val += Cint * val_fh * val_gh;
                    // std::cout << "val = " << val << std::endl;
                }
            }
        }
    }

    double val_receive = 0;
    MPIcf::AllReduce(val, val_receive, MPI_SUM);
    // return sqrt(val_receive);
    return val_receive;
}

#endif
