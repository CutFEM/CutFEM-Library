#ifndef _GENERAL_NORM_
#define _GENERAL_NORM_

#include "../FESpace/FESpace.hpp"
#include "itemVF.hpp"


template<typename Mesh>
double integralCut(const ListItemVF<Mesh::Rd::d>& VF, const FunFEM<Mesh>& f,const FunFEM<Mesh>& g, const CutGFESpace<Mesh>& Vh) {

  // TYPEDEF NAMES OF TEMPLATE CLASSES
  typedef CutGFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement  FElement;
  typedef typename TypeCutData<Mesh::Rd::d>::CutData CutData;
  typedef typename Mesh::Partition Partition;
  typedef typename FElement::QF QF;

  const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));

  double val = 0.;

  for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {

    // COMPUTE THE PARTITION OF THE ELEMENT
    const FElement& FK(Vh[k]);
    const int domain = FK.whichDomain();
    const int kb = Vh.Th(FK.T);
    CutData cutData(Vh.getInterface(0).getCutData(kb));     // get the cut data
    const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
    ElementSignEnum domPart = cutK.what_part(domain);


    // LOOP OVER PART IN THE DOMAIN
    for(typename Partition::const_element_iterator it = cutK.element_begin(domPart);
    it != cutK.element_end(domPart); ++it){

      const R meas = cutK.mesure(it);

      // LOOP OVER TERMS IN THE WEAK FORMULATION
      for(int l=0; l<VF.size();++l) {

        R coef = VF[l].computeCoef(domain,cutK);

        for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

          typename QF::QuadraturePoint ip(qf[ipq]);
          typename Mesh::Rd mip = cutK.toK(it, ip);

          double Cint = meas * ip.getWeight() * coef;
          double val_fh = f.evalOnBackMesh(kb, mip, VF[l].cu, VF[l].du, domain);
          double val_gh = g.evalOnBackMesh(kb, mip, VF[l].cv, VF[l].dv, domain);

          // VF[l].fxu_backMesh(kb, domain, mip);
          val += Cint * val_fh * val_gh * VF[l].c;
        }
      }
    }
  }

  double val_receive = 0;
  MPIcf::AllReduce(val, val_receive, MPI_SUM);
  return sqrt(val_receive);
}

























#endif
