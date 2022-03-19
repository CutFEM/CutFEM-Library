#include "baseCutProblem.hpp"

//
// template<>
// void BaseCutProblem<Mesh3>::addElementMatEdge(const ListItemVF<Rd::d>& VF, const int k, const int ifac) {
//
//   // typedef typename Mesh::Partition Partition;
//   typedef typename QFB::QuadraturePoint QuadraturePoint;
//   typedef typename FElement::RdHatBord RdHatBord;
//   typedef typename Mesh::Element Element;
//   typedef typename Element::Face Face;
//
//   assert(VF[0].fespaceU && VF[0].fespaceV);
//   const FESpace& Vhu = *VF[0].fespaceU;
//   const FESpace& Vhv = *VF[0].fespaceV;
//   const FElement& FK(Vhu[k]);
//   const int kb = Vhu.idxElementInBackMesh(k);
//
//   CutData cutData(Vh->getInterface(0).getCutData(kb));
//
//   if(!Vhu.isCutSpace() || cutData.face_isnot_cut(ifac)){
//     return;
//     BaseProblem<Mesh3>::addElementMatEdge(VF, k, ifac);
//     return;
//   }
//
//   double ss = 0.;
//   // for(int the_domain=0;the_domain<2;++the_domain){
//   int the_domain = FK.whichDomain();
//   int ifacn = ifac;
//   int kn = Vhu.getNeighborElement(k, ifacn, the_domain);
//
//
//   if(kn == -1) return;         // border edge
//   if(k > kn) return;           // only compute one time
//
//   // const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
//   //
//   // Rd normal = FK.T.N(ifac);
//   // const R meas_cut = cutK.mesure(the_domain);
//   // // const R measK = cutK.mesure(the_domain);
//   //
//   // const R h = FK.T.lenEdge(ifac);
//   // const R meas = cutK.mesureEdge(ifac, the_domain);
//   // ss += meas;
//   //
//   // const FElement& FKn(Vhu[kn]);
//   // const int kbn = Vhu.idxElementInBackMesh(kn);
//   // CutData cutDatan(Vh->getInterface(0).getCutData(kbn));
//   // const Partition& cutKn =  Partition(FKn.T, cutDatan);  // build the cut
//   // double measK = cutK.mesure(the_domain) + cutKn.mesure(the_domain);
//
//
//   std::cout << " Start build partition CutFace " << std::endl;
//
//   Face face;
//   FK.T.set_face(ifac, face);
//   byte ls_sign[Face::nv];
//   for(int i=0;i<Face::nv;++i) ls_sign[i] = cutData.sign[Element::nvface[ifac][i]];
//   Partition<Face> face_partition(face, ls_sign);
//
//
//   getchar();
//   return;
//
//
//
// //   for(int l=0; l<VF.size();++l) {
// //
// //     assert(VF[l].fespaceU == VF[l].fespaceV);
// //     int lastop = getLastop(VF[l].du, VF[l].dv);
// //
// //     // R coef = BaseProblem<M>::computeCoef(VF[l],h,meas, measK,the_domain) * VF[l].getCoef(normal);
// //     R coef = BaseCutProblem<Mesh3>::computeCoefEdge(VF[l],h,meas, measK,meas_cut, the_domain) * VF[l].getCoef(normal);
// //
// //     const int ku = (VF[l].domu == 0)? k : kn;
// //     const int kv = (VF[l].domv == 0)? k : kn;
// //     const FElement& FKu(Vhu[ku]);
// //     const FElement& FKv(Vhv[kv]);
// //     this->initIndex(FKu, FKv);
// //     int kuback = Vh->idxElementInBackMesh(ku);
// //     int kvback = Vh->idxElementInBackMesh(kv);
// //
// //     bool same = (ku == kv );
// //     What_d Fop = Fwhatd(lastop);
// //     RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
// //     RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
// //
// //
// //     for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
// //
// //       QuadraturePoint ip(qfb[ipq]); // integration point
// //       const Rd mip = cutK.toEdge(ifac, (RdHatBord)ip, the_domain); // mip is in the global edge
// //       const R Cint = meas * ip.getWeight();
// //
// //       FKu.BF(Fop,FKu.T.toKref(mip), fu); // need point in local reference element
// //       if(!same)FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element
// //
// //       double val = VF[l].fx_backMesh_U(kuback,the_domain, mip, normal)
// //                   *VF[l].fx_backMesh_V(kvback,the_domain, mip, normal);
// //       double Cst = Cint * VF[l].c  * val * coef;
// //
// //       for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
// //         for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
// //           this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j))  +=  Cst * fv(i,VF[l].cv,VF[l].dv) * fu(j,VF[l].cu,VF[l].du);
// //         }
// //       }
// //     }
// //   }
// // // }
// //   this->resetIndex();
// //   this->addLocalContribution();
//
// }
//
// template<>
// void BaseCutProblem<Mesh3>::addElementMatEdge(const ListItemVF<Rd::d>& VF, const int ifac) {
//
//   // typedef typename Mesh::Partition Partition;
//   typedef typename QFB::QuadraturePoint QuadraturePoint;
//   typedef typename FElement::RdHatBord RdHatBord;
//   typedef typename Mesh::Element Element;
//   typedef typename Element::Face Face;
//
//
//
//   assert(VF[0].fespaceU && VF[0].fespaceV);
//   const FESpace& Vhu = *VF[0].fespaceU;
//   const FESpace& Vhv = *VF[0].fespaceV;
//   assert(&Vhu.Th == &Vhv.Th);
//
//   std::cout << " Start build partition CutFace " << std::endl;
//
//
//   getchar();
//   return;
//
//
//   // const FElement& FK(Vhu[k]);
//   // const int kb = Vhu.idxElementInBackMesh(k);
//
// //   CutData cutData(Vh->getInterface(0).getCutData(kb));
// //
// //   if(!Vhu.isCutSpace() || cutData.face_isnot_cut(ifac)){
// //     return;
// //     BaseProblem<Mesh3>::addElementMatEdge(VF, ifac);
// //     return;
// //   }
// //
// //   double ss = 0.;
// //   // for(int the_domain=0;the_domain<2;++the_domain){
// //   int the_domain = FK.whichDomain();
// //   int ifacn = ifac;
// //   int kn = Vhu.getNeighborElement(k, ifacn, the_domain);
// //
// //
// //   if(kn == -1) return;         // border edge
// //   if(k > kn) return;           // only compute one time
// //
// //   const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
// //
// //   Rd normal = FK.T.N(ifac);
// //   const R meas_cut = cutK.mesure(the_domain);
// //   // const R measK = cutK.mesure(the_domain);
// //
// //   const R h = FK.T.lenEdge(ifac);
// //   const R meas = cutK.mesureEdge(ifac, the_domain);
// //   ss += meas;
// //
// //   const FElement& FKn(Vhu[kn]);
// //   const int kbn = Vhu.idxElementInBackMesh(kn);
// //   CutData cutDatan(Vh->getInterface(0).getCutData(kbn));
// //   const Partition& cutKn =  Partition(FKn.T, cutDatan);  // build the cut
// //   double measK = cutK.mesure(the_domain) + cutKn.mesure(the_domain);
// //
// //
//
// //
// //
// //
// //   for(int l=0; l<VF.size();++l) {
// //
// //     assert(VF[l].fespaceU == VF[l].fespaceV);
// //     int lastop = getLastop(VF[l].du, VF[l].dv);
// //
// //     // R coef = BaseProblem<M>::computeCoef(VF[l],h,meas, measK,the_domain) * VF[l].getCoef(normal);
// //     R coef = BaseCutProblem<Mesh3>::computeCoefEdge(VF[l],h,meas, measK,meas_cut, the_domain) * VF[l].getCoef(normal);
// //
// //     const int ku = (VF[l].domu == 0)? k : kn;
// //     const int kv = (VF[l].domv == 0)? k : kn;
// //     const FElement& FKu(Vhu[ku]);
// //     const FElement& FKv(Vhv[kv]);
// //     this->initIndex(FKu, FKv);
// //     int kuback = Vh->idxElementInBackMesh(ku);
// //     int kvback = Vh->idxElementInBackMesh(kv);
// //
// //     bool same = (ku == kv );
// //     What_d Fop = Fwhatd(lastop);
// //     RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
// //     RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
// //
// //
// //     for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
// //
// //       QuadraturePoint ip(qfb[ipq]); // integration point
// //       const Rd mip = cutK.toEdge(ifac, (RdHatBord)ip, the_domain); // mip is in the global edge
// //       const R Cint = meas * ip.getWeight();
// //
// //       FKu.BF(Fop,FKu.T.toKref(mip), fu); // need point in local reference element
// //       if(!same)FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element
// //
// //       double val = VF[l].fx_backMesh_U(kuback,the_domain, mip, normal)
// //                   *VF[l].fx_backMesh_V(kvback,the_domain, mip, normal);
// //       double Cst = Cint * VF[l].c  * val * coef;
// //
// //       for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
// //         for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
// //           this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j))  +=  Cst * fv(i,VF[l].cv,VF[l].dv) * fu(j,VF[l].cu,VF[l].du);
// //         }
// //       }
// //     }
// //   }
// // // }
// //   this->resetIndex();
// //   this->addLocalContribution();
//
// }
//

/*
Add integral on cut Edges
*/

// template<>
// void BaseCutProblem<Mesh2>::addElementMatEdge(const ListItemVF<Rd::d>& VF, const int ifac) {
//
//   typedef typename Mesh::Partition Partition;
//   typedef typename QFB::QuadraturePoint QuadraturePoint;
//   typedef typename FElement::RdHatBord RdHatBord;
//   typedef typename Mesh::Element Element;
//
//
//   std::cout << " Start build partition CutFace " << std::endl;
//
//
//   getchar();
//   return;
//
//
//
//
//   // assert(VF[0].fespaceU && VF[0].fespaceV);
//   // const FESpace& Vhu = *VF[0].fespaceU;
//   // const FESpace& Vhv = *VF[0].fespaceV;
//   // const FElement& FK(Vhu[k]);
//   // const int kb = Vhu.idxElementInBackMesh(k);
//
// //   CutData cutData(Vh->getInterface(0).getCutData(kb));
// //
// //   if(!Vhu.isCutSpace() || cutData.face_isnot_cut(ifac)){
// //     return;
// //     BaseProblem<Mesh3>::addElementMatEdge(VF, ifac);
// //     return;
// //   }
// //
// //   double ss = 0.;
// //   // for(int the_domain=0;the_domain<2;++the_domain){
// //   int the_domain = FK.whichDomain();
// //   int ifacn = ifac;
// //   int kn = Vhu.getNeighborElement(k, ifacn, the_domain);
// //
// //
// //   if(kn == -1) return;         // border edge
// //   if(k > kn) return;           // only compute one time
// //
// //   const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
// //
// //   Rd normal = FK.T.N(ifac);
// //   const R meas_cut = cutK.mesure(the_domain);
// //   // const R measK = cutK.mesure(the_domain);
// //
// //   const R h = FK.T.lenEdge(ifac);
// //   const R meas = cutK.mesureEdge(ifac, the_domain);
// //   ss += meas;
// //
// //   const FElement& FKn(Vhu[kn]);
// //   const int kbn = Vhu.idxElementInBackMesh(kn);
// //   CutData cutDatan(Vh->getInterface(0).getCutData(kbn));
// //   const Partition& cutKn =  Partition(FKn.T, cutDatan);  // build the cut
// //   double measK = cutK.mesure(the_domain) + cutKn.mesure(the_domain);
// //
// //
//
// //
// //
// //
// //   for(int l=0; l<VF.size();++l) {
// //
// //     assert(VF[l].fespaceU == VF[l].fespaceV);
// //     int lastop = getLastop(VF[l].du, VF[l].dv);
// //
// //     // R coef = BaseProblem<M>::computeCoef(VF[l],h,meas, measK,the_domain) * VF[l].getCoef(normal);
// //     R coef = BaseCutProblem<Mesh3>::computeCoefEdge(VF[l],h,meas, measK,meas_cut, the_domain) * VF[l].getCoef(normal);
// //
// //     const int ku = (VF[l].domu == 0)? k : kn;
// //     const int kv = (VF[l].domv == 0)? k : kn;
// //     const FElement& FKu(Vhu[ku]);
// //     const FElement& FKv(Vhv[kv]);
// //     this->initIndex(FKu, FKv);
// //     int kuback = Vh->idxElementInBackMesh(ku);
// //     int kvback = Vh->idxElementInBackMesh(kv);
// //
// //     bool same = (ku == kv );
// //     What_d Fop = Fwhatd(lastop);
// //     RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
// //     RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
// //
// //
// //     for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
// //
// //       QuadraturePoint ip(qfb[ipq]); // integration point
// //       const Rd mip = cutK.toEdge(ifac, (RdHatBord)ip, the_domain); // mip is in the global edge
// //       const R Cint = meas * ip.getWeight();
// //
// //       FKu.BF(Fop,FKu.T.toKref(mip), fu); // need point in local reference element
// //       if(!same)FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element
// //
// //       double val = VF[l].fx_backMesh_U(kuback,the_domain, mip, normal)
// //                   *VF[l].fx_backMesh_V(kvback,the_domain, mip, normal);
// //       double Cst = Cint * VF[l].c  * val * coef;
// //
// //       for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
// //         for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
// //           this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j))  +=  Cst * fv(i,VF[l].cv,VF[l].dv) * fu(j,VF[l].cu,VF[l].du);
// //         }
// //       }
// //     }
// //   }
// // // }
// //   this->resetIndex();
// //   this->addLocalContribution();
//
// }
