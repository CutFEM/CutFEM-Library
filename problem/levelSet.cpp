#include "levelSet.hpp"
// 
//
// void LevelSet2::assembly(const Fun& up, const Fun& Betap, const Fun& Beta, double dt) {
//   typedef typename QF::QuadraturePoint QuadraturePoint;
//
//   const int d = Rd::d;
//   const FESpace2 & Lh(*this->Vh);
//   RNMK_ fu(this->databf,Lh[0].NbDoF(),1,op_all);
//
//   for(int k=Lh.first_element(); k<Lh.last_element(); k+= Lh.next_element()) {
//
//     const FElement& FK(Lh[k]);
//     double h = FK.T.hElement();
//     double meas = FK.getMeasure();
//
//     for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//
//       QuadraturePoint ip(qf[ipq]); // integration point
//       Rd mip = FK.map(ip);
//       double Cint = meas * ip.getWeight();
//
//       FK.BF(ip, fu); // need point in local reference element
//       double Bx = Beta.eval(k, mip, 0, op_id);
//       double By = Beta.eval(k, mip, 1, op_id);
//       double Bpx = Betap.eval(k, mip, 0, op_id);
//       double Bpy = Betap.eval(k, mip, 1, op_id);
//       double Up = up.eval(k, mip, 0, op_id);
//       double dxup = up.eval(k, mip, 0, op_dx);
//       double dyup = up.eval(k, mip, 0, op_dy);
//       double normB = Bx*Bx + By*By;
//       double Tsd =  2. / (sqrt(1./dt/dt + normB * 1. / h / h));
//
//
//       for(int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
//         for(int j = FK.dfcbegin(0); j < FK.dfcend(0); ++j) {
//           (*this)(FK.loc2glb(i),FK.loc2glb(j)) +=  Cint * (
//             ((1./dt)*fu(j,0,op_id) + 0.5*(Bx*fu(j,0,op_dx)+ By*fu(j,0,op_dy)))
//             *
//             (fu(i,0,op_id) + Tsd * (Bx*fu(i,0,op_dx) + By*fu(i,0,op_dy)))
//           );
//         }
//         (*this)(FK.loc2glb(i)) +=
//         Cint*(
//           ((1./dt)*Up - 0.5*(Bpx*dxup + Bpy*dyup))
//           *
//           (fu(i,0,op_id) + Tsd * (Bx*fu(i,0,op_dx) + By*fu(i,0,op_dy)))
//         )
//         ;
//       }
//
//     }
//     if(label_strongBC.size() > 0){
//       ExpressionFunFEM<Mesh2> Up(up, 0, op_id);
//       this->addStrongBC(Up, label_strongBC);
//     }
//   }
// }
// void LevelSet3::assembly(const Fun& up, const Fun& Betap, const Fun& Beta, double dt) {
//   typedef typename QF::QuadraturePoint QuadraturePoint;
//
//   const int d = Rd::d;
//   const FESpace3 & Lh(*this->Vh);
//   RNMK_ fu(this->databf,Lh[0].NbDoF(),1,op_all);
//
//   for(int k=Lh.first_element(); k<Lh.last_element(); k+= Lh.next_element()) {
//
//     const FElement& FK(Lh[k]);
//     double h = FK.T.hElement();
//     double meas = FK.getMeasure();
//
//     for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//
//       QuadraturePoint ip(qf[ipq]); // integration point
//       Rd mip = FK.map(ip);
//       double Cint = meas * ip.getWeight();
//
//       FK.BF(ip, fu); // need point in local reference element
//       double Bx = Beta.eval(k, mip, 0, op_id);
//       double By = Beta.eval(k, mip, 1, op_id);
//       double Bz = Beta.eval(k, mip, 2, op_id);
//       double Bpx = Betap.eval(k, mip, 0, op_id);
//       double Bpy = Betap.eval(k, mip, 1, op_id);
//       double Bpz = Betap.eval(k, mip, 2, op_id);
//       double Up = up.eval(k, mip, 0, op_id);
//       double dxup = up.eval(k, mip, 0, op_dx);
//       double dyup = up.eval(k, mip, 0, op_dy);
//       double dzup = up.eval(k, mip, 0, op_dz);
//       double normB = Bx*Bx + By*By + Bz*Bz;
//       double Tsd =  2. / (sqrt(1./dt/dt + normB * 1. / h / h ));
//
//
//       for(int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
//         for(int j = FK.dfcbegin(0); j < FK.dfcend(0); ++j) {
//           (*this)(FK.loc2glb(i),FK.loc2glb(j)) +=  Cint * (
//             ((1./dt)*fu(j,0,op_id) + 0.5*(Bx*fu(j,0,op_dx)+ By*fu(j,0,op_dy)+ Bz*fu(j,0,op_dz)))
//             *
//             (fu(i,0,op_id) + Tsd * (Bx*fu(i,0,op_dx) + By*fu(i,0,op_dy)+ Bz*fu(i,0,op_dz)))
//             // + epsilon_diffusion*fu(i,0,op_dx)*fu(j,0,op_dx)
//             // + epsilon_diffusion*fu(i,0,op_dy)*fu(j,0,op_dy)
//             // + epsilon_diffusion*fu(i,0,op_dz)*fu(j,0,op_dz)
//           );
//
//         }
//         (*this)(FK.loc2glb(i)) +=
//         Cint*(
//           ((1./dt)*Up - 0.5*(Bpx*dxup + Bpy*dyup+ Bpz*dzup))
//           *
//           (fu(i,0,op_id) + Tsd * (Bx*fu(i,0,op_dx) + By*fu(i,0,op_dy)+ Bz*fu(i,0,op_dz)))
//         )
//         ;
//       }
//     }
//     if(label_strongBC.size() > 0){
//       ExpressionFunFEM<Mesh3> Up(up, 0, op_id);
//       this->addStrongBC(Up, label_strongBC);
//     }
//   }
// }
