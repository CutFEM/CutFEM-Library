/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/
#include "reinitialization.hpp"

template <>
void Reinitialization<Mesh2>::assembly(const Fun &up, const Fun &u0) {
   typedef Mesh2 Mesh;
   typedef typename Mesh::Partition Partition;
   typedef typename QF::QuadraturePoint QuadraturePoint;

   RNMK_ fu(this->databf, ndofK, 1, op_all); //  the value for basic fonction

   for (int k = Vh.first_element(); k < Vh.last_element();
        k += Vh.next_element()) {

      const FElement &FK(Vh[k]);
      const int kb = Vh.Th(FK.T);
      CutData cutData(interface.getCutData(kb));        // get the cut data
      const Partition &cutK = Partition(FK.T, cutData); // build the cut

      const R h = FK.T.hElement();

      for (typename Partition::const_element_iterator it = cutK.element_begin();
           it != cutK.element_end(); ++it) {

         const R meas        = cutK.mesure(it);
         const double signU0 = cutK.whatSign(it);
         R Tsd               = 1. / (sqrt(1. / dt / dt + 1. / h / h));

         for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

            QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip       = cutK.toK(it, ip);
            const R Cint = meas * ip.getWeight();

            FK.BF(FK.T.toKref(mip),
                  fu); // need point in local reference element

            R ux = up.eval(k, mip, 0, op_id);
            Rd dup;
            dup.x      = up.eval(k, mip, 0, op_dx);
            dup.y      = up.eval(k, mip, 0, op_dy);
            R normDup  = (fabs(dup.norm()) < 1e-14) ? 1e-9 : dup.norm();
            const Rd w = dup / normDup * signU0;
            const R ff = signU0;

            if (ODE_method == "RK2") {
               Rd dphiHalf;
               for (int i = 0; i < ndofK; ++i) {
                  Rd ddu;
                  ddu.x     = up.eval(k, FK.Pt(i), 0, op_dx);
                  ddu.y     = up.eval(k, FK.Pt(i), 0, op_dy);
                  R uu      = up.eval(k, FK.Pt(i), 0, op_id);
                  R uu0     = u0.eval(k, FK.Pt(i), 0, op_id);
                  R signUU0 = (uu0 > 0) ? 1 : -1;
                  R uHalf   = uu + dt / 2 * (-signUU0 + ddu.norm());
                  dphiHalf.x += uHalf * fu(i, 0, op_dx);
                  dphiHalf.y += uHalf * fu(i, 0, op_dy);
               }
               R normDupHalf =
                   (fabs(dphiHalf.norm()) < 1e-14) ? 1e-9 : dphiHalf.norm();

               // RK2
               for (int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
                  for (int j = FK.dfcbegin(0); j < FK.dfcend(0); ++j) {
                     (*this)(FK.loc2glb(i), FK.loc2glb(j)) +=
                         Cint * ((1. / dt) * fu(j, 0, op_id) *
                                     (fu(i, 0, op_id) +
                                      Tsd * (w.x * fu(i, 0, op_dx) +
                                             w.y * fu(i, 0, op_dy))) +
                                 epsilon_diffusion * fu(i, 0, op_dx) *
                                     fu(j, 0, op_dx) +
                                 epsilon_diffusion * fu(i, 0, op_dy) *
                                     fu(j, 0, op_dy));
                  }
                  (*this)(FK.loc2glb(i)) +=
                      Cint * ((1. / dt) * ux - signU0 * normDupHalf + ff) *
                      (fu(i, 0, op_id) +
                       Tsd * (w.x * fu(i, 0, op_dx) + w.y * fu(i, 0, op_dy)));
               }
            }

            // // explicit Euler
            else {
               for (int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
                  for (int j = FK.dfcbegin(0); j < FK.dfcend(0); ++j) {
                     (*this)(FK.loc2glb(i), FK.loc2glb(j)) +=
                         Cint * ((1. / dt) * fu(j, 0, op_id) *
                                     (fu(i, 0, op_id) +
                                      Tsd * (w.x * fu(i, 0, op_dx) +
                                             w.y * fu(i, 0, op_dy))) +
                                 epsilon_diffusion * fu(i, 0, op_dx) *
                                     fu(j, 0, op_dx) +
                                 epsilon_diffusion * fu(i, 0, op_dy) *
                                     fu(j, 0, op_dy));
                  }
                  (*this)(FK.loc2glb(i)) +=
                      Cint *
                      ((1. / dt) * ux - (w, dup) + ff)
                      // Cint*((1./dt)*ux  - signU0*normDup + ff)
                      * (fu(i, 0, op_id) +
                         Tsd * (w.x * fu(i, 0, op_dx) + w.y * fu(i, 0, op_dy)));
               }
            }
         }
      }
   }
}

template <>
void Reinitialization<Mesh3>::assembly(const Fun &up, const Fun &u0) {
   typedef Mesh3 Mesh;
   typedef typename Mesh::Partition Partition;
   typedef typename QF::QuadraturePoint QuadraturePoint;

   RNMK_ fu(this->databf, ndofK, 1, op_all); //  the value for basic fonction

   for (int k = Vh.first_element(); k < Vh.last_element();
        k += Vh.next_element()) {

      const FElement &FK(Vh[k]);
      const int kb = Vh.Th(FK.T);
      CutData cutData(interface.getCutData(kb));        // get the cut data
      const Partition &cutK = Partition(FK.T, cutData); // build the cut

      const R h = FK.T.hElement();

      for (typename Partition::const_element_iterator it = cutK.element_begin();
           it != cutK.element_end(); ++it) {

         const R meas        = cutK.mesure(it);
         const double signU0 = cutK.whatSign(it);
         R Tsd               = 1. / (sqrt(1. / dt / dt + 1. / h / h));

         for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

            QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip       = cutK.toK(it, ip);
            const R Cint = meas * ip.getWeight();

            FK.BF(FK.T.toKref(mip),
                  fu); // need point in local reference element

            R ux = up.eval(k, mip, 0, op_id);
            Rd dup;
            dup.x      = up.eval(k, mip, 0, op_dx);
            dup.y      = up.eval(k, mip, 0, op_dy);
            dup.z      = up.eval(k, mip, 0, op_dz);
            R normDup  = (fabs(dup.norm()) < 1e-14) ? 1e-9 : dup.norm();
            const Rd w = dup / normDup * signU0;
            const R ff = signU0;

            // // explicit Euler
            for (int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
               for (int j = FK.dfcbegin(0); j < FK.dfcend(0); ++j) {
                  (*this)(FK.loc2glb(i), FK.loc2glb(j)) +=
                      Cint *
                      ((1. / dt) * fu(j, 0, op_id) *
                           (fu(i, 0, op_id) + Tsd * (w.x * fu(i, 0, op_dx) +
                                                     w.y * fu(i, 0, op_dy) +
                                                     w.z * fu(i, 0, op_dz))) +
                       epsilon_diffusion * fu(i, 0, op_dx) * fu(j, 0, op_dx) +
                       epsilon_diffusion * fu(i, 0, op_dy) * fu(j, 0, op_dy) +
                       epsilon_diffusion * fu(i, 0, op_dz) * fu(j, 0, op_dz));
               }
               (*this)(FK.loc2glb(i)) +=
                   Cint * ((1. / dt) * ux - (w, dup) + ff) *
                   (fu(i, 0, op_id) +
                    Tsd * (w.x * fu(i, 0, op_dx) + w.y * fu(i, 0, op_dy) +
                           w.z * fu(i, 0, op_dz)))
                   //   Cint*((1./dt)*ux * fu(i,0,op_id))
                   // + Cint*(- (w,dup) + ff)*(fu(i,0,op_id) + Tsd *
                   // (w.x*fu(i,0,op_dx) + w.y*fu(i,0,op_dy)+
                   // w.z*fu(i,0,op_dz)))
                   ;
            }
         }
      }
   }
}
