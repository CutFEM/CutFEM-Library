#include "limiter.hpp"



// void Limiter::KXRCF_indicator(const Fun2_h& uh, const Fun2_h& flux){
//   const FESpace2& Vh(*uh.Vh);
//   indicator.init(Vh.NbElement()); indicator = 0.;
//   Rn data_send(Vh.NbElement()); data_send = 1e+300;
//   // const QFB& qfb(*QF_Simplex<R1>(3));
//   // const QF& qf(*QF_Simplex<R2>(3));
//
//
//   for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
//
//     const int kb = Vh.idxElementInBackMesh(k);
//     const FElement& FK(Vh[k]);
//     int the_domain = FK.whichDomain();
//
//     const R2 x_mid = FK.T(R2(1./3, 1./3));
//     R2 local_flux;
//     local_flux.x = flux.evalOnBackMesh(kb, x_mid, 0, op_id, the_domain);
//     local_flux.y = flux.evalOnBackMesh(kb, x_mid, 1, op_id, the_domain);
//
//     const R h = FK.T.hMax() / 2; // lenght of the radius of the circomscribed triangle
//     double measB = 1e-15;
//     double maxQj = 1e-15;
//
//     // compute the local norm
//     for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//       typename QF::QuadraturePoint ip(qf[ipq]); // integration point
//       const R2 mip = FK.T((R2)ip);
//       maxQj = max(maxQj, fabs(uh.eval(k, mip, 0, op_id)));
//     }
//
//     for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces
//
//       const R2 normal = FK.T.N(ifac);
//       if( (normal, local_flux) > 0 ) continue;
//
//       int ifacn = ifac;
//       int kn = Vh.getNeighborElement(k, ifacn, the_domain);
//       if(kn == -1) continue;         // border edge
//       const FElement& FKn(Vh[kn]);
//
//       const R meas    = FK.T.mesureBord(ifac);
//       measB += meas;
//
//       for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
//
//         typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
//         const R2 mip = FK.T(FK.T.toKref((R1)ip, ifac));
//         const R Cint = meas * ip.getWeight();
//
//         indicator(k) += Cint * (uh.eval(k, mip, 0, op_id) - uh.eval(kn, mip, 0, op_id));
//       }
//     }
//     data_send(k) = fabs(indicator(k))/ (h*measB*maxQj);
//     // indicator(k) = fabs(indicator(k))/ (h*measB*maxQj);
//   }
//   MPIcf::AllReduce(data_send, indicator, MPI_MIN);
// }
//
// void Limiter::min_and_max_average_neighbor(const FElement& FK, const Fun2_h&uh, double&m, double&M){
//   m=1e300, M=-1e300;
//   int the_domain = FK.whichDomain();
//   for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces
//     int ifacn = ifac;
//     int kn = FK.Vh.getNeighborElement(FK.index(), ifacn, the_domain);
//     if(kn == -1) continue;         // border edge
//     const FElement& FKn(FK.Vh[kn]);
//     double val = 0.;
//     for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//       typename QF::QuadraturePoint ip(qf[ipq]); // integration point
//       const R2 mip = FKn.T((R2)ip);
//       val += ip.getWeight() * uh.eval(kn, mip, 0, op_id);
//     }
//     m = min(m, val);
//     M = max(M, val);
//   }
// }
//
// double Limiter::average_on_Kj(const FElement& FK, const Fun2_h&uh){
//   double Uj = 0.;
//   for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//     typename QF::QuadraturePoint ip(qf[ipq]); // integration point
//     const R2 mip = FK.T((R2)ip);
//     Uj += ip.getWeight() * uh.eval(FK.index(), mip, 0, op_id);
//   }
//   return Uj;
// }
//
// double Limiter::compute_alpha(const FElement& FK, const Fun2_h&uh, double Uj, double m, double M) {
//   double min_alpha_i = 1e300;
//   int the_domain = FK.whichDomain();
//
//   for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces
//     int ifacn = ifac;
//     int kn = FK.Vh.getNeighborElement(FK.index(), ifacn, the_domain);
//     if(kn == -1) continue;
//     double alpha_i = 1.;
//     for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
//       typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
//       const R2 mip = FK.T(FK.T.toKref((R1)ip, ifac));
//       // compute U_j(x_i^*)
//       double ux = uh.eval(FK.index(), mip, 0, op_id);
//       if(ux - M > 0){
//         alpha_i = (M-Uj)/(ux-Uj);
//       }
//       else if(ux - m < 0){
//         alpha_i = (m-Uj)/(ux-Uj);
//       }
//       // std::cout << alpha_i << std::endl;
//       alpha_i = max(alpha_i, 0.);
//       min_alpha_i = min(min_alpha_i, alpha_i);
//       assert(min_alpha_i < 1e300);
//     }
//   }
//   return min_alpha_i;
// }
//
// void Limiter::limiting(const Fun2_h& uh, Rn& lh){
//
//   const FESpace2& Vh(*uh.Vh);
//   assert(Vh.NbElement() == indicator.size());
//   lh.init(uh.v);
//   Rn data_send(lh.size()); data_send = 1e+300;
//
//   // for(int k=0; k<indicator.size(); ++k){
//   for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
//
//     const FElement& FK(Vh[k]);
//     if(indicator(k) < 1) {
//       for(int df=FK.dfcbegin(0);df<FK.dfcend(0);++df){
//         data_send(FK.loc2glb(df)) = uh.v[FK.loc2glb(df)];
//       }
//       continue;
//     }
//
//     //compute m and M , min max average of neighbor
//     double m, M;
//     min_and_max_average_neighbor(FK, uh, m, M);
//
//     // compute solution average on K_j
//     double Uj = average_on_Kj(FK, uh);
//
//     // compute the min max (alpha_i , 0)
//     double min_alpha_i = compute_alpha(FK, uh, Uj, m, M);
//
//     for(int df=FK.dfcbegin(0);df<FK.dfcend(0);++df){
//       // lh(FK.loc2glb(df)) *= min_alpha_i;
//       // if(fabs(lh(FK.loc2glb(df))) < 1e-16){ lh(FK.loc2glb(df)) = 1e-16;}
//       int alpha = (df == 0)? 1 : min_alpha_i;
//       data_send(FK.loc2glb(df)) = uh.v[FK.loc2glb(df)]*alpha;
//     }
//   }
//   // std::cout << uh.v << std::endl;
//   MPIcf::AllReduce(data_send, lh, MPI_MIN);
// }
//
//
// void Filter_Box::cut_triangle(const typename Mesh2::Element& K) {
//   // initialize with the full triangle
//   for(int i=0; i<3;++i) {
//     vertices_.push_back((R2)(K[i]));
//   }
//   triangles_.push_back(TriaIdx( 0, 1, 2));
//   std::vector<TriaIdx> temp_triangles_;
//   int iii = 0;
//   for(auto it=sides_.begin(); it != sides_.end();++it) {
//     temp_triangles_ = triangles_;
//     triangles_.clear();
//     for(int k=0; k<temp_triangles_.size();++k){
//
//       // check if triIdx is cut
//       bool cutTriangle = false;
//       int idx_edge_cut[2];
//       int iee = 0;
//       for(int ie=0;ie<3;++ie) {
//         int i0 = temp_triangles_[k][Element::nvedge[ie][0]];
//         int i1 = temp_triangles_[k][Element::nvedge[ie][1]];
//         if(it->cut_edge(vertices_[i0], vertices_[i1])) {
//           cutTriangle = true;
//           idx_edge_cut[iee++] = ie;
//           std::cout << "side " << iii << "  find a cut edge" << std::endl;
//         }
//       }
//       // triangle is cut so need to create new one
//       if(cutTriangle){
//         assert(iee == 2); // found 2 cuts
//         int cut[2];  // max 2
//         int idx_nodes[2]; // max 2
//         // find cut nodes.
//         for(int ie=0;ie<2;++ie) {
//           int ie_cut = idx_edge_cut[ie];
//           int i0 = temp_triangles_[k][Element::nvedge[ie_cut][0]];
//           int i1 = temp_triangles_[k][Element::nvedge[ie_cut][1]];
//           cut[ie] = vertices_.size();
//           vertices_.push_back(it->get_cut_node(vertices_[i0], vertices_[i1]));
//         }
//         // find nodes in domain ( 1 or 2)
//         int nb_node = 0;
//         for(int i=0;i<3;++i) {
//           int i0 = temp_triangles_[k][i];
//           if(it->is_negative(vertices_[i0])) {
//             idx_nodes[nb_node] = i0;
//             nb_node++;
//           }
//         }
//         const int v = Element::commonVertOfEdges[idx_edge_cut[0]][idx_edge_cut[1]];
//         if(nb_node == 1){
//           // ADD TRIANGLE
//           triangles_.push_back(TriaIdx(idx_nodes[0], cut[0], cut[1]));
//         }
//         else {
//           // ADD 2 TRIANGLES
//           const int i0 = temp_triangles_[k][Element::oppVertOfEdge( idx_edge_cut[0], v)];
//           const int i1 = cut[0];
//           const int i2 = temp_triangles_[k][Element::oppVertOfEdge( idx_edge_cut[1], v)];
//           const int i3 = cut[1];
//           triangles_.push_back(TriaIdx(i0, i1, i2));
//           triangles_.push_back(TriaIdx(i1, i3, i2));
//         }
//
//       }
//       else{  // not cut triangle => put it back in array
//         triangles_.push_back(temp_triangles_[k]);
//       }
//     }
//     iii++;
//   }
// }
//
// void Filter_Box::print(std::string filename ) const {
//   std::ofstream plot;
//   plot.open(filename.c_str(), std::ofstream::out);
//   for(int k=0; k<triangles_.size();++k) {
//     for(int i=0;i<3;++i) {
//       plot << vertices_[triangles_[k][i]] << std::endl;
//     }
//     plot << vertices_[triangles_[k][0]] << std::endl;
//     plot << std::endl;
//     plot << std::endl;
//   }
// plot.close();
// }
