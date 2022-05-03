#include "marker.hpp"
#include "../util/util.hpp"
#include "../FESpace/expression.hpp"


//
// Marker::Marker(const Mesh2& Thh) : Th_(Thh) { }
Marker::Marker(const Mesh2& Thh, R2(*fparam)(double t), double x_begin, double x_end, int npoint) : Th_(Thh) {
  T_.reserve(npoint);
  X_.reserve(npoint);
  Y_.reserve(npoint);

  double h = (x_end - x_begin) / (npoint - 1);
  for(int i=0;i<npoint;++i) {
    double t = x_begin + i*h;
    R2 val = fparam(t);
    this->add(t, val);
  }
  // check periodicity
  int l = X_.size()-1;
  R2 AB(X_[l] - X_[0], Y_[l]-Y_[0]);;
  periodic_ = (AB.norm() < 1e-12);

  // Find element containing markers.
  element_of_marker_.reserve(this->size());
  int k = 0;
  for(int e=0;e<this->size();++e) {
    R2 P = this->get_marker(e);
    int k = geometry::find_triangle_contenant_p(Th_, P , k);
    element_of_marker_.push_back(k);
  }
  assert(element_of_marker_.size() == this->size());
  if(periodic_) assert(element_of_marker_[0] == element_of_marker_[this->size()-1]);

}

void Marker::add(double t, R2 val) {
  T_.push_back(t);
  X_.push_back(val.x);
  Y_.push_back(val.y);
}

void Marker::move(const FunFEMVirtual& uh, double dt) {
  int k = 0;
  for(auto k=0; k<this->size(); ++k) {
    int kb = element_of_marker_[k];  // get index in backMesh
    R2 val(0.,0.);
    R2 P = this->get_marker(k);
    for(int j=0; j<2; ++j) {
      val[j] = uh.evalOnBackMesh(kb, 0, P, j, op_id);
    }
    P += dt*val;
    this->set(k, P);
  }
}




// Marker::Marker(const Mesh2& Thh, std::string path) : Th_(Thh) {
//   add(path);
// }

//
// R2 Marker::make_normal (int i, int j) {
//   Rd a = vertices_[i];
//   Rd b = vertices_[j];
//   Rd normal_ls(a.y-b.y, b.x-a.x);
//   normal_ls /= normal_ls.norm();
//   return -normal_ls;
// }
//
// void Marker::build_sign_array() {
//
//   vector<byte> element_seen(Th.nt, -1);
//
//   // set all element as Omega2
//   ls_sign.resize(Th.nv); ls_sign = -1;
//
//   // start from a border element to find Omega1 elements
//   int ifaceK;
//   int k = Th.BoundaryElement(0, ifaceK);
//
//   visit_element_sign(k, element_seen);
//
//   element_seen.clear();
//
// };
//
// void Marker::visit_element_sign(int k, vector<byte>& elementSeen ) {
//
//   elementSeen[k] = 1;
//   for(int i=0;i<3;++i) {
//     ls_sign(Th(k, i)) = 1;
//   }
//
//   for(int e=0;e<3;++e) {
//
//     int ib=e;
//     int k_next = Th.ElementAdj(k,ib);
//     if(k_next == -1 || isCut(k_next) || elementSeen[k_next] == 1) continue;
//
//     visit_element_sign(k_next,elementSeen);
//   }
// }
//

//
// R2 Marker::get_intersect_edge(R2 A, R2 B, int previousK, int k, int& k_next) {
//   assert(0);
//   // const Element& K(Th[k]);
//   //
//   // for(int i=0;i<Element::ne;++i) {
//   //   int jn = i;
//   //   k_next = Th.ElementAdj(k, jn);
//   //   if(k_next == -1) continue;
//   //   if(k_next == previousK) continue;
//   //
//   //   R2 C = K[Element::nvedge[i][0]];
//   //   R2 D = K[Element::nvedge[i][1]];
//   //   R2 AB(A,B), CD(C,D), AC(A,C);
//   //   R t;
//   //   if(check_intersect(AB, CD, AC, t)) {
//   //     R2 P = A+t*AB;
//   //     return P;
//   //   }
//   // }
//   // assert(0);
//   return R2();
// }
//
// int Marker::find_next_element(R2 A, R2 B, int previousK, int k) {
//   // const Element& K(Th[k]);
//   //
//   // for(int i=0;i<Element::ne;++i) {
//   //   int jn = i;
//   //   int k_next = Th.ElementAdj(k, jn);
//   //   if(k_next == -1) continue;
//   //   if(k_next == previousK) continue;
//   //
//   //   R2 C = K[Element::nvedge[i][0]];
//   //   R2 D = K[Element::nvedge[i][1]];
//   //   R2 AB(A,B), CD(C,D), AC(A,C);
//   //   R t;
//   //   if(check_intersect(AB, CD, AC, t)) {
//   //     return k_next;
//   //   }
//   // }
//   assert(0);
//   return -1;
// }
//
// void Marker::find_marker_limit(int k, int& i){
//   assert(0);
//   // const Element& K(Th[k]);
//   // while(true){
//   //   R2 node = markers[i];
//   //   if(point_inside_tri(node, K)) {
//   //     elementOfMarker.push_back(k);
//   //     i++;
//   //   }
//   //   else break;
//   // }
// }
//
// R2 Marker::find_intersection(R2 A, R2 B, int k, int kn, int& ie){
//   // const Element& K(Th[k]);
//   // for(int i=0;i<Element::ne;++i) {
//   //   ie = i;
//   //   int jn = i;
//   //   int k_next = Th.ElementAdj(k, jn);
//   //   if(k_next != kn) continue;
//   //
//   //   R2 C = K[Element::nvedge[i][0]];
//   //   R2 D = K[Element::nvedge[i][1]];
//   //   R2 AB(A,B), CD(C,D), AC(A,C);
//   //   R t;
//   //   assert(check_intersect(AB, CD, AC, t));
//   //   return A+t*AB;
//   // }
//   assert(0); return R2();
// }
//
// int Marker::find_edge(int k, int kn) {
//   const Element& K(Th[k]);
//   int ie;
//   for(int i=0;i<Element::ne;++i) {
//     ie = i;
//     int jn = i;
//     int k_next = Th.ElementAdj(k, jn);
//     if(k_next == kn) return ie;
//   }
// }
//
// void Marker::find_vertices() {
//   assert(0);
//   // // Find the nodes that intersect the edges
//   // // nstartMarker.push_back(edges_node.size());
//   // // nstartFace.push_back(faces.size());
//   // int nbeg = nMarker_begin[0];
//   // int nlast = nMarker_end[0];
//   //
//   //
//   // // get the starting triangle
//   // R2 firstNode = markers[nbeg];
//   // int i1 = vertices_.size(), i2 = 0;
//   // int k0 = find_triangle_belong_point(Th, firstNode);
//   // elementOfMarker.push_back(k0);
//   // assert(k0 != -1);
//   //
//   // int idxK = k0;
//   // int nextK = -1;
//   // int previousK = -1;
//   // int i=nbeg+1;
//   // int nloop = 0;
//   // while(i<nlast){
//   //
//   //   find_marker_limit(idxK, i);
//   //   assert(i<nlast);
//   //   nextK = find_next_element(markers[i-1], markers[i], previousK, idxK);
//   //   previousK = idxK;
//   //   if(nextK == k0) break;
//   //
//   //   idxK = nextK;
//   //   int ed1;
//   //   R2 P1 = find_intersection(markers[i-1], markers[i], idxK, previousK, ed1);
//   //
//   //   int j = i;
//   //   find_marker_limit(idxK, j);
//   //   if(j>= nlast) break;
//   //   assert(j<nlast);
//   //
//   //   nextK = find_next_element(markers[j-1], markers[j], previousK, idxK);
//   //   int ed2;
//   //   R2 P2 = find_intersection(markers[j-1], markers[j], idxK, nextK, ed2);
//   //
//   //
//   //
//   //   if (nloop == 0) {
//   //     vertices_.push_back(P1);
//   //   }
//   //   vertices_.push_back(P2);
//   //   if(ed1 < ed2) {
//   //     faces_.push_back(Face(vertices_.size()-2,vertices_.size()-1, 0));
//   //   }
//   //   else {
//   //     faces_.push_back(Face(vertices_.size()-1,vertices_.size()-2, 0));
//   //   }
//   //   // if(ed1 < ed2) {
//   //   //   vertices_.push_back(P1);
//   //   //   vertices_.push_back(P2);
//   //   // } else {
//   //   //   vertices_.push_back(P2);
//   //   //   vertices_.push_back(P1);
//   //   // }
//   //
//   //
//   //   element_of_face_.push_back(idxK);
//   //   face_of_element_[idxK] = element_of_face_.size()-1;
//   //
//   //   // faces_.push_back(Face(vertices_.size()-2,vertices_.size()-1, 0));
//   //   outward_normal_.push_back(make_normal(vertices_.size()-2,vertices_.size()-1));
//   //   nloop += 1;
//   // }
//   //
//   // // Check boundary
//   //
//   //
//   // // Find last face
//   // if(!periodic) return;
//   // int ed1,ed2;
//   // R2 P1 = find_intersection(markers[i-1], markers[i], k0, previousK, ed1);
//   // int j=0;
//   // find_marker_limit(k0, j);
//   // nextK = element_of_face_[0];
//   // R2 P2 = find_intersection(markers[j-1], markers[j], k0, nextK, ed2);
//   //
//   //
//   //
//   // // if(ed1 < ed2) {
//   // //   vertices_.push_back(P1);
//   // //   vertices_.push_back(P2);
//   // // } else {
//   // //   vertices_.push_back(P2);
//   // //   vertices_.push_back(P1);
//   // // }
//   // vertices_.push_back(P2);
//   // if(ed1 < ed2) {
//   //   faces_.push_back(Face(vertices_.size()-2,vertices_.size()-1, 0));
//   // }
//   // else {
//   //   faces_.push_back(Face(vertices_.size()-1,vertices_.size()-2, 0));
//   // }
//   // element_of_face_.push_back(k0);
//   // face_of_element_[k0] = element_of_face_.size()-1;
//   // // faces_.push_back(Face(vertices_.size()-2,vertices_.size()-1, 0));
//   // outward_normal_.push_back(make_normal(vertices_.size()-2,vertices_.size()-1));
//
//
// }
//
// void Marker::make_interface(bool leftB, bool rightB){
//
//   vertices_.resize(0);
//   element_of_face_.resize(0);
//   face_of_element_.clear();
//   faces_.resize( 0);                          // reinitialize arrays
//   outward_normal_.resize(0);
//   is_boundary_face_.resize( 0);
//
//   // get vertices & element_of_face_ & face_of_element_ & faces_ & outward_normal
//   find_vertices();
//
//   if(leftB && (fabs(vertices_[0].x - markers[0].x) > 1e-12
//            ||  fabs(vertices_[0].y - markers[0].y) > 1e-12) ) {
//
//     int k0 = elementOfMarker[0];
//     int ed1 = find_edge(k0, -1);
//     R2 P1 = markers[0];
//     int ed2 = find_edge(k0, element_of_face_[0]);
//     R2 P2 = markers[markers.size()-1];
//     if(ed1 < ed2) {
//       vertices_.push_back(P1);
//       vertices_.push_back(P2);
//     } else {
//       vertices_.push_back(P2);
//       vertices_.push_back(P1);
//     }
//     element_of_face_.push_back(k0);
//     face_of_element_[k0] = element_of_face_.size()-1;
//     faces_.push_back(Face(vertices_.size()-2,vertices_.size()-1, 0));
//     outward_normal_.push_back(make_normal(vertices_.size()-2,vertices_.size()-1));
//
//   }
//
//
//
//   // need to get a sign array
//   build_sign_array();
// };
//
// void Marker::add(std::string path) {
//
//   nMarker_begin.push_back(markers.size());
//
//   int nloc;
//   ifstream f(path.c_str());
//   R2 P;
//   if(!f) {cerr << " cannot open marker file " << path << endl; exit(1);}
//   cout << " Read On file \"" << path <<"\""<<  endl;
//   f >> nloc;
//   for (int i=0;i<nloc;i++) {
//     f >> P;
//     markers.push_back(P);
//     assert(f.good());
//   }
//
//   nMarker_end.push_back(markers.size());
// }
//























//
//
// void Marker::add_interface(int ll) {
//
//   int nbeg = nstart[nbBoundary-1];
//   int nlast = nend[nbBoundary-1];
//
//   // get the starting triangle
//   R2 firstNode = vertices[nbeg];
//   int i1 = edges_node.size(), i2 = 0;
//
//   edges_node.push_back(firstNode);
//   int k0 = find_triangle_belong_point(Th, firstNode);
//   cut_element.push_back(k0);
//   element_seen[k0] = cut_element.size()-1;
//
//   for(int i=nbeg+1; i<nlast;++i) {
//
//     // get two point that change triangle
//     R2 nextNode = vertices[i];
//
//     {
//       const Element& K(Th[k0]);
//
//       if(point_inside_tri(nextNode, K)) {
//         firstNode = nextNode;
//         continue;
//       }
//     }
//
//     while(true) {
//
//       get_intersect_edge(firstNode, nextNode, Th, k0);
//
//       i2 = edges_node.size()-1;
//       faces.push_back(FaceMarker(*this, cut_element[cut_element.size()-2], i1,i2, ll));
//       i1 = i2;
//       const Element& K(Th[k0]);
//
//       if( point_inside_tri(nextNode, K) ) {
//         firstNode = nextNode;
//         break;
//       }
//     }
//
//   }
//   R2 lastNode = vertices[nlast-1];
//   R2 P = *(edges_node.end()-1);
//
//   if(Norme2(lastNode-P) > 1e-12) {
//     R2 lastNode = vertices[nlast-1];
//     i1 = edges_node.size()-1;
//     i2 = edges_node.size();
//     edges_node.push_back(lastNode);
//     // cut_element.push_back(k0);
//     // element_seen(k0) = cut_element.size()-1;
//     std::cout << k0 << std::endl;
//     faces.push_back(FaceMarker(*this, cut_element[cut_element.size()-1], i1,i2, ll));
//   }
//
// }
//
//
// void Marker::getNode(int k, R2 loc_node[2]) const {
//
//   if(isCut(k)) {
//     int iface = element_seen(k);
//     Face face = faces[iface];
//
//     int i1 = face[0];
//     int i2 = face[1];
//
//     loc_node[0] = edges_node[i1];
//     loc_node[1] = edges_node[i2];
//
//   }
// }
//
//
// void Marker::getEdges(int k, int array_edge[3]) const {
//
//   if(isCut(k)) {
//     R2 loc_node[2];
//     const Element& K(Th[k]);
//     int iface = element_seen(k);
//     Face face = faces[iface];
//     int i1 = face[0];
//     int i2 = face[1];
//     loc_node[0] = edges_node[i1];
//     loc_node[1] = edges_node[i2];
//
//     for(int i=0;i<3;++i) {
//       array_edge[i] = -1;
//     }
//
//     for(int i=0;i<2;++i) {
//       R2 Pref = K.toKref(loc_node[i]);
//       if( fabs(Pref.x) < 100*Epsilon)
//       array_edge[1] = i;
//       else if ( fabs(Pref.y) < 100*Epsilon)
//       array_edge[2] = i;
//       else
//       array_edge[0] = i;
//     }
//
//   }
// }
//
//
// void Marker::getSign(int k, byte loc_sign[3]) const {
//
//   for(int i=0;i<Element::nv;++i) {
//     loc_sign[i] = domain_sign(Th(k, i));
//   }
// }
//
//
// CutData2 Marker::getCutData(int k) const {
//   byte loc_sign[3];
//   for(int i=0;i<Element::nv;++i) {
//     loc_sign[i] = domain_sign(Th(k, i));
//   }
//
//   R2 point[2];
//   getNode(k, point);
//
//   int ar_edge[3];
//   getEdges(k, ar_edge);
//
//   if(isCut(k))
//     return CutData2(loc_sign, point, ar_edge);
//   else
//     return CutData2(loc_sign);
// }
//
//
// void Marker::setDomainSign(const KN<double>& lls) {
//
//   copy_levelset_sign( lls, domain_sign);
//   // domain_sign.init(lls);
//
// //   domain_sign = 1;// all in Omega1
// //   list_element1.init(Th.nt);
// //   list_element1 = 0;
// //   list_element2.init(Th.nt);
// //   list_element2 = 0;
//   // std::ofstream plot;
//   // plot.open("domain2.dat" , std::ofstream::out);
//
//
// //   // Find k0 as starting element
// //   R2 node1(0, 2000);
// //   int k0 = find_triangle_belong_point(Th, node1);
// //   addNeighbor(k0, 2);
//
// //   // node1 = R2(1.6e6, -42);
// //   // k0 = find_triangle_belong_point(Th, node1);
// //   // addNeighbor(k0, 2);
//
// //   // Find k0 as starting element
// //   // R2 node1(0, -1000);
// //   // k0 = find_triangle_belong_point(Th, node1);
// //   // addNeighbor(k0,1);
//
// //   // R2 node2(0, 5000);
// //   // k0 = find_triangle_belong_point(Th, node2);
// //   // addNeighbor(k0,1);
//
//
//
//   // for(int i=0;i<Th.nt;++i) {
//   //   const Element& K(Th[i]);
//
//   //   if(domain_sign(Th(i,0)) == -1
//   //      && domain_sign(Th(i,1)) == -1
//   //      && domain_sign(Th(i,2)) == -1
//   //      ) {
//   //     plot << K[0] << "\n"
//   // 	   << K[1] << "\n"
//   // 	   << K[2] << "\n"
//   // 	   << K[0] << "\n \n";
//   //   }
//
//   // }
//
//
// //   // const Element& K(Th[k0]);
// //   // plot << K[0] << "\n"
// //   //      << K[1] << "\n"
// //   //      << K[2] << "\n"
// //   //      << K[0] << "\n \n";
//
//
//
// //   // const Element& K1(Th[k0]);
// //   // plot << K1[0] << "\n"
// //   //      << K1[1] << "\n"
// //   //      << K1[2] << "\n"
// //   //      << K1[0] << "\n \n" ;
//
//   // plot.close();
//
// }
//
// // void Marker::addNeighbor(int k, int dom) {
//
// //   if(dom == 1)list_element1(k) = dom;  // add itself to the list
// //   else list_element2(k) = dom;  // add itself to the list
// //   for(int i=0;i<Element::nv;++i) {
// //     domain_sign(Th(k,i)) = 1;
// //   }
//
// //   for(int i=0;i<Element::ne;++i) {
// //     int jn = i;
// //     int k_next = Th.ElementAdj(k, jn);
// //     std::cout << k << "\t" << k_next << std::endl;
//
//
// //     if(k_next == -1) continue;               // boundary
// //     // std::cout << list_element(k_next) << "\t" << element_seen(k_next) << std::endl;
// //     if(dom == list_element1(k_next) || dom == list_element2(k_next) ) continue;  // already seen
// //     if(element_seen(k_next) != -1) continue; // element cut
//
// //     addNeighbor(k_next, dom);
//
// //   }
// // }
