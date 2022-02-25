#include "marker.hpp"

#include "../util/util.hpp"
#include "GenericInterface.hpp"

#include "../FESpace/expression.hpp"


static bool point_inside_tri(R2 s, R2 a, R2 b, R2 c) {

  int as_x = s.x-a.x;
  int as_y = s.y-a.y;

  bool s_ab = (b.x-a.x)*as_y-(b.y-a.y)*as_x >= 0;
  if(((c.x-a.x)*as_y-(c.y-a.y)*as_x >= 0) == s_ab) return false;

  if(((c.x-b.x)*(s.y-b.y)-(c.y-b.y)*(s.x-b.x) >= 0) != s_ab) return false;

  return true;
}

static bool point_inside_tri(R2 s,  const typename Mesh2::Element& K) {
  R2 Pref = K.toKref(s);
  // std::cout << Pref << std::endl;
  if(Pref.x < -Epsilon || Pref.x > 1+Epsilon) return false;
  if(Pref.y < -Epsilon || Pref.y > 1+Epsilon) return false;
  if(Pref.y-Epsilon > 1 - Pref.x) return false;

  return true;

}

static int find_triangle_belong_point(const Mesh2& Th, R2 s) {
  typedef typename Mesh2::Element Element;

  for( int k=0; k<Th.nbElmts();++k) {
    const Element& K(Th[k]);
    // if(point_inside_tri(s, K[0], K[1], K[2])) {
    if(point_inside_tri(s, K)) {

      return k;
    }
  }
  return -1;
}

static bool check_intersect(R2 u, R2 v, R2 w, R& t) {

  R det = -u.x*v.y + u.y*v.x;
  if(fabs(det) < 1e-14) return false;

  t = -w[0]*v[1] + w[1]*v[0];
  t /= det;
  R t2 = -w[0]*u[1] + w[1]*u[0];
  t2 /= det;
  if(t >= 0 && t <= 1 && t2 >=0 && t2 <= 1) return true;
  return false;
}


Marker::Marker(const Mesh2& Thh) : Interface2(Thh), Th(Thh) { }
Marker::Marker(const Mesh2& Thh, R2(*fparam)(double t), double x_begin, double x_end, int npoint) :Marker(Thh) {
  {
    nMarker_begin.push_back(0);


    double h = (x_end - x_begin) / (npoint - 1);
    for(int i=0;i<npoint;++i) {
      double t = x_begin + i*h;
      R2 val = fparam(t);
      T.push_back(t);
      X.push_back(val.x);
      Y.push_back(val.y);
      markers.push_back(val);
    }
    nMarker_end.push_back(markers.size());
  }
  // check periodicity
  R2 AB = markers[markers.size()-1] - markers[0];
  periodic = (AB.norm() < 1e-12);
}
Marker::Marker(const Mesh2& Thh, std::string path) : Interface2(Thh), Th(Thh) {
  add(path);
}


R2 Marker::make_normal (int i, int j) {
  Rd a = vertices_[i];
  Rd b = vertices_[j];
  Rd normal_ls(a.y-b.y, b.x-a.x);
  normal_ls /= normal_ls.norm();
  return -normal_ls;
}

void Marker::build_sign_array() {

  vector<byte> element_seen(Th.nt, -1);

  // set all element as Omega2
  ls_sign.resize(Th.nv); ls_sign = -1;

  // start from a border element to find Omega1 elements
  int ifaceK;
  int k = Th.BoundaryElement(0, ifaceK);

  visit_element_sign(k, element_seen);

  element_seen.clear();

};

void Marker::visit_element_sign(int k, vector<byte>& elementSeen ) {

  elementSeen[k] = 1;
  for(int i=0;i<3;++i) {
    ls_sign(Th(k, i)) = 1;
  }

  for(int e=0;e<3;++e) {

    int ib=e;
    int k_next = Th.ElementAdj(k,ib);
    if(k_next == -1 || isCut(k_next) || elementSeen[k_next] == 1) continue;

    visit_element_sign(k_next,elementSeen);
  }
}

void Marker::move(const FunFEMVirtual& uh, double dt) {
  int i = 0;
  for(auto it = markers.begin(); it!= markers.end();++it, ++i) {
    int idxK_bm = elementOfMarker[i];  // get index in backMesh
    int idxK = uh.idxElementFromBackMesh(idxK_bm);  // index for function evaluation
    R2 val(0.,0.);
    for(int i=0; i<2; ++i) {
      val[i] = uh.eval(idxK, *it, i, op_id);
    }
    *it += dt*val;
  }
}

R2 Marker::get_intersect_edge(R2 A, R2 B, int previousK, int k, int& k_next) {
  const Element& K(Th[k]);

  for(int i=0;i<Element::ne;++i) {
    int jn = i;
    k_next = Th.ElementAdj(k, jn);
    if(k_next == -1) continue;
    if(k_next == previousK) continue;

    R2 C = K[Element::nvedge[i][0]];
    R2 D = K[Element::nvedge[i][1]];
    R2 AB(A,B), CD(C,D), AC(A,C);
    R t;
    if(check_intersect(AB, CD, AC, t)) {
      R2 P = A+t*AB;
      return P;
    }
  }
  assert(0);
  return R2();
}

int Marker::find_next_element(R2 A, R2 B, int previousK, int k) {
  const Element& K(Th[k]);

  for(int i=0;i<Element::ne;++i) {
    int jn = i;
    int k_next = Th.ElementAdj(k, jn);
    if(k_next == -1) continue;
    if(k_next == previousK) continue;

    R2 C = K[Element::nvedge[i][0]];
    R2 D = K[Element::nvedge[i][1]];
    R2 AB(A,B), CD(C,D), AC(A,C);
    R t;
    if(check_intersect(AB, CD, AC, t)) {
      return k_next;
    }
  }
  assert(0);
  return -1;
}

void Marker::find_marker_limit(int k, int& i){
  const Element& K(Th[k]);
  while(true){
    R2 node = markers[i];
    if(point_inside_tri(node, K)) {
      elementOfMarker.push_back(k);
      i++;
    }
    else break;
  }
}

R2 Marker::find_intersection(R2 A, R2 B, int k, int kn, int& ie){
  const Element& K(Th[k]);
  for(int i=0;i<Element::ne;++i) {
    ie = i;
    int jn = i;
    int k_next = Th.ElementAdj(k, jn);
    if(k_next != kn) continue;

    R2 C = K[Element::nvedge[i][0]];
    R2 D = K[Element::nvedge[i][1]];
    R2 AB(A,B), CD(C,D), AC(A,C);
    R t;
    assert(check_intersect(AB, CD, AC, t));
    return A+t*AB;
  }
}

int Marker::find_edge(int k, int kn) {
  const Element& K(Th[k]);
  int ie;
  for(int i=0;i<Element::ne;++i) {
    ie = i;
    int jn = i;
    int k_next = Th.ElementAdj(k, jn);
    if(k_next == kn) return ie;
  }
}

void Marker::find_vertices() {
  // Find the nodes that intersect the edges
  // nstartMarker.push_back(edges_node.size());
  // nstartFace.push_back(faces.size());
  int nbeg = nMarker_begin[0];
  int nlast = nMarker_end[0];


  // get the starting triangle
  R2 firstNode = markers[nbeg];
  int i1 = vertices_.size(), i2 = 0;
  int k0 = find_triangle_belong_point(Th, firstNode);
  elementOfMarker.push_back(k0);
  assert(k0 != -1);

  int idxK = k0;
  int nextK = -1;
  int previousK = -1;
  int i=nbeg+1;
  int nloop = 0;
  while(i<nlast){

    find_marker_limit(idxK, i);
    assert(i<nlast);
    nextK = find_next_element(markers[i-1], markers[i], previousK, idxK);
    previousK = idxK;
    if(nextK == k0) break;

    idxK = nextK;
    int ed1;
    R2 P1 = find_intersection(markers[i-1], markers[i], idxK, previousK, ed1);

    int j = i;
    find_marker_limit(idxK, j);
    if(j>= nlast) break;
    assert(j<nlast);

    nextK = find_next_element(markers[j-1], markers[j], previousK, idxK);
    int ed2;
    R2 P2 = find_intersection(markers[j-1], markers[j], idxK, nextK, ed2);



    if (nloop == 0) {
      vertices_.push_back(P1);
    }
    vertices_.push_back(P2);
    if(ed1 < ed2) {
      faces_.push_back(Face(vertices_.size()-2,vertices_.size()-1, 0));
    }
    else {
      faces_.push_back(Face(vertices_.size()-1,vertices_.size()-2, 0));
    }
    // if(ed1 < ed2) {
    //   vertices_.push_back(P1);
    //   vertices_.push_back(P2);
    // } else {
    //   vertices_.push_back(P2);
    //   vertices_.push_back(P1);
    // }


    element_of_face_.push_back(idxK);
    face_of_element_[idxK] = element_of_face_.size()-1;

    // faces_.push_back(Face(vertices_.size()-2,vertices_.size()-1, 0));
    outward_normal_.push_back(make_normal(vertices_.size()-2,vertices_.size()-1));
    nloop += 1;
  }

  // Check boundary


  // Find last face
  if(!periodic) return;
  int ed1,ed2;
  R2 P1 = find_intersection(markers[i-1], markers[i], k0, previousK, ed1);
  int j=0;
  find_marker_limit(k0, j);
  nextK = element_of_face_[0];
  R2 P2 = find_intersection(markers[j-1], markers[j], k0, nextK, ed2);



  // if(ed1 < ed2) {
  //   vertices_.push_back(P1);
  //   vertices_.push_back(P2);
  // } else {
  //   vertices_.push_back(P2);
  //   vertices_.push_back(P1);
  // }
  vertices_.push_back(P2);
  if(ed1 < ed2) {
    faces_.push_back(Face(vertices_.size()-2,vertices_.size()-1, 0));
  }
  else {
    faces_.push_back(Face(vertices_.size()-1,vertices_.size()-2, 0));
  }
  element_of_face_.push_back(k0);
  face_of_element_[k0] = element_of_face_.size()-1;
  // faces_.push_back(Face(vertices_.size()-2,vertices_.size()-1, 0));
  outward_normal_.push_back(make_normal(vertices_.size()-2,vertices_.size()-1));


}

void Marker::make_interface(bool leftB, bool rightB){

  vertices_.resize(0);
  element_of_face_.resize(0);
  face_of_element_.clear();
  faces_.resize( 0);                          // reinitialize arrays
  outward_normal_.resize(0);
  is_boundary_face_.resize( 0);

  // get vertices & element_of_face_ & face_of_element_ & faces_ & outward_normal
  find_vertices();

  if(leftB && (fabs(vertices_[0].x - markers[0].x) > 1e-12
           ||  fabs(vertices_[0].y - markers[0].y) > 1e-12) ) {

    int k0 = elementOfMarker[0];
    int ed1 = find_edge(k0, -1);
    R2 P1 = markers[0];
    int ed2 = find_edge(k0, element_of_face_[0]);
    R2 P2 = markers[markers.size()-1];
    if(ed1 < ed2) {
      vertices_.push_back(P1);
      vertices_.push_back(P2);
    } else {
      vertices_.push_back(P2);
      vertices_.push_back(P1);
    }
    element_of_face_.push_back(k0);
    face_of_element_[k0] = element_of_face_.size()-1;
    faces_.push_back(Face(vertices_.size()-2,vertices_.size()-1, 0));
    outward_normal_.push_back(make_normal(vertices_.size()-2,vertices_.size()-1));

  }



  // need to get a sign array
  build_sign_array();
};

void Marker::add(std::string path) {

  nMarker_begin.push_back(markers.size());

  int nloc;
  ifstream f(path.c_str());
  R2 P;
  if(!f) {cerr << " cannot open marker file " << path << endl; exit(1);}
  cout << " Read On file \"" << path <<"\""<<  endl;
  f >> nloc;
  for (int i=0;i<nloc;i++) {
    f >> P;
    markers.push_back(P);
    assert(f.good());
  }

  nMarker_end.push_back(markers.size());
}
























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
