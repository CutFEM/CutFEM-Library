#include "cut_method.hpp"

template<> int SignElement<Hexa>::nb_node_positif() const{
  int s = util::fsign(sum_);
  return 0.5*(8+sum_) ;
}

template<> RefPatch<Triangle2>::InitializerCL::InitializerCL (){
  if (init_count_++ > 0) return;

  byte ls[3];
  for (ls[0]= -1; ls[0] < 2; ++ls[0]) {
    for (ls[1]= -1; ls[1] < 2; ++ls[1]) {
      for (ls[2]= -1; ls[2] < 2; ++ls[2]) {
        if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0) continue;
        RefPatch<Triangle2>::instance( ls);
      }
    }
	}
}

template<> RefPatch<Quad2>::InitializerCL::InitializerCL (){
  if (init_count_++ > 0) return;

  byte ls[4];
  for (ls[0]= -1; ls[0] < 2; ++ls[0]){
    for (ls[1]= -1; ls[1] < 2; ++ls[1]){
      for (ls[2]= -1; ls[2] < 2; ++ls[2]){
        for (ls[3]= -1; ls[3] < 2; ++ls[3]) {
          if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0) continue;
          RefPatch<Quad2>::instance( ls);
        }
      }
    }
	}
}

template<> RefPatch<Tet>::InitializerCL::InitializerCL () {
  if (init_count_++ > 0) return;

  byte ls[4];
  for (ls[0]= -1; ls[0] < 2; ++ls[0]){
    for (ls[1]= -1; ls[1] < 2; ++ls[1]){
      for (ls[2]= -1; ls[2] < 2; ++ls[2]){
        for (ls[3]= -1; ls[3] < 2; ++ls[3]) {
          if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0) continue;
          RefPatch<Tet>::instance( ls);
        }
      }
    }
  }
}

template<> RefPatch<Hexa>::InitializerCL::InitializerCL () {
  init_count_++;
  if (init_count_++ > 0) return;
}


template<> bool RefPatch<Triangle2>::assign (const SignPattern<Triangle2>& cut) {

  for (size_= 0; size_ < num_elements( cut); ++size_) {
    face_[size_] = FaceIdx( cut(size_), cut(size_ + 1));
  }
  return empty();
}

template<> bool RefPatch<Quad2>::assign (const SignPattern<Quad2>& cut) {

  for (size_= 0; size_ < num_elements( cut); ++size_) {
    face_[size_] = FaceIdx( cut(size_), cut(size_ + 1));
  }
  return empty();
}

template<> bool RefPatch<Tet>::assign (const SignPattern<Tet>& cut){
  for (size_= 0; size_ < num_elements( cut); ++size_) {
    face_[size_] = FaceIdx( cut(size_), cut(size_ + 1), cut(size_ + 2));
  }
  return empty();
}

template<> Ubyte RefPatch<Hexa>::num_elements (const SignPattern<Hexa>& cut) const {
  switch (cut.num_cut_simplexes()) {
    case 3 : return 1;
    case 4 : return 2;
    case 5 : return 3;
    case 6 : return 4;
    default : return 0;
  }
}

template<> bool RefPatch<Hexa>::assign (const SignPattern<Hexa>& cut) {
  switch (static_cast<int>(cut.num_cut_simplexes())) {
    case 3 :
    face_[size_++] = FaceIdx(cut(0), cut(1), cut(2));
    break;
    case 4 :
    face_[size_++] = FaceIdx(cut(0), cut(1), cut(2));
    face_[size_++] = FaceIdx(cut(1), cut(2), cut(3));
    break;
    case 5 :
    face_[size_++] = FaceIdx(cut(0), cut(1), cut(2));
    face_[size_++] = FaceIdx(cut(0), cut(2), cut(3));
    face_[size_++] = FaceIdx(cut(0), cut(3), cut(4));
    break;
    case 6 :
    face_[size_++] = FaceIdx(cut(0), cut(1), cut(2));
    face_[size_++] = FaceIdx(cut(0), cut(2), cut(3));
    face_[size_++] = FaceIdx(cut(0), cut(3), cut(5));
    face_[size_++] = FaceIdx(cut(3), cut(4), cut(5));
    break;
  }
  return empty();
}


// RefPartition
template<> RefPartition<Triangle2>::InitializerCL::InitializerCL (){
  // std::cout << " initialize instances Triangle" << std::endl;

  if (init_count_++ > 0) return;

  byte ls[3];
  for (ls[0]= -1; ls[0] < 2; ++ls[0]) {
    for (ls[1]= -1; ls[1] < 2; ++ls[1]) {
      for (ls[2]= -1; ls[2] < 2; ++ls[2]) {
        if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0) continue;
        RefPartition<Triangle2>::instance( ls);
      }
    }
	}
}
template<> RefPartition<Quad2>::InitializerCL::InitializerCL (){
  // std::cout << " initialize instances Quad" << std::endl;

  if (init_count_++ > 0) return;

  byte ls[4];
  for (ls[0]= -1; ls[0] < 2; ++ls[0]){
    for (ls[1]= -1; ls[1] < 2; ++ls[1]){
      for (ls[2]= -1; ls[2] < 2; ++ls[2]){
        for (ls[3]= -1; ls[3] < 2; ++ls[3]) {

          if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0) continue;
          RefPartition<Quad2>::instance( ls);
        }
      }
    }
	}
}
template<> RefPartition<Tet>::InitializerCL::InitializerCL () {

  if (init_count_++ > 0) return;

  byte ls[4];
  for (ls[0]= -1; ls[0] < 2; ++ls[0]){
    for (ls[1]= -1; ls[1] < 2; ++ls[1]){
      for (ls[2]= -1; ls[2] < 2; ++ls[2]){
        for (ls[3]= -1; ls[3] < 2; ++ls[3]) {
          if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0) continue;
          RefPartition<Tet>::instance( ls);
        }
      }
    }
  }
}
template<> RefPartition<Hexa>::InitializerCL::InitializerCL () {
  init_count_++;
  if (init_count_++ > 0) return;
}


template<> bool RefPartition<Triangle2>::assign (const SignPattern<Triangle2>& cut) {
  end_= begin_= elements_ + start_array;

  if (cut.empty()) { // Most common case: no cut.
    Ubyte list_v[] = {0,1,2};
    // AddTriangle( 0, 1, 2, cut.sign( 0));
    AddElement( list_v, cut.sign( 0));

  }
  else if (cut.no_zero_vertex()) { // next common case: cuts without vertices on the zero level

    const Ubyte v = Triangle2::commonVertOfEdges[cut[0]][cut[1]];
    Ubyte list_v[] = {v, cut(0), cut(1)};
    AddElement( list_v, cut.sign( v));
    // AddTriangle( v, cut(0), cut(1), cut.sign( v));
    AddQuadrilateral( Triangle2::oppVertOfEdge( cut[0], v), cut(0),
              Triangle2::oppVertOfEdge( cut[1], v), cut(1),
	      -cut.sign( v));

  }
  // else if (cut.num_cut_simplexes() > cut.num_zero_vertexes()) { // next common case: there are cut edges, and also 1 or 2 vertices of the tetra with value 0 (the latter as we are in the else-part of cut.no_zero_vertex())
  //     if (cut.num_zero_vertexes() == 1) { // triangular cut through a vertex: a tetra and a remaining pyramid with quadrilateral base
  //         const Ubyte e= cut[1], f= cut[2];
  //         const Ubyte v= VertByEdge( e, f);
  //         AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));
  //         const Ubyte opp_v_in_e= v == VertOfEdge( e, 0) ? VertOfEdge( e, 1) : VertOfEdge( e, 0);
  //         const Ubyte opp_v_in_f= v == VertOfEdge( f, 0) ? VertOfEdge( f, 1) : VertOfEdge( f, 0);
  //         // the pyramid
  //         AddTetra( cut(0), cut(1), opp_v_in_f, opp_v_in_e, -cut.sign( v));
  //         AddTetra( cut(0), cut(1), opp_v_in_f, cut(2), -cut.sign( v));
  //     }
  //     else if (cut.num_zero_vertexes() == 2) { // triangular cut through 2 vertexes: two tetras
  //         const Ubyte e= OppEdge( EdgeByVert( cut[0], cut[1]));
  //         const Ubyte v0= VertOfEdge( e, 0), v1= VertOfEdge( e, 1);
  //         AddTetra( cut(0), cut(1), v0, cut(2), cut.sign( v0));
  //         AddTetra( cut(0), cut(1), v1, cut(2), cut.sign( v1));
  //     }
  // }
  // else // remaining cases: 1, 2 or 3 cuts, which are vertices of the tetra
  //     AddTetra( 0, 1, 2, 3, cut.sign( some_non_zero_vertex( cut)));
  return is_uncut();
}

template<> bool RefPartition<Quad2>::assign (const SignPattern<Quad2>& cut) {
  end_= begin_= elements_ + start_array;
  if (cut.empty()) { // Most common case: no cut.
    AddQuadrilateral( 0, 1, 2, 3, cut.sign( 0));
  }
  else if (cut.no_zero_vertex()) { // next common case: cuts without vertices on the zero level
    // case when interface cut two adjacent edges
    if(cut.sign_element.sum() != 0) {
      const Ubyte v = Quad2::commonVertOfEdges[cut[0]][cut[1]];
      Ubyte list_v[] = {v, cut(0), cut(1)};
      AddElement (list_v, cut.sign( v));
      // AddTriangle (v, cut(0), cut(1), cut.sign( v));
      const Ubyte v_op = (v+2)%4;
      AddPentagone( cut(0), Quad2::oppVertOfEdge( cut[0], v),
                    cut(1), Quad2::oppVertOfEdge( cut[1], v),
                    v_op, -cut.sign( v));
    }
    else {
      const Ubyte v0 = Quad2::nvedge[cut[0]][0];
      const Ubyte u0 = Quad2::nvedge[cut[0]][1];
      const Ubyte v1 = (u0+2)%4;
      const Ubyte u1 = (v0+2)%4;
      AddQuadrilateral(v0, cut(0), v1, cut(1),  cut.sign( v0));
      AddQuadrilateral(u1, cut(1), u0, cut(0), -cut.sign( v0));
    }
  }

  return is_uncut();
}

template<> bool RefPartition<Tet>::assign (const SignPattern<Tet>& cut) {
  end_= begin_= elements_ + start_array;

  if (cut.empty()) { // Most common case: no cut.
    Ubyte list_v[] = {0,1,2,3};
    AddElement( list_v, cut.sign( 0));
    // AddTetra( 0, 1, 2, 3, cut.sign( 0));
  }
  else if (cut.no_zero_vertex()) { // next common case: cuts without vertices on the zero level
    if (cut.num_cut_simplexes() == 3) { // triangular cut: a tetra and a remaining prism

      const Ubyte v = Tet::commonVertOfEdges[cut[0]][cut[1]];
      Ubyte list_v[] = {v, cut(0), cut(1), cut(2)};
      AddElement( list_v, cut.sign( v));
      // AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));

      AddPrism( Tet::oppVertOfEdge( cut[0], v), cut(0),
      Tet::oppVertOfEdge( cut[1], v), cut(1),
      Tet::oppVertOfEdge( cut[2], v), cut(2),
      -cut.sign( v));
    }
    else if (cut.num_cut_simplexes() == 4) { // quadrilateral cut: two prisms

      const Ubyte e = first_uncut_edge( cut);
      const Ubyte f = Tet::oppEdgeOfEdge[e];
      AddPrism( Tet::nvedge[e][0], Tet::nvedge[e][1],
        cut(0), cut(2),
        cut(1), cut(3),
        cut.sign(Tet::nvedge[e][0]));
        AddPrism( Tet::nvedge[f][0], Tet::nvedge[f][1],
          cut(0), cut(1),
          cut(2), cut(3),
          cut.sign(Tet::nvedge[f][0]));
        }
      }
      // else if (cut.num_cut_simplexes() > cut.num_zero_vertexes()) { // next common case: there are cut edges, and also 1 or 2 vertices of the tetra with value 0 (the latter as we are in the else-part of cut.no_zero_vertex())
      //     if (cut.num_zero_vertexes() == 1) { // triangular cut through a vertex: a tetra and a remaining pyramid with quadrilateral base
      //         const Ubyte e= cut[1], f= cut[2];
      //         const Ubyte v= VertByEdge( e, f);
      //         AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));
      //         const Ubyte opp_v_in_e= v == VertOfEdge( e, 0) ? VertOfEdge( e, 1) : VertOfEdge( e, 0);
      //         const Ubyte opp_v_in_f= v == VertOfEdge( f, 0) ? VertOfEdge( f, 1) : VertOfEdge( f, 0);
      //         // the pyramid
      //         AddTetra( cut(0), cut(1), opp_v_in_f, opp_v_in_e, -cut.sign( v));
      //         AddTetra( cut(0), cut(1), opp_v_in_f, cut(2), -cut.sign( v));
      //     }
      //     else if (cut.num_zero_vertexes() == 2) { // triangular cut through 2 vertexes: two tetras
      //         const Ubyte e= OppEdge( EdgeByVert( cut[0], cut[1]));
      //         const Ubyte v0= VertOfEdge( e, 0), v1= VertOfEdge( e, 1);
      //         AddTetra( cut(0), cut(1), v0, cut(2), cut.sign( v0));
      //         AddTetra( cut(0), cut(1), v1, cut(2), cut.sign( v1));
      //     }
      // }
      // else // remaining cases: 1, 2 or 3 cuts, which are vertices of the tetra
      //     AddTetra( 0, 1, 2, 3, cut.sign( some_non_zero_vertex( cut)));
      return is_uncut();
    }

template<> bool RefPartition<Hexa>::assign (const SignPattern<Hexa>& cut) {
  end_= begin_= elements_ + start_array;
  if (cut.empty()) { // Most common case: no cut.
    Ubyte list_v[] = {0,1,2,3,4,5,6,7};
    // addHexa(0,1,2,3,4,5,6,7,cut.sign( 0));
  }
  else {
    assert(cut.no_zero_vertex());
    int nb_pos = cut.sign_element.nb_node_positif();

    switch (static_cast<int>(cut.num_cut_simplexes())) {
      case 3 :
      std::cout << " cut not yet implemented " << std::endl;
      assert(0);
      break;
      case 4 :
      // add 1 prisme in part with 2 nodes
      // add 1 prisme && 1 hexa in part with 6 nodes
      if(nb_pos == 2 || nb_pos == 6){
        // Prisme => add 3 edges connecting the 2 triangular faces
          int e0 = 0;
          int e1 = (Hexa::commonVertOfEdges[cut[0]][cut[1]]==-1)? 1 : 3;
          int f0 = (Hexa::commonVertOfEdges[cut[0]][cut[1]]==-1)? 3 : 1;
          int f1 = 2;
          int g0 = Hexa::commonVertOfEdges[cut[e0]][cut[f0]];
          int g1 = Hexa::commonVertOfEdges[cut[e1]][cut[f1]];
          AddPrism( cut(e0), cut(e1), cut(f0), cut(f1), g0, g1, cut.sign( g0));

          // // add Hexa
          // // send a cartesian index way
          // int h0 = Hexa::oppVertOfEdge( cut[e0], g0);
          // int h1 = Hexa::oppVertOfEdge( cut[e1], g1);
          // int l0 = Hexa::oppVertOfEdge( cut[f0], g0);
          // int l1 = Hexa::oppVertOfEdge( cut[f1], g1);
          // AddHexa(cut(e0), cut(e1), cut(f0), cut(f1),h0,h1,l0,l1,-cut.sign( g0));
          //
          // int ed_op00 = Hexa::oppEdgeOfEdge[cut[e1]];
          // int ed_op01 = Hexa::oppEdgeOfEdge[cut[f1]];
          // int ed_op10 = Hexa::oppEdgeOfEdge[cut[e0]];
          // int ed_op11 = Hexa::oppEdgeOfEdge[cut[f0]];
          // int k0 = Hexa::commonVertOfEdges[ed_op00][ed_op01];
          // int k1 = Hexa::commonVertOfEdges[ed_op10][ed_op11];
          // AddPrism( h0,h1,l0,l1, k0, k1, -cut.sign( g0));
      }
      else {


      }



      break;
      case 5 :
      std::cout << " cut not yet implemented " << std::endl;
      assert(0);
      break;
      case 6 :
      std::cout << " cut not yet implemented " << std::endl;
      assert(0);
      break;
    };

  }
  return is_uncut();
}

template<> Ubyte RefPartition<Hexa>::first_uncut_edge (const SignPattern<Hexa>& cut) const {
  assert(0);
}

template<> void Partition<Triangle2>::get_list_node(vector<R2>& node, int s) const {
  node.resize(0);
  const SignPattern<Triangle2> cut(ls);
  if (cut.empty()) {
    for(int i=0;i<3;++i) node.push_back(T[i]);
  }
  else {
    const Ubyte v = Triangle2::commonVertOfEdges[cut[0]][cut[1]];
    if(cut.sign( v) == s) { // add triangle v, cut(0), cut(1)
      node.push_back(T[v]);
      node.push_back(get_vertex(cut(0)));
      node.push_back(get_vertex(cut(1)));
    }
    else {
      node.push_back(get_vertex(Triangle2::oppVertOfEdge( cut[0], v)));
      node.push_back(get_vertex(cut(0)));
      node.push_back(get_vertex(cut(1)));
      node.push_back(get_vertex(Triangle2::oppVertOfEdge( cut[1], v)));
    }
  }

}

template<> void Partition<Quad2>::get_list_node(vector<R2>& node, int s) const {
  node.resize(0);
  const SignPattern<Quad2> cut(ls);
  if (cut.empty()) {
    for(int i=0;i<4;++i) node.push_back(T[i]);
  }
  else {
    if(cut.sign_element.sum() != 0) {
      const Ubyte v = Quad2::commonVertOfEdges[cut[0]][cut[1]];
      if(cut.sign( v) == s) { // add triangle v, cut(0), cut(1)
        node.push_back(T[v]);
        node.push_back(get_vertex(cut(0)));
        node.push_back(get_vertex(cut(1)));
      }
      else{ // add pentagone
        const Ubyte v_op = (v+2)%4;
        node.push_back(get_vertex(Quad2::oppVertOfEdge( cut[0], v)));
        node.push_back(get_vertex(cut(0)));
        node.push_back(get_vertex(cut(1)));
        node.push_back(get_vertex(Quad2::oppVertOfEdge( cut[1], v)));
        node.push_back(get_vertex(v_op));
      }
    }else{

      const Ubyte v0 = Quad2::nvedge[cut[0]][0];
      const Ubyte u0 = Quad2::nvedge[cut[0]][1];
      const Ubyte v1 = (u0+2)%4;
      const Ubyte u1 = (v0+2)%4;
      if(cut.sign( v0) == s) {
        node.push_back(get_vertex(v0));
        node.push_back(get_vertex(cut(0)));
        node.push_back(get_vertex(cut(1)));
        node.push_back(get_vertex(v1));
      }
      else{
        node.push_back(get_vertex(u0));
        node.push_back(get_vertex(cut(0)));
        node.push_back(get_vertex(cut(1)));
        node.push_back(get_vertex(u1));
      }
    }
  }
}
