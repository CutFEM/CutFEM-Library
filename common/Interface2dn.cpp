#include "Interface2dn.hpp"


std::ostream& operator<< (std::ostream& out, const SignPatternTrait2& c)
{
  out << static_cast<int>( c.num_root_vert_) << ' ' << static_cast<int>( c.num_root_) << '\n';
  for (int i= 0; i < 3; ++i)
    out << static_cast<int>( c.sign_[i]) << ' ';
  out << '\n';
  for (int i= 0; i < 2; ++i)
    out << static_cast<int>( c.cut_simplex_[i]) << ' ';
  out << '\n';
  for (int i= 0; i < 2; ++i)
    out << static_cast<int>( c.cut_simplex_rep_[i]) << ' ';
  return out << '\n';
}

RefPatch2 RefPatch2::instance_array_[27];

int RefPatch2::InitializerCL::init_count_= 0;

RefPatch2::InitializerCL::InitializerCL ()
{
  if (init_count_++ > 0)
    return;

  byte ls[3];
  for (ls[0]= -1; ls[0] < 2; ++ls[0])
    for (ls[1]= -1; ls[1] < 2; ++ls[1])
      for (ls[2]= -1; ls[2] < 2; ++ls[2]) {
	  if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0)
	    continue;
	  RefPatch2::instance( ls);
	}
}

bool RefPatch2::assign (const SignPatternTrait2& cut)
{
  for (size_= 0; size_ < num_elements( cut); ++size_)
    face_[size_]= MakeTriangle( cut(size_), cut(size_ + 1));
  return empty();
}





RefPartition2 RefPartition2::instance_array_[27];

int RefPartition2::InitializerCL::init_count_= 0;

RefPartition2::InitializerCL::InitializerCL ()
{
  if (init_count_++ > 0)
    return;

  byte ls[3];
  for (ls[0]= -1; ls[0] < 2; ++ls[0])
    for (ls[1]= -1; ls[1] < 2; ++ls[1])
      for (ls[2]= -1; ls[2] < 2; ++ls[2]) {
	if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0)
	  continue;
	RefPatch2::instance( ls);
	}
}


bool
RefPartition2::assign (const SignPatternTrait2& cut)
{
  end_= begin_= triangles_ + 2;

  if (cut.empty()) { // Most common case: no cut.
    AddTriangle( 0, 1, 2, cut.sign( 0));
  }
  else if (cut.no_zero_vertex()) { // next common case: cuts without vertices on the zero level

    const Ubyte v = Element::commonVertOfEdges[cut[0]][cut[1]];
    AddTriangle( v, cut(0), cut(1), cut.sign( v));
    AddPrism( Element::oppVertOfEdge( cut[0], v), cut(0),
	      Element::oppVertOfEdge( cut[1], v), cut(1),
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
