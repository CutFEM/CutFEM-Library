#include "Interface3dn.hpp"
#include <iostream>
#include <iterator>     // std::ostream_iterator


std::ostream& operator<< (std::ostream& out, const SignPatternTrait3& c)
{
  out << static_cast<int>( c.num_root_vert_) << ' ' << static_cast<int>( c.num_root_) << '\n';
  for (int i= 0; i < 4; ++i)
    out << static_cast<int>( c.sign_[i]) << ' ';
  out << '\n';
  for (int i= 0; i < 4; ++i)
    out << static_cast<int>( c.cut_simplex_[i]) << ' ';
  out << '\n';
  for (int i= 0; i < 4; ++i)
    out << static_cast<int>( c.cut_simplex_rep_[i]) << ' ';
  return out << '\n';
}

RefPatch3 RefPatch3::instance_array_[81];

int RefPatch3::InitializerCL::init_count_= 0;

RefPatch3::InitializerCL::InitializerCL ()
{
  if (init_count_++ > 0)
    return;

  byte ls[4];
  for (ls[0]= -1; ls[0] < 2; ++ls[0])
    for (ls[1]= -1; ls[1] < 2; ++ls[1])
      for (ls[2]= -1; ls[2] < 2; ++ls[2]) 
	for (ls[3]= -1; ls[3] < 2; ++ls[3]) {
	  if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0)
	    continue;
	  RefPatch3::instance( ls);
	}
}

bool RefPatch3::assign (const SignPatternTrait3& cut)
{
  for (size_= 0; size_ < num_elements( cut); ++size_) {
    face_[size_] = MakeTriangle( cut(size_), cut(size_ + 1), cut(size_ + 2));
  }

  
  return empty();
}



void write_paraview_vtu (std::ostream& file_, const Interface3& t)
{
  file_ << "<?xml version=\"1.0\"?>"  << '\n'
	<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"   << '\n'
	<< "<UnstructuredGrid>"   << '\n';
  
  file_<< "<Piece NumberOfPoints=\""<< t.vertices_.size() <<"\" NumberOfCells=\""<< t.faces_.size() << "\">";
  file_<< "\n\t<Points>"
       << "\n\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"" << "ascii\">\n\t\t";
  for(Interface3::const_vertex_iterator it= t.vertex_begin(), end= t.vertex_end();
      it != end; ++it) {
    file_ << it[0][0] << ' ' << it[0][1] << ' ' << it[0][2] << ' ';
  }
  file_<< "\n\t\t</DataArray> \n"
       << "\t</Points>\n";

  file_   << "\t<Cells>\n"
	  << "\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\""
	  <<"ascii\">\n\t\t";
  std::copy( t.faces_.begin(), t.faces_.end(), std::ostream_iterator<Interface3::FaceIdx>( file_));
  file_ << "\n\t\t</DataArray>\n";
  file_ << "\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t";
  for(Uint i= 1; i <= t.faces_.size(); ++i) {
    file_ << i*3 << " ";
  }
  file_ << "\n\t\t</DataArray>";
  file_ << "\n\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n\t\t";
  const int Type= 5; // Triangles
  for(Uint i= 1; i <= t.faces_.size(); ++i) {
    file_ << Type<<" ";
  }
  file_<<"\n\t\t</DataArray>"
       <<"\n\t</Cells>";

  file_ <<"\n</Piece>"
	<<"\n</UnstructuredGrid>"
	<<"\n</VTKFile>";
}




RefPartition3 RefPartition3::instance_array_[81];

int RefPartition3::InitializerCL::init_count_= 0;

RefPartition3::InitializerCL::InitializerCL ()
{
  if (init_count_++ > 0)
    return;

  byte ls[4];
  for (ls[0]= -1; ls[0] < 2; ++ls[0])
    for (ls[1]= -1; ls[1] < 2; ++ls[1])
      for (ls[2]= -1; ls[2] < 2; ++ls[2]) 
	for (ls[3]= -1; ls[3] < 2; ++ls[3]) {
	  if ( ls[0] == 0 && ls[1] == 0 && ls[2] == 0 && ls[3] == 0)
	    continue;
	  RefPatch3::instance( ls);
	}
}


bool
RefPartition3::assign (const SignPatternTrait3& cut)
{
  end_= begin_= tetras_ + 3;
  
  if (cut.empty()) { // Most common case: no cut.
    AddTetra( 0, 1, 2, 3, cut.sign( 0));
  }
  else if (cut.no_zero_vertex()) { // next common case: cuts without vertices on the zero level
    if (cut.num_cut_simplexes() == 3) { // triangular cut: a tetra and a remaining prism
      
      const Ubyte v = Element::commonVertOfEdges[cut[0]][cut[1]];
      AddTetra( v, cut(0), cut(1), cut(2), cut.sign( v));
      AddPrism( Element::oppVertOfEdge( cut[0], v), cut(0),
		Element::oppVertOfEdge( cut[1], v), cut(1),
		Element::oppVertOfEdge( cut[2], v), cut(2),
		-cut.sign( v));
    }
    else if (cut.num_cut_simplexes() == 4) { // quadrilateral cut: two prisms

      const Ubyte e = first_uncut_edge( cut);
      const Ubyte f = Element::oppEdgeOfEdge[e];
      AddPrism( Element::nvedge[e][0], Element::nvedge[e][1],
		cut(0), cut(2),
		cut(1), cut(3),
		cut.sign(Element::nvedge[e][0]));
      AddPrism( Element::nvedge[f][0], Element::nvedge[f][1],
		cut(0), cut(1),
		cut(2), cut(3),
		cut.sign(Element::nvedge[f][0]));
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



// Element RefPartition3::
