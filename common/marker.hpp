#ifndef _MARKER_HPP
#define _MARKER_HPP

#include <cstring>
#include "R3.hpp"
#include "Interface2dn.hpp"




class Marker;
class FunFEMVirtual ;

struct FaceMarker : public  SortArray<Uint, 2>, public Label {
  typedef SortArray<Uint, 2> FaceIdx;
  const Marker* marker = nullptr;
  int k;

  FaceMarker(const Marker& mark, int kk,
	      const Uint& a0,const Uint &a1, int l=0) : FaceIdx(a0,a1), Label(l),
							marker(&mark), k(kk){}
  FaceMarker(const Marker& mark, int kk, int l=0) : FaceMarker(mark,kk,0,0,l) {}
  FaceMarker(int kk, int l=0) : Label(l), k(kk){}

  // R2 operator()(int i) const;
  // double mesure() const;
  // R2 normal() const;
  // R2 mapToFace( R1 ip) const;

};


class Marker  : public Interface2 {
public:
  typedef typename Mesh2::Element Element;
  // typedef FaceMarker Face;
  typedef FaceIdx Face;

  typedef CutData2 CutData;

  const Mesh2& Th;
  bool periodic = false;


  std::vector<R2>  markers ;             // initial points
  std::vector<int> elementOfMarker;

  vector<int> nMarker_begin, nMarker_end;  // same but for node on edges

public :

  Marker(const Mesh2& Thh);
  Marker(const Mesh2& Thh, R2(*fparam)(double t), double x_begin, double x_end, int npoint);
  Marker(const Mesh2& Thh, std::string);

  void add(std::string path);
  void make_interface(bool leftB = false, bool rightB = false);
//   Marker(const Mesh2& Thh);
//   // void add(std::string, int ll = 0);
//   // void setDomainSign(const KN<double>& lls);
//
//   int NbNode() const {return vertices.size();}
//   const FaceMarker& getFace(int i) const {return faces[i];}
//   const FaceMarker& operator[](int i) const {return faces[i];}
//
  void move(const FunFEMVirtual&, double);
//
//   // void getNode(int, R2 loc_node[2]) const ;
//   // void getEdges(int, int ed[2]) const;
//   // void getSign(int, byte loc_sign[3]) const ;
    // CutData getCutData(int) const {return CutData2();};
//
//
//   bool isCut(int k) const    { return (element_seen.find(k) != element_seen.end());}
//   bool isNotCut(int k) const { return (element_seen.find(k) == element_seen.end());}
//   //
//   // bool isCut(int k) const { return (element_seen(k) != -1);}
//   // bool isNotCut(int k) const { return (element_seen(k) == -1);}
//   #ifdef USE_MPI
//   int first_element() const { return MPIcf::my_rank();}
//   int next_element() const {return MPIcf::size();}
//   #else
//   int first_element() const { return 0;}
//   int next_element() const {return 1;}
//   #endif
//   int last_element() const { return faces.size();}
//
  // int idxElementOfFace(int i) const {return 0;}//cut_element[i];}
  // R2 computeDx(const Face& face) const {
  //   return R2(0.,0.);//R2(edges_node[face[0]], edges_node[face[1]]);
  // }
  // // R2 normal(int i) const { return R2(0.,0.);}//return faces[i].normal();}
  // R2 mapToFace(const Face& f, const typename Element::RdHatBord x ) const{
  //   return R2(0.,0.);//return f.mapToFace( x);
  // };
  R2 make_normal (int i, int j);

 private:
//   // void add_interface(int ll);
  void find_vertices();
  void build_sign_array();
  void visit_element_sign(int k, vector<byte>& elementSeen);
  R2 get_intersect_edge(R2 firstNode, R2 nextNode, int previousK, int k, int& knext );

  void find_marker_limit(int k, int& i);
  int find_next_element(R2 A, R2 B, int previousK, int k);
  R2 find_intersection(R2 A, R2 B, int k, int kn, int& ie);
  int find_edge(int k, int kn);
  // bool check_intersect(R2 u, R2 v, R2 w, R& t);
//   // void addNeighbor(int , int);
//


};
































#endif
