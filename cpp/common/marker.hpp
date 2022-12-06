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
#ifndef COMMON_MARKER_HPP
#define COMMON_MARKER_HPP

#include <cstring>
#include "R2.hpp"
#include "Mesh2dn.hpp"
#include <vector>

//
// template<typename Mesh>
// class Marker  {
//
//
//   const Mesh& Th_;
//   bool periodic_ = false;
//   std::vector<double>  T_ ;             // initial points
//   std::vector<double>  X_ ;             // initial points
//   std::vector<double>  Y_ ;             // initial points
//
//   std::vector<int> element_of_marker_;
//
// public:
//   Marker(const Mesh& Thh, R2(*fparam)(double t), double t_begin, double
//   t_end, int npoint);
//
//   void add(double t, R2 x);
//   int nb_marker() const { return X_.size();}
//
// };
//

// #include "Interface2dn.hpp"

// class Marker;
class FunFEMVirtual;
//
// struct FaceMarker : public  SortArray<Uint, 2>, public Label {
//   typedef SortArray<Uint, 2> FaceIdx;
//   const Marker* marker = nullptr;
//   int k;
//
//   FaceMarker(const Marker& mark, int kk,
// 	      const Uint& a0,const Uint &a1, int l=0) : FaceIdx(a0,a1),
// Label(l), 							marker(&mark),
// k(kk){}
//   FaceMarker(const Marker& mark, int kk, int l=0) : FaceMarker(mark,kk,0,0,l)
//   {} FaceMarker(int kk, int l=0) : Label(l), k(kk){}
//
//   // R2 operator()(int i) const;
//   // double mesure() const;
//   // R2 normal() const;
//   // R2 mapToFace( R1 ip) const;
//
// };

class Marker {
 public:
   // typedef typename Mesh2::Element Element;

   const Mesh2 &Th_;
   bool periodic_ = false;

   std::vector<double> T_; // initial points
   std::vector<double> X_; // initial points
   std::vector<double> Y_; // initial points

   std::vector<int> element_of_marker_;

   vector<int> nMarker_begin_, nMarker_end_;

 public:
   Marker(const Mesh2 &Thh);
   Marker(const Mesh2 &Thh, R2 (*fparam)(double t), double x_begin,
          double x_end, int npoint);
   Marker(const Mesh2 &Thh, std::string);

   int size() const { return T_.size(); }
   R2 operator[](int i) const {
      assert(i < T_.size());
      return R2(X_[i], Y_[i]);
   }
   R2 get_marker(int i) const {
      assert(i < T_.size());
      return R2(X_[i], Y_[i]);
   }
   double get_parameter(int i) const {
      assert(i < T_.size());
      return T_[i];
   }
   void add(double t, R2 val);
   void set(int i, R2 val) {
      X_[i] = val.x;
      Y_[i] = val.y;
   }
   void move(const FunFEMVirtual &, double);

   // void add(std::string path);
   // void make_interface(bool leftB = false, bool rightB = false);
   //   Marker(const Mesh2& Thh);
   //   // void add(std::string, int ll = 0);
   //   // void setDomainSign(const KN<double>& lls);
   //
   //   int NbNode() const {return vertices.size();}
   //   const FaceMarker& getFace(int i) const {return faces[i];}
   //   const FaceMarker& operator[](int i) const {return faces[i];}
   //
   //
   //   // void getNode(int, R2 loc_node[2]) const ;
   //   // void getEdges(int, int ed[2]) const;
   //   // void getSign(int, byte loc_sign[3]) const ;
   // CutData getCutData(int) const {return CutData2();};
   //
   //
   //   bool isCut(int k) const    { return (element_seen.find(k) !=
   //   element_seen.end());} bool isNotCut(int k) const { return
   //   (element_seen.find(k) == element_seen.end());}
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
   // R2 make_normal (int i, int j);

 private:
   //   // void add_interface(int ll);
   // void find_vertices();
   // void build_sign_array();
   // void visit_element_sign(int k, vector<byte>& elementSeen);
   // R2 get_intersect_edge(R2 firstNode, R2 nextNode, int previousK, int k,
   // int& knext );
   //
   // void find_marker_limit(int k, int& i);
   // int find_next_element(R2 A, R2 B, int previousK, int k);
   // R2 find_intersection(R2 A, R2 B, int k, int kn, int& ie);
   // int find_edge(int k, int kn);
   // bool check_intersect(R2 u, R2 v, R2 w, R& t);
   //   // void addNeighbor(int , int);
   //
};

#endif
