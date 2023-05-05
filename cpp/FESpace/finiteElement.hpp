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
#ifndef FINITE_ELEMENT_HPP
#define FINITE_ELEMENT_HPP

#include "GTypeOfFE_Sum.hpp"
#include "../concept/function.hpp"

enum class ContinuityType { continuous, discontinuous };

template <typeMesh mesh_t, ContinuityType C> class BaseFE_Array : public GTypeOfFESum<mesh_t> {
  public:
    static std::vector<std::vector<const GTypeOfFE<mesh_t> *>> FE_;
    BaseFE_Array(int k);
};

class Lagrange2 : public BaseFE_Array<Mesh2, ContinuityType::continuous> {
  public:
    Lagrange2(int k);
};
class LagrangeQuad2 : public BaseFE_Array<MeshQuad2, ContinuityType::continuous> {
  public:
    LagrangeQuad2(int k);
};

class Lagrange3 : public BaseFE_Array<Mesh3, ContinuityType::continuous> {
  public:
    Lagrange3(int k);
};

class LagrangeDC2 : public BaseFE_Array<Mesh2, ContinuityType::discontinuous> {
  public:
    LagrangeDC2(int k);
};

// class TaylorHood2 : public GTypeOfFESum<Mesh2> {
//     typedef KN<const GTypeOfFE<Mesh2> *> FEarray;
//     static const GTypeOfFE<Mesh2> *FE_[3];

//   public:
//     TaylorHood2() : GTypeOfFESum<Mesh2>(FEarray(3, FE_)) {}
// };

// class TaylorHood3 : public GTypeOfFESum<Mesh3> {
//     typedef KN<const GTypeOfFE<Mesh3> *> FEarray;
//     static const GTypeOfFE<Mesh3> *FE_[4];

//   public:
//     TaylorHood3() : GTypeOfFESum<Mesh3>(FEarray(4, FE_)) {}
// };

#endif
