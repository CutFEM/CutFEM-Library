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
#include "finiteElement.hpp"

template <>
std::vector<std::vector<const GTypeOfFE<Mesh2> *>> BaseFE_Array<Mesh2, ContinuityType::continuous>::FE_ = {
    {&DataFE<Mesh2>::P0, &DataFE<Mesh2>::P0},
    {&DataFE<Mesh2>::P1, &DataFE<Mesh2>::P1},
    {&DataFE<Mesh2>::P2, &DataFE<Mesh2>::P2},
    {&DataFE<Mesh2>::P3, &DataFE<Mesh2>::P3},
    {&DataFE<Mesh2>::P4, &DataFE<Mesh2>::P4}};

template <>
BaseFE_Array<Mesh2, ContinuityType::continuous>::BaseFE_Array(int k) : GTypeOfFESum<Mesh2>(this->FE_.at(k)) {}
Lagrange2::Lagrange2(int k) : BaseFE_Array<Mesh2, ContinuityType::continuous>(k){};

template <>
std::vector<std::vector<const GTypeOfFE<MeshQuad2> *>> BaseFE_Array<MeshQuad2, ContinuityType::continuous>::FE_ = {
    {&DataFE<MeshQuad2>::P0, &DataFE<MeshQuad2>::P0}, {&DataFE<MeshQuad2>::P1, &DataFE<MeshQuad2>::P1}};

template <>
BaseFE_Array<MeshQuad2, ContinuityType::continuous>::BaseFE_Array(int k) : GTypeOfFESum<MeshQuad2>(this->FE_.at(k)) {}
LagrangeQuad2::LagrangeQuad2(int k) : BaseFE_Array<MeshQuad2, ContinuityType::continuous>(k){};

template <>
std::vector<std::vector<const GTypeOfFE<Mesh3> *>> BaseFE_Array<Mesh3, ContinuityType::continuous>::FE_ = {
    {&DataFE<Mesh3>::P0, &DataFE<Mesh3>::P0, &DataFE<Mesh3>::P0},
    {&DataFE<Mesh3>::P1, &DataFE<Mesh3>::P1, &DataFE<Mesh3>::P1},
    {&DataFE<Mesh3>::P2, &DataFE<Mesh3>::P2, &DataFE<Mesh3>::P2}};

template <>
BaseFE_Array<Mesh3, ContinuityType::continuous>::BaseFE_Array(int k) : GTypeOfFESum<Mesh3>(this->FE_.at(k)) {}
Lagrange3::Lagrange3(int k) : BaseFE_Array<Mesh3, ContinuityType::continuous>(k){};

template <>
std::vector<std::vector<const GTypeOfFE<Mesh2> *>> BaseFE_Array<Mesh2, ContinuityType::discontinuous>::FE_ = {
    {&DataFE<Mesh2>::P0, &DataFE<Mesh2>::P0},
    {&DataFE<Mesh2>::P1dc, &DataFE<Mesh2>::P1dc},
    {&DataFE<Mesh2>::P2dc, &DataFE<Mesh2>::P2dc},
    {&DataFE<Mesh2>::P3dc, &DataFE<Mesh2>::P3dc}};

template <>
BaseFE_Array<Mesh2, ContinuityType::discontinuous>::BaseFE_Array(int k) : GTypeOfFESum<Mesh2>(this->FE_.at(k)) {}
LagrangeDC2::LagrangeDC2(int k) : BaseFE_Array<Mesh2, ContinuityType::discontinuous>(k){};

// const GTypeOfFE<Mesh2> *TaylorHood2::FE_[3] = {&DataFE<Mesh2>::P2, &DataFE<Mesh2>::P2,
//                                                &DataFE<Mesh2>::P1}; //&DataFE<Mesh2>::P2;
//                                                                     // const GTypeOfFE<Mesh2> *Lagrange2::FE_[5][2] =
//                                                                     {

// const GTypeOfFE<Mesh3> *TaylorHood3::FE_[4] = {&DataFE<Mesh3>::P2, &DataFE<Mesh3>::P2, &DataFE<Mesh3>::P2,
//                                                &DataFE<Mesh3>::P1}; //&DataFE<Mesh2>::P2;

// template<>
// const GTypeOfFE<MeshHexa>* Lagrange3<MeshHexa>::FE_[3][3] =
// {{&DataFE<MeshHexa>::P0,&DataFE<MeshHexa>::P0,&DataFE<MeshHexa>::P0},
// {&DataFE<MeshHexa>::P1,&DataFE<MeshHexa>::P1,&DataFE<MeshHexa>::P1}//,
// // {&DataFE<MeshHexa>::P2,&DataFE<MeshHexa>::P2,&DataFE<MeshHexa>::P2}
// // ,{&DataFE<Mesh3>::P3,&DataFE<Mesh3>::P3,&DataFE<Mesh3>::P3}
//                                                                                             };
