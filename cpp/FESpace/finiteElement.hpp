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

class Lagrange2 : public GTypeOfFESum<Mesh2> {
   typedef KN<const GTypeOfFE<Mesh2> *> FEarray;
   static const GTypeOfFE<Mesh2> *FE_[5][2];

 public:
   Lagrange2(int k = 1) : GTypeOfFESum<Mesh2>(FEarray(2, FE_[k])) {}
};

class LagrangeQuad2 : public GTypeOfFESum<MeshQuad2> {
   typedef KN<const GTypeOfFE<MeshQuad2> *> FEarray;
   static const GTypeOfFE<MeshQuad2> *FE_[1][2];

 public:
   LagrangeQuad2(int k = 1) : GTypeOfFESum<MeshQuad2>(FEarray(2, FE_[k])) {}
};

class LagrangeDC2 : public GTypeOfFESum<Mesh2> {
   typedef KN<const GTypeOfFE<Mesh2> *> FEarray;
   static const GTypeOfFE<Mesh2> *FE_[4][2];

 public:
   LagrangeDC2(int k = 1) : GTypeOfFESum<Mesh2>(FEarray(2, FE_[k])) {}
};

class TaylorHood2 : public GTypeOfFESum<Mesh2> {
   typedef KN<const GTypeOfFE<Mesh2> *> FEarray;
   static const GTypeOfFE<Mesh2> *FE_[3];

 public:
   TaylorHood2() : GTypeOfFESum<Mesh2>(FEarray(3, FE_)) {}
};

class TaylorHood3 : public GTypeOfFESum<Mesh3> {
   typedef KN<const GTypeOfFE<Mesh3> *> FEarray;
   static const GTypeOfFE<Mesh3> *FE_[4];

 public:
   TaylorHood3() : GTypeOfFESum<Mesh3>(FEarray(4, FE_)) {}
};

class Lagrange3 : public GTypeOfFESum<Mesh3> {
   typedef KN<const GTypeOfFE<Mesh3> *> FEarray;
   static const GTypeOfFE<Mesh3> *FE_[3][3];

 public:
   Lagrange3(int k = 1) : GTypeOfFESum<Mesh3>(FEarray(3, FE_[k])) {}
};

#endif
