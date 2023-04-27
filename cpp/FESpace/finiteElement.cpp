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

const GTypeOfFE<Mesh2> *TaylorHood2::FE_[3] = {&DataFE<Mesh2>::P2, &DataFE<Mesh2>::P2,
                                               &DataFE<Mesh2>::P1}; //&DataFE<Mesh2>::P2;

const GTypeOfFE<Mesh2> *Lagrange2::FE_[5][2] = {{&DataFE<Mesh2>::P0, &DataFE<Mesh2>::P0},
                                                {&DataFE<Mesh2>::P1, &DataFE<Mesh2>::P1},
                                                {&DataFE<Mesh2>::P2, &DataFE<Mesh2>::P2},
                                                {&DataFE<Mesh2>::P3, &DataFE<Mesh2>::P3},
                                                {&DataFE<Mesh2>::P4, &DataFE<Mesh2>::P4}};

const GTypeOfFE<MeshQuad2> *LagrangeQuad2::FE_[1][2] = {{&DataFE<MeshQuad2>::P1, &DataFE<MeshQuad2>::P1}};

const GTypeOfFE<Mesh2> *LagrangeDC2::FE_[4][2] = {{&DataFE<Mesh2>::P0, &DataFE<Mesh2>::P0},
                                                  {&DataFE<Mesh2>::P1dc, &DataFE<Mesh2>::P1dc},
                                                  {&DataFE<Mesh2>::P2dc, &DataFE<Mesh2>::P2dc},
                                                  {&DataFE<Mesh2>::P3dc, &DataFE<Mesh2>::P3dc}};

const GTypeOfFE<Mesh3> *TaylorHood3::FE_[4] = {&DataFE<Mesh3>::P2, &DataFE<Mesh3>::P2, &DataFE<Mesh3>::P2,
                                               &DataFE<Mesh3>::P1}; //&DataFE<Mesh2>::P2;

const GTypeOfFE<Mesh3> *Lagrange3::FE_[3][3] = {
    {&DataFE<Mesh3>::P0, &DataFE<Mesh3>::P0, &DataFE<Mesh3>::P0},
    {&DataFE<Mesh3>::P1, &DataFE<Mesh3>::P1, &DataFE<Mesh3>::P1},
    {&DataFE<Mesh3>::P2, &DataFE<Mesh3>::P2, &DataFE<Mesh3>::P2}
    // ,{&DataFE<Mesh3>::P3,&DataFE<Mesh3>::P3,&DataFE<Mesh3>::P3}
};
// template<>
// const GTypeOfFE<MeshHexa>* Lagrange3<MeshHexa>::FE_[3][3] =
// {{&DataFE<MeshHexa>::P0,&DataFE<MeshHexa>::P0,&DataFE<MeshHexa>::P0},
// {&DataFE<MeshHexa>::P1,&DataFE<MeshHexa>::P1,&DataFE<MeshHexa>::P1}//,
// // {&DataFE<MeshHexa>::P2,&DataFE<MeshHexa>::P2,&DataFE<MeshHexa>::P2}
// // ,{&DataFE<Mesh3>::P3,&DataFE<Mesh3>::P3,&DataFE<Mesh3>::P3}
//                                                                                             };
