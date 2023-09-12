// NEW AND NON-WORKING FOR MAXWELL3D

// /*
// This file is part of CutFEM-Library.

// CutFEM-Library is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

// CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
// */
// /*

//  This file is part of Freefem++

//  Freefem++ is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation; either version 2.1 of the License, or
//  (at your option) any later version.

//  Freefem++  is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.

//  You should have received a copy of the GNU Lesser General Public License
//  along with Freefem++; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//  */

// #ifndef _OPERATIONTYPEFE_HPP
// #define _OPERATIONTYPEFE_HPP

// enum operatortype {
//     op_id    = 0,
//     op_dx    = 1,
//     op_dy    = 2,
//     op_dz    = 3,
//     op_dxx   = 4,
//     op_dyy   = 5,
//     op_dzz   = 6,
//     op_dyx   = 7,
//     op_dxy   = 7,
//     op_dzx   = 8,
//     op_dxz   = 8,
//     op_dzy   = 9,
//     op_dyz   = 9,
//     op_all   = 10,
//     op_Dall  = 10,
//     // need to change those value, do -1, start from 10
//     // op_dxxx  = 11,
//     // op_dxxy  = 14,
//     // op_dxxz  = 15,
//     // op_dxyx  = 14,
//     // op_dxyy  = 16,
//     // op_dxyz  = 17,
//     // op_dxzx  = 15,
//     // op_dxzy  = 17,
//     // op_dxzz  = 18,
//     // op_dyxx  = 14,
//     // op_dyxy  = 16,
//     // op_dyxz  = 17,
//     // op_dyyx  = 16,
//     // op_dyyy  = 12,
//     // op_dyyz  = 19,
//     // op_dyzx  = 17,
//     // op_dyzy  = 19,
//     // op_dyzz  = 20,
//     // op_dzxx  = 15,
//     // op_dzxy  = 17,
//     // op_dzxz  = 18,
//     // op_dzyx  = 17,
//     // op_dzyy  = 19,
//     // op_dzyz  = 20,
//     // op_dzzx  = 18,
//     // op_dzzy  = 20,
//     // op_dzzz  = 13,
//     // op_DDall = 21,
//     // op_All   = 21

//     op_dxxx  = 10,
//     op_dxxy  = 13,
//     op_dxxz  = 14,
//     op_dxyx  = 13,
//     op_dxyy  = 15,
//     op_dxyz  = 16,
//     op_dxzx  = 14,
//     op_dxzy  = 16,
//     op_dxzz  = 17,
//     op_dyxx  = 13,
//     op_dyxy  = 15,
//     op_dyxz  = 16,
//     op_dyyx  = 15,
//     op_dyyy  = 11,
//     op_dyyz  = 18,
//     op_dyzx  = 16,
//     op_dyzy  = 18,
//     op_dyzz  = 19,
//     op_dzxx  = 14,
//     op_dzxy  = 16,
//     op_dzxz  = 17,
//     op_dzyx  = 16,
//     op_dzyy  = 18,
//     op_dzyz  = 19,
//     op_dzzx  = 17,
//     op_dzzy  = 19,
//     op_dzzz  = 12,
//     op_DDall = 21,
//     op_All   = 21
// };

// typedef unsigned int What_d;

// const unsigned int Fop_id = 1 << op_id;

// const unsigned int Fop_dx = 1 << op_dx;
// const unsigned int Fop_dy = 1 << op_dy;
// const unsigned int Fop_dz = 1 << op_dz;

// const unsigned int Fop_dxx = 1 << op_dxx;
// const unsigned int Fop_dxy = 1 << op_dxy;
// const unsigned int Fop_dxz = 1 << op_dxz;

// const unsigned int Fop_dyx = 1 << op_dyx;
// const unsigned int Fop_dyy = 1 << op_dyy;
// const unsigned int Fop_dyz = 1 << op_dyz;

// const unsigned int Fop_dzx = 1 << op_dzx;
// const unsigned int Fop_dzy = 1 << op_dzy;
// const unsigned int Fop_dzz = 1 << op_dzz;

// // const unsigned int Fop_dxxx = 1 << op_dxxx;
// // const unsigned int Fop_dxxy = 1 << op_dxxy;
// // const unsigned int Fop_dxxz = 1 << op_dxxz;

// // const unsigned int Fop_dxyy = 1 << op_dxyy;
// // const unsigned int Fop_dxyz = 1 << op_dxyz;
// // const unsigned int Fop_dxzz = 1 << op_dxzz;

// // const unsigned int Fop_dyyy = 1 << op_dyyy;
// // const unsigned int Fop_dyyz = 1 << op_dyyz;
// // const unsigned int Fop_dyzz = 1 << op_dyzz;

// // const unsigned int Fop_dzzz = 1 << op_dzzz;

// const unsigned int Fop_dxxx = 1 << op_dxxx;
// const unsigned int Fop_dxxy = 1 << op_dxxy;
// const unsigned int Fop_dxxz = 1 << op_dxxz;
// const unsigned int Fop_dxyx = 1 << op_dxyx;
// const unsigned int Fop_dxyy = 1 << op_dxyy;
// const unsigned int Fop_dxyz = 1 << op_dxyz;
// const unsigned int Fop_dxzx = 1 << op_dxzx;
// const unsigned int Fop_dxzy = 1 << op_dxzy;
// const unsigned int Fop_dxzz = 1 << op_dxzz;
// const unsigned int Fop_dyxx = 1 << op_dyxx;
// const unsigned int Fop_dyxy = 1 << op_dyxy;
// const unsigned int Fop_dyxz = 1 << op_dyxz;
// const unsigned int Fop_dyyx = 1 << op_dyyx;
// const unsigned int Fop_dyyy = 1 << op_dyyy;
// const unsigned int Fop_dyyz = 1 << op_dyyz;
// const unsigned int Fop_dyzx = 1 << op_dyzx;
// const unsigned int Fop_dyzy = 1 << op_dyzy;
// const unsigned int Fop_dyzz = 1 << op_dyzz;
// const unsigned int Fop_dzxx = 1 << op_dzxx;
// const unsigned int Fop_dzxy = 1 << op_dzxy;
// const unsigned int Fop_dzxz = 1 << op_dzxz;
// const unsigned int Fop_dzyx = 1 << op_dzyx;
// const unsigned int Fop_dzyy = 1 << op_dzyy;
// const unsigned int Fop_dzyz = 1 << op_dzyz;
// const unsigned int Fop_dzzx = 1 << op_dzzx;
// const unsigned int Fop_dzzy = 1 << op_dzzy;
// const unsigned int Fop_dzzz = 1 << op_dzzz;

// const unsigned int Fop_D0 = Fop_id;
// const unsigned int Fop_D1 = Fop_dx | Fop_dy | Fop_dz;
// const unsigned int Fop_D2 = Fop_dxx | Fop_dyy | Fop_dzz | Fop_dxy | Fop_dxz | Fop_dyz;
// const unsigned int Fop_D3 =
//     Fop_dxxx | Fop_dxxy | Fop_dxxz | Fop_dxyy | Fop_dxyz | Fop_dxzz | Fop_dyyy | Fop_dyyz | Fop_dyzz | Fop_dzzz;

// const unsigned int Fop_Dall = Fop_D0 | Fop_D1 | Fop_D2 | Fop_D3;
// const unsigned int Fop_DAll = 20;

// inline What_d Fwhatd(const operatortype op) { return 1 << op; }
// inline What_d Fwhatd(const int op) {
//     if (op == 1)
//         return Fop_D0;
//     else if (op < 5)
//         return Fop_D0 | Fop_D1;
//     else if (op < 10)
//         return Fop_D0 | Fop_D1 | Fop_D2;
//     else
//         return Fop_Dall;
// }

// inline int getLastop(int i, int j) {
//     int lastop = std::max(i, j);
//     if (lastop == 0)
//         return lastop + 1;
//     if (lastop <= op_dz)
//         return op_dz+1;
//     if (lastop <= op_dyz)
//         return op_dyz;
//     if (lastop <= op_dzzy)
//         return op_DDall;
//     return 0;
// }

// const int last_operatortype                     = 20;
// const bool operatortypeValue[last_operatortype] = {
//     true,  false, false, false, false, false, false, false, false, false,
//     false, false, false, false, false, false, false, false, false, false,
// };

// inline void initwhatd(bool *whatd, int k) {
//     for (int i = 0; i < last_operatortype; i++)
//         whatd[i] = false;
//     whatd[k] = true;
// }

// #endif



// OLD AND WORKING
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
/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef _OPERATIONTYPEFE_HPP
#define _OPERATIONTYPEFE_HPP

enum operatortype {
    op_id    = 0,
    op_dx    = 1,
    op_dy    = 2,
    op_dz    = 3,
    op_dxx   = 4,
    op_dyy   = 5,
    op_dzz   = 6,
    op_dyx   = 7,
    op_dxy   = 7,
    op_dzx   = 8,
    op_dxz   = 8,
    op_dzy   = 9,
    op_dyz   = 9,
    op_all   = 10,
    op_Dall  = 10,
    // need to change those value, do -1, start from 10
    op_dxxx  = 11,
    op_dxxy  = 14,
    op_dxxz  = 15,
    op_dxyx  = 14,
    op_dxyy  = 16,
    op_dxyz  = 17,
    op_dxzx  = 15,
    op_dxzy  = 17,
    op_dxzz  = 18,
    op_dyxx  = 14,
    op_dyxy  = 16,
    op_dyxz  = 17,
    op_dyyx  = 16,
    op_dyyy  = 12,
    op_dyyz  = 19,
    op_dyzx  = 17,
    op_dyzy  = 19,
    op_dyzz  = 20,
    op_dzxx  = 15,
    op_dzxy  = 17,
    op_dzxz  = 18,
    op_dzyx  = 17,
    op_dzyy  = 19,
    op_dzyz  = 20,
    op_dzzx  = 18,
    op_dzzy  = 20,
    op_dzzz  = 13,
    op_DDall = 21,
    op_All   = 21
};

typedef unsigned int What_d;

const unsigned int Fop_id = 1 << op_id;

const unsigned int Fop_dx = 1 << op_dx;
const unsigned int Fop_dy = 1 << op_dy;
const unsigned int Fop_dz = 1 << op_dz;

const unsigned int Fop_dxx = 1 << op_dxx;
const unsigned int Fop_dxy = 1 << op_dxy;
const unsigned int Fop_dxz = 1 << op_dxz;

const unsigned int Fop_dyx = 1 << op_dyx;
const unsigned int Fop_dyy = 1 << op_dyy;
const unsigned int Fop_dyz = 1 << op_dyz;

const unsigned int Fop_dzx = 1 << op_dzx;
const unsigned int Fop_dzy = 1 << op_dzy;
const unsigned int Fop_dzz = 1 << op_dzz;

const unsigned int Fop_D0 = Fop_id;
const unsigned int Fop_D1 = Fop_dx | Fop_dy | Fop_dz;
const unsigned int Fop_D2 =
    Fop_dxx | Fop_dyy | Fop_dzz | Fop_dxy | Fop_dxz | Fop_dyz;
const unsigned int Fop_Dall = Fop_D0 | Fop_D1 | Fop_D2;
const unsigned int Fop_DAll = 20;

inline What_d Fwhatd(const operatortype op) { return 1 << op; }
inline What_d Fwhatd(const int op) {
    if (op == 1)
        return Fop_D0;
    else if (op < 5)
        return Fop_D0 | Fop_D1;
    else
        return Fop_Dall;
}

inline int getLastop(int i, int j) {
    int lastop = std::max(i, j);
    if (lastop == 0)
        return lastop + 1;
    if (lastop <= op_dz)
        return op_dz + 1;
    if (lastop <= op_dyz)
        return op_Dall;
    return 0;
}

const int last_operatortype                     = 20;
const bool operatortypeValue[last_operatortype] = {
    true,  false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false,
};

inline void initwhatd(bool *whatd, int k) {
    for (int i = 0; i < last_operatortype; i++)
        whatd[i] = false;
    whatd[k] = true;
}

#endif
