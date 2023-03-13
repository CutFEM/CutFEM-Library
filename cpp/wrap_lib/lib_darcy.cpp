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
#include "lib_darcy.hpp"

extern "C" {

Darcy2 *Darcy2_new(Darcy2 *darcy) { return new Darcy2(); }

void Darcy2_add_natural_BC(Darcy2 *darcy, double (*f)(double *, int, int)) { darcy->add_natural_BC(f); }
}
