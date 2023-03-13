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
#include "lib_stokes.hpp"

extern "C" {

FictitiousStokesRT_2 *FictitiousStokes_new(FictitiousStokesRT_2 *stokes) { return new FictitiousStokesRT_2(); }
FictitiousStokesVorticity_2 *FictitiousStokesVorticity_new(FictitiousStokesVorticity_2 *stokes) {
    return new FictitiousStokesVorticity_2();
}
}
