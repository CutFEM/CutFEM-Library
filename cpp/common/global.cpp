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
#include "global.hpp"

int globalVariable::verbose        = 1;
double globalVariable::UnSetMesure = -1e+200;
double globalVariable::Epsilon    = 10 * std::numeric_limits<double>::epsilon();
double globalVariable::DoubleEpsC = 1e-9;
double globalVariable::Pi         = 3.14159265359;
double globalVariable::pi         = 3.14159265359;
