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
#include "parametrization.hpp"

LinearParametrization::LinearParametrization(const std::vector<double> &x,
                                             const std::vector<double> &y) {
   this->init(x, y);
}

void LinearParametrization::init(const std::vector<double> &x,
                                 const std::vector<double> &y) {

   int n       = x.size() - 1;
   nb_interval = n;

   std::vector<double> b(n);
   for (int i = 1; i < n; ++i)
      b[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);

   myLineSet = new LineSet[n];
   for (int i = 0; i < n; ++i) {
      mySplineSet[i].a = y[i];
      mySplineSet[i].b = b[i];
      mySplineSet[i].x = x[i];
   }
   return;
}
std::vector
