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
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "../tool.hpp"

TEST_CASE("Test Time interface class in 2D", "[TimeInterface]") {

   Mesh2 Th(13, 13, -1., -1., 2, 2.);
   FESpace2 Vh(Th, DataFE<Mesh2>::P1);

   auto f = [](double *x, int ci, double t) -> double {
      return x[0] - (t - 2. / 3);
   };

   SECTION("Test with Lobatto 3 points") {
      const auto *qTime = Lobatto(3);
      auto n            = qTime->n;
      double dt         = 1e-3;
      TimeInterface<Mesh2> interface(qTime);
      REQUIRE(interface.interface().size() == 3);
      double t = 0;
      for (int iter = 0; iter < 10; ++iter) {
         for (int i = 0; i < n; ++i) {
            FunFEM<Mesh2> ls(Vh, f, t + i * dt);
            interface.init(i, Th, ls);
            const auto *gamma = interface[i];

            REQUIRE(isEqual((*gamma)(3)[0], t + i * dt - 2. / 3));
         }
         t += (n - 1) * dt;
      }
   }
}