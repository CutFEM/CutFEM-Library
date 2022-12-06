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

TEST_CASE("Test Expression class for scalar functions", "[Expression]") {

   Mesh2 mesh(2, 2, -1., -1., 2., 2);
   FESpace2 Vh(mesh, DataFE<Mesh2>::P1);
   //    FESpace2 Rh(mesh, DataFE<Mesh2>::RT0);

   auto f1 = [](double *x) -> double { return x[0] + x[1]; };
   auto f2 = [](double *x) -> double { return x[0] - 2; };

   FunFEM<Mesh2> f1h(Vh, f1);
   FunFEM<Mesh2> f2h(Vh, f2);

   auto expr1 = f1h.expr();
   auto expr2 = f2h.expr();
   R2 x1(-1. / Pi, 3. / Pi);
   R2 x2(0, 0.3333);

   {
      REQUIRE(isEqual(f1h.eval(0, x1, 0, op_id), f1(x1)));
      REQUIRE(isEqual(f1h.eval(0, x1), f1(x1)));
      REQUIRE(isEqual(f1h.eval(0, x1), expr1->eval(0, x1, nullptr)));

      REQUIRE(isEqual(f1h.eval(0, x2, 0, op_id), f1(x2)));
      REQUIRE(isEqual(f1h.eval(0, x2), f1(x2)));
      REQUIRE(isEqual(f1h.eval(0, x2), expr1->eval(0, x2, nullptr)));
   }
   {
      const auto p = expr1 * expr2;
      REQUIRE(isEqual(p->eval(0, x1, nullptr), f1(x1) * f2(x1)));
      REQUIRE(isEqual(p->eval(0, x2, nullptr), f1(x2) * f2(x2)));
   }
   {
      const auto p = (expr2 ^ 3);
      REQUIRE(isEqual(p->eval(0, x1, nullptr), pow(f2(x1), 3), 1e-12));
      REQUIRE(isEqual(p->eval(0, x2, nullptr), pow(f2(x2), 3), 1e-12));
   }
   {
      const auto p = pow(expr2, 3); //(expr2 ^ 3);
      REQUIRE(isEqual(p->eval(0, x1, nullptr), pow(f2(x1), 3), 1e-12));
      REQUIRE(isEqual(p->eval(0, x2, nullptr), pow(f2(x2), 3), 1e-12));

      const auto q = fabs(p);
      REQUIRE(
          isEqual(q->eval(0, x1, nullptr), std::fabs(pow(f2(x1), 3)), 1e-12));
      REQUIRE(
          isEqual(q->eval(0, x2, nullptr), std::fabs(pow(f2(x2), 3)), 1e-12));

      const auto r = sqrt(q);
      REQUIRE(isEqual(r->eval(0, x1, nullptr),
                      std::sqrt(std::fabs(pow(f2(x1), 3))), 1e-12));
      REQUIRE(isEqual(r->eval(0, x2, nullptr),
                      std::sqrt(std::fabs(pow(f2(x2), 3))), 1e-12));
   }
   {
      const auto p = dx(expr2);
      REQUIRE(p->cu == 0);
      REQUIRE(p->op == op_dx);
      REQUIRE(p->opt == op_id);
      REQUIRE(isEqual(p->eval(0, x1, nullptr), 1.));

      const auto q = dy(expr2);
      REQUIRE(q->cu == 0);
      REQUIRE(q->op == op_dy);
      REQUIRE(isEqual(q->eval(0, x1, nullptr), 0.));

      const auto r = dt(expr2);
      REQUIRE(r->cu == 0);
      REQUIRE(r->opt == op_dx);
   }
   {
      const auto p = expr1 + expr2;
      REQUIRE(isEqual(p->eval(0, x1, nullptr), f1(x1) + f2(x1)));
      REQUIRE(isEqual(p->eval(0, x2, nullptr), f1(x2) + f2(x2)));

      const auto q = p / expr2;
      REQUIRE(isEqual(q->eval(0, x1, nullptr), (f1(x1) + f2(x1)) / f2(x1)));
      REQUIRE(isEqual(q->eval(0, x2, nullptr), (f1(x2) + f2(x2)) / f2(x2)));

      const auto r = p - expr2;
      REQUIRE(isEqual(r->eval(0, x1, nullptr), f1(x1)));
      REQUIRE(isEqual(r->eval(0, x2, nullptr), f1(x2)));
   }
}
