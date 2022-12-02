#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "../calculus/RNM.hpp"
#include "../mesh/mesh.hpp"
#include "../fem/basisFct.hpp"
#include "../fem/fespace.hpp"
#include "../fem/function_expression.hpp"
#include "../operator/test_function.hpp"
#include "/usr/local/include/catch2/catch_all.hpp"

TEST_CASE("Test Expression class", "[Expression]") {
   typedef typename Mesh<1>::v_t v_t;
   Mesh<1> mesh(2, -1., 1.);
   P1_1D shapeFct;
   FESpace Vh(mesh, shapeFct);

   auto f1 = [](const v_t x) -> double { return x[0]; };
   auto f2 = [](const v_t x) -> double { return x[0] + 2; };

   Function1 f1h(Vh, f1);
   Function1 f2h(Vh, f2);

   auto expr1 = f1h.expr();
   auto expr2 = f2h.expr();

   SECTION("Evaluation") {
      REQUIRE(std::fabs(f1h.eval(0, -1. / Pi) + 1. / Pi) < Epsilon);
      REQUIRE(std::fabs(f1h.eval(0, 0.5) - 0.5) < Epsilon);
      REQUIRE(std::fabs(f1h.eval(0, 0.25) - 0.25) < Epsilon);
      REQUIRE(std::fabs(f1h.eval(0, -1. / Pi) - expr1->eval(0, -1. / Pi)) <
              Epsilon);
      REQUIRE(std::fabs(f1h.eval(0, 0.5) - expr1->eval(0, 0.5)) < Epsilon);
      REQUIRE(std::fabs(f1h.eval(0, 0.25) - expr1->eval(0, 0.25)) < Epsilon);
   }

   SECTION("Product of expression") {
      const auto p = expr1 * expr2;
      REQUIRE(std::fabs(p->eval(0, -1. / Pi) - (-1. / Pi) * (-1. / Pi + 2)) <
              Epsilon);
      REQUIRE(std::fabs(p->eval(0, 0.5) - 0.5 * (0.5 + 2)) < Epsilon);
      REQUIRE(std::fabs(p->eval(0, 0.25) - 0.25 * (0.25 + 2)) < Epsilon);
   }

   SECTION("operator ^ for expression") {
      const auto p = (expr2 ^ 3);
      REQUIRE(std::fabs(p->eval(0, -1. / Pi) - pow(-1. / Pi + 2, 3)) < Epsilon);
      REQUIRE(std::fabs(p->eval(0, 0.5) - pow(0.5 + 2, 3)) < Epsilon);
      REQUIRE(std::fabs(p->eval(0, 0.25) - pow(0.25 + 2, 3)) < Epsilon);
   }
   SECTION("Pow function on expression") {
      const auto p = pow(expr2, 3);
      REQUIRE(std::fabs(p->eval(0, -1. / Pi) - pow(-1. / Pi + 2, 3)) < Epsilon);
      REQUIRE(std::fabs(p->eval(0, 0.5) - pow(0.5 + 2, 3)) < Epsilon);
      REQUIRE(std::fabs(p->eval(0, 0.25) - pow(0.25 + 2, 3)) < Epsilon);
   }

   SECTION("Product of expression with constant") {
      const double c = Pi;
      const auto p   = c * expr1 * expr2;
      REQUIRE(std::fabs(p->eval(0, -1. / Pi) -
                        c * (-1. / Pi) * (-1. / Pi + 2)) < Epsilon);
      REQUIRE(std::fabs(p->eval(0, 0.5) - c * 0.5 * (0.5 + 2)) < Epsilon);
      REQUIRE(std::fabs(p->eval(0, 0.25) - c * 0.25 * (0.25 + 2)) < Epsilon);
   }

   SECTION("Test when wrap in a test function") {
      TestFunction<1> u(Vh, 1, 0);
      const auto u1 = expr1 * u;
      const auto &item(u1.getItem(0, 0));
      REQUIRE(item.function_p.get() != nullptr);
      REQUIRE(std::fabs(item.function_p->eval(0, 0.5) - 0.5) < Epsilon);

      const auto u2 = expr2 * u1;
      const auto &item2(u2.getItem(0, 0));
      REQUIRE(item2.function_p.get() != nullptr);
      REQUIRE(std::fabs(item2.function_p->eval(0, 0.5) - 0.5 * (0.5 + 2)) <
              Epsilon);

      const auto u3 = u2 * expr1;
      const auto &item3(u3.getItem(0, 0));
      REQUIRE(item3.function_p.get() != nullptr);
      REQUIRE(std::fabs(item3.function_p->eval(0, 0.5) -
                        0.5 * (0.5 + 2) * 0.5) < Epsilon);
   }
}
