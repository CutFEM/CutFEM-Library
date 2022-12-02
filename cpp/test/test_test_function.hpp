#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "../tool.hpp"
#include "/usr/local/include/catch2/catch_all.hpp"

TEST_CASE("Test Test Function class for scalar case", "[TestFunction]") {

   Mesh2 mesh(3, 3, 0., 0., 1., 1);
   FESpace2 Vh(mesh, DataFE<Mesh2>::P1);
   FESpace2 Rh(mesh, DataFE<Mesh2>::RT0);

   SECTION("Constructor") {
      TestFunction<2> u(Vh, 1);

      REQUIRE(u.isScalar());
      const auto item(u.getItem({0, 0}, 0));
      REQUIRE(item.du == op_id);
      REQUIRE(item.c == 1.);
      REQUIRE(item.cu == 0);
   }

   SECTION("Test dx, dy, dz, dt operator") {
      TestFunction<2> u(Vh, 1);
      {
         TestFunction<2> dxu(dx(u));
         const auto item(dxu.getItem({0, 0}, 0));
         REQUIRE(dxu.isScalar());
         REQUIRE(dxu.sizeItemList(0, 0) == 1);
         REQUIRE(item.du == op_dx);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
      }
      {
         TestFunction<2> dyu(dy(u));
         const auto item(dyu.getItem({0, 0}, 0));
         REQUIRE(dyu.isScalar());
         REQUIRE(dyu.sizeItemList(0, 0) == 1);
         REQUIRE(item.du == op_dy);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
      }
      {
         TestFunction<2> dzu(dz(u));
         const auto item(dzu.getItem({0, 0}, 0));
         REQUIRE(dzu.isScalar());
         REQUIRE(dzu.sizeItemList(0, 0) == 1);
         REQUIRE(item.du == op_dz);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
      }
      {
         TestFunction<2> dtu(dt(u));
         const auto item(dtu.getItem({0, 0}, 0));
         REQUIRE(dtu.isScalar());
         REQUIRE(dtu.sizeItemList(0, 0) == 1);
         REQUIRE(item.dtu == op_dx);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
      }
   }

   SECTION("Test operator +, -") {
      TestFunction<2> u(Vh, 1);
      TestFunction<2> v(Vh, 1);

      {
         TestFunction<2> uv(u + dx(v));

         REQUIRE(uv.isScalar());
         REQUIRE(uv.sizeItemList(0, 0) == 2);
         const auto item(uv.getItem({0, 0}, 0));
         REQUIRE(item.du == op_id);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         const auto item2(uv.getItem({0, 0}, 1));
         REQUIRE(item2.du == op_dx);
         REQUIRE(item2.c == 1.);
         REQUIRE(item2.cu == 0);
      }
      {
         TestFunction<2> vu(dx(u) - dz(v));
         REQUIRE(vu.isScalar());
         REQUIRE(vu.sizeItemList(0, 0) == 2);
         const auto item(vu.getItem({0, 0}, 0));
         REQUIRE(item.du == op_dx);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         const auto item2(vu.getItem({0, 0}, 1));
         REQUIRE(item2.du == op_dz);
         REQUIRE(item2.c == -1.);
         REQUIRE(item2.cu == 0);
      }
   }

   SECTION("Test operator *") {
      TestFunction<2> u(Vh, 1);
      TestFunction<2> v(Vh, 1);
      TestFunction<2> uu(5. * u * 3.);

      REQUIRE(uu.isScalar());
      REQUIRE(uu.sizeItemList(0, 0) == 1);

      const auto &item(uu.getItem({0, 0}, 0));
      REQUIRE(item.du == op_id);
      REQUIRE(item.c == 15.);
      REQUIRE(item.cu == 0);
      REQUIRE(item.face_side_ == -1);
      REQUIRE(item.fespace == &Vh);
      REQUIRE(item.root_fun_p == &u);
      REQUIRE(item.expru == nullptr);
   }

   SECTION("Test simplify") {
      TestFunction<2> u(Vh, 1);
      TestFunction<2> v(Vh, 1);
      TestFunction<2> uv(u + u + v);

      REQUIRE(uv.isScalar());
      REQUIRE(uv.sizeItemList(0, 0) == 2);

      const auto &item0(uv.getItem({0, 0}, 0));
      const auto &item1(uv.getItem({0, 0}, 1));

      REQUIRE(item0.du == op_id);
      REQUIRE(item0.c == 2.);
      REQUIRE(item0.cu == 0);
      REQUIRE(item0.face_side_ == -1);
      REQUIRE(item1.du == op_id);
      REQUIRE(item1.c == 1.);
      REQUIRE(item1.cu == 0);
      REQUIRE(item1.face_side_ == -1);
   }

   SECTION("Test grad") {
      TestFunction<2> u(Vh, 1);
      TestFunction<2> grad_u(grad(u));

      REQUIRE(!grad_u.isScalar());
      REQUIRE(grad_u.sizeItemList(0, 0) == 1);
      REQUIRE(grad_u.sizeItemList(1, 0) == 1);

      {
         const auto &item(grad_u.getItem({0, 0}, 0));
         REQUIRE(item.du == op_dx);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
      }
      {
         const auto &item(grad_u.getItem({1, 0}, 0));
         REQUIRE(item.du == op_dy);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
      }
   }

   SECTION("Test grad^2") {
      TestFunction<2> u(Vh, 1);
      TestFunction<2> grad_u(grad(grad(u)));

      REQUIRE(!grad_u.isScalar());
      REQUIRE(grad_u.sizeItemList(0, 0) == 1);
      REQUIRE(grad_u.sizeItemList(1, 0) == 1);
      REQUIRE(grad_u.sizeItemList(0, 1) == 1);
      REQUIRE(grad_u.sizeItemList(1, 1) == 1);

      {
         const auto &item(grad_u.getItem({0, 0}, 0));
         REQUIRE(item.du == op_dxx);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
      }
      {
         const auto &item(grad_u.getItem({1, 0}, 0));
         REQUIRE(item.du == op_dxy);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
      }
      {
         const auto &item(grad_u.getItem({0, 1}, 0));
         REQUIRE(item.du == op_dxy);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
      }
      {
         const auto &item(grad_u.getItem({1, 1}, 0));
         REQUIRE(item.du == op_dyy);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
      }
   }

   SECTION("Test Average") {
      TestFunction<2> u(Vh, 1);
      TestFunction<2> jump_u(average(u, 1.5, -1.5));

      REQUIRE(jump_u.isScalar());
      REQUIRE(jump_u.sizeItemList(0, 0) == 2);

      {
         const auto &item(jump_u.getItem({0, 0}, 0));
         REQUIRE(item.du == op_id);
         REQUIRE(item.c == 1.5);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == 0);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
      }
      {
         const auto &item(jump_u.getItem({0, 0}, 1));
         REQUIRE(item.du == op_id);
         REQUIRE(item.c == -1.5);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == 1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
      }
   }

   SECTION("Test Normal multiplication with scalar") {
      Normal N;
      TestFunction<2> u(Vh, 1);
      TestFunction<2> un(u * N);

      REQUIRE(!un.isScalar());
      REQUIRE(un.nbRow() == 2);
      REQUIRE(un.nbCol() == 1);
      {
         const auto &item(un.getItem({0, 0}, 0));
         REQUIRE(item.du == op_id);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
         REQUIRE(item.ar_nu.size() == 1);
         REQUIRE(item.ar_nu[0] == 0);
      }
      {
         const auto &item(un.getItem({1, 0}, 0));
         REQUIRE(item.du == op_id);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
         REQUIRE(item.ar_nu.size() == 1);
         REQUIRE(item.ar_nu[0] == 1);
      }
   }

   SECTION("Test Normal multiplication with vector") {
      Normal N;
      TestFunction<2> u(Vh, 1);
      TestFunction<2> un(grad(u) * N);

      REQUIRE(un.isScalar());
      REQUIRE(un.sizeItemList(0, 0) == 2);

      {
         const auto &item(un.getItem({0, 0}, 0));
         REQUIRE(item.du == op_dx);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
         REQUIRE(item.ar_nu.size() == 1);
         REQUIRE(item.ar_nu[0] == 0);
      }
      {
         const auto &item(un.getItem({0, 0}, 1));
         REQUIRE(item.du == op_dy);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
         REQUIRE(item.ar_nu.size() == 1);
         REQUIRE(item.ar_nu[0] == 1);
      }
   }

   SECTION("Test Normal multiplication with matrix") {
      Normal N;
      TestFunction<2> u(Vh, 1);
      TestFunction<2> un(grad(grad(u)) * N);

      REQUIRE(!un.isScalar());
      REQUIRE(un.sizeItemList(0, 0) == 2);
      REQUIRE(un.sizeItemList(1, 0) == 2);

      {
         const auto &item(un.getItem({0, 0}, 0));
         REQUIRE(item.du == op_dxx);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
         REQUIRE(item.ar_nu.size() == 1);
         REQUIRE(item.ar_nu[0] == 0);
      }
      {
         const auto &item(un.getItem({0, 0}, 1));
         REQUIRE(item.du == op_dxy);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
         REQUIRE(item.ar_nu.size() == 1);
         REQUIRE(item.ar_nu[0] == 1);
      }
      {
         const auto &item(un.getItem({1, 0}, 0));
         REQUIRE(item.du == op_dxy);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
         REQUIRE(item.ar_nu.size() == 1);
         REQUIRE(item.ar_nu[0] == 0);
      }
      {
         const auto &item(un.getItem({1, 0}, 1));
         REQUIRE(item.du == op_dyy);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
         REQUIRE(item.ar_nu.size() == 1);
         REQUIRE(item.ar_nu[0] == 1);
      }
   }

   SECTION("Test Tangent multiplication with scalar") {
      Tangent T;
      TestFunction<2> u(Vh, 1);
      TestFunction<2> un(u * T);

      REQUIRE(!un.isScalar());
      REQUIRE(un.nbRow() == 2);
      REQUIRE(un.nbCol() == 1);
      {
         const auto &item(un.getItem({0, 0}, 0));
         REQUIRE(item.du == op_id);
         REQUIRE(item.c == -1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
         REQUIRE(item.ar_nu.size() == 1);
         REQUIRE(item.ar_nu[0] == 1);
      }
      {
         const auto &item(un.getItem({1, 0}, 0));
         REQUIRE(item.du == op_id);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Vh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru == nullptr);
         REQUIRE(item.ar_nu.size() == 1);
         REQUIRE(item.ar_nu[0] == 0);
      }
   }

   // SECTION("Test Expression multiplication with scalar") {

   //    auto g = [](double *x) -> double { return 1; };
   //    FunFEM<Mesh2> gh(Vh, g);
   //    // ExpressionFunFEM<Mesh2> expr(gh, 0, 0);
   //    auto expr = gh.expr();
   //    TestFunction<2> u(Vh, 1);
   //    TestFunction<2> unn(gh.expr() * u);
   //    TestFunction<2> un(expr * expr * u);

   //    REQUIRE(un.isScalar());
   //    REQUIRE(un.nbRow() == 1);
   //    REQUIRE(un.nbCol() == 1);
   //    {
   //       const auto &item(un.getItem({0, 0}, 0));
   //       REQUIRE(item.du == op_id);
   //       REQUIRE(item.c == 1.);
   //       REQUIRE(item.cu == 0);
   //       REQUIRE(item.face_side_ == -1);
   //       REQUIRE(item.fespace == &Vh);
   //       REQUIRE(item.root_fun_p == &u);
   //       REQUIRE(item.expru.get() != nullptr);
   //       REQUIRE(item.ar_nu.size() == 0);
   //    }
   // }

   SECTION("Test Eps operator") {

      TestFunction<2> u(Rh, 2);
      TestFunction<2> epsU(Eps(u));

      REQUIRE(!epsU.isScalar());
      REQUIRE(epsU.nbRow() == 2);
      REQUIRE(epsU.nbCol() == 2);
      {
         REQUIRE(epsU.sizeItemList(0, 0) == 1);
         const auto &item(epsU.getItem({0, 0}, 0));
         REQUIRE(item.du == op_dx);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Rh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru.get() == nullptr);
         REQUIRE(item.ar_nu.size() == 0);
      }
      {
         REQUIRE(epsU.sizeItemList(1, 0) == 2);
         const auto &item(epsU.getItem({1, 0}, 0));
         REQUIRE(item.du == op_dx);
         REQUIRE(item.c == 0.5);
         REQUIRE(item.cu == 1);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Rh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru.get() == nullptr);
         REQUIRE(item.ar_nu.size() == 0);
         const auto &item1(epsU.getItem({1, 0}, 1));
         REQUIRE(item1.du == op_dy);
         REQUIRE(item1.c == 0.5);
         REQUIRE(item1.cu == 0);
         REQUIRE(item1.face_side_ == -1);
         REQUIRE(item1.fespace == &Rh);
         REQUIRE(item1.root_fun_p == &u);
         REQUIRE(item1.expru.get() == nullptr);
         REQUIRE(item1.ar_nu.size() == 0);
      }
      {
         REQUIRE(epsU.sizeItemList(0, 1) == 2);
         const auto &item(epsU.getItem({0, 1}, 0));
         REQUIRE(item.du == op_dy);
         REQUIRE(item.c == 0.5);
         REQUIRE(item.cu == 0);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Rh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru.get() == nullptr);
         REQUIRE(item.ar_nu.size() == 0);
         const auto &item1(epsU.getItem({0, 1}, 1));
         REQUIRE(item1.du == op_dx);
         REQUIRE(item1.c == 0.5);
         REQUIRE(item1.cu == 1);
         REQUIRE(item1.face_side_ == -1);
         REQUIRE(item1.fespace == &Rh);
         REQUIRE(item1.root_fun_p == &u);
         REQUIRE(item1.expru.get() == nullptr);
         REQUIRE(item1.ar_nu.size() == 0);
      }
      {
         REQUIRE(epsU.sizeItemList(1, 1) == 1);
         const auto &item(epsU.getItem({1, 1}, 0));
         REQUIRE(item.du == op_dy);
         REQUIRE(item.c == 1.);
         REQUIRE(item.cu == 1);
         REQUIRE(item.face_side_ == -1);
         REQUIRE(item.fespace == &Rh);
         REQUIRE(item.root_fun_p == &u);
         REQUIRE(item.expru.get() == nullptr);
         REQUIRE(item.ar_nu.size() == 0);
      }
   }
   // SECTION("Test average operator") {
   //   TestFunction<1> u(Vh, 1);
   //   TestFunction<1> avgu(average(dx(u), 1., -1.));

   //   REQUIRE(avgu.isScalar());
   //   REQUIRE(avgu.sizeItemList(0) == 2);

   //   const auto &item0(avgu.getItem(0, 0));
   //   const auto &item1(avgu.getItem(0, 1));

   //   REQUIRE(item0.oper == op_dx);
   //   REQUIRE(item0.cst == 1.);
   //   REQUIRE(item0.comp == 0);
   //   REQUIRE(item0.side_edge == 0);
   //   REQUIRE(item0.fespace_p == &Vh);
   //   REQUIRE(item0.root_fun_p == &u);
   //   REQUIRE(item0.function_p == nullptr);

   //   REQUIRE(item1.oper == op_dx);
   //   REQUIRE(item1.cst == -1.);
   //   REQUIRE(item1.comp == 0);
   //   REQUIRE(item1.side_edge == 1);
   //   REQUIRE(item1.fespace_p == &Vh);
   //   REQUIRE(item1.root_fun_p == &u);
   //   REQUIRE(item1.function_p == nullptr);
   // }
}