

using namespace globalVariable;

TEST_CASE("Test quadrature Lobatto 1D", "[Quadrature Lobatto]") {
   for (int order = 1; order <= 7; ++order) {

      const auto *lobatto = Lobatto(order);
      for (int i = 0; i < lobatto->getNbrOfQuads(); ++i) {
         //  REQUIRE(std::fabs(qf.p(i)[0] - ip[i][0]) < Epsilon);
         //  REQUIRE(std::fabs(qf.w(i) - w[i]) < Epsilon);
      }
   }
}

// TEST_CASE("Test quadrature Legendre 1D", "[Quadrature]")
// {

//     for (int order = 1; order < 15; ++order)
//     {
//         const QuadratureFormular<1> qf(order);
//         const auto &ip(qf.p());
//         const auto &w(qf.w());
//         REQUIRE(ip.size() == qf.size());
//         for (int i = 0; i < ip.size(); ++i)
//         {
//             REQUIRE(std::fabs(qf.p(i)[0] - ip[i][0]) < Epsilon);
//             REQUIRE(std::fabs(qf.w(i) - w[i]) < Epsilon);
//         }
//         Real val = 0.;
//         for (int i = 0; i < ip.size(); ++i)
//         {
//             val += w[i] * pow(ip[i][0], order);
//         }
//         REQUIRE(std::fabs(val - 1. / (order + 1.)) < 10*Epsilon);
//     }
// }
