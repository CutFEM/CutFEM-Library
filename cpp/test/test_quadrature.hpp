

using namespace globalVariable;

TEST_CASE("Test quadrature Lobatto 1D", "[Quadrature Lobatto]") {
   for (int order = 1; order <= 7; ++order) {

      const auto *lobatto = Lobatto(order);
      double val          = 0.;
      for (int i = 0; i < lobatto->getNbrOfQuads(); ++i) {
         const auto &ip(lobatto->at(i));
         val += ip.a * pow(ip.x, order);
      }
      REQUIRE(std::fabs(val - 1. / (order + 1.)) < 10 * Epsilon);
   }
}

TEST_CASE("Test quadrature Gauss-Legendre1D", "[Quadrature Gauss-Legendre]") {
   for (int order = 1; order <= 17; ++order) {

      const auto *qf = QF_Simplex<R1>(order);
      double val     = 0.;
      for (int i = 0; i < qf->getNbrOfQuads(); ++i) {
         const auto &ip(qf->at(i));
         val += ip.a * pow(ip.x, order);
      }
      REQUIRE(std::fabs(val - 1. / (order + 1.)) < 10 * Epsilon);
   }
}

// TEST_CASE("Test quadrature Gauss-Legendre 2D for quad",
//           "[Quadrature Gauss-Legendre 2D]") {
//    for (int order = 1; order <= 17; ++order) {

//       const auto *qf = QF_Simplex<R1>(order);
//       double val     = 0.;
//       for (int i = 0; i < qf->getNbrOfQuads(); ++i) {
//          const auto &ip(qf->at(i));
//          val += ip.a * pow(ip.x, order);
//       }
//       REQUIRE(std::fabs(val - 1. / (order + 1.)) < 10 * Epsilon);
//    }
// }
