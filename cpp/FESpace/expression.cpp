#include "expression.hpp"
#include "../problem/CutFEM_parameter.hpp"

std::shared_ptr<ExpressionMultConst>
operator*(const std::shared_ptr<ExpressionVirtual> &f1, double cc) {
   return std::make_shared<ExpressionMultConst>(f1, cc);
}
std::shared_ptr<ExpressionMultConst>
operator*(double cc, const std::shared_ptr<ExpressionVirtual> &f1) {
   return std::make_shared<ExpressionMultConst>(f1, cc);
}
std::shared_ptr<ExpressionMultConst>
operator*(const std::shared_ptr<ExpressionVirtual> &f1,
          const Normal_Component_X &cc) {
   return std::make_shared<ExpressionMultConst>(f1, cc);
}
std::shared_ptr<ExpressionMultConst>
operator*(const std::shared_ptr<ExpressionVirtual> &f1,
          const Normal_Component_Y &cc) {
   return std::make_shared<ExpressionMultConst>(f1, cc);
}
std::shared_ptr<ExpressionMultConst>
operator*(const std::shared_ptr<ExpressionVirtual> &f1,
          const Normal_Component_Z &cc) {
   return std::make_shared<ExpressionMultConst>(f1, cc);
}

std::shared_ptr<ExpressionAbs>
fabs(const std::shared_ptr<ExpressionVirtual> &f1) {
   return std::make_shared<ExpressionAbs>(f1);
}

std::shared_ptr<ExpressionProduct>
operator*(const std::shared_ptr<ExpressionVirtual> &f1,
          const std::shared_ptr<ExpressionVirtual> &f2) {
   return std::make_shared<ExpressionProduct>(f1, f2);
}

std::shared_ptr<ExpressionPow> pow(const std::shared_ptr<ExpressionVirtual> &f1,
                                   const double nn) {
   return std::make_shared<ExpressionPow>(f1, nn);
}
std::shared_ptr<ExpressionPow>
operator^(const std::shared_ptr<ExpressionVirtual> &f1, const double nn) {
   return std::make_shared<ExpressionPow>(f1, nn);
}
std::shared_ptr<ExpressionPow>
sqrt(const std::shared_ptr<ExpressionVirtual> &f1) {
   return pow(f1, 1. / 2);
}

std::shared_ptr<ExpressionDivision>
operator/(const std::shared_ptr<ExpressionVirtual> &f1,
          const std::shared_ptr<ExpressionVirtual> &f2) {
   return std::make_shared<ExpressionDivision>(f1, f2);
}

std::shared_ptr<ExpressionSum>
operator+(const std::shared_ptr<ExpressionVirtual> &f1,
          const std::shared_ptr<ExpressionVirtual> &f2) {
   return std::make_shared<ExpressionSum>(f1, f2);
}

std::shared_ptr<ExpressionSum>
operator-(const std::shared_ptr<ExpressionVirtual> &f1,
          const std::shared_ptr<ExpressionVirtual> &f2) {
   return f1 + (-1. * f2);
}

ExpressionNormal2 operator*(const FunFEM<Mesh2> &f1, const Normal &n) {
   return ExpressionNormal2(f1, n);
}
ExpressionNormal2 operator*(const FunFEM<Mesh2> &f1, const Tangent &n) {
   return ExpressionNormal2(f1, n);
}
ExpressionNormal2 operator*(const FunFEM<Mesh2> &f1, const Conormal &n) {
   return ExpressionNormal2(f1, n);
}
ExpressionNormal3 operator*(const FunFEM<Mesh3> &f1, const Normal &n) {
   return ExpressionNormal3(f1);
}

ExpressionAverage average(const std::shared_ptr<ExpressionVirtual> &f1,
                          const double kk1, const double kk2) {
   return ExpressionAverage(f1, kk1, kk2);
}
ExpressionAverage jump(const std::shared_ptr<ExpressionVirtual> &f1,
                       const double kk1, const double kk2) {
   return ExpressionAverage(f1, 1, -1);
}
ExpressionAverage operator*(double c, const ExpressionAverage &fh) {
   return ExpressionAverage(fh.fun1, c * fh.k1, c * fh.k2);
}
ExpressionAverage operator*(const ExpressionAverage &fh, double c) {
   return ExpressionAverage(fh.fun1, c * fh.k1, c * fh.k2);
}

ExpressionBurgerFlux burgerFlux(const ExpressionVirtual &f1) {
   return ExpressionBurgerFlux(f1);
}
ExpressionNormalBurgerFlux burgerFlux(const ExpressionVirtual &f1,
                                      const Normal &n) {
   return ExpressionNormalBurgerFlux(f1);
}

ExpressionDSx2 dxS(const FunFEM<Mesh2> &f1) { return ExpressionDSx2(f1); }
ExpressionDSy2 dyS(const FunFEM<Mesh2> &f1) { return ExpressionDSy2(f1); }
ExpressionDivS2 divS(const FunFEM<Mesh2> &f1) { return ExpressionDivS2(f1); }
ExpressionDSx3 dxS(const FunFEM<Mesh3> &f1) { return ExpressionDSx3(f1); }
ExpressionDSy3 dyS(const FunFEM<Mesh3> &f1) { return ExpressionDSy3(f1); }
ExpressionDSz3 dzS(const FunFEM<Mesh3> &f1) { return ExpressionDSz3(f1); }
ExpressionDivS3 divS(const FunFEM<Mesh3> &f1) { return ExpressionDivS3(f1); }
