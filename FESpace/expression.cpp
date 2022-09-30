#include "expression.hpp"
#include "../problem/CutFEM_parameter.hpp"


ExpressionMultConst operator*(const ExpressionVirtual& f1, double cc) {
  return ExpressionMultConst(f1, cc);
}
ExpressionMultConst operator*(double cc, const ExpressionVirtual& f1) {
  return ExpressionMultConst(f1, cc);
}
ExpressionMultConst operator*(const ExpressionVirtual& f1, const Normal_Component_X& cc){
  return ExpressionMultConst(f1, cc);
}
ExpressionMultConst operator*(const ExpressionVirtual& f1, const Normal_Component_Y& cc){
  return ExpressionMultConst(f1, cc);
}
ExpressionMultConst operator*(const ExpressionVirtual& f1, const Normal_Component_Z& cc){
  return ExpressionMultConst(f1, cc);
}
// ExpressionMultConst operator*(const CutFEM_Parameter& v, const ExpressionVirtual& f1){
//   return ExpressionMultConst(f1, R2(v.val1, v.val2));
// }

ExpressionAbs fabs(const ExpressionVirtual& f1) {
  return ExpressionAbs(f1);
}

ExpressionSqrt sqrt(const ExpressionVirtual& f1) {
  return ExpressionSqrt(f1);
}

ExpressionProduct operator*(const ExpressionVirtual& f1, const ExpressionVirtual& f2) {
  return ExpressionProduct(f1, f2);
}
ExpressionDivision operator/(const ExpressionVirtual& f1, const ExpressionVirtual& f2) {
  return ExpressionDivision(f1, f2);
}

ExpressionSum operator+(const ExpressionVirtual& f1, const ExpressionVirtual& f2) {
  return ExpressionSum(f1, f2);
}


ExpressionDif operator-(const ExpressionVirtual& f1, const ExpressionVirtual& f2) {
  return ExpressionDif(f1, f2);
}
ExpressionNormal2 operator*(const FunFEM<Mesh2>& f1, const Normal& n){
  return ExpressionNormal2(f1,n);
}
ExpressionNormal2 operator*(const FunFEM<Mesh2>& f1, const Tangent& n){
  return ExpressionNormal2(f1,n);
}
ExpressionNormal2 operator*(const FunFEM<Mesh2>& f1, const Conormal& n){
  return ExpressionNormal2(f1,n);
}
ExpressionNormal3 operator*(const FunFEM<Mesh3>& f1, const Normal& n){
  return ExpressionNormal3(f1);
}

// ExpressionAverage average(const ExpressionVirtual & f1, const double kk1, const double kk2){
//   return ExpressionAverage(f1,kk1,kk2);
// }
// ExpressionAverage jump(const ExpressionVirtual & f1, const double kk1, const double kk2){
//   return ExpressionAverage(f1,1,-1);
// }
// ExpressionAverage operator* (double c, const ExpressionAverage& fh){
//   return ExpressionAverage(fh.fun1,c*fh.k1, c*fh.k2);
// }
// ExpressionAverage operator* (const ExpressionAverage& fh, double c){
//   return ExpressionAverage(fh.fun1,c*fh.k1, c*fh.k2);
// }

ExpressionBurgerFlux burgerFlux(const ExpressionVirtual & f1) {
  return ExpressionBurgerFlux(f1);
}
ExpressionNormalBurgerFlux burgerFlux(const ExpressionVirtual & f1, const Normal& n ) {
  return ExpressionNormalBurgerFlux(f1);
}

ExpressionDSx2  dxS (const FunFEM<Mesh2>& f1){return ExpressionDSx2(f1);}
ExpressionDSy2  dyS (const FunFEM<Mesh2>& f1){return ExpressionDSy2(f1);}
ExpressionDivS2 divS(const FunFEM<Mesh2>& f1){return ExpressionDivS2(f1);}
ExpressionDSx3  dxS (const FunFEM<Mesh3>& f1){return ExpressionDSx3(f1);}
ExpressionDSy3  dyS (const FunFEM<Mesh3>& f1){return ExpressionDSy3(f1);}
ExpressionDSz3  dzS (const FunFEM<Mesh3>& f1){return ExpressionDSz3(f1);}
ExpressionDivS3 divS(const FunFEM<Mesh3>& f1){return ExpressionDivS3(f1);}
