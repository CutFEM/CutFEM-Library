#include "expression.hpp"


ExpressionMultConst operator*(const ExpressionVirtual& f1, const double& cc) {
  return ExpressionMultConst(f1, cc);
}
ExpressionMultConst operator*(const double& cc, const ExpressionVirtual& f1) {
  return ExpressionMultConst(f1, cc);
}
ExpressionAbs fabs(const ExpressionVirtual& f1) {
  return ExpressionAbs(f1);
}

ExpressionProduct operator*(const ExpressionVirtual& f1, const ExpressionVirtual& f2) {
  return ExpressionProduct(f1, f2);
}

ExpressionSum operator+(const ExpressionVirtual& f1, const ExpressionVirtual& f2) {
  return ExpressionSum(f1, f2);
}


ExpressionDif operator-(const ExpressionVirtual& f1, const ExpressionVirtual& f2) {
  return ExpressionDif(f1, f2);
}
const int Normal::idx[3] = {0,1,2};
const int Tangent::idx[3] = {1,0,2}; // only in 2D


ExpressionDSx2  dxS (const FunFEM<Mesh2>& f1){return ExpressionDSx2(f1);}
ExpressionDSy2  dyS (const FunFEM<Mesh2>& f1){return ExpressionDSy2(f1);}
ExpressionDivS2 divS(const FunFEM<Mesh2>& f1){return ExpressionDivS2(f1);}
ExpressionDSx3  dxS (const FunFEM<Mesh3>& f1){return ExpressionDSx3(f1);}
ExpressionDSy3  dyS (const FunFEM<Mesh3>& f1){return ExpressionDSy3(f1);}
ExpressionDSz3  dzS (const FunFEM<Mesh3>& f1){return ExpressionDSz3(f1);}
ExpressionDivS3 divS(const FunFEM<Mesh3>& f1){return ExpressionDivS3(f1);}
