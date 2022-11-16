#include "testFunction.hpp"

void f_id(RNMK_ &x, int cu, int du) {}
void f_ln(RNMK_ &x, int cu, int du) {
   for (int i = 0; i < x.N(); ++i) {
      x(i, cu, du) = std::log(x(i, cu, du));
   }
}

std::string whichOperator(int op) {
   std::string s;
   if (op == 0)
      s = "t ";
   else if (op == 1)
      s = "dt ";
   else if (op == -1)
      s = " ";
   else
      s = " ";

   return s;
}

std::string whichOperator(int op, int cu) {
   std::string s;
   if (op == 0)
      s = "u" + std::to_string(cu);
   else if (op == 1)
      s = "dx(u" + std::to_string(cu) + ")";
   else if (op == 2)
      s = "dy(u" + std::to_string(cu) + ")";
   else if (op == 3)
      s = "dz(u" + std::to_string(cu) + ")";
   else if (op == op_dxx)
      s = "dxx(u" + std::to_string(cu) + ")";
   else if (op == op_dxy)
      s = "dxy(u" + std::to_string(cu) + ")";
   else if (op == op_dxz)
      s = "dxz(u" + std::to_string(cu) + ")";
   else if (op == op_dyy)
      s = "dyy(u" + std::to_string(cu) + ")";
   else if (op == op_dyz)
      s = "dyz(u" + std::to_string(cu) + ")";
   else if (op == op_dzz)
      s = "dzz(u" + std::to_string(cu) + ")";
   else
      s = "no Op";

   return s;
}
std::string whichOperatorV(int op, int cu) {
   std::string s;
   if (op == 0)
      s = "v" + std::to_string(cu);
   else if (op == 1)
      s = "dx(v" + std::to_string(cu) + ")";
   else if (op == 2)
      s = "dy(v" + std::to_string(cu) + ")";
   else if (op == 3)
      s = "dz(v" + std::to_string(cu) + ")";
   else if (op == op_dxx)
      s = "dxx(v" + std::to_string(cu) + ")";
   else if (op == op_dxy)
      s = "dxy(v" + std::to_string(cu) + ")";
   else if (op == op_dxz)
      s = "dxz(v" + std::to_string(cu) + ")";
   else if (op == op_dyy)
      s = "dyy(v" + std::to_string(cu) + ")";
   else if (op == op_dyz)
      s = "dyz(v" + std::to_string(cu) + ")";
   else if (op == op_dzz)
      s = "dzz(v" + std::to_string(cu) + ")";
   else
      s = "no Op";

   return s;
}
