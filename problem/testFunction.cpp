#include "testFunction.hpp"


string whichOperator(int op) {
  string s;
  if(op == 0)
  s="t ";
  else if(op == 1)
  s= "dt " ;
  else if(op == -1)
  s = " ";
  else
  s = " ";

  return s;
}

string whichOperator (int op, int cu) {
  string s;
  if(op == 0)
  s="u"+to_string(cu);
  else if(op == 1)
  s= "dx(u"+to_string(cu)+")" ;
  else if(op == 2)
  s= "dy(u"+to_string(cu)+")" ;
  else if(op == 3)
  s= "dz(u"+to_string(cu)+")" ;
  else if(op == op_dxx)
  s= "dxx(u"+to_string(cu)+")" ;
  else if(op == op_dxy)
  s= "dxy(u"+to_string(cu)+")" ;
  else if(op == op_dxz)
  s= "dxz(u"+to_string(cu)+")" ;
  else if(op == op_dyy)
  s= "dyy(u"+to_string(cu)+")" ;
  else if(op == op_dyz)
  s= "dyz(u"+to_string(cu)+")" ;
  else if(op == op_dzz)
  s= "dzz(u"+to_string(cu)+")" ;
  else
  s = "no Op";

  return s;
}
string whichOperatorV(int op, int cu) {
  string s;
  if(op == 0)
  s="v"+to_string(cu);
  else if(op == 1)
  s= "dx(v"+to_string(cu)+")" ;
  else if(op == 2)
  s= "dy(v"+to_string(cu)+")" ;
  else if(op == 3)
  s= "dz(v"+to_string(cu)+")" ;
  else if(op == op_dxx)
  s= "dxx(v"+to_string(cu)+")" ;
  else if(op == op_dxy)
  s= "dxy(v"+to_string(cu)+")" ;
  else if(op == op_dxz)
  s= "dxz(v"+to_string(cu)+")" ;
  else if(op == op_dyy)
  s= "dyy(v"+to_string(cu)+")" ;
  else if(op == op_dyz)
  s= "dyz(v"+to_string(cu)+")" ;
  else if(op == op_dzz)
  s= "dzz(v"+to_string(cu)+")" ;
  else
  s = "no Op";

  return s;
}
