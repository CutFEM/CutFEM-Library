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
#include "CutFEM_parameter.hpp"

const MeshParameter Parameter::h               = MeshParameter(1, 0, 0, 0);
const MeshParameter Parameter::measureIntegral = MeshParameter(0, 1, 0, 0);
const MeshParameter Parameter::measureElement  = MeshParameter(0, 0, 1, 0);
const MeshParameter Parameter::measureCut      = MeshParameter(0, 0, 0, 1);

Mul_Cst_Parameter operator*(const double a, const VirtualParameter &B) {
   return Mul_Cst_Parameter(a, B);
}
Mul_Cst_Parameter operator*(const VirtualParameter &B, const double a) {
   return Mul_Cst_Parameter(a, B);
}

Pow_Parameter pow(const VirtualParameter &A, int n) {
   return Pow_Parameter(A, n);
}
Pow_Parameter operator^(const VirtualParameter &A, int n) {
   return Pow_Parameter(A, n);
}

SumDiff_Parameter operator+(const VirtualParameter &A,
                            const VirtualParameter &B) {
   return SumDiff_Parameter(1, A, 1, B);
}
SumDiff_Parameter operator-(const VirtualParameter &A,
                            const VirtualParameter &B) {
   return SumDiff_Parameter(1, A, -1, B);
}

Mult_Parameter operator*(const VirtualParameter &A, const VirtualParameter &B) {
   return Mult_Parameter(A, B);
}
Inverse_Parameter operator/(double a, const VirtualParameter &A) {
   return Inverse_Parameter(a, A);
}
Inverse_Parameter inv(const VirtualParameter &A) {
   return Inverse_Parameter(1., A);
}

///  OLD STUFF
// --------------------------------------------------------------

//
// CutFEMParameter Parameter::kappa1   = CutFEMParameter("kappa1", fun_kappa1);
// CutFEMParameter Parameter::kappa2   = CutFEMParameter("kappa2", fun_kappa2);
// CutFEMParameter Parameter::lambdaG  = CutFEMParameter("lambdaG",
// fun_lambdaG); CutFEMParameter Parameter::lambdaB  =
// CutFEMParameter("lambdaB", fun_lambdaB);

// CutFEMParameter Parameter::lambdaB3 = CutFEMParameter("lambdaB3",
// fun_lambdaB3); CutFEMParameter Parameter::lambdaG3 =
// CutFEMParameter("lambdaG3", fun_lambdaG3);
//

// static double fun_kappa1(int i, double hh, double meas, double measK, double
// meas_Cut) {
//   if(CutFEM_ParameterList::find("mu") ) {
//     double mu1 = CutFEM_ParameterList::listParameter["mu"]->expression(0);
//     double mu2 = CutFEM_ParameterList::listParameter["mu"]->expression(1);
//     double alpha1 = 1;//meas_Cut / measK;
//     double alpha2 = 1;//(1 - alpha1);
//     return mu2*alpha1/(alpha2*mu1+alpha1*mu2);
//   }else{
//     return 0.5;
//   }
// }
// static double fun_kappa2(int i, double hh, double meas, double measK, double
// meas_Cut) {
//   if(CutFEM_ParameterList::find("mu") ) {
//     double mu1 = CutFEM_ParameterList::listParameter["mu"]->expression(0);
//     double mu2 = CutFEM_ParameterList::listParameter["mu"]->expression(1);
//     double alpha2 = 1;//meas_Cut / measK;
//     double alpha1 = 1;//(1 - alpha2);
//     return alpha2*mu1/(alpha2*mu1+alpha1*mu2);
//   }else{
//     return 0.5;
//   }
// }
// static double fun_lambdaG(int i, double hh, double meas, double measK, double
// meas_Cut) {
//   double mu1 =
//   (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(0)
//   : 1; double mu2 =
//   (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(1)
//   : 1; double kappa1 =
//   (CutFEM_ParameterList::find("kappa1"))?CutFEM_ParameterList::listParameter["kappa1"]->expression(i,hh,meas,meas_Cut)
//   : 0.5; double kappa2 =
//   (CutFEM_ParameterList::find("kappa2"))?CutFEM_ParameterList::listParameter["kappa2"]->expression(i,hh,meas,meas_Cut)
//   : 0.5; double gamma = meas / hh; double alphaK = hh;//meas_Cut/hh/hh;
//   // return 20 / hh / hh;//
//   return  (kappa1*mu1 + kappa2*mu2)*(100 + 10*gamma)/(alphaK);
// }
// static double fun_lambdaG3(int i, double hh, double meas, double measK,
// double meas_Cut) {
//   double mu1 =
//   (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(0)
//   : 1; double mu2 =
//   (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(1)
//   : 1; double kappa1 =
//   (CutFEM_ParameterList::find("kappa1"))?CutFEM_ParameterList::listParameter["kappa1"]->expression(i,hh,meas,meas_Cut)
//   : 0.5; double kappa2 =
//   (CutFEM_ParameterList::find("kappa2"))?CutFEM_ParameterList::listParameter["kappa2"]->expression(i,hh,meas,meas_Cut)
//   : 0.5; double gamma = meas / hh / hh; double alphaK =
//   hh/hh;//meas_Cut/hh/hh/hh;
//
//   // return 20 / hh / hh;//
//   return  (kappa1*mu1 + kappa2*mu2)*(100 + 10*gamma)/(alphaK);
// }
// static double fun_lambdaB(int i, double hh, double meas, double measK, double
// meas_Cut) {
//   double mu = (CutFEM_ParameterList::find("mu"))?
//   CutFEM_ParameterList::listParameter["mu"]->expression(0) : 1; double gammaK
//   = meas/hh; double alphaK = measK/hh/hh; return
//   mu/hh*(10+1e2*gammaK/alphaK);
//    // return 20 / hh / hh;
// }
// static double fun_lambdaB3(int i, double hh, double meas, double measK,
// double meas_Cut) {
//   double mu = (CutFEM_ParameterList::find("mu"))?
//   CutFEM_ParameterList::listParameter["mu"]->expression(0) : 1; double gammaK
//   = meas/hh/hh; double alphaK = measK/hh/hh/hh;
//
//   return mu/hh*(10+1e2*gammaK/alphaK);
//   // return 20 / hh / hh;
// }
//
//

//
//
