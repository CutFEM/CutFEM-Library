#include "CutFEM_parameter.hpp"
// std::map<std::string,ParameterCutFEM*> CutFEM_ParameterList::listParameter = std::map<std::string,ParameterCutFEM*>();

// static double fun_meas(int i, double hh, double meas, double measK, double meas_Cut) {return meas;}
// static double fun_invmeas(int i, double hh, double meas, double measK, double meas_Cut) {return 1./meas;}
// static double fun_h(int i, double hh, double meas, double measK, double meas_Cut) {return hh;}
// static double fun_invh(int i, double hh, double meas, double measK, double meas_Cut) {return 1./hh;}
//
// ParameterCutFEM Parameter::h        = ParameterCutFEM("h", fun_h);
// ParameterCutFEM Parameter::invh     = ParameterCutFEM("invh", fun_invh);
// ParameterCutFEM Parameter::meas     = ParameterCutFEM("meas", fun_meas);
// ParameterCutFEM Parameter::invmeas  = ParameterCutFEM("invmeas", fun_invmeas);




Mul_Cst_Parameter operator*(const double a, const Virtual_Parameter& B) {
  return Mul_Cst_Parameter(a, B);
}
Mul_Cst_Parameter operator*(const Virtual_Parameter& B, const double a){
  return Mul_Cst_Parameter(a, B);
}

Pow_Parameter pow(const Virtual_Parameter& A, int n) {
  return Pow_Parameter(A,n);
}
Pow_Parameter operator^(const Virtual_Parameter& A, int n){
  return Pow_Parameter(A,n);
}

SumDiff_Parameter operator+(const Virtual_Parameter& A, const Virtual_Parameter& B){
  return SumDiff_Parameter(1,A,1,B);
}
SumDiff_Parameter operator-(const Virtual_Parameter& A, const Virtual_Parameter& B){
  return SumDiff_Parameter(1,A,-1,B);
}

Mult_Parameter operator*(const Virtual_Parameter& A, const Virtual_Parameter& B){
  return Mult_Parameter(A,B);
}
Inverse_Parameter operator/(double a, const Virtual_Parameter& A){
  return Inverse_Parameter(a,A);
}
Inverse_Parameter inv(const Virtual_Parameter& A){
  return Inverse_Parameter(1.,A);
}


///  OLD STUFF
// --------------------------------------------------------------

//
// ParameterCutFEM Parameter::kappa1   = ParameterCutFEM("kappa1", fun_kappa1);
// ParameterCutFEM Parameter::kappa2   = ParameterCutFEM("kappa2", fun_kappa2);
// ParameterCutFEM Parameter::lambdaG  = ParameterCutFEM("lambdaG", fun_lambdaG);
// ParameterCutFEM Parameter::lambdaB  = ParameterCutFEM("lambdaB", fun_lambdaB);

// ParameterCutFEM Parameter::lambdaB3 = ParameterCutFEM("lambdaB3", fun_lambdaB3);
// ParameterCutFEM Parameter::lambdaG3 = ParameterCutFEM("lambdaG3", fun_lambdaG3);
//

// static double fun_kappa1(int i, double hh, double meas, double measK, double meas_Cut) {
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
// static double fun_kappa2(int i, double hh, double meas, double measK, double meas_Cut) {
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
// static double fun_lambdaG(int i, double hh, double meas, double measK, double meas_Cut) {
//   double mu1 = (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(0) : 1;
//   double mu2 = (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(1) : 1;
//   double kappa1 = (CutFEM_ParameterList::find("kappa1"))?CutFEM_ParameterList::listParameter["kappa1"]->expression(i,hh,meas,meas_Cut) : 0.5;
//   double kappa2 = (CutFEM_ParameterList::find("kappa2"))?CutFEM_ParameterList::listParameter["kappa2"]->expression(i,hh,meas,meas_Cut) : 0.5;
//   double gamma = meas / hh;
//   double alphaK = hh;//meas_Cut/hh/hh;
//   // return 20 / hh / hh;//
//   return  (kappa1*mu1 + kappa2*mu2)*(100 + 10*gamma)/(alphaK);
// }
// static double fun_lambdaG3(int i, double hh, double meas, double measK, double meas_Cut) {
//   double mu1 = (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(0) : 1;
//   double mu2 = (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(1) : 1;
//   double kappa1 = (CutFEM_ParameterList::find("kappa1"))?CutFEM_ParameterList::listParameter["kappa1"]->expression(i,hh,meas,meas_Cut) : 0.5;
//   double kappa2 = (CutFEM_ParameterList::find("kappa2"))?CutFEM_ParameterList::listParameter["kappa2"]->expression(i,hh,meas,meas_Cut) : 0.5;
//   double gamma = meas / hh / hh;
//   double alphaK = hh/hh;//meas_Cut/hh/hh/hh;
//
//   // return 20 / hh / hh;//
//   return  (kappa1*mu1 + kappa2*mu2)*(100 + 10*gamma)/(alphaK);
// }
// static double fun_lambdaB(int i, double hh, double meas, double measK, double meas_Cut) {
//   double mu = (CutFEM_ParameterList::find("mu"))? CutFEM_ParameterList::listParameter["mu"]->expression(0) : 1;
//   double gammaK = meas/hh;
//   double alphaK = measK/hh/hh;
//   return mu/hh*(10+1e2*gammaK/alphaK);
//    // return 20 / hh / hh;
// }
// static double fun_lambdaB3(int i, double hh, double meas, double measK, double meas_Cut) {
//   double mu = (CutFEM_ParameterList::find("mu"))? CutFEM_ParameterList::listParameter["mu"]->expression(0) : 1;
//   double gammaK = meas/hh/hh;
//   double alphaK = measK/hh/hh/hh;
//
//   return mu/hh*(10+1e2*gammaK/alphaK);
//   // return 20 / hh / hh;
// }
//
//

//
//
