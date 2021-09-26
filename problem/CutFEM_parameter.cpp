#include "CutFEM_parameter.hpp"
std::map<std::string,CutFEM_Parameter*> CutFEM_ParameterList::listParameter = std::map<std::string,CutFEM_Parameter*>();

static double fun_meas(int i, double hh, double meas, double measK, double meas_Cut) {return meas;}
static double fun_h(int i, double hh, double meas, double measK, double meas_Cut) {return hh;}
static double fun_invh(int i, double hh, double meas, double measK, double meas_Cut) {return 1./hh;}
static double fun_kappa1(int i, double hh, double meas, double measK, double meas_Cut) {
  if(CutFEM_ParameterList::find("mu") ) {
    double mu1 = CutFEM_ParameterList::listParameter["mu"]->expression(0);
    double mu2 = CutFEM_ParameterList::listParameter["mu"]->expression(1);
    double alpha1 = 1;//meas_Cut / measK;
    double alpha2 = 1;//(1 - alpha1);
    return mu2*alpha1/(alpha2*mu1+alpha1*mu2);
  }else{
    return 0.5;
  }
}
static double fun_kappa2(int i, double hh, double meas, double measK, double meas_Cut) {
  if(CutFEM_ParameterList::find("mu") ) {
    double mu1 = CutFEM_ParameterList::listParameter["mu"]->expression(0);
    double mu2 = CutFEM_ParameterList::listParameter["mu"]->expression(1);
    double alpha2 = 1;//meas_Cut / measK;
    double alpha1 = 1;//(1 - alpha2);
    return alpha2*mu1/(alpha2*mu1+alpha1*mu2);
  }else{
    return 0.5;
  }
}

static double fun_lambdaG(int i, double hh, double meas, double measK, double meas_Cut) {
  double mu1 = (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(0) : 1;
  double mu2 = (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(1) : 1;
  double kappa1 = (CutFEM_ParameterList::find("kappa1"))?CutFEM_ParameterList::listParameter["kappa1"]->expression(i,hh,meas,meas_Cut) : 0.5;
  double kappa2 = (CutFEM_ParameterList::find("kappa2"))?CutFEM_ParameterList::listParameter["kappa2"]->expression(i,hh,meas,meas_Cut) : 0.5;
  double gamma = meas / hh;
  double alphaK = hh;//meas_Cut/hh/hh;
  // return 20 / hh / hh;//
  return  (kappa1*mu1 + kappa2*mu2)*(100 + 10*gamma)/(alphaK);
}
static double fun_lambdaG3(int i, double hh, double meas, double measK, double meas_Cut) {
  double mu1 = (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(0) : 1;
  double mu2 = (CutFEM_ParameterList::find("mu"))?CutFEM_ParameterList::listParameter["mu"]->expression(1) : 1;
  double kappa1 = (CutFEM_ParameterList::find("kappa1"))?CutFEM_ParameterList::listParameter["kappa1"]->expression(i,hh,meas,meas_Cut) : 0.5;
  double kappa2 = (CutFEM_ParameterList::find("kappa2"))?CutFEM_ParameterList::listParameter["kappa2"]->expression(i,hh,meas,meas_Cut) : 0.5;
  double gamma = meas / hh / hh;
  double alphaK = hh/hh;//meas_Cut/hh/hh/hh;

  // return 20 / hh / hh;//
  return  (kappa1*mu1 + kappa2*mu2)*(100 + 10*gamma)/(alphaK);
}

static double fun_lambdaB(int i, double hh, double meas, double measK, double meas_Cut) {
  double mu = (CutFEM_ParameterList::find("mu"))? CutFEM_ParameterList::listParameter["mu"]->expression(0) : 1;
  double gammaK = meas/hh;
  double alphaK = measK/hh/hh;
  return mu/hh*(10+1e2*gammaK/alphaK);
   // return 20 / hh / hh;
}
static double fun_lambdaB3(int i, double hh, double meas, double measK, double meas_Cut) {
  double mu = (CutFEM_ParameterList::find("mu"))? CutFEM_ParameterList::listParameter["mu"]->expression(0) : 1;
  double gammaK = meas/hh/hh;
  double alphaK = measK/hh/hh/hh;

  return mu/hh*(10+1e2*gammaK/alphaK);
  // return 20 / hh / hh;
}

CutFEM_Parameter Parameter::h = CutFEM_Parameter("h", fun_h);
CutFEM_Parameter Parameter::invh = CutFEM_Parameter("invh", fun_invh);
CutFEM_Parameter Parameter::kappa1 = CutFEM_Parameter("kappa1", fun_kappa1);
CutFEM_Parameter Parameter::kappa2 = CutFEM_Parameter("kappa2", fun_kappa2);
CutFEM_Parameter Parameter::lambdaG = CutFEM_Parameter("lambdaG", fun_lambdaG);
CutFEM_Parameter Parameter::lambdaB = CutFEM_Parameter("lambdaB", fun_lambdaB);
CutFEM_Parameter Parameter::meas = CutFEM_Parameter("meas", fun_meas);
CutFEM_Parameter Parameter::lambdaB3 = CutFEM_Parameter("lambdaB3", fun_lambdaB3);
CutFEM_Parameter Parameter::lambdaG3 = CutFEM_Parameter("lambdaG3", fun_lambdaG3);




CutFEM_Parameter& CutFEM_Parameter::operator=(const CutFEM_Parameter& F) {
  name = F.name;
  val1 = F.val1;
  val2 = F.val2;
  fun_expression = F.fun_expression;
  return *this;
}

void CutFEM_Parameter::addToList() {
  if(checkAlreadyExist()){
    *CutFEM_ParameterList::listParameter[name] = *this;
    // std::cout << name << "already exist" << std::endl;
  }
  CutFEM_ParameterList::listParameter[name] = this;
}
bool CutFEM_Parameter::checkAlreadyExist() const {
  if(CutFEM_ParameterList::listParameter.find(name) != CutFEM_ParameterList::listParameter.end()){
    return true;
  }
  else return false;
}



Pow_Par operator*(const CutFEM_Parameter& A, const CutFEM_Parameter& B){
  return Pow_Par(A,B);
}
Pow_Par operator*(const double a, const CutFEM_Parameter& B) {
  return Pow_Par(a,B);
}
Pow_Par operator*(const CutFEM_Parameter& B, const double a){
  return Pow_Par(a,B);
}
Pow_Par operator*(const Pow_Par& B, const double a){
  return Pow_Par(a,B);
}
Pow_Par operator*(const double a, const Pow_Par& B) {
  return Pow_Par(a,B);
}
Pow_Par operator*(const CutFEM_Parameter& A, const Pow_Par& B) {
  return Pow_Par(A,B);
}
Pow_Par operator*(const Pow_Par& B, const CutFEM_Parameter& A) {
  return Pow_Par(A,B);
}
Pow_Par operator*(const Pow_Par& A, const Pow_Par& B) {
  return Pow_Par(A,B);
}
Pow_Par pow(const CutFEM_Parameter& A, int n) {
  return Pow_Par(A,n);
}
Pow_Par operator^(const CutFEM_Parameter& A, int n){
  return Pow_Par(A,n);
}
