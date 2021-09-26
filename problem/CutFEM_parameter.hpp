#ifndef CutFEM_PARAMETER_HPP_
#define CutFEM_PARAMETER_HPP_
#include <iostream>
#include <map>
#include <string>
#include <cassert>
#include <math.h>
#include <vector>
#include "../common/R3.hpp"
struct CutFEM_Parameter;

class CutFEM_ParameterList {
public:
  static std::map<std::string,CutFEM_Parameter*> listParameter;
  static bool find(std::string s)  {
    return (listParameter.find(s) != listParameter.end());
  }
};

// Parameter defined for the FEM/CutFEM problem
class Parameter { public:
  static CutFEM_Parameter h;
  static CutFEM_Parameter invh;
  static CutFEM_Parameter kappa1;
  static CutFEM_Parameter kappa2;
  static CutFEM_Parameter lambdaG;
  static CutFEM_Parameter lambdaG3;
  static CutFEM_Parameter lambdaB;
  static CutFEM_Parameter lambdaB3;
  static CutFEM_Parameter meas;
};

// struct CutFEM_double {
//   double v1, v2;
//   CutFEM_double(const double a, const double b) : v1(a), v2(b) {}
// };



// defined parameter
class CutFEM_Parameter {
  public :
  typedef double (pfun)(int, double, double, double, double);
  std::string name;
  double val1, val2;
  double (*fun_expression)(int i, double hh, double meas, double measK, double meas_Cut) = nullptr;
public:
  CutFEM_Parameter(std::string n, pfun f) :
  name(n), fun_expression(f), val1(0), val2(0) {
    addToList();
  }
  CutFEM_Parameter(std::string n, double v1) :
  name(n), val1(v1), val2(v1), fun_expression(nullptr)  {
    addToList();
  }
  CutFEM_Parameter(std::string n, double v1, double v2) :
  name(n), val1(v1), val2(v2) , fun_expression(nullptr) {
    addToList();
  }

  void set(double v1, double v2 = 0) {
    val1 = v1;
    val2 = v2;
    fun_expression = nullptr;
  }


  double operator()(int i, double hh, double meas, double measK, double meas_Cut=0) const {
    if(fun_expression != nullptr) return fun_expression(i,hh,meas,measK, meas_Cut);
    else return (i==0)?val1 : val2;
  };
  virtual double expression(int i=0, double hh=0, double meas = 0, double measK=0, double meas_Cut=0) const {
    if(fun_expression != nullptr) return fun_expression(i,hh,meas,measK,meas_Cut);
    else return (i==0)?val1 : val2;
  };


  CutFEM_Parameter& operator=(const CutFEM_Parameter& F);
  std::string getName() const {return name;}

private:
  void addToList() ;
  bool checkAlreadyExist() const ;

  // CutFEM_Parameter(const CutFEM_Parameter&);
};

struct CutFEM_R2 {
  CutFEM_Parameter p1, p2;
  std::string name;
  CutFEM_R2(std::string n, const R2 a, const R2 b) : p1(n+"_x", a.x, b.x), p2(n+"_y", a.y, b.y) {}
  std::string getName(int i) const {
    return (i == 0)? p1.getName() : p2.getName();}
};

// Class for optimizing constant multiplication parameters
//-----------------------------------------------------------------------------------
class Pow_Par {
public:
  std::vector<std::string> list_name;
  std::vector<double> list_cst;
  Pow_Par(const Pow_Par& B) {
    for(int ll=0;ll<B.list_name.size();++ll){
      list_name.push_back(B.list_name[ll]);
    }
    for(int ll=0;ll<B.list_cst.size();++ll){
      list_cst.push_back(B.list_cst[ll]);
    }
  }
  Pow_Par(const Pow_Par& A, const Pow_Par& B) : Pow_Par(A) {
    for(int ll=0;ll<B.list_name.size();++ll){
      list_name.push_back(B.list_name[ll]);
    }
    for(int ll=0;ll<B.list_cst.size();++ll){
      list_cst.push_back(B.list_cst[ll]);
    }
  }
  Pow_Par(const CutFEM_Parameter& A,  const CutFEM_Parameter& B){
    list_name.push_back(A.getName());
    list_name.push_back(B.getName());
  }
  Pow_Par(const double a,  const CutFEM_Parameter& B){
    list_cst.push_back(a);
    list_name.push_back(B.getName());
  }
  Pow_Par(const double a,  const Pow_Par& B) : Pow_Par(B){
    list_cst.push_back(a);
  }
  Pow_Par(const CutFEM_Parameter& A,  const Pow_Par& B) : Pow_Par(B){
    list_name.push_back(A.getName());
  }

  Pow_Par(const CutFEM_Parameter& A, int n) {
    for(int i=0;i<n;++i) list_name.push_back(A.getName());
  }

};
Pow_Par operator*(const CutFEM_Parameter& A, const CutFEM_Parameter& B);
Pow_Par operator*(const double a, const CutFEM_Parameter& B);
Pow_Par operator*(const CutFEM_Parameter& B, const double a);
Pow_Par operator*(const Pow_Par& B, const double a);
Pow_Par operator*(const double a, const Pow_Par& B);
Pow_Par operator*(const CutFEM_Parameter& A, const Pow_Par& B);
Pow_Par operator*(const Pow_Par& B, const CutFEM_Parameter& A);
Pow_Par operator*(const Pow_Par& A, const Pow_Par& B);
Pow_Par pow(const CutFEM_Parameter& A, int n);
Pow_Par operator^(const CutFEM_Parameter& A, int n);

#endif
