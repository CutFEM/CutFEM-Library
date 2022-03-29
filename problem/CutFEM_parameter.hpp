#ifndef CutFEM_PARAMETER_HPP_
#define CutFEM_PARAMETER_HPP_
#include <iostream>
#include <map>
#include <string>
#include <cassert>
#include <math.h>
#include <vector>
#include "../common/R3.hpp"
struct ParameterCutFEM;



class CutFEM_ParameterList {
public:
  static std::map<std::string,ParameterCutFEM*> listParameter;
  static bool find(std::string s)  {
    return (listParameter.find(s) != listParameter.end());
  }
};

// Parameter defined for the FEM/CutFEM problem
class Parameter { public:
  static ParameterCutFEM h;
  static ParameterCutFEM invh;
  static ParameterCutFEM meas;
  static ParameterCutFEM invmeas;
};

// VIRTUAL CLASS FOR HANDLING OPERATION ON PARAMETER
//------------------------------------------------------------------------------
class Virtual_Parameter {

public:
  virtual double evaluate(int domain, double h, double meas, double measK, double meas_Cut) const = 0;

};


template<int N>
struct CutFEM_Rd {
  typedef typename typeRd<N>::Rd Rd;

  std::vector<ParameterCutFEM> p;
  std::string name_;

  CutFEM_Rd(std::vector<Rd> a) : p(3) {
    int n = a.size();
    std::vector<double> v(n);
    for(int i=0;i<Rd::d;++i){
      for(int j=0;j<n;++j) v[j] = a[j][i];
      p[i].set(v);
    }
  }
  CutFEM_Rd(std::string l, std::vector<Rd> a) : p(3) {
    int n = a.size();
    std::string head[3] = {"_x", "_y", "_z"};
    std::vector<double> v(n);
    for(int i=0;i<Rd::d;++i){
      for(int j=0;j<n;++j) {v[j] = a[j][i];}
      p[i].set(l+head[i], v);
    }
  }

  const ParameterCutFEM& get_parameter(int i) const {
    return p[i];
  }
};

typedef CutFEM_Rd<2> CutFEM_R2;
typedef CutFEM_Rd<3> CutFEM_R3;

// THE CLASSIC CUTFEM_PARAMETER CLASS
//------------------------------------------------------------------------------
class ParameterCutFEM : public Virtual_Parameter {

  typedef double (pfun)(int, double, double, double, double);
  std::string name_;
  std::vector<double> val_;
  double (*fun_expression_)(int domain, double hh, double meas, double measK, double meas_Cut) = nullptr;

public:
  ParameterCutFEM(): name_(), val_(0), fun_expression_(nullptr){}
  ParameterCutFEM(const std::vector<double>& x) : name_(), val_(x), fun_expression_(nullptr) {}
  ParameterCutFEM(std::string n,const std::vector<double>& x) : name_(n), val_(x), fun_expression_(nullptr) {
    addToList();
  }
  ParameterCutFEM(pfun f) : name_(), val_(0), fun_expression_(f){}
  ParameterCutFEM(std::string n, pfun f) : name_(n), fun_expression_(f), val_(0) {
    addToList();
  }

  ParameterCutFEM(double v1) : name_(), val_(v1), fun_expression_(nullptr)  {
    val_.push_back(v1);
  }
  ParameterCutFEM(std::string n, double v1) : name_(n), val_(v1), fun_expression_(nullptr)  {
    val_.push_back(v1);
    addToList();
  }
  ParameterCutFEM(double v1, double v2) : ParameterCutFEM(v1)  {
    val_.push_back(v2);
  }
  ParameterCutFEM(std::string n, double v1, double v2) : ParameterCutFEM(n,v1)  {
    val_.push_back(v2);
  }
  ParameterCutFEM(const ParameterCutFEM& F):name_(F.name_), val_(F.val_), fun_expression_(F.fun_expression_) {}

  void set(std::string n, const std::vector<double>& x) {
    name_=n;
    val_ = x;
    fun_expression_ = nullptr;
    addToList();
  }
  void set(const std::vector<double>& x) {
    name_= "";
    val_ = x;
    fun_expression_ = nullptr;
  }

  double evaluate(int domain, double h, double meas, double measK, double meas_Cut) const {
    if(fun_expression_ != nullptr) return fun_expression_(domain ,h, meas, measK, meas_Cut);
    else {assert(domain >=0 && domain< val_.size()); return val_[domain];}
  };


private:
  void addToList() ;
  // bool checkAlreadyExist() const ;
  ParameterCutFEM& operator=(const ParameterCutFEM& F) ;
};


// CLASS TO BE ABLE MULTIPLY BY A CONSTANT
//------------------------------------------------------------------------------
class Mul_Cst_Parameter : public Virtual_Parameter{

  double cst_;
  const Virtual_Parameter* parameter_;

public:

  Mul_Cst_Parameter(const double a,  const Virtual_Parameter& B)
  : cst_(a), parameter_(&B) { }

  double evaluate(int domain, double h, double meas, double measK, double meas_cut) const {
    double val = parameter_->evaluate(domain,h,meas,measK,meas_cut);
    return cst_*val;
  }
};
Mul_Cst_Parameter operator*(const double a, const Virtual_Parameter& B);
Mul_Cst_Parameter operator*(const Virtual_Parameter& B, const double a);


// CLASS TO USE POWER OF A PARAMETER
//------------------------------------------------------------------------------
class Pow_Parameter : public Virtual_Parameter{

  // std::vector<double> list_cst_;
  const Virtual_Parameter* parameter_;
  int n;

public:
  Pow_Parameter(const Virtual_Parameter& A, int nn): parameter_(&A) , n(nn) {}
  double evaluate(int domain, double h, double meas, double measK, double meas_cut) const {
    double val = parameter_->evaluate(domain,h,meas,measK,meas_cut);
    return pow(val, n);
  }
};
Pow_Parameter pow(const Virtual_Parameter& A, int n) ;
Pow_Parameter operator^(const Virtual_Parameter& A, int n) ;


// TO HANDLE LINEAR COMBINAISON OF PARAMETER
//------------------------------------------------------------------------------
class SumDiff_Parameter : public Virtual_Parameter{

  double a_,b_;
  const Virtual_Parameter* A_;
  const Virtual_Parameter* B_;

public:

  SumDiff_Parameter(double a, const Virtual_Parameter& A,  double b, const Virtual_Parameter& B) :
  a_(a), A_(&A), b_(b), B_(&B){}


  double evaluate(int domain, double h, double meas, double measK, double meas_cut) const {
    double val1 = A_->evaluate(domain,h,meas,measK,meas_cut);
    double val2 = B_->evaluate(domain,h,meas,measK,meas_cut);

    return a_*val1 + b_*val2;
  }
};
SumDiff_Parameter operator+(const Virtual_Parameter& A, const Virtual_Parameter& B);
SumDiff_Parameter operator-(const Virtual_Parameter& A, const Virtual_Parameter& B);

// TO MULTIPLY PARAMETERS
//------------------------------------------------------------------------------
class Mult_Parameter : public Virtual_Parameter{

  const Virtual_Parameter* A_;
  const Virtual_Parameter* B_;

public:

  Mult_Parameter(const Virtual_Parameter& A, const Virtual_Parameter& B) : A_(&A),B_(&B){}


  double evaluate(int domain, double h, double meas, double measK, double meas_cut) const {
    double val1 = A_->evaluate(domain,h,meas,measK,meas_cut);
    double val2 = B_->evaluate(domain,h,meas,measK,meas_cut);
    return val1*val2;
  }
};
Mult_Parameter operator*(const Virtual_Parameter& A, const ParameterCutFEM& B);


// TO INVERSE PARAMETERS
//------------------------------------------------------------------------------
class Inverse_Parameter : public Virtual_Parameter{

  double a_;
  const Virtual_Parameter* A_;

public:

  Inverse_Parameter(double a, const Virtual_Parameter& A) : a_(a), A_(&A){}


  double evaluate(int domain, double h, double meas, double measK, double meas_cut) const {
    double val1 = A_->evaluate(domain,h,meas,measK,meas_cut);
    return a_/val1;
  }
};
Inverse_Parameter operator/(double a, const Virtual_Parameter& A);
Inverse_Parameter inv(const Virtual_Parameter& A);




#endif
