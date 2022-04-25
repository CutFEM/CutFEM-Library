#ifndef CutFEM_PARAMETER_HPP_
#define CutFEM_PARAMETER_HPP_
#include <iostream>
#include <map>
#include <string>
#include <cassert>
#include <math.h>
#include <vector>
#include "../common/R3.hpp"




// VIRTUAL CLASS FOR HANDLING OPERATION ON PARAMETER
//------------------------------------------------------------------------------
class VirtualParameter {

public:
  virtual double evaluate(int domain, double h, double meas, double measK, double meas_Cut) const = 0;

};

struct MeshParameter : public VirtualParameter {
  int a1, a2, a3, a4;
  MeshParameter(int b1, int b2, int b3, int b4) : a1(b1), a2(b2), a3(b3), a4(b4) {}
  double evaluate(int domain, double h, double meas, double measK, double meas_Cut) const {
    return a1*h + a2*meas + a3*measK + a4*meas_Cut;
  }
};
// Parameter defined for the FEM/CutFEM problem
class Parameter { public:
  static const MeshParameter h;
  static const MeshParameter measureIntegral;
  static const MeshParameter measureElement;
  static const MeshParameter measureCut;
};



// THE CLASSIC CUTFEM_PARAMETER CLASS
//------------------------------------------------------------------------------
class CutFEMParameter : public VirtualParameter {

  typedef double (pfun)(int, double, double, double, double);
  std::vector<double> val_;
public:
  CutFEMParameter():val_(0){}
  CutFEMParameter(const std::vector<double>& x) : val_(x){}
  CutFEMParameter(double v1) {
    val_.push_back(v1);
  }
  CutFEMParameter(double v1, double v2) : CutFEMParameter(v1)  {
    val_.push_back(v2);
  }
  CutFEMParameter(const CutFEMParameter& F) : val_(F.val_){}

  void set(const std::vector<double>& x) {
    val_ = x;
  }

  double evaluate(int domain, double h, double meas, double measK, double meas_Cut) const {
    assert(domain >=0 && domain< val_.size());
    return val_[domain];
  }

private:
  CutFEMParameter& operator=(const CutFEMParameter& F) ;
};


template<int N>
struct CutFEM_Rd {
  typedef typename typeRd<N>::Rd Rd;

  std::vector<CutFEMParameter> p;
  std::string name_;

  CutFEM_Rd(std::vector<Rd> a) : p(3) {
    int n = a.size();
    std::vector<double> v(n);
    for(int i=0;i<Rd::d;++i){
      for(int j=0;j<n;++j) v[j] = a[j][i];
      p[i].set(v);
    }
  }

  const CutFEMParameter& get_parameter(int i) const {
    return p[i];
  }
};

typedef CutFEM_Rd<2> CutFEM_R2;
typedef CutFEM_Rd<3> CutFEM_R3;

// CLASS TO BE ABLE MULTIPLY BY A CONSTANT
//------------------------------------------------------------------------------
class Mul_Cst_Parameter : public VirtualParameter{

  double cst_;
  const VirtualParameter* parameter_;

public:

  Mul_Cst_Parameter(const double a,  const VirtualParameter& B)
  : cst_(a), parameter_(&B) { }

  double evaluate(int domain, double h, double meas, double measK, double meas_cut) const {
    double val = parameter_->evaluate(domain,h,meas,measK,meas_cut);
    return cst_*val;
  }
};
Mul_Cst_Parameter operator*(const double a, const VirtualParameter& B);
Mul_Cst_Parameter operator*(const VirtualParameter& B, const double a);


// CLASS TO USE POWER OF A PARAMETER
//------------------------------------------------------------------------------
class Pow_Parameter : public VirtualParameter{

  // std::vector<double> list_cst_;
  const VirtualParameter* parameter_;
  int n;

public:
  Pow_Parameter(const VirtualParameter& A, int nn): parameter_(&A) , n(nn) {}
  double evaluate(int domain, double h, double meas, double measK, double meas_cut) const {
    double val = parameter_->evaluate(domain,h,meas,measK,meas_cut);
    return pow(val, n);
  }
};
Pow_Parameter pow(const VirtualParameter& A, int n) ;
Pow_Parameter operator^(const VirtualParameter& A, int n) ;


// TO HANDLE LINEAR COMBINAISON OF PARAMETER
//------------------------------------------------------------------------------
class SumDiff_Parameter : public VirtualParameter{

  double a_,b_;
  const VirtualParameter* A_;
  const VirtualParameter* B_;

public:

  SumDiff_Parameter(double a, const VirtualParameter& A,  double b, const VirtualParameter& B) :
  a_(a), A_(&A), b_(b), B_(&B){}


  double evaluate(int domain, double h, double meas, double measK, double meas_cut) const {
    double val1 = A_->evaluate(domain,h,meas,measK,meas_cut);
    double val2 = B_->evaluate(domain,h,meas,measK,meas_cut);

    return a_*val1 + b_*val2;
  }
};
SumDiff_Parameter operator+(const VirtualParameter& A, const VirtualParameter& B);
SumDiff_Parameter operator-(const VirtualParameter& A, const VirtualParameter& B);

// TO MULTIPLY PARAMETERS
//------------------------------------------------------------------------------
class Mult_Parameter : public VirtualParameter{

  const VirtualParameter* A_;
  const VirtualParameter* B_;

public:

  Mult_Parameter(const VirtualParameter& A, const VirtualParameter& B) : A_(&A),B_(&B){}


  double evaluate(int domain, double h, double meas, double measK, double meas_cut) const {
    double val1 = A_->evaluate(domain,h,meas,measK,meas_cut);
    double val2 = B_->evaluate(domain,h,meas,measK,meas_cut);
    return val1*val2;
  }
};
Mult_Parameter operator*(const VirtualParameter& A, const CutFEMParameter& B);


// TO INVERSE PARAMETERS
//------------------------------------------------------------------------------
class Inverse_Parameter : public VirtualParameter{

  double a_;
  const VirtualParameter* A_;

public:

  Inverse_Parameter(double a, const VirtualParameter& A) : a_(a), A_(&A){}


  double evaluate(int domain, double h, double meas, double measK, double meas_cut) const {
    double val1 = A_->evaluate(domain,h,meas,measK,meas_cut);
    return a_/val1;
  }
};
Inverse_Parameter operator/(double a, const VirtualParameter& A);
Inverse_Parameter inv(const VirtualParameter& A);




#endif
