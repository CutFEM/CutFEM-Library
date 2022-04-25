#ifndef FINITE_ELEMENT_HPP
#define FINITE_ELEMENT_HPP

#include "GTypeOfFE_Sum.hpp"



class Lagrange2 :  public GTypeOfFESum<Mesh2>{
  typedef KN<const GTypeOfFE<Mesh2>*> FEarray;
  static const GTypeOfFE<Mesh2>* FE_[4][2];
public:
  Lagrange2(int k=1)
  : GTypeOfFESum<Mesh2>(FEarray(2,FE_[k])){}
};

class LagrangeDC2 :  public GTypeOfFESum<Mesh2>{
  typedef KN<const GTypeOfFE<Mesh2>*> FEarray;
  static const GTypeOfFE<Mesh2>* FE_[4][2];
public:
  LagrangeDC2(int k=1)
  : GTypeOfFESum<Mesh2>(FEarray(2,FE_[k])){}
};

class TaylorHood2 : public GTypeOfFESum<Mesh2>{
  typedef KN<const GTypeOfFE<Mesh2>*> FEarray;
  static const GTypeOfFE<Mesh2>* FE_[3];
public:
  TaylorHood2()
  : GTypeOfFESum<Mesh2>(FEarray(3,FE_)){}
};


class TaylorHood3 : public GTypeOfFESum<Mesh3>{
  typedef KN<const GTypeOfFE<Mesh3>*> FEarray;
  static const GTypeOfFE<Mesh3>* FE_[4];
public:
  TaylorHood3()
  : GTypeOfFESum<Mesh3>(FEarray(4,FE_)){}
};

template<typename Mesh>
class Lagrange3 :  public GTypeOfFESum<Mesh>{
  typedef KN<const GTypeOfFE<Mesh>*> FEarray;
  static const GTypeOfFE<Mesh>* FE_[3][3];
public:
  Lagrange3(int k=1)
  : GTypeOfFESum<Mesh>(FEarray(3,FE_[k])){}
};






#endif
