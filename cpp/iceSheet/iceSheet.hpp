#ifndef ICE_SHEET_HPP_
#define ICE_SHEET_HPP_

#include "../problem/problem.hpp"
#include "../util/cputime.h"

enum ElementIceSheet { Ice,  Water};


class IceSheet : public ShapeOfLinProblem, public Solver {

protected:

  typedef Mesh2 Mesh;
  typedef CutFESpace2 CutFESpace;
  typedef FESpace2 FESpace;
  typedef Interface2 Interface;
  typedef CutData2 CutData;

  typedef typename Mesh::Partition Partition;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::BorderElement BorderElement;
  typedef typename Mesh::Rd Rd;

  typedef typename FESpace::FElement FElement;
  typedef typename FElement::QF QF;
  typedef typename FElement::QFB QFB;
  typedef typename FElement::RdHat RdHat;
  typedef typename FElement::RdHatBord RdHatBord;

  // typedef GenericInterface<Mesh> Interface;
  typedef typename Interface::FaceIdx Face;
  typedef FEMFunction<FESpace> Fun;


public:
  const CutFESpace* Vh;

public:
  R sigma = 1.;
  R lambdaG;
  R Edu = 1e2, Edp = 1e2;

  const double yearinsec = 365.25*24*60*60;
  const R rhoi = 900.0 /(1.0e6*yearinsec*yearinsec);
  const R rhow = 1000.0/(1.0e6*yearinsec*yearinsec);
  const R A = (4.6416e-25)*yearinsec*1.0e18;
  const R mu_ = 1.0/pow(2.0*A,1.0/3) ;
  const R C = 7.624e6/(1.0e6*yearinsec) ;
  const R gravity = 9.8*yearinsec*yearinsec;

  IceSheet(const CutFESpace& vh);
  void init(const CutFESpace& vh);
  void solve(const Marker& );


private:
  void assembly(const Marker& );
  void assembly_full(const Marker& );
  void assembly_surface(const Marker&);
  void assembly_boundary(const Marker&);
  void stabilization(const Marker&);

public:
  R Pwater(R2 P) const {return (P.y<0)?-rhow*gravity*P.y : 0;}
  R2 frhs(R2 P) const {return R2(0, -gravity);}

  R mu() const {return  mu_;}
  R rho(enum ElementIceSheet a) const {return (a == Ice)? rhoi : rhow;}


};






#endif
