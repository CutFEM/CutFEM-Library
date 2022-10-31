#ifndef CUTLEVEL_SET_HPP
#define CUTLEVEL_SET_HPP

#include "baseCutProblem.hpp"

template<typename M>
class GCutLevelSet : public BaseCutProblem<M> {//public ShapeOfLinProblem , public Solver  {
public:

  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename Mesh::Element Element;
  typedef typename FElement::Rd Rd;

  typedef GenericInterface<Mesh> Interface;
  typedef GLevelSet<M> LevelSet;
  typedef FunFEM<Mesh> Fun;

  list<int> label_strongBC;


  GCutLevelSet(const FESpace &);
  GCutLevelSet(const Fun& up, const Fun& Betap, const Fun& Beta, double dt);

  void solve(const Fun& up, const Fun& Betap, const Fun& Beta, double dt);
  void assembly(const Fun&, const Fun&, const Fun&, double);



  void setStrongBC(  std::list<int> lab) {
    label_strongBC = lab;
  }


  };

  template<typename M>
  GCutLevelSet<M>::GCutLevelSet(const FESpace & vh) : BaseCutProblem<M>(vh) {};

  template<typename M>
  GCutLevelSet<M>::GCutLevelSet(const Fun& up, const Fun& Betap, const Fun& Beta, double dt) :
  BaseCutProblem<M>(*up.Vh) {
    solve(up, Betap, Beta, dt);
  };


  template<typename M>
  void GCutLevelSet<M>::assembly(const Fun& up, const Fun& Betap, const Fun& Beta, double dt) {

    const int d = Rd::d;
    TestFunction<Rd::d> u(*this->Vh,1), v(*this->Vh,1);
    // Normal n;
    ExpressionFunFEM<M> Bx(Beta, 0, op_id), By(Beta,1, op_id);
    ExpressionFunFEM<M> Bpx(Betap, 0, op_id), Bpy(Betap,1, op_id);
    ExpressionFunFEM<M> dxU(up, 0, op_dx), dyU(up,0, op_dy);
    ExpressionFunFEM<M> Up(up, 0, op_id);
    TestFunction<Rd::d> Bu = Bx*dx(u) + By*dy(u);

    double h = sqrt(this->Vh->Th.mesure()/this->Vh->Th.nbVertices() );
    const KN_<double> dataBeta(Beta.getArray());
    double normB = dataBeta.norm() / dataBeta.size();
    double Tsd =  2. / (sqrt(dt*dt + normB * 1. / h / h));

    //a(u,v)_Gamma
    this->addBilinearFormDomain(
      innerProduct((1./dt)*u + 0.5*Bu , v + Tsd*Bu )
      // + innerProduct(0.01*grad(u), grad(v))
    );
    // l(v)_Omega
    this->addLinearFormDomain(
      innerProduct((1./dt)*Up - 0.5*(Bpx*dxU + Bpy*dyU) , v + Tsd*Bu )
    );

    // if(label_strongBC.size() > 0)
    //   this->addStrongBC(Up, label_strongBC);

    setBoundary(up);

  }


// class CutLevelSet2 : public GCutLevelSet<Mesh2> {
//
// public :
//   CutLevelSet2(const FESpace & vh, R f0(const Rd)) : CutLevelSet(vh,f0) {};
// };
//
//
//
//
// class  CutLevelSet3 : public GCutLevelSet<Mesh3> {
//
// public :
//
//   CutLevelSet3(const FESpace & vh, R f0(const Rd))  : CutLevelSet(vh,f0) {};
// };


#endif
