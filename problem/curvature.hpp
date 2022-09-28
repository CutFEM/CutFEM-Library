#ifndef CURVATURE_HPP_
#define CURVATURE_HPP_

#include "baseProblem.hpp"


template<typename M>
class Curvature {
  public :
  typedef M Mesh;
  typedef CutFESpace<Mesh> CutSpace;
  typedef typename CutSpace::FElement FElement;
  typedef typename FElement::Rd Rd;
  // typedef GenericMapping<Mesh> Mapping;
  typedef FunFEM<Mesh> Fun_h;
  typedef TestFunction<Rd::d> FunTest;

  static const int D = Rd::d;
  // const Mapping& mapping = DataMapping<Mesh>::Id;

  // using BaseProblem<M>::Vh;
  // using BaseProblemSurface<M>::interface;
  const CutSpace& Vh;
  const ActiveMesh<Mesh>& Kh;
  const Interface<Mesh>& interface;

  // GCurvature(const FESpace& vh, const Interface& inter, const Mapping& mapp = DataMapping<Mesh>::Id) :
  // BaseProblem<M>(vh), interface(inter), mapping(mapp) {
  //   this->solve();
  // }
  Curvature(const CutSpace& vh, const Interface<Mesh>* inter) :
  Vh(vh), Kh(vh.cutTh), interface(*inter) {

  }

  // GCurvature(const FESpace* vh, const Interface* inter, const Mapping& mapp = DataMapping<Mesh>::Id) :
  // BaseProblem<M>(*vh), interface(*inter), mapping(mapp) {
  //   this->solve();
  // }
public:
  Rn solve() {


    CutFEM<Mesh2> problem(Vh);
    FunTest H(Vh, D), v(Vh,D);

    // TestFunction<Rd::d> gradSu = gradS(-1*uu);
    Normal n;
    // TestFunction<Rd::d> gradun = grad(uu)*n;
    // double h = (*Vh)[0].T.lenEdge(0);
    // int deg = ((*Vh)[0].NbNode() > Rd::d+1) ? 2 : 1;
    // Projection Pg;
    //
    Rnm Id(D,D); Id = 0.;
    for(int i=0;i<D;++i) Id(i,i) = 1.;
    //
    // //a(u,v)_Gamma
    // double t0 = CPUtime();
    problem.addBilinear(
      (H,v) //+ (gradun,gradun)*1e-2
      , interface
      // , {}
      // , mapping
    );
    // l(v)_Omega
    problem.addLinear(
      -contractProduct(Id,gradS(v))
      , interface
      // , {}
      // , mapping
    );
    //
    // ListItemVF<Rd::d> Sh = (jump(gradun),jump(gradun))*1e-2;
    // this->addEdgeIntegral(Sh);
    problem.addBilinear(
      (jump(grad(H)*n),jump(grad(v)*n))*1e-2
      , Kh
      , INTEGRAL_INNER_FACET
    );


    // if(deg == 2) {
    //   TestFunction<Rd::d> grad2un = grad(gradun)*n;
    //   this->addEdgeIntegral(innerProduct(1e-2*h*h*jump(grad2un), jump(grad2un)));
    // }
    problem.solve();
    return problem.rhs_;
  }
};












#endif
