#ifndef CURVATURE_HPP_
#define CURVATURE_HPP_

#include "baseProblem.hpp"


template<typename M>
class GCurvature : public BaseProblem<M> {
  public :
  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::Rd Rd;
  typedef GenericInterface<Mesh> Interface;
  typedef GenericMapping<Mesh> Mapping;
  typedef FunFEM<Mesh> Fun_h;

  const Mapping& mapping = DataMapping<Mesh>::Id;

  using BaseProblem<M>::Vh;
  // using BaseProblemSurface<M>::interface;
  const Interface& interface;

  GCurvature(const FESpace& vh, const Interface& inter, const Mapping& mapp = DataMapping<Mesh>::Id) :
  BaseProblem<M>(vh), interface(inter), mapping(mapp) {
    this->solve();
  }
  GCurvature(const FESpace& vh, const Interface* inter, const Mapping& mapp = DataMapping<Mesh>::Id) :
  BaseProblem<M>(vh), interface(*inter), mapping(mapp) {
    this->solve();
  }
  GCurvature(const FESpace* vh, const Interface* inter, const Mapping& mapp = DataMapping<Mesh>::Id) :
  BaseProblem<M>(*vh), interface(*inter), mapping(mapp) {
    this->solve();
  }
private:
  void assembly() {

    const int d = Rd::d;

    TestFunction<Rd::d> uu(*Vh, d);
    TestFunction<Rd::d> gradSu = gradS(-1*uu);
    Normal n;
    TestFunction<Rd::d> gradun = grad(uu)*n;
    double h = (*Vh)[0].T.lenEdge(0);
    int deg = ((*Vh)[0].NbNode() > Rd::d+1) ? 2 : 1;
    Projection Pg;

    Rnm Id(d,d); Id = 0.;
    for(int i=0;i<d;++i) Id(i,i) = 1.;

    // std::cout << contractProduct(Id,gradSu) << std::endl;
    // std::cout << grad(uu)  << std::endl;
    // std::cout << gradun    << std::endl;
    // std::cout << gradS(uu) << std::endl;
    // std::cout << contractProduct(Pg,gradS(uu)) << std::endl;
    // getchar();

    //a(u,v)_Gamma
    double t0 = CPUtime();
    this->addBilinear(
      (uu,uu) //+ (gradun,gradun)*1e-2
      , interface
      , mapping
    );
    // l(v)_Omega
    this->addLinear(
      contractProduct(Id,gradSu)
      , interface
      , mapping
    );

    ListItemVF<Rd::d> Sh = (jump(gradun),jump(gradun))*1e-2;
    this->addEdgeIntegral(Sh);

    if(deg == 2) {
      TestFunction<Rd::d> grad2un = grad(gradun)*n;
      this->addEdgeIntegral(innerProduct(1e-2*h*h*jump(grad2un), jump(grad2un)));
    }
  }
};












#endif
