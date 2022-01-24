#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "../util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif

#include "finiteElement.hpp"
#include "baseCutProblem.hpp"

#include "../num/gnuplot.hpp"


#define PROBLEM_STOKES_DIV

#ifdef PROBLEM_STOKES_DIV

/* REMAINING ISSUES & CURRENT KNOWLEDGE:
= Normal must point outward from domain
= Normal must point outward from + element in jump term
= Adding terms on the border should not be done when sol is zero there
= The BC is ignored when sigma >> penparam
= Big sigma seems to steer the velocity SE to incorrect sols - eg sigma*100 gives vel/100
== Get correct sol (u,p) for sigma=1, pen>1e3/meshsize, with BC normal pen & BC tangent natural
== Same sol if include all natural BC (ie avg and jumps)
== Same sol if remove BC tangent natural, but then the pressure gets the NE SW corner issue
" For example, in the H(div) method, the normal component of the velocity
is set as an essential boundary condition and is strongly enforced, but the tangential
component of the velocity is treated as a natural boundary condition and is weakly
imposed. "


- Velocity is incorrect (weird symmetry along y=-x+c)

+ Pressure gets corrected by adding interior edge terms also on the boundary
*/
namespace Erik_Data_StokesDiv {

  // R fun_rhs(const R2 P, int i) {
  //   return 0;
  // }
  // R fun_exact(const R2 P, int i) {
  //   if(i==0) return 0;
  //   else if(i==1) return P.x;
  //   else return 0;
  // }

  R fun_rhs(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    if(i==0) return -4*(6*x*x - 6*x + 1)*y*(2*y*y - 3*y + 1) - 12*(x-1)*(x-1)*x*x*(2*y - 1);
    else if(i==1) return 4*x*(2*x*x - 3 *x + 1)*(6*y*y - 6*y + 1) + 12*(2*x - 1)*(y - 1)*(y-1)*y*y;
    else return 0;
  }
  R fun_exact(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    if(i==0) return 2*x*y*(x-1)*(y-1)*x*(1-x)*(2*y-1);
    else if(i==1) return 2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x-1);
    else return 0;
  }
  // R fun_rhs(const R2 P, int i) {
  //   R mu=1;
  //   R x = P.x;
  //   R y = P.y;
  //   if(i==0) return -2*mu*(pow(x,4)*(6*y-3)+pow(x,3)*(6-12*y)+3*pow(x,2)*(4*pow(y,3)-6*pow(y,2)+4*y-1)-6*x*y*(2*pow(y,2)-3*y+1)+y*(2*pow(y,2)-3*y+1)) + 10*(3*pow(x-0.5,2)*pow(y,2)-3*pow(1-x,2)*pow(y-0.5,3));
  //   else if(i==1) return 2*mu*(2*pow(x,3)*(6*pow(y,2)-6*y+1)-3*pow(x,2)*(6*pow(y,2)-6*y+1)+x*(6*pow(y,4)-12*pow(y,3)+12*pow(y,2)-6*y+1)-3*pow(y-1,2)*pow(y,2)) + 10*(2*pow(x-0.5,3)*y+3*pow(1-x,3)*pow(y-0.5,2));
  //   else return 0;
  // }
  // R fun_exact(const R2 P, int i) {
  //   R x = P.x;
  //   R y = P.y;
  //   if(i==0) return pow(x,2)*pow(1-x,2)*y*(1-y)*(1-2*y);
  //   else if(i==1) return (-x)*(1-x)*(1-2*x)*pow(y,2)*pow(1-y,2);
  //   else return 10*(pow(x-0.5,3)*pow(y,2)+pow(1-x,3)*pow(y-0.5,3));
  // }
}

using namespace Erik_Data_StokesDiv;

int main(int argc, char** argv )
{
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  const double cpubegin = CPUtime();
  MPIcf cfMPI(argc,argv);

  int nx = 10;
  int ny = 10;
  // int d = 2;

  vector<double> ul2,pl2,h, convu, convp;

  for(int i=0;i<4;++i) { // i<3

    std::cout << "\n ------------------------------------- " << std::endl;
    Mesh2 Th(nx, ny, 0., 0., 1., 1.);
    // Mesh2 Th("../mesh/RTmesh100.msh");


    KN<const GTypeOfFE<Mesh2>* > arrayFE(2);
    arrayFE(0) = &DataFE<Mesh2>::RT1;
    arrayFE(1) = &DataFE<Mesh2>::P1dc;
    GTypeOfFESum<Mesh2> mixedFE(arrayFE);
    FESpace2 mixedSpace(Th, mixedFE); // (Vh,Qh)
    mixedSpace.info();

    FESpace2 Vh(Th, DataFE<Mesh2>::RT1); // REMOVE THESE TWO LATER
    FESpace2 Qh(Th, DataFE<Mesh2>::P0); // FOR MIXEDSPACE

    const R meshSize = Th[0].lenEdge(0);
    const R penaltyParam = 1e1/meshSize;
    const R sigma = 10; // HIGHER 1e2,1e3,etc WORSENS THE SOL [??]

    Fun_h fh(Vh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
    Fun_h gh(Vh, fun_exact);
    // Fun_h gh(Qh, fun_div); // THIS IS IF div u != 0

    FEM<Mesh2> stokesDiv(mixedSpace);

    Normal n;
    Tangent t;
    /* Syntax:
    FunTest (fem space, #components, place in space)
    */
    FunTest u(mixedSpace,2,0), p(mixedSpace,1,2), v(mixedSpace,2,0), q(mixedSpace,1,2);

    R mu = 1;
    const CutFEM_Parameter& invh(Parameter::invh);
    const CutFEM_Parameter& invlEdge(Parameter::invmeas);

    // a_h(u,v)
    stokesDiv.addBilinear(
      contractProduct(mu*grad(u),grad(v))
      - innerProduct(p,div(v))
      + innerProduct(div(u),q)
    );
    stokesDiv.addBilinear(
      - innerProduct(average(grad(u*t)*n,0.5,0.5), jump(v*t))
      + innerProduct(jump(u*t), average(grad(v*t)*n,0.5,0.5))
      + innerProduct(invlEdge*(sigma*jump(u*t)), jump(v*t))
      , innerEdge
    );

    stokesDiv.addBilinear(
      - innerProduct(grad(u*t)*n, v*t)
      + innerProduct(u*t, grad(v*t)*n)
      + innerProduct(invlEdge*sigma*(u*t), v*t)
      + innerProduct(p, v*n)
      - innerProduct(u*n, q)
      , boundary
    );
    // stokesDiv.addBilinearFormBorder(
    //   innerProduct(penaltyParam*u.t()*n,v.t()*n)
    // );
    // stokesDiv.addBilinearFormBorder(
    //   innerProduct(penaltyParam*u, v)
    // );
    // b(v,p), b(u,q)
    // stokesDiv.addBilinear(
    //   - innerProduct(p,div(v))
    //   - innerProduct(div(u),q)
    // );


    // stokesDiv.addLinear(
    //     innerProduct(gh.expression(2), invh*v)
    //   - innerProduct(gh.expression(2), grad(v)*n)
    //   - innerProduct(gh.expression(2), p*n)
    //   , boundary
    // );

    // l(v)_Omega
    stokesDiv.addLinear(
      innerProduct(fh.expression(2),v)
    );

    // Sets uniqueness of the pressure
    stokesDiv.addLagrangeMultiplier(
      innerProduct(1.,p), 0.
    );

    // stokesDiv.imposeStrongBC(
    //   +0
    // );

    stokesDiv.solve();

    Fun_h sol(mixedSpace, stokesDiv.rhs);

    Paraview2 writerS(mixedSpace, "stokesDivN"+to_string(i)+".vtk");
    writerS.add(sol, "uh", 0, 2);
    writerS.add(sol, "ph", 2, 1);

    Rn solExVec(mixedSpace.NbDoF());
    interpolate(mixedSpace, solExVec, fun_exact);
    R errU = L2norm(sol,fun_exact,0,2);
    R errP = L2norm(sol,fun_exact,2,1);


    // solExVec -= stokesDiv.rhs;
    for(int i=0;i<solExVec.size();++i){
      solExVec(i) = fabs(solExVec(i)-stokesDiv.rhs(i));
    }

    Fun_h solEx(mixedSpace, solExVec);
    writerS.add(solEx,"uh_err", 0, 2);

    // writerS.add(solEx,"uh_ex", 0, 2);
    // writerS.add(solEx,"ph_ex", 2, 1);



    ul2.push_back(errU);
    pl2.push_back(errP);
    h.push_back(1./nx);
    if(i==0) {convu.push_back(0); convp.push_back(0);}
    else {
      convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
      convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
    }

    nx *= 2;
    ny *= 2;
  }
  std::cout << "\n\n h \t err u \t\t convU \t\t err p \t\t convP" << std::endl;
  for(int i=0;i<ul2.size();++i) {
    std::cout << h[i] << "\t" << ul2[i] << "\t" << convu[i] << "\t\t"
  	      << pl2[i] << "\t" << convp[i] << std::endl;
  }

}
#endif
