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
#include "baseProblem.hpp"
#include "paraview.hpp"

// #include "../num/gnuplot.hpp"


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
     if(i==0) return -2*y;
     else return 2*x;
    // if(i==0) return -4*(6*x*x - 6*x + 1)*y*(2*y*y - 3*y + 1) - 12*(x-1)*(x-1)*x*x*(2*y - 1);
    // else if(i==1) return 4*x*(2*x*x - 3 *x + 1)*(6*y*y - 6*y + 1) + 12*(2*x - 1)*(y - 1)*(y-1)*y*y;
    // else return 0;
  }
  // R fun_exact(const R2 P, int i) {
  //   R x = P.x;
  //   R y = P.y;
  //   if(i==0) return 2*x*y*(x-1)*(y-1)*x*(1-x)*(2*y-1);
  //   else if(i==1) return 2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x-1);
  //   else return 0;
  // }
  R fun_exact_u(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    // if(i==0) return      2*x*y*(x-1)*(y-1)*x*(1-x)*(2*y-1);
    // else if(i==1) return 2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x-1);
    if(i==0)      return  x*x*y;
    else if(i==1) return -x*y*y;
  }
  R fun_exact_p(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    return 0;
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
  typedef Mesh2 Mesh;
  typedef ActiveMeshT2 CutMesh;
  typedef FESpace2   Space;
  typedef CutFESpaceT2 CutSpace;

  const double cpubegin = CPUtime();
  MPIcf cfMPI(argc,argv);

  int nx = 10;
  int ny = 10;
  // int d = 2;

  vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

  for(int i=0;i<4;++i) { // i<3

    std::cout << "\n ------------------------------------- " << std::endl;
    Mesh Th(nx, ny, 0., 0., 1., 1.);


    Lagrange2 FEvelocity(2);
    Space Uh(Th, FEvelocity);
    // Space Qh(Th, DataFE<Mesh>::P1);

    Space Vh(Th, DataFE<Mesh2>::RT1); // REMOVE THESE TWO LATER
    Space Qh(Th, DataFE<Mesh2>::P1dc); // FOR MIXEDSPACE

    const R hi = Th[0].lenEdge(0);
    const R penaltyParam = 1e0/hi;
    const R sigma = 1; // HIGHER 1e2,1e3,etc WORSENS THE SOL [??]

    Fun_h fh(Uh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
    Fun_h gh(Uh, fun_exact_u);
    // Fun_h gh(Qh, fun_div); // THIS IS IF div u != 0

    FEM<Mesh2> stokes(Vh); stokes.add(Qh);

    Normal n;
    Tangent t;
    /* Syntax:
    FunTest (fem space, #components, place in space)
    */
    FunTest u(Vh,2,0), p(Qh,1,0), v(Vh,2,0), q(Qh,1,0);

    R mu = 1;
    // const MeshParameter& h(Parameter::h);
    // const MeshParameter& invlEdge(Parameter::invmeas);

    // a_h(u,v)
    stokes.addBilinear(
      contractProduct(mu*grad(u),grad(v))
      - innerProduct(p,div(v))
      + innerProduct(div(u),q)
      , Th
    );

    stokes.addBilinear(
      - innerProduct(average(grad(u*t)*n,0.5,0.5), jump(v*t))
      + innerProduct(jump(u*t), average(grad(v*t)*n,0.5,0.5))
      + innerProduct(1./hi*(sigma*jump(u*t)), jump(v*t))
      , Th
      , innerFacet
    );

    // stokes.addBilinear(
    //   - innerProduct(grad(u*t)*n, v*t) // natural
    //   - innerProduct(u, grad(v)*n)
    //   + innerProduct(p, v*n)
    //   - innerProduct(u*n, q)
    //   + innerProduct(penaltyParam*u*t, v*t)
    //   , Th
    //   , boundary
    // );
    stokes.addBilinear(
      // + innerProduct(penaltyParam*u*t, v*t)
      // - innerProduct(u*t, grad(v*t)*n)
      - innerProduct(grad(u)*n, v)  //natural
      - innerProduct(u, grad(v)*n)
      + innerProduct(p, v*n)  // natural
      - innerProduct(u*n, q)  // essential
      + innerProduct(penaltyParam*u, v)
      , Th
      , boundary
    );
    stokes.addLinear(
      // + innerProduct(gh*t, penaltyParam*v*t)
      // - innerProduct(gh*t, grad(v*t)*n)
      + innerProduct(gh.expression(2), penaltyParam*v)
      - innerProduct(gh.expression(2), grad(v)*n)
      - innerProduct(gh.expression(2), q*n)

      , Th
      , boundary
    );



    // stokesDiv.addBilinearFormBorder(
    //   innerProduct(penaltyParam*u.t()*n,v.t()*n)
    // );
    // stokesDiv.addBilinear(
    //   + innerProduct(penaltyParam*u, v)
    //
    // );
    // b(v,p), b(u,q)
    // stokesDiv.addBilinear(
    //   - innerProduct(p,div(v))
    //   - innerProduct(div(u),q)
    // );




    // l(v)_Omega
    stokes.addLinear(
      innerProduct(fh.expression(2),v)
      , Th
    );

    // Sets uniqueness of the pressure
    stokes.addLagrangeMultiplier(
      innerProduct(1.,p), 0., Th
    );

    // stokesDiv.imposeStrongBC(
    //   +0
    // );

    stokes.solve();

    // EXTRACT SOLUTION
    int idx0_s = Vh.get_nb_dof();
    Rn_ data_uh = stokes.rhs_(SubArray(Vh.get_nb_dof(),0));
    Rn_ data_ph = stokes.rhs_(SubArray(Qh.get_nb_dof(),idx0_s));
    Fun_h uh(Vh, data_uh);
    Fun_h ph(Qh, data_ph);
    ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
    ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);


    // {
    //   Fun_h soluErr(Vh, fun_exact_u);
    //   Fun_h soluh(Vh, fun_exact_u);
    //   soluErr.v -= uh.v;
    //   soluErr.v.map(fabs);
    //   // Fun_h divSolh(Wh, fun_div);
    //   // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);
    //
    //   Paraview<Mesh> writer(Th, "stokes_"+to_string(i)+".vtk");
    //   writer.add(uh, "velocity" , 0, 2);
    //   writer.add(ph, "pressure" , 0, 1);
    //   writer.add(femSol_0dx+femSol_1dy, "divergence");
    //   writer.add(soluh, "velocityExact" , 0, 2);
    //   writer.add(soluErr, "velocityError" , 0, 2);
    //   // writer.add(solh, "velocityError" , 0, 2);
    //
    //   // writer.add(fabs(femDiv, "divergenceError");
    // }

    // Rn solExVec(mixedSpace.NbDoF());
    // interpolate(mixedSpace, solExVec, fun_exact);
    // R errU = L2norm(sol,fun_exact,0,2);
    // R errP = L2norm(sol,fun_exact,2,1);

    R errU      = L2norm(uh,fun_exact_u,0,2);
    R errP      = L2norm(ph,fun_exact_p,0,1);
    R errDiv    = L2norm (femSol_0dx+femSol_1dy,Th);
    R maxErrDiv = maxNorm(femSol_0dx+femSol_1dy,Th);


    // // solExVec -= stokesDiv.rhs;
    // for(int i=0;i<solExVec.size();++i){
    //   solExVec(i) = fabs(solExVec(i)-stokesDiv.rhs(i));
    // }
    //
    // Fun_h solEx(mixedSpace, solExVec);
    // writerS.add(solEx,"uh_err", 0, 2);

    // writerS.add(solEx,"uh_ex", 0, 2);
    // writerS.add(solEx,"ph_ex", 2, 1);



    ul2.push_back(errU);
    pl2.push_back(errP);
    divl2.push_back(errDiv);
    divmax.push_back(maxErrDiv);
    h.push_back(1./nx);
    if(i==0) {convu.push_back(0); convp.push_back(0);}
    else {
      convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
      convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
    }

    nx *= 2;
    ny *= 2;
  }
  std::cout << "\n" << std::left
  << std::setw(10) << std::setfill(' ') << "h"
  << std::setw(15) << std::setfill(' ') << "err_p u"
  // << std::setw(15) << std::setfill(' ') << "conv p"
  << std::setw(15) << std::setfill(' ') << "err u"
  // << std::setw(15) << std::setfill(' ') << "conv u"
  << std::setw(15) << std::setfill(' ') << "err divu"
  // << std::setw(15) << std::setfill(' ') << "conv divu"
  // << std::setw(15) << std::setfill(' ') << "err_new divu"
  // << std::setw(15) << std::setfill(' ') << "convLoc divu"
  << std::setw(15) << std::setfill(' ') << "err maxdivu"
  // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
  << "\n" << std::endl;
  for(int i=0;i<h.size();++i) {
    std::cout << std::left
    << std::setw(10) << std::setfill(' ') << h[i]
    << std::setw(15) << std::setfill(' ') << pl2[i]
    // << std::setw(15) << std::setfill(' ') << convpPr[i]
    << std::setw(15) << std::setfill(' ') << ul2[i]
    // << std::setw(15) << std::setfill(' ') << convuPr[i]
    << std::setw(15) << std::setfill(' ') << divl2[i]
    // << std::setw(15) << std::setfill(' ') << convdivPr[i]
    // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
    // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
    << std::setw(15) << std::setfill(' ') << divmax[i]
    // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
    << std::endl;
  }

}
#endif
