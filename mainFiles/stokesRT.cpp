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
#include "../num/matlab.hpp"

// #include "../num/gnuplot.hpp"

// #define PROBLEM_STOKES_DIV
#define DYNAMIC_DROP_EXAMPLE

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

int main(int argc, char** argv ) {
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
  ProblemOption optionStokes;
  optionStokes.order_space_element_quadrature_ = 9;
  vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

  for(int i=0;i<3;++i) { // i<3

    std::cout << "\n ------------------------------------- " << std::endl;
    Mesh Th(nx, ny, 0., 0., 1., 1.);


    Lagrange2 FEvelocity(4);
    Space Uh(Th, FEvelocity);
    // Space Qh(Th, DataFE<Mesh>::P1);

    Space Vh(Th, FEvelocity);
    // Space Vh(Th, DataFE<Mesh2>::RT1); // REMOVE THESE TWO LATER
    Space Qh(Th, DataFE<Mesh2>::P3dc); // FOR MIXEDSPACE

    const R hi = Th[0].lenEdge(0);
    const R penaltyParam = 1e1/pow(hi,3);
    const R sigma = 1; // HIGHER 1e2,1e3,etc WORSENS THE SOL [??]

    Fun_h fh(Uh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
    Fun_h gh(Uh, fun_exact_u);
    // Fun_h gh(Qh, fun_div); // THIS IS IF div u != 0

    FEM<Mesh2> stokes(Vh, optionStokes); stokes.add(Qh);

    Normal n;
    Tangent t;
    /* Syntax:
    FunTest (fem space, #components, place in space)
    */
    FunTest u(Vh,2,0), p(Qh,1,0), v(Vh,2,0), q(Qh,1,0);

    R mu = 1;
    // const MeshParameter& h(Parameter::h);
    // const MeshParameter& invlEdge(Parameter::invmeas);

    // stokes.addBilinear(
    //   (grad(div(u)),grad(div(v)))
    //   , Th
    // );

    // a_h(u,v)
    stokes.addBilinear(
      contractProduct(mu*grad(u),grad(v))
      - innerProduct(p,div(v))
      + innerProduct(div(u),q)
      , Th
    );
    // stokes.addBilinear(
    //   - innerProduct(average(grad(u*t)*n,0.5,0.5), jump(v*t))
    //   + innerProduct(jump(u*t), average(grad(v*t)*n,0.5,0.5))
    //   + innerProduct(1./hi*(sigma*jump(u*t)), jump(v*t))
    //   , Th
    //   , innerFacet
    // );
    stokes.addBilinear(
      - innerProduct(grad(u)*n, v)  //natural
      + innerProduct(u, grad(v)*n)
      + innerProduct(p, v*n)  // natural
      // - innerProduct(u*n, q)  // essential
      + innerProduct(penaltyParam*u, v)
      , Th
      , boundary
    );

    // l(v)_Omega
    stokes.addLinear(
      + innerProduct(gh.expression(2), penaltyParam*v)
      + innerProduct(gh.expression(2), grad(v)*n)
      // - innerProduct(gh.expression(2), q*n)
      , Th
      , boundary
    );
    stokes.addLinear(
      innerProduct(fh.expression(2),v)
      , Th
    );


    stokes.addLagrangeMultiplier(
      innerProduct(1.,p), 0., Th
    );



    stokes.solve();

    // EXTRACT SOLUTION
    int idx0_s = Vh.get_nb_dof();
    Rn_ data_uh = stokes.rhs_(SubArray(Vh.get_nb_dof(),0));
    Rn_ data_ph = stokes.rhs_(SubArray(Qh.get_nb_dof(),idx0_s));
    Fun_h uh(Vh, data_uh);
    Fun_h ph(Qh, data_ph);
    ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
    ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);


    {
      Fun_h soluErr(Vh, fun_exact_u);
      Fun_h soluh(Vh, fun_exact_u);
      soluErr.v -= uh.v;
      soluErr.v.map(fabs);
      // Fun_h divSolh(Wh, fun_div);
      // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

      Paraview<Mesh> writer(Th, "stokes_"+to_string(i)+".vtk");
      writer.add(uh, "velocity" , 0, 2);
      writer.add(ph, "pressure" , 0, 1);
      writer.add(femSol_0dx+femSol_1dy, "divergence");
      writer.add(soluh, "velocityExact" , 0, 2);
      writer.add(soluErr, "velocityError" , 0, 2);
      // writer.add(solh, "velocityError" , 0, 2);

      // writer.add(fabs(femDiv, "divergenceError");
    }

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
      convu.push_back(log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
      convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
    }

    nx *= 2;
    ny *= 2;
  }
  std::cout << "\n" << std::left
  << std::setw(10) << std::setfill(' ') << "h"
  << std::setw(15) << std::setfill(' ') << "err_p u"
  << std::setw(15) << std::setfill(' ') << "conv p"
  << std::setw(15) << std::setfill(' ') << "err u"
  << std::setw(15) << std::setfill(' ') << "conv u"
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
    << std::setw(15) << std::setfill(' ') << convp[i]
    << std::setw(15) << std::setfill(' ') << ul2[i]
    << std::setw(15) << std::setfill(' ') << convu[i]
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

#ifdef DYNAMIC_DROP_EXAMPLE
namespace DynamicDropData {
R2 shift(0,0);
double sigma = 0;
double mu1 = 1;
double mu2 = 1;

  R fun_levelSet(const R2 P, int i) {
    return sqrt((P.x-shift.x)*(P.x-shift.x) + (P.y-shift.y)*(P.y-shift.y)) - 2./3;
  }

  static R falpha1(const R r) {
    // const R MU1 = 1e-1;
    // const R MU2 = 1e-3;
    const R MU1 = mu1;
    const R MU2 = mu2;
    const R r0 = 2./3;
    return 1./MU1 + (1./MU2-1./MU1)*exp(r*r-r0*r0);
  }
  static R falpha2(const R r) {
    R MU2 = mu2;//1e-3;
    return 1./MU2;
  }

  static R falpha(const R r) {
    const R r0 = 2./3;
    return (r < r0)? falpha2(r) : falpha1(r);
  }

  static R fun_rhs(const R2 P, int i) {
    const R r2 = Norme2_2(P-shift);
    const R s = exp(-r2);
    R2 R(4*s*(r2-2)*P.y + 3*P.x*P.x, -4*s*(r2-2)*P.x );
    return (i<2)?R[i] : 0;
  }
  static R fun_boundary(const R2 P,const int i) {
    const R r = Norme2(P-shift);
    R2 R(-P.y,P.x); R = falpha(r)*exp(-r*r)*R;
    return (i<2)?R[i] : 0;
  }

  static R2 fun_velocity1(const R2 P) {
    R r = Norme2(P-shift);
    R2 R(-P.y,P.x); R = falpha1(r)*exp(-r*r)*R; return R;
  }
  static R2 fun_velocity2(const R2 P) {
    R r = Norme2(P-shift);
    R2 R(-P.y,P.x); R = falpha2(r)*exp(-r*r)*R; return R;
  }
  static R fun_pressure1(const R2 P) { return pow(P.x,3) ;}
  static R fun_pressure2(const R2 P) {
    //R sigma = 1;//24.5;//700;
    return pow(P.x,3) + sigma*3./2.;
  }

  // static R fun_solution(const R2 P, int ci, int domain) {
  //   if(domain == 0) return (ci < 2)? fun_velocity1(P)[ci] : fun_pressure1(P);
  //   else return (ci < 2)? fun_velocity2(P)[ci] : fun_pressure2(P);
  // }
  static R fun_exact_u(const R2 P, int ci, int domain) {
    if(domain == 0) return fun_velocity1(P)[ci];
    else            return fun_velocity2(P)[ci];
  }
  static R fun_exact_p(const R2 P, int ci, int domain) {
    if(domain == 0) return fun_pressure1(P);
    else            return fun_pressure2(P);
  }

  R2 fparam(double t){ return R2(2./3*cos(t+1./3), 2./3*sin(t+1./3));}

}

namespace Erik_Data_StokesDiv {

  // R fun_rhs(const R2 P, int i) {
  //   return 0;
  // }
  // R fun_exact(const R2 P, int i) {
  //   if(i==0) return 0;
  //   else if(i==1) return P.x;
  //   else return 0;
  // }
  R shift = 0.5;
  R interfaceRad = 0.250001;//2./3; // not exactly 1/4 to avoid interface cutting exaclty a vertex
  R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift)) - interfaceRad;
  }

  R fun_rhs(const R2 P, int i, int dom) {
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
  R fun_exact_u(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    // if(i==0) return      2*x*y*(x-1)*(y-1)*x*(1-x)*(2*y-1);
    // else if(i==1) return 2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x-1);
    if(i==0)      return  x*x*y;
    else if(i==1) return -x*y*y;
  }
  R fun_exact_p(const R2 P, int i, int dom) {
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
using namespace DynamicDropData;

// class LambdaGamma : public VirtualParameter { public :
// const CutFEMParameter& mu_;
//
// LambdaBoundary(const CutFEM_Parameter& mu, double G, double H) : mu_(mu) , G_(G) , H_(H) {} double evaluate(int domain, double h, double meas, double measK, double measCut) const {
// double gamma = meas / h ;
// double alpha = measK / h / h;
// double val = mu_.evaluate(domain,h,meas,measK,measCut) return val/h(G_+H_*gamma/alpha) ;}
// };

int main(int argc, char** argv ) {
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

  ProblemOption optionStokes;
  optionStokes.order_space_element_quadrature_ = 9;

  vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

  for(int i=0;i<3;++i) { // i<3

    std::cout << "\n ------------------------------------- " << std::endl;
    Mesh Kh(nx, ny, -1., -1., 2., 2.);
    const R hi = 2./(nx-1);
    const R penaltyParam = 1e1/hi;
    // const R sigma = 1;
    double uPenParam = 1e0;
    double pPenParam = 1e0;
    double jumpParam = 1e0;

    // CutFEMParameter mu(1e-1,1e-3);
    // CutFEMParameter invmu(1e1,1e3);
    CutFEMParameter mu(mu1,mu2);
    // double sigma = 0;//24.5;//700;
    double kappa1 = 0.5;//mu(1)/(mu(0)+mu(1));
    double kappa2 = 0.5;//mu(0)/(mu(0)+mu(1));


    Space Lh(Kh, DataFE<Mesh2>::P1);
    Fun_h levelSet(Lh, fun_levelSet);
    InterfaceLevelSet<Mesh> interface(Kh, levelSet);


    Lagrange2 FEvelocity(2);
    Space UUh(Kh, FEvelocity);

    // Space Qh(Kh, DataFE<Mesh>::P1);


    // Space Wh(Kh, DataFE<Mesh2>::RT1); // REMOVE KhESE TWO LATER
    Space Wh(Kh, FEvelocity);
    Space Qh(Kh, DataFE<Mesh2>::P1); // FOR MIXEDSPACE

    ActiveMesh<Mesh> Khi(Kh, interface);
    // Khi.truncate(interface, 1);
    CutSpace Vh(Khi, Wh);
    CutSpace Ph(Khi, Qh);
    CutSpace Uh(Khi, UUh);



    Fun_h fh(Uh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
    Fun_h gh(Uh, fun_exact_u);
    // Fun_h gh(Qh, fun_div); // THIS IS IF div u != 0

    CutFEM<Mesh2> stokes(Vh,optionStokes); stokes.add(Ph);

    Normal n;
    Tangent t;
    FunTest u(Vh,2,0), p(Ph,1,0), v(Vh,2,0), q(Ph,1,0),p1(Ph,1,0,0);

    stokes.addBilinear(
      // contractProduct(mu*grad(u),grad(v))
      contractProduct(2*mu*Eps(u),Eps(v))
      - innerProduct(p,div(v))
      + innerProduct(div(u),q)
      , Khi
    );

    // stokes.addBilinear(
    //   - innerProduct(average(2*mu*Eps(u)*t*n,kappa1,kappa2), jump(v*t))
    //   + innerProduct(jump(u*t), average(2*mu*Eps(v)*t*n,kappa1,kappa2))
    //   + innerProduct(1./hi*(jump(u*t)), jump(v*t))
    //   , Khi
    //   , innerFacet
    // );
    stokes.addBilinear(
      innerProduct(jump(u), -2*mu*average(Eps(v)*n, kappa1,kappa2))
      + innerProduct(-2*mu*average(Eps(u)*n,kappa1,kappa2), jump(v))
      + innerProduct(1./hi/hi*jump(u), jump(v))
      + innerProduct(average(p,kappa1,kappa2), jump(v*n))
      // - innerProduct(jump(u*n), average(q,0.5,0.5))
      , interface
    );
    // stokes.addBilinear(
    //   innerProduct(jump(u), -2*mu*average(grad(v)*n, 0.5,0.5))
    //   + innerProduct(-2*mu*average(grad(u)*n,0.5,0.5), jump(v))
    //   + innerProduct(10*jump(u), jump(v))
    //   + innerProduct(average(p,0.5,0.5), jump(v*n))
    //   - innerProduct(jump(u*n), average(q,0.5,0.5))
    //   , interface
    // );
    stokes.addLinear(
      innerProduct(3./2, average(v*n, kappa2,kappa1))*sigma
      , interface
    );

    stokes.addBilinear(
      - innerProduct(2*mu*Eps(u)*n, v)  //natural
      + innerProduct(u, 2*mu*Eps(v)*n)
      + innerProduct(p, v*n)  // natural
      + innerProduct(penaltyParam*u, v)
      // - innerProduct(u*n, q)  // essential
      , Khi
      , boundary
    );
    stokes.addLinear(
      + innerProduct(gh.expression(2), penaltyParam*v)
      + innerProduct(gh.expression(2), 2*mu*Eps(v)*n)
      // - innerProduct(gh.expression(2), q*n)
      , Khi
      , boundary
    );

    // l(v)_Omega
    stokes.addLinear(
      innerProduct(fh.expression(2),v)
      , Khi
    );

  //   FunTest grad2un = grad(grad(u)*n)*n;
  //   stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
  //   //  innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v)) // [Method 1: Remove jump in vel]
  //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
  //   // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
  //   // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
  //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))
  //
  //   // +innerProduct(uPenParam*pow(hi,1)*jump(u), mu*jump(v)) // [Method 1: Remove jump in vel]
  //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), mu*jump(grad(v)*n))
  //   // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
  //
  //   // -innerProduct(pPenParam*hi*jump(p), (1./mu)*jump(div(v)))
  //   // +innerProduct(pPenParam*hi*jump(div(u)), (1./mu)*jump(q))
  //   // -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), (1./mu)*jump(grad(div(v))))
  //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(v))) , (1./mu)*jump(grad(q)))
  //
  //   // +innerProduct(uPenParam*pow(hi,3)*jump(div(u)), mu*jump(div(v))) // [Method 1: Remove jump in vel]
  //   // +innerProduct(uPenParam*pow(hi,5)*jump(grad(div(u))), mu*jump(grad(div(v))))
  //   , Khi
  //   // , macro
  // );

  // Sets uniqueness of the pressure
  // double tt0 = MPIcf::Wtime();
  // int N = stokes.get_nb_dof();
  // std::map<int, double> df2rm;
  // R2 P = Qh[0].Pt(0);
  // double val = fun_exact_p(P, 0, 0);
  // df2rm.insert({Vh.get_nb_dof(), val});
  // eraseRow(N, stokes.mat_, stokes.rhs_, df2rm);
  // std::cout << "Time setting p val " << MPIcf::Wtime() - tt0 << std::endl;
  stokes.addLagrangeMultiplier(
    innerProduct(1.,p1), 0., Khi
  );

  matlab::Export(stokes.mat_, "matB"+to_string(i)+".dat");

    stokes.solve();

    // EXTRACT SOLUTION
    int idx0_s = Vh.get_nb_dof();
    Rn_ data_uh = stokes.rhs_(SubArray(Vh.get_nb_dof(),0));
    Rn_ data_ph = stokes.rhs_(SubArray(Ph.get_nb_dof(),idx0_s));
    Fun_h uh(Vh, data_uh);
    Fun_h ph(Ph, data_ph);
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
    //   Paraview<Mesh> writer(Khi, "stokes_"+to_string(i)+".vtk");
    //   writer.add(uh, "velocity" , 0, 2);
    //   writer.add(ph, "pressure" , 0, 1);
    //   writer.add(fabs(femSol_0dx+femSol_1dy), "divergence");
    //   writer.add(soluh, "velocityExact" , 0, 2);
    //   writer.add(soluErr, "velocityError" , 0, 2);
    //   // writer.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceError");
    //
    //   // writer.add(solh, "velocityError" , 0, 2);
    //
    //   // writer.add(fabs(femDiv, "divergenceError");
    // }

    R errU      = L2normCut(uh,fun_exact_u,0,2);
    R errP      = L2normCut(ph,fun_exact_p,0,1);
    R errDiv1    = L2normCut (femSol_0dx+femSol_1dy,Khi, 1);
    R errDiv0    = L2normCut (femSol_0dx+femSol_1dy,Khi, 0);
    R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy,Khi);


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
    divl2.push_back(errDiv0);
    divmax.push_back(errDiv1);

    // divmax.push_back(maxErrDiv);
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
  << std::setw(15) << std::setfill(' ') << "conv p"
  << std::setw(15) << std::setfill(' ') << "err u"
  << std::setw(15) << std::setfill(' ') << "conv u"
  << std::setw(15) << std::setfill(' ') << "err divu0"
  << std::setw(15) << std::setfill(' ') << "err divu1"

  // << std::setw(15) << std::setfill(' ') << "conv divu"
  // << std::setw(15) << std::setfill(' ') << "err_new divu"
  // << std::setw(15) << std::setfill(' ') << "convLoc divu"
  // << std::setw(15) << std::setfill(' ') << "err maxdivu"
  // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
  << "\n" << std::endl;
  for(int i=0;i<h.size();++i) {
    std::cout << std::left
    << std::setw(10) << std::setfill(' ') << h[i]
    << std::setw(15) << std::setfill(' ') << pl2[i]
    << std::setw(15) << std::setfill(' ') << convp[i]
    << std::setw(15) << std::setfill(' ') << ul2[i]
    << std::setw(15) << std::setfill(' ') << convu[i]
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
