#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "../util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif

// #include "finiteElement.hpp"
#include "baseProblem.hpp"
#include "../num/matlab.hpp"
#include "paraview.hpp"
// #include "projection.hpp"
// #include "../common/fracture.hpp"

// #define TEST_PIOLA
// #define DARCY_FEM
#define DARCY_EXAMPLE_SCOTTI
// #define DARCY_EXAMPLE_PUPPI
// #define DARCY_FRACTURE
// #define DARCY_MULTI_FRACTURE
// #define PROBLEM_RT_PROJECTION
// #define PROBLEM_MIXED_DARCY


void Scotty_diagonal_preconditioner(int N, std::map<std::pair<int,int>,R>& P){

  // SparseMatrixRC<double> B(N ,N ,P );
  // P.clear();
  Rn v(N);

  // create the diagonal Matrix
  for(int i=0;i<N;++i){
    v[i] = P[std::make_pair(i,i)] ;
    // for(int k=B.p[i];k<B.p[i+1];++k){
      // P[std::make_pair(i,i)] += B.a[k];
    // }
  }
  P.clear();
  for(int i=0;i<N;++i){
    P[std::make_pair(i,i)] = 1./sqrt(v[i]);//P[std::make_pair(i,i)];
  }

}

#ifdef TEST_PIOLA


namespace Data_p{
  R d_x = 2.;
  R d_y = 2.;

  R fun_neumann(const R2 P, int compInd) {
    return 0;
  }
  R fun_div(const R2 P, int compInd, int dom) {// is also exact divergence
    return 0;
  }


}
using namespace Data_p;


int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;
  typedef Mesh2 Mesh;
  typedef ActiveMeshT2 CutMesh;
  typedef FESpace2   Space;
  typedef CutFESpaceT2 CutSpace;


  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();
  bool set_p = false;
  int nx = 2; // 6
  int ny = 2; // 6

  Mesh Kh(nx, ny, 0., 0., 1., 1.);
  Kh.info();
  double h = 1./(nx-1);

  Space Vh(Kh, DataFE<Mesh2>::RT0);
  Space Vhm(Kh, DataFE<Mesh2>::RT0m);

  typename Space::FElement FK(Vh[0]);
  typename Space::FElement FKm(Vhm[0]);

  KNMK<double> f (FK.NbDoF(),2,Fop_D1);
  KNMK<double> fm(FKm.NbDoF(),2,Fop_D1);


}

#endif

#ifdef DARCY_FEM


namespace Data_p{
  R d_x = 2.;
  R d_y = 2.;

  R fun_neumann(const R2 P, int compInd) {
    return 0;
  }
  R fun_div(const R2 P, int compInd, int dom) {// is also exact divergence
    return 0;
  }


}
using namespace Data_p;


int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;
  typedef Mesh2 Mesh;
  typedef ActiveMeshT2 CutMesh;
  typedef FESpace2   Space;
  typedef CutFESpaceT2 CutSpace;


  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();
  bool set_p = false;
  int nx = 2; // 6
  int ny = 2; // 6

  vector<double> uPrint,pPrint,divPrint,divPrintLoc,maxDivPrint,h,convuPr,convpPr,convdivPr,convdivPrLoc,convmaxdivPr;
  vector<double> ratioCut1, ratioCut2;
  int iters =5;

  for(int i=0;i<iters;++i) {
    Mesh Kh(nx, ny, 0., 0., 1., 1.);
    Kh.info();
    double h = 1./(nx-1);


    if(set_p){
      Space Vh(Kh, DataFE<Mesh2>::RT0);
      FEM<Mesh> darcy(Vh);
      FunTest u(Vh,2), v(Vh,2);
      // [ASSEMBLY]
      darcy.addBilinear(
        innerProduct(u, v) + innerProduct(div(u), div(v))
        , Kh
      );
      matlab::Export(darcy.mat_, "matB"+to_string(i)+".dat");
      nx = 2*nx -1;
      ny = nx;
      continue;

    }
    else{
      Space Vh(Kh, DataFE<Mesh>::RT0);
      Space Ph(Kh, DataFE<Mesh>::P0sc);
      FEM<Mesh> darcy(Vh); darcy.add(Ph);
      FunTest u(Vh,2), v(Vh,2),p(Ph,1), q(Ph,1);
      Vh.info();
      Ph.info();
      // [ASSEMBLY]
      darcy.addBilinear(
        innerProduct(u, v)
        + innerProduct(div(u), q) - innerProduct(p, div(v))
        , Kh
      );

      matlab::Export(darcy.mat_, "matB"+to_string(i)+".dat");
      nx = 2*nx-1;
      ny = nx;
      continue;

    }

  }

}

#endif

#ifdef DARCY_EXAMPLE_SCOTTI

namespace Data_CutMixedDarcy {
  bool solHasJump = true;
  R shift = 0.5;
  R interfaceRad = 0.250001; // not exactly 1/4 to avoid interface cutting exaclty a vertex
  R mu_G = 2./3*interfaceRad; // xi0*mu_G = 1/8*2/3*1/4

  R fun_radius2(const R2 P){
    return (P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift);
  }
  R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift)) - interfaceRad;
  }

  R fun_dirichlet(const R2 P, int compInd) {
    if (compInd == 0)
    // return 2*sin(2*pi*P[0])*sin(4*pi*P[1]);
    return 0;
    else
    return 0;
  }
  R fun_neumann(const R2 P, int compInd, int dom) {
    R r2 = fun_radius2(P);
    R radius2 = interfaceRad*interfaceRad;
    // if (dom==0) // r2>radius2
    return r2/(2*radius2)+3./2.;
    // else
    //   return r2/(radius2);
  }
  R fun_force(const R2 P, int compInd) {
    if (compInd == 0)
    // return 40*pow(pi,2)*sin(2*pi*P.x)*sin(4*pi*P.y);
    return 0; //1
    else
    return 0;
  }
  R fun_div(const R2 P, int compInd, int dom) {// is also exact divergence
    R radius2 = interfaceRad*interfaceRad;
    if (dom==0) // r2>radius2
    return -2./radius2;
    else
    return -4./radius2;
  }
  R fun_exact_u(const R2 P, int compInd, int dom) {
    Diff<R,2> X(P.x,0), Y(P.y,1);
    Diff<R,2> r2 = (X-shift)*(X-shift) + (Y-shift)*(Y-shift);
    R radius2 = interfaceRad*interfaceRad;
    R cst = (dom==0)*3./2;
    R mul = (dom==0)*2 + (dom==1)*1;
    Diff<R, 2> val = r2/(mul*radius2) + cst;
    return -val.d[compInd];
  }
  R fun_exact_p(const R2 P, int compInd, int dom) {
    Diff<R,2> X(P.x,0), Y(P.y,1);
    Diff<R,2> r2 = (X-shift)*(X-shift) + (Y-shift)*(Y-shift);
    R radius2 = interfaceRad*interfaceRad;
    R cst = (dom==0)*3./2;
    R mul = (dom==0)*2 + (dom==1)*1;
    Diff<R, 2> val = r2/(mul*radius2) + cst;
    return val.val;
  }
  R fun_interfacePr(const R2 P, int compInd) {
    return 19./12;
  }
}
using namespace Data_CutMixedDarcy;
int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;
  typedef Mesh2 Mesh;
  typedef ActiveMeshT2 CutMesh;
  typedef FESpace2   Space;
  typedef CutFESpaceT2 CutSpace;


  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 51; // 6
  int ny = 51; // 6

  vector<double> uPrint,pPrint,divPrint,divPrintLoc,maxDivPrint,h,convuPr,convpPr,convdivPr,convdivPrLoc,convmaxdivPr;
  vector<double> ratioCut1, ratioCut2;
  int iters =1;

  for(int i=0;i<iters;++i) {
    Mesh Kh(nx, ny, 0., 0., 1., 1.);
    Kh.info();


    Space Lh(Kh, DataFE<Mesh2>::P1);
    Fun_h levelSet(Lh, fun_levelSet);
    InterfaceLevelSet<Mesh> interface(Kh, levelSet);

    Space Vh(Kh, DataFE<Mesh2>::RT0);
    Space Qh(Kh, DataFE<Mesh2>::P0);


    ActiveMesh<Mesh> Kh_i(Kh, interface);
    Kh_i.info();
    CutSpace Wh(Kh_i, Vh);
    Wh.info();
    CutSpace Ph(Kh_i, Qh);
    Ph.info();

    CutFEM<Mesh2> darcy(Wh); darcy.add(Ph);
    const R h_i = 1./(nx-1);
    const R invh = 1./h_i;
    MacroElement<Mesh> macro(Kh_i, 0.1);
    // ratioCut1.push_back((double)macro.nb_element_0 / (double)interface.nbElement());
    // ratioCut2.push_back((double)macro.nb_element_1 / (double)interface.nbElement());
    R xi = 3./4;
    R xi0 = (xi-0.5)/2.;
    // CutFEM_Parameter mu(1.,1.);
    // const CutFEM_Parameter& invh(Parameter::invh);
    // const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
    // const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
    // const R hei = Th[0].lenEdge(0);

    // We define fh on the cutSpace
    Fun_h fq(Ph, fun_div);
    Fun_h p0(Lh, fun_neumann);
    Fun_h phat(Ph, fun_interfacePr);


    Normal n;
    Tangent t;
    FunTest p(Ph,1), q(Ph,1), u(Wh,2), v(Wh,2);

    double uPenParam = 1e0;//1e-2; //cont 1e-1`
    double pPenParam = 1e0;//1e1; // cont 1e2
    double jumpParam = 1e0; // [anything<1e0 doesn't give full u convergence]

    darcy.addBilinear(
      innerProduct(u, v)
      -innerProduct(p, div(v))
      +innerProduct(div(u), q)
      , Kh_i
    );

    darcy.addBilinear(
      innerProduct(mu_G*average(u*n), average(v*n))
      +innerProduct(xi0*mu_G*jump(u*n), jump(v*n)) // b(p,v)-b(q,u) bdry terms
      ,interface
    );


    // matlab::Export(darcy.mat_, "matNew.dat");
    // nx = 2*nx-1;
    // ny = 2*ny-1;
    // continue;

    // l(v)_Omega
    darcy.addLinear(
      innerProduct(fq.expression(), q) , Kh_i
    );


    darcy.addLinear(
      -innerProduct(p0.expression(), v*n) // Only on Gamma_N (pressure)
      // + innerProduct(u0, invh*(v.t()*n)*lambdaB) // Only on Gamma_D (vel normal comp)
      // - innerProduct(u0, q) // Only on Gamma_D (vel normal comp)
      , Kh_i, boundary
    );
   //  darcy.addBilinear(
   //  +innerProduct(p, v*n)
   //  +innerProduct(u*n, 1*invh*v*n) // invh
   //  -innerProduct(u*n,q)
   // ,Kh_i, boundary);

    // std::cout << darcy.rhs_ << std::endl;
    // return 0;

    darcy.addLinear(
      -innerProduct(phat.expression(), jump(v*n))

      // -innerProduct(phat.expression(), jump(v.t()*n))
      ,interface
    );
    // matlab::Export(darcy.rhs_, "rhsNew.dat");

    int N = Wh.get_nb_dof()+Ph.get_nb_dof();



//     matlab::Export(darcy.mat_, "matB"+to_string(i)+".dat");
// return 0;
  // FunTest grad2un = grad(grad(u)*n)*n;
  darcy.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
  //  innerProduct(uPenParam*pow(h_i,0)*jump(u), jump(v)) // [Method 1: Remove jump in vel]
  // +innerProduct(uPenParam*pow(h_i,2)*jump(grad(u)*n), jump(grad(v)*n))
  // // +innerProduct(uPenParam*pow(h_i,5)*jump(grad2un), jump(grad2un))
  // +innerProduct(pPenParam*pow(h_i,0)*jump(p), jump(q))
  // +innerProduct(pPenParam*pow(h_i,2)*jump(grad(p)), jump(grad(q)))

   innerProduct(uPenParam*h_i*jump(u*n), jump(v*n)) // [Method 1: Remove jump in vel]
  +innerProduct(uPenParam*pow(h_i,3)*jump(grad(u)*n), jump(grad(v)*n))
  // +innerProduct(uPenParam*pow(h_i,5)*jump(grad2un), jump(grad2un))
  -innerProduct(pPenParam*h_i*jump(p), jump(div(v)))
  +innerProduct(pPenParam*h_i*jump(div(u)), jump(q))
  -innerProduct(pPenParam*pow(h_i,3)*jump(grad(p)), jump(grad(div(v))))
  +innerProduct(pPenParam*pow(h_i,3)*jump(grad(div(v))) , jump(grad(q)))
  , Kh_i
  , macro
);

// matlab::Export(darcy.mat_, "matB"+to_string(i)+".dat");
// nx = 2*nx-1;
// ny = 2*ny-1;
// continue;

    //
    // {  // Diagonal Preconditionning
    //   std::map<std::pair<int,int>,R> P;
    //   darcy.set_map(P);
    //   darcy.addBilinear(
    //     innerProduct(u, v)
    //     + innerProduct(p, q)
    //     , Kh_i
    //   );
    //   darcy.addBilinear(
    //     innerProduct(mu_G*average(u*n), average(v*n))
    //     +innerProduct(xi0*mu_G*jump(u*n), jump(v*n)) // b(p,v)-b(q,u) bdry terms
    //     ,interface
    //   );
    //
    //   Scotty_diagonal_preconditioner(N, P);
    //   darcy.set_map();
    //
    //
    //   // darcy.leftPreconditioning(P);
    //   darcy.applyPreconditioning(P);
    //   //
    //   matlab::Export(darcy.mat_, "matB"+to_string(i)+".dat");
    //   // // matlab::Export(P, "matA"+to_string(i)+".dat");
    //   // darcy.solve("umfpack");
    //
    //   nx = 2*nx-1;
    //   ny = 2*ny-1;
    //   continue;
    //
    //   darcy.recoverSolution(P);
    //
    //
    // }
    //


  // Test 1 - on all triangle
  // It destroy the pressure
  // darcyCutMix.addBilinear(
  //   innerProduct(1e-2*div(u), div(v))
  // );

  // Test 2 - on triangle in macro elements
  // It destroy the pressure and not better condition number
  // darcyCutMix.addBilinear(
  //   innerProduct(1e-2*div(u), div(v))
  //   , macro
  // );

  // Test 3 - on inner edges in macro elements
  // darcyCutMix.addEdgeIntegral(
  //   innerProduct(1e-2*hei*jump(div(u)), jump(div(v)))
  //   , macro
  // );


  // matlab::Export(darcy.mat_, "matB"+to_string(i)+".dat");
  // matlab::Export(darcy.rhs_, "rhsA"+to_string(i)+".dat");

  darcy.solve("umfpack");

  // darcyCutMix.solve();
  // EXTRACT SOLUTION
  int idx0_s = Wh.get_nb_dof();
  Rn_ data_uh = darcy.rhs_(SubArray(Wh.get_nb_dof(),0));
  Rn_ data_ph = darcy.rhs_(SubArray(Ph.get_nb_dof(),idx0_s));
  Fun_h uh(Wh, data_uh);
  Fun_h ph(Ph, data_ph);
  ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
  ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);



  // L2 norm vel
  R errU      = L2normCut(uh,fun_exact_u,0,2);
  R errP      = L2normCut(ph,fun_exact_p,0,1);
  R errDiv    = L2normCut (femSol_0dx+femSol_1dy,fun_div,Kh_i);
  R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy,fun_div,Kh_i);
  // std::cout << " hey" << std::endl;
  // getchar();
  // [PLOTTING]
  {
    // Fun_h solh(Wh, fun_exact);
    // solh.v -= uh.v;
    // solh.v.map(fabs);
    Fun_h divSolh(Ph, fun_div);
    ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

    // Paraview<Mesh> writer(Kh_i, "darcyNew_"+to_string(i)+".vtk");
    Paraview<Mesh> writer(Kh_i, "darcyRT0scotti.vtk");

    writer.add(uh, "velocity" , 0, 2);
    writer.add(ph, "pressure" , 0, 1);
    writer.add(femSol_0dx+femSol_1dy, "divergence");
    // writer.add(solh, "velocityError" , 0, 2);
    writer.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceError");
  }

  pPrint.push_back(errP);
  uPrint.push_back(errU);
  divPrint.push_back(errDiv);
  maxDivPrint.push_back(maxErrDiv);
  h.push_back(h_i);

  if(i==0) {convpPr.push_back(0);convuPr.push_back(0);convdivPr.push_back(0);convdivPrLoc.push_back(0);convmaxdivPr.push_back(0);}
  else {
    convpPr.push_back( log(pPrint[i]/pPrint[i-1])/log(h[i]/h[i-1]));
    convuPr.push_back( log(uPrint[i]/uPrint[i-1])/log(h[i]/h[i-1]));
    convdivPr.push_back( log(divPrint[i]/divPrint[i-1])/log(h[i]/h[i-1]));
    convmaxdivPr.push_back( log(maxDivPrint[i]/maxDivPrint[i-1])/log(h[i]/h[i-1]));
  }
  nx = 2*nx-1;
  ny = 2*ny-1;
}

std::cout << "\n" << std::left
<< std::setw(10) << std::setfill(' ') << "h"
<< std::setw(15) << std::setfill(' ') << "err_p"
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
for(int i=0;i<uPrint.size();++i) {
  std::cout << std::left
  << std::setprecision(5)
  << std::setw(10) << std::setfill(' ') << h[i]
  << std::setw(15) << std::setfill(' ') << pPrint[i]
  // << std::setw(15) << std::setfill(' ') << convpPr[i]
  << std::setw(15) << std::setfill(' ') << uPrint[i]
  // << std::setw(15) << std::setfill(' ') << convuPr[i]
  << std::setw(15) << std::setfill(' ') << divPrint[i]
  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
  << std::setw(15) << std::setfill(' ') << maxDivPrint[i]
  // << std::setw(15) << std::setfill(' ') << ratioCut1[i]
  // << std::setw(15) << std::setfill(' ') << ratioCut2[i]
  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
  << std::endl;
}

}

#endif

#ifdef DARCY_EXAMPLE_PUPPI
namespace Data_CutMixedDarcyPUPPI_quarter_circle {
R d_x = 2.;
R d_y = 2.;
R shift = 0;
R interfaceRad = 0.501;
R mu_G = 2./3*interfaceRad; // xi0*mu_G = 1/8*2/3*1/4

R fun_radius2(const R2 P){
  return (P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift);
}
R fun_levelSet(const R2 P, const int i) {
  return sqrt((P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift)) - interfaceRad;
}

// R fun_neumann_east(const R2 P, int compInd) {
//   return -sin(P.x)*sinh(P.y) - (cos(1) - 1)*(cosh(1) - 1);
// }
R fun_natural(const R2 P, int compInd) {
  return -sin(P.x)*sinh(P.y) - (cos(1) - 1)*(cosh(1) - 1);
}

R fun_enforced(const R2 P, int compInd) {
  return (compInd == 0)? cos(P.x)*sinh(P.y) : sin(P.x)*cosh(P.y);
}


R fun_force(const R2 P, int compInd) {
  return 0;
}
R fun_div(const R2 P, int compInd, int dom) {// is also exact divergence
  return 0;
}
R fun_exact(const R2 P, int compInd, int dom) {
  // Diff<R,2> X(P.x,0), Y(P.y,1);
  // Diff<R, 2> val = - sin(X)*sinh(Y) - (cos(1) - 1)*(cosh(1) - 1);
  // if (compInd==2) {
  //   return val.val;
  // }
  // else {
  //   return -val.d[compInd];
  // }
  if (compInd==2) {
    return -sin(P.x)*sinh(P.y) - (cos(1) - 1)*(cosh(1) - 1);
  } else if (compInd==0) {
    return cos(P.x)*sinh(P.y);
  } else {
    return sin(P.x)*cosh(P.y);
  }
}
R fun_interfacePr(const R2 P, int compInd) {
  if (compInd == 2)
  return 19./12;
  else
  return 0;
}
}
namespace Data_CutMixedDarcyPUPPI_circle {
  bool solHasJump = true;
  R d_x = 1.;
  R d_y = 1.;
  R shift = 0.5;
  R interfaceRad = 0.4501; // not exactly 1/4 to avoid interface cutting exaclty a vertex
  R mu_G = 2./3*interfaceRad; // xi0*mu_G = 1/8*2/3*1/4
  R pie = 3.14159265359;

  R fun_levelSet(const R2 P, const int i) {
    return interfaceRad - sqrt((P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift));
  }

  // R fun_levelSet(const R2 P, const int i) {
  //   return 1.5001-P.y;
  // }

  R fun_natural(const R2 P, int compInd) {
    return -sin(2*pie*P.x)*cos(2*pie*P.y);
  }
  R fun_enforced(const R2 P, int compInd) { // [these need to be enforced!]
    if (compInd == 0)
      return 2*pie*cos(2*pie*P.x)*cos(2*pie*P.y);
    else if (compInd == 1)
      return -2*pie*sin(2*pie*P.x)*sin(2*pie*P.y);
    else return 0;
  }
  R fun_force(const R2 P, int compInd) {
    return 0;
  }
  R fun_div(const R2 P, int compInd, int dom) {// is also exact divergence
    return -8*pie*pie*sin(2*pie*P.x)*cos(2*pie*P.y);
  }
  R fun_exact_u(const R2 P, int compInd, int dom) {
    if (compInd==0) {
      return 2*pie*cos(2*pie*P.x)*cos(2*pie*P.y);
    } else {
      return -2*pie*sin(2*pie*P.x)*sin(2*pie*P.y);
    }
  }
  R fun_exact_p(const R2 P, int compInd, int dom) {
    return -sin(2*pie*P.x)*cos(2*pie*P.y);
  }

  // R fun_natural(const R2 P, int compInd) {
  //   return -sin(P.x)*sinh(P.y) - (cos(1) - 1)*(cosh(1) - 1);
  // }
  // R fun_enforced(const R2 P, int compInd) { // [these need to be enforced!]
  //   if (compInd == 0)
  //     return cos(P.x)*sinh(P.y);
  //   else if (compInd == 1)
  //     return sin(P.x)*cosh(P.y);
  //   else return 0;
  // }
  // R fun_force(const R2 P, int compInd) {
  //   return 0;
  // }
  // R fun_div(const R2 P, int compInd, int dom) {// is also exact divergence
  //   return 0;
  // }
  // R fun_exact(const R2 P, int compInd, int dom) {
  //   if (compInd==2) {
  //     return -sin(P.x)*sinh(P.y) - (cos(1) - 1)*(cosh(1) - 1);
  //   } else if (compInd==0) {
  //     return cos(P.x)*sinh(P.y);
  //   } else {
  //     return sin(P.x)*cosh(P.y);
  //   }
  // }
}
using namespace Data_CutMixedDarcyPUPPI_circle;

int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;
  typedef Mesh2 Mesh;
  typedef ActiveMeshT2 CutMesh;
  typedef FESpace2   Space;
  typedef CutFESpaceT2 CutSpace;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx =10; // 6
  int ny =10; // 6

  vector<double> uPrint,pPrint,divPrint,divPrintLoc,maxDivPrint,h,convuPr,convpPr,convdivPr,convdivPrLoc,convmaxdivPr;

  int iters = 4;

  for(int i=0;i<iters;++i) {
    Mesh Kh(nx, ny, 0., 0., d_x, d_y);

    Space Lh(Kh, DataFE<Mesh2>::P1);
    Fun_h levelSet(Lh, fun_levelSet);
    InterfaceLevelSet<Mesh> interface(Kh, levelSet);

    Space Vh(Kh, DataFE<Mesh>::RT0);
    Space Qh(Kh, DataFE<Mesh>::P0); // for the RHS
    Space Q2h(Kh, DataFE<Mesh>::P2); // for the RHS

    // [CutFEM]
    ActiveMesh<Mesh> Kh_i(Kh);
    Kh_i.truncate(interface, -1);

    // Kh_i.info();
    CutSpace Wh(Kh_i, Vh);
    Wh.info();
    CutSpace Ph(Kh_i, Qh);
    CutSpace P2h(Kh_i, Q2h);

    MacroElement<Mesh> macro(Kh_i, 0.25);


    CutFEM<Mesh2> darcy(Wh); darcy.add(Ph);
    const R h_i = 1./(nx-1);
    const R invh = 1./h_i;

    R xi = 3./4; // > 1/2
    R xi0 = (xi-0.5)/2.;

    // We define fh on the cutSpace
    Fun_h fv(Wh, fun_force);
    Fun_h fq(P2h, fun_div);
    // Fun_h p0(Ph, fun_natural);
    Fun_h u0(Wh, fun_enforced);
    // Fun_h phat(P2h, fun_interfacePr);

    Normal n;
    Tangent t;
    FunTest p(Ph,1), q(Ph,1), u(Wh,2), v(Wh,2);

    // u=1e-2 p=1e1 REALLY GOOD for jump problem!
    double uPenParam = 1e0;//1e-2; //cont 1e-1`
    double pPenParam = 1e0;//1e1; // cont 1e2
    double penParam = 1e0; // [anything<1e0 doesn't give full u convergence]
    // [ASSEMBLY]
    //:: a(u,v)_Omega
    darcy.addBilinear(
      innerProduct(u, v)
      -innerProduct(p, div(v))
      +innerProduct(div(u), q)
      , Kh_i
    );
    // l(v)_Omega
    darcy.addLinear(
      innerProduct(fq.expression(), q) , Kh_i
    );

    // darcy.cleanMatrix();
    darcy.addLagrangeMultiplier(
      innerProduct(1.,p), 0. , Kh_i
    );
    // matlab::Export(darcy.mat_, "matBthomas.dat");
    // return 0;

    darcy.addBilinear(
      +innerProduct(p, v*n)
      +innerProduct(penParam*u*n, 1./h_i*v*n)
      ,interface
    );


    darcy.addLinear(
      +innerProduct(u0*n, penParam*1./h_i*v*n)
      ,interface
    );
    FunTest grad2un = grad(grad(u)*n)*n;
  //   darcy.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
  //   //  innerProduct(uPenParam*pow(h_i,1)*jump(u), jump(v)) // [Method 1: Remove jump in vel]
  //   // +innerProduct(uPenParam*pow(h_i,3)*jump(grad(u)*n), jump(grad(v)*n))
  //   // +innerProduct(uPenParam*pow(h_i,5)*jump(grad2un), jump(grad2un))
  //   // +innerProduct(pPenParam*pow(h_i,1)*jump(p), jump(q))
  //   // +innerProduct(pPenParam*pow(h_i,3)*jump(grad(p)), jump(grad(q)))
  //    innerProduct(uPenParam*h_i*jump(u), jump(v)) // [Method 1: Remove jump in vel]
  //   +innerProduct(uPenParam*pow(h_i,3)*jump(grad(u)*n), jump(grad(v)*n))
  //   +innerProduct(uPenParam*pow(h_i,5)*jump(grad2un), jump(grad2un))
  //   -innerProduct(pPenParam*pow(h_i,1)*jump(p), jump(div(v)))
  //   +innerProduct(pPenParam*pow(h_i,1)*jump(div(u)), jump(q))
  //   -innerProduct(pPenParam*pow(h_i,3)*jump(grad(p)), jump(grad(div(v))))
  //   +innerProduct(pPenParam*pow(h_i,3)*jump(grad(div(v))) , jump(grad(q)))
  //   , Kh_i
  //   // , macro
  // );

    // matlab::Export(darcy.mat_, "matB"+to_string(i)+".dat");



    darcy.solve();

    // EXTRACT SOLUTION
    int idx0_s = Wh.get_nb_dof();
    Rn_ data_uh = darcy.rhs_(SubArray(Wh.get_nb_dof(),0));
    Rn_ data_ph = darcy.rhs_(SubArray(Ph.get_nb_dof(),idx0_s));
    Fun_h uh(Wh, data_uh);
    Fun_h ph(Ph, data_ph);
    ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
    ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);



    // L2 norm vel
    R errU      = L2normCut(uh,fun_exact_u,0,2);
    R errP      = L2normCut(ph,fun_exact_p,0,1);
    R errDiv    = L2normCut (femSol_0dx+femSol_1dy,fun_div,Kh_i);
    R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy,fun_div,Kh_i);

    // [PLOTTING]
    {
      // Fun_h solh(Wh, fun_exact);
      // solh.v -= uh.v;
      // solh.v.map(fabs);
      Fun_h divSolh(Wh, fun_div);
      ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

      Paraview<Mesh> writer(Kh_i, "darcyPuppi_"+to_string(i)+".vtk");
      writer.add(uh, "velocity" , 0, 2);
      writer.add(ph, "pressure" , 0, 1);
      writer.add(femSol_0dx+femSol_1dy, "divergence");
      // writer.add(solh, "velocityError" , 0, 2);
      writer.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceError");
    }


    pPrint.push_back(errP);
    uPrint.push_back(errU);
    // continue;
    divPrint.push_back(errDiv);
    // divPrintLoc.push_back(errDivLoc);

    maxDivPrint.push_back(maxErrDiv);
    h.push_back(h_i);

    if(i==0) {convpPr.push_back(0);convuPr.push_back(0);convdivPr.push_back(0);convdivPrLoc.push_back(0);convmaxdivPr.push_back(0);}
    else {
      convpPr.push_back( log(pPrint[i]/pPrint[i-1])/log(h[i]/h[i-1]));
      convuPr.push_back( log(uPrint[i]/uPrint[i-1])/log(h[i]/h[i-1]));
      convdivPr.push_back( log(divPrint[i]/divPrint[i-1])/log(h[i]/h[i-1]));
      // convdivPrLoc.push_back( log(divPrintLoc[i]/divPrintLoc[i-1])/log(h[i]/h[i-1]));

      convmaxdivPr.push_back( log(maxDivPrint[i]/maxDivPrint[i-1])/log(h[i]/h[i-1]));
    }

    nx = 2*nx-1;
    ny = 2*ny-1;
    // nx += 1;
    // ny += 1;
    // nx = (int)round( (1+0.2*i)*nx/2 )*2; // Makes a nonuniform refinement to an EVEN integer
    // ny = (int)round( (1+0.2*i)*ny/2 )*2;
    // std::cout << nx << std::endl;
    // shift = 0.5+(i+1)*h_i/iters; // moves one grid cell over entire span
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
  for(int i=0;i<uPrint.size();++i) {
    std::cout << std::left
    << std::setw(10) << std::setfill(' ') << h[i]
    << std::setw(15) << std::setfill(' ') << pPrint[i]
    // << std::setw(15) << std::setfill(' ') << convpPr[i]
    << std::setw(15) << std::setfill(' ') << uPrint[i]
    // << std::setw(15) << std::setfill(' ') << convuPr[i]
    << std::setw(15) << std::setfill(' ') << divPrint[i]
    // << std::setw(15) << std::setfill(' ') << convdivPr[i]
    // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
    // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
    << std::setw(15) << std::setfill(' ') << maxDivPrint[i]
    // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
    << std::endl;
  }
}



#endif

#ifdef DARCY_FRACTURE
R fun_levelSet(const R2 P, const int i) {
  return -(P.y+2*P.x) + 1.401;
  // return -(2*P.x) + 1.4;

}
R fun_neumann(const R2 P, int compInd, int dom) {
  return P.y;
}
R fun_1(const R2 P, int compInd, int dom) {
  return 1;
}
int main(int argc, char** argv ) {

  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx =50;
  int ny =50;
  std::vector<int> idx0;

  int iters = 1;

  for(int i=0;i<iters;++i) {

    Mesh2 Th(nx, ny, 0., 0., 1., 1.);


    FESpace2 Lh(Th, DataFE<Mesh2>::P1);
    Fun_h levelSet(Lh, fun_levelSet);

    Interface2 interface(Th, levelSet.v);

    KN<const GTypeOfFE<Mesh2>* > arrayFE(2);
    arrayFE(0) = &DataFE<Mesh2>::RT0;
    arrayFE(1) = &DataFE<Mesh2>::P0;
    GTypeOfFESum<Mesh2> mixedFE(arrayFE);


    FESpace2 Vh(Th, mixedFE);
    FESpace2 Ph (Th, DataFE<Mesh2>::P1);

    // CutFEM stuff
    CutFESpace2 Wh(Vh , interface, {1,-1});

    Mesh2 cutTh(interface);
    FESpace2 Sh(cutTh, mixedFE);
    Sh.backSpace = &Vh;




    // OBJECTS NEEDED FOR THE PROBLEM
    // =====================================================
    CutFEM<Mesh2> darcy(Wh);
    darcy.add(Sh);


    MacroElement macro(Wh, 0.05);
    MacroElementSurface macroInterface(interface, 0.5);
    // getchar();
    // Extension extension(darcy);
    // extension.tag_extension_edges(macro);

    // extension.tag_exhaust_edges(macro);
    //
    // gnuplot::save(Th);
    // gnuplot::save(interface);
    // gnuplot::save(macro, extension);
    // gnuplot::save(macroInterface);
    // gnuplot::save(Wh,0,"Vh1.dat");
    // gnuplot::save(Wh,1,"Vh2.dat");



    int idx_s0 = Wh.NbDoF();
    idx0.push_back(idx_s0);
    KN_<double> data_uh(darcy.rhs(SubArray(Wh.NbDoF(),0)));
    KN_<double> data_sh(darcy.rhs(SubArray(Sh.NbDoF(),idx_s0)));
    Fun_h uh(Wh, data_uh);  // fun representing bulk variables  [v1, v2]
    Fun_h us(Sh, data_sh);  // fun representing surface variable v0

    Normal n;
    Tangent t;
    FunTest p(Wh,1,2), q(Wh,1,2), u(Wh,2), v(Wh,2);
    FunTest p_hat(Sh,1,2), q_hat(Sh,1,2), u_hat(Sh,2), v_hat(Sh,2);

    Fun_h p0(Ph, fun_neumann);
    Fun_h uuu(Sh, fun_1);



    double nu = 1;
    double nu_G = 1;
    double nu_hat = 1;
    double gamma_D = 10;
    const CutFEM_Parameter& invh(Parameter::invh);
    const CutFEM_Parameter& h(Parameter::h);
    double xi = 3./4;
    double xi0 = (xi-0.5)/2.;
    double l_G = 0.01;
    double uPenParam = 1e-1;
    double pPenParam = 1e-1;
    double sPenParam = 1e0;
    double hei = 1. /(nx-1);

    // BULK PROBLEM
    // =====================================================
    darcy.addBilinear(
      innerProduct(nu*u, v)
      -innerProduct(p, div(v))
      +innerProduct(div(u), q)
    );

    darcy.addBilinear(
      innerProduct(gamma_D*invh*(u*n), (v*n))
      + innerProduct(p, v*n) - innerProduct(u*n, q)
      , boundary
      , {2,4}
    );

    darcy.addBilinear(
      innerProduct(     nu_G*average(u*n), average(v*n))
      +innerProduct(xi0*nu_G*jump(u*n)   , jump(v*n))
      ,interface
    );


    // l(v)_Omega
    darcy.addLinear(
      // innerProduct(fv, v)  // fv = ?????
      innerProduct(4, q)  // fq = 4
    );
    // FunTest fphat(phat, 2, 1);
    darcy.addBilinear(
      innerProduct(p_hat, jump(v*n))
      ,interface
    );
    // darcy.addLinear(
    //   -innerProduct(p0.expression(), jump(v*n))
    //   ,interface
    // );
    darcy.addLinear(
      -innerProduct(p0.expression(), v*n)
      , boundary
      , {1,3}
    );

    // SURFACE PROBLEM
    // =====================================================
    darcy.addBilinear(
      innerProduct(nu_hat*u_hat, v_hat)
      -innerProduct(p_hat, divS(v_hat))
      +innerProduct(divS(u_hat), q_hat)
      , interface
    );

    darcy.addBilinear(
      - innerProduct(jump(u*n), q_hat)
      , interface
    );

    darcy.addLinear(
      innerProduct(l_G*4, q_hat)
      , interface
    );

    darcy.addBilinear(
      - innerProduct(average(p_hat, 0.5, 0.5),jump(v_hat*t, 1,1))
      - innerProduct(jump(u_hat*t,1,1), average(q_hat, 0.5,0.5))
      + innerProduct(1e2*jump(u_hat*t,1,1), jump(v_hat*t,1,1))
      , interface
      , innerEdge
    );

    darcy.addLinear(
      -innerProduct(p0.expression(), v_hat*t)
      , interface
      , boundary
      , {1,3}
    );

    // Full/macro Stabilization
    // darcy.addFaceStabilization( // [h^(2k+1) h^(2k)]
    //     innerProduct(uPenParam*pow(h,3)*jump(grad(u)), jump(grad(v)))
    //    +innerProduct(uPenParam*h*jump(u*t), jump(v*t)) // [Method 1: Remove jump in vel]
    //    +innerProduct(pPenParam*pow(h,2)*jump(grad(p)), jump(grad(q)))
    //    +innerProduct(pPenParam*jump(p), jump(q))
    //    , macro
    //  );


     // Full stabilization (3.22)
     // darcy.addBilinear(
     //           innerProduct(sPenParam*hei*jump(u_hat*t), jump(v_hat*t))
     //         + innerProduct(sPenParam*jump(grad(u_hat)), jump(grad(v_hat)))
     //         , innerEdge
     //         // , macroInterface
     // );
     // darcy.addBilinear(
     //         innerProduct(sPenParam*pow(h,2)*grad(u_hat)*n, grad(v_hat)*n)
     //         , interface
     //         , macroInterface
     // );


    // extension.do_extension();
    // extension.solve();
    // matlab::Export(darcy.mat, "matB"+to_string(i)+".dat");

    darcy.solve();
    nx += 2;
    ny += 2;
    // continue;
    // //

    // continue;
    // Rn sigma(N);
    // darcy.solve();


    // std::cout << darcy.rhs << std::endl;

    // int N = darcy.nDoF;
    // Rn tmp(N);
    // // Rn sigma(N);
    // int M = N-extension.dof2rm.size(); // = darcyCutMix.rhs.size()
    // multiply(N,M, extension.Pt, darcy.rhs, tmp);
    // darcy.rhs.resize(N);
    // multiply(N,N, extension.S, tmp, darcy.rhs);


    // const CutFEM_Parameter mu("mu", 2, 2);
    R e_divU = integralCut(
      innerProduct(div(u), div(v))
      , uh, uh
      , Wh
    );
    // const CutFEM_Parameter mu("mu", 2, 2);
    R es_divU = integral(
      innerProduct(div(u_hat), div(v_hat))
      , us, us
      , interface
    );
    std::cout << e_divU<< std::endl;
    std::cout << es_divU<< std::endl;

    // PRINT THE SOLUTION TO PARAVIEW
    // =====================================================
    if(MPIcf::IamMaster()){
      Paraview2 writer(Wh, levelSet, "fracture_darcy"+to_string(i)+".vtk");
      writer.add(uh, "uh", 0, 2);
      writer.add(uh, "ph", 2, 1);


      // Fun_h uexsurf (Sh, fun_uSurface);
      Paraview2 writerU0(Sh, "fracture_surface_darcy"+to_string(i)+".vtk");
      writerU0.add(us, "us", 0, 2);
      writerU0.add(us, "ps", 2, 1);
      writerU0.add(levelSet, "levelSet", 0, 1);

    }
  }
  for(int i=0; i<idx0.size();++i){
    std::cout << idx0[i] << "\t";
  }

}
#endif

#ifdef DARCY_MULTI_FRACTURE
R fun_levelSet(const R2 P) {
  return -(P.y+2*P.x) + 1.4;
  // R2 shift(0.5,-1);
  // return sqrt((P.x-shift.x)*(P.x-shift.x) + (P.y-shift.y)*(P.y-shift.y)) - 1.5;
}
// R fun_neumann(const R2 P, int compInd, int dom) {
//   return P.y;
// }

int main(int argc, char** argv ) {

  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx =10;
  int ny =10;
  Mesh2 Th(nx, ny, 0., 0., 1., 1.);


  Cut_LevelSet L1(fun_levelSet, 1);
  Cut_Line L2(1,-1,0.05, 2);
  Cut_Line L3(0,1,-0.48, 3);

  Fracture fracture(Th);
  fracture.add(L1);
  fracture.add(L2);
  fracture.add(L3);


  gnuplot::save(Th, fracture, "Th_fractured.dat");
  gnuplot::save(fracture, "fracture.dat");



  // FESpace2 Lh(Th, DataFE<Mesh2>::P1);

/*
  Fun_h levelSet(Lh, fun_levelSet);

  Interface2 interface(Th, levelSet.v);

  KN<const GTypeOfFE<Mesh2>* > arrayFE(2);
  arrayFE(0) = &DataFE<Mesh2>::RT0;
  arrayFE(1) = &DataFE<Mesh2>::P0;
  GTypeOfFESum<Mesh2> mixedFE(arrayFE);
  FESpace2 Vh(Th, mixedFE);

  FESpace2 Ph (Th, DataFE<Mesh2>::P2);

  // CutFEM stuff
  CutFESpace2 Wh(Vh , interface, {1,-1});
  CutFESpace2 CPh(Ph, interface, {1,-1});

  Mesh2 cutTh(interface);
  FESpace2 Sh(cutTh, mixedFE);
  Sh.backSpace = &Vh;


  // OBJECTS NEEDED FOR THE PROBLEM
  // =====================================================
  CutFEM<Mesh2> darcy(Wh);
  darcy.add(Sh);

  int idx_s0 = Wh.NbDoF();
  std::cout << "idx_s0 = " << idx_s0 << std::endl;


  KN_<double> data_uh(darcy.rhs(SubArray(Wh.NbDoF(),0)));
  KN_<double> data_sh(darcy.rhs(SubArray(Sh.NbDoF(),idx_s0)));
  Fun_h uh(Wh, data_uh);  // fun representing bulk variables  [v1, v2]
  Fun_h us(Sh, data_sh);  // fun representing surface variable v0

  Normal n;
  Tangent t;
  FunTest p(Wh,1,2), q(Wh,1,2), u(Wh,2), v(Wh,2);
  FunTest p_hat(Sh,1,2), q_hat(Sh,1,2), u_hat(Sh,2), v_hat(Sh,2);

  // We define fh on the cutSpace
  // Fun_h fv(Wh, fun_force);
  // Fun_h fun_v_hat(Wh, fun_force);
  // Fun_h fq(CPh, fun_div);
  // Fun_h u0(Wh, fun_dirichlet);
  Fun_h p0(Ph, fun_neumann);
  // Fun_h phat(Wh, fun_interfacePr);

  double nu = 1;
  double nu_G = 0.01;
  double nu_hat = 100;
  double gamma_D = 10;
  const CutFEM_Parameter& invh(Parameter::invh);
  double xi = 3./4;
  double xi0 = (xi-0.5)/2.;
  double l_G = 0.01;


  // BULK PROBLEM
  // =====================================================
  darcy.addBilinear(
     innerProduct(nu*u, v)
    -innerProduct(p, div(v))
    +innerProduct(div(u), q)
  );

  darcy.addBilinear(
    innerProduct( gamma_D*invh*(u*n), (v*n))
    + innerProduct(p, v*n) - innerProduct(u*n, q)
    , boundary
    , {2,4}
  );

  darcy.addBilinear(
    innerProduct(     nu_G*average(u.t()*n), average(v*n))
    +innerProduct(xi0*nu_G*jump(u*n)   , jump(v*n))
    ,interface
  );

  // l(v)_Omega
  darcy.addLinear(
    // innerProduct(fv, v)  // fv = ?????
     innerProduct(4, q)  // fq = 4
  );
  // FunTest fphat(phat, 2, 1);
  // darcy.addBilinear(
  //   innerProduct(p_hat, jump(v*n))
  //   ,interface
  // );
  darcy.addLinear(
    -innerProduct(p0.expression(), jump(v*n))
    ,interface
  );
  darcy.addLinear(
    -innerProduct(p0.expression(), v*n)
    , boundary
    , {1,3}
  );
  // darcy.addLinear(
  //   -innerProduct(u0.expression(), q)
  //   +innerProduct(gamma_D*invh*u0.expression(), (v*n))
  //   , boundary
  //   , {2,4}
  // );

// std::cout << divS(v_hat) << std::endl;
// std::cout << divT(v_hat) << std::endl;
// getchar();
  // SURFACE PROBLEM
  // =====================================================
  darcy.addBilinear(
     innerProduct(nu_hat*u_hat, v_hat)
     -innerProduct(p_hat, divS(v_hat))
     +innerProduct(divS(u_hat), q_hat)
     , interface
   );
   darcy.addBilinear(
     - innerProduct(jump(u*n), q_hat)
     , interface
   );
   // FunTest fv_hat(fun_v_hat, 0, 2);
   // FunTest fq_tild();
   darcy.addLinear(
     // innerProduct(fv_hat, v_hat)
       // innerProduct(l_G*fq_tild, q_hat)
       innerProduct(l_G*4, q_hat)
     , interface
   );

   darcy.addLinear(
     -innerProduct(p0.expression(), v_hat*n)
     , interface
     , boundary
     , {1,3}
   );

   // matlab::Export(darcy.mat, "coupledMat.dat");

   darcy.solve();



  // PRINT THE SOLUTION TO PARAVIEW
  // =====================================================
  if(MPIcf::IamMaster()){
      Paraview2 writer(Wh, levelSet, "fracture_darcy.vtk");
      writer.add(uh, "uh", 0, 2);
      writer.add(uh, "ph", 2, 1);


      // Fun_h uexsurf (Sh, fun_uSurface);
      Paraview2 writerU0(Sh, "fracture_surface_darcy.vtk");
      writerU0.add(us, "us", 0, 2);
      writerU0.add(us, "ps", 2, 1);
      writerU0.add(levelSet, "levelSet", 0, 1);

    }
*/
}
#endif
