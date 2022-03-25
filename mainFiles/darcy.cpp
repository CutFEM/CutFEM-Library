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
#include "projection.hpp"
#include "../common/fracture.hpp"
#include "extension.hpp"
#include "../num/gnuplot.hpp"
#include "generalNorm.hpp"


#define PROBLEM_CUT_MIXED_DARCY
// #define PROBLEM_CUT_MIXED_DARCY_BLOCK
// #define DARCY_EXAMPLE_PUPPI
// #define DARCY_FRACTURE
// #define DARCY_MULTI_FRACTURE
// #define PROBLEM_RT_PROJECTION
// #define PROBLEM_MIXED_DARCY


void Scotty_diagonal_preconditioner(int N, std::map<std::pair<int,int>,R>& P){

  SparseMatrixRC<double> B (N ,N ,P );
  P.clear();

  // create the diagonal Matrix
  for(int i=0;i<B.n;++i){
    for(int k=B.p[i];k<B.p[i+1];++k){
      P[std::make_pair(i,i)] += B.a[k];
    }
  }
  for(int i=0;i<B.n;++i){
    P[std::make_pair(i,i)] = 1./P[std::make_pair(i,i)];
  }

}



#ifdef PROBLEM_CUT_MIXED_DARCY

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
  R fun_exact(const R2 P, int compInd, int dom) {
    Diff<R,2> X(P.x,0), Y(P.y,1);
    Diff<R,2> r2 = (X-shift)*(X-shift) + (Y-shift)*(Y-shift);
    R radius2 = interfaceRad*interfaceRad;
    R cst = (dom==0)*3./2;
    R mul = (dom==0)*2 + (dom==1)*1;
    Diff<R, 2> val = r2/(mul*radius2) + cst;
    if (compInd==2) {
      return val.val;
    }
    else {
      return -val.d[compInd];
    }
  }
  R fun_interfacePr(const R2 P, int compInd) {
    if (compInd == 2)
    return 19./12;
    else
    return 0;
  }
}
namespace Data_CutMixedDarcyNOJUMP {
  bool solHasJump = false;

  R interfaceRad = 0.250001; // not exactly 1/4 to avoid bad cuts
  R shift = 0.5;

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
    if (compInd == 2) {
      R r2 = fun_radius2(P);
      R radius2 = interfaceRad*interfaceRad;
      return r2/(2*radius2)+3./2.;
    }
    else
    return 0;
  }

  R fun_force(const R2 P, int compInd) {
    return 0;
  }
  R fun_div(const R2 P, int compInd, int dom) { // is also exact divergence
    R radius2 = interfaceRad*interfaceRad;
    return -2./radius2;
  }
  R fun_exact(const R2 P, int compInd, int dom) {
    Diff<R,2> X(P.x,0), Y(P.y,1);
    Diff<R,2> r2 = (X-shift)*(X-shift) + (Y-shift)*(Y-shift);
    R radius2 = interfaceRad*interfaceRad;
    Diff<R, 2> val = r2/(2*radius2) + 3./2;
    if (compInd==2) {
      return val.val;
    }
    else {
      return -val.d[compInd];
    }
  }
  R fun_interfacePr(const R2 P, int compInd) {
    return fun_exact(P,2,0);
  }
}

namespace Data_Puppy{
  R d_x = 2.;
  R d_y = 2.;
  R shift = 0;
  R interfaceRad = 0.501;
  R mu_G = 2./3*interfaceRad; // xi0*mu_G = 1/8*2/3*1/4

  R fun_radius2(const R2 P){
    return (P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift);
  }
  R fun_levelSet(const R2 P, const int i) {
    return P.y-1.051;
  }

  // R fun_neumann_east(const R2 P, int compInd) {
  //   return -sin(P.x)*sinh(P.y) - (cos(1) - 1)*(cosh(1) - 1);
  // }
  R fun_neumann(const R2 P, int compInd) {
    return -sin(P.x)*sinh(P.y) - (cos(1) - 1)*(cosh(1) - 1);
  }



  R fun_force(const R2 P, int compInd) {
    return 0;
  }
  R fun_div(const R2 P, int compInd, int dom) {// is also exact divergence
    return 0;
  }
  R fun_exact(const R2 P, int compInd, int dom) {
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
    return 0;//19./12;
    else
    return 0;
  }

}
using namespace Data_CutMixedDarcy;

static void exactRHSintegration(const CutFESpace2& Vh, Rn& rhs, R (*f)(const R2)){
  typedef typename Mesh2::BorderElement BorderElement;
  typedef typename FESpace2::FElement FElement;
  typedef typename Mesh2::Element Element;
  typedef typename FElement::QFB QFB;
  typedef typename QFB::QuadraturePoint QuadraturePoint;


  KNMK<double> fv(Vh[0].NbDoF(),3,1); //  the value for basic fonction
  What_d Fop = Fwhatd(1);
  const QFB &qfb(*QF_Simplex<R1>(5));

  for( int ifac = Vh.first_boundary_element(); ifac < Vh.last_boundary_element(); ifac+=Vh.next_boundary_element()) {

    const BorderElement & BE(Vh.Th.be(ifac)); // The face
    int ifaceK; // index of face of triangle corresp to edge (0,1,2)
    const int kb = Vh.Th.BoundaryElement(ifac, ifaceK); // index of element (triangle), ifaceK gets modified inside
    const int k = Vh.idxElementFromBackMesh(kb, 0);

    const FElement& FK((Vh)[k]);
    R2 normal = FK.T.N(ifaceK);
    const R meas = BE.mesure();

    int nb_face_onB = 0;
    for(int i=0;i<Element::nea;++i){
      int ib = i;
      if(Vh.Th.ElementAdj(kb,ib) == -1) nb_face_onB += 1;
    }
    assert(nb_face_onB > 0);
    double measOnB = nb_face_onB*meas;
    const int kv = Vh.idxElementFromBackMesh(kb, 0);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const R2 mip = BE((R1)ip); // mip is in the global edge
      const R Cint = meas * ip.getWeight();

      FK.BF(Fop,FK.T.toKref(mip), fv);
      double val_fh = f(mip);

      for(int d=0;d<2;++d){
        for(int i = FK.dfcbegin(d); i < FK.dfcend(d); ++i) {
          rhs(FK.loc2glb(i)) +=  -Cint * val_fh * fv(i,d,op_id)*normal[d];
        }
      }
    }
  }
}

int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 10; // 6
  int ny = 10; // 6

  vector<double> uPrint,pPrint,divPrint,divPrintLoc,maxDivPrint,h,convuPr,convpPr,convdivPr,convdivPrLoc,convmaxdivPr;
  vector<double> ratioCut1, ratioCut2;
  int iters =5;

  for(int i=0;i<iters;++i) {
    Mesh2 Th(nx, ny, 0., 0., 1., 1.);
    // Mesh2 Th(nx, ny, 0., 0., 2., 2.);
    Th.info();
    std::cout << "nx = \t" << nx << std::endl;
    // Mesh2 Th("../mesh/RTmesh20.msh");
    FESpace2 Lh(Th, DataFE<Mesh2>::P1);
    Fun_h levelSet(Lh, fun_levelSet);

    Interface2 interface(Th, levelSet.v);

    KN<const GTypeOfFE<Mesh2>* > arrayFE(2);
    arrayFE(0) = &DataFE<Mesh2>::RT1;
    arrayFE(1) = &DataFE<Mesh2>::P1dc;
    GTypeOfFESum<Mesh2> mixedFE(arrayFE);
    FESpace2 Vh(Th, mixedFE);
    FESpace2 P0h(Th, DataFE<Mesh2>::P1); // for the RHS
    FESpace2 Ph (Th, DataFE<Mesh2>::P1); // for the RHS

    // CutFEM stuff
    CutFESpace2 Wh(Vh , interface, {1,-1});
    CutFESpace2 CPh(Ph, interface, {1,-1});

    // gnuplot::save(Wh,0,"Vh1.dat");
    // gnuplot::save(Wh,1,"Vh2.dat");
    CutFEM<Mesh2> darcyCutMix(Wh);
    const R hei = 1./(nx-1);//Th[0].lenEdge(0);

    MacroElement macro(Wh, 1e-2);
    // Extension extension(darcyCutMix);
    // extension.tag_extension_edges(macro);
    // // extension.tag_extension_edges(macro, innerEdge);
    // // extension.tag_extension_edges(macro, innerEdge, boundary);
    // extension.tag_exhaust_edges(macro);
    ratioCut1.push_back((double)macro.nb_element_0 / (double)interface.nbElement());
    ratioCut2.push_back((double)macro.nb_element_1 / (double)interface.nbElement());
    // gnuplot::save(macro, extension);
    // gnuplot::save(interface);
    // return 0;

    R xi = 3./4; // > 1/2
    R xi0 = (xi-0.5)/2.;
    CutFEM_Parameter mu("mu",1.,1.);
    // const CutFEM_Parameter& invh(Parameter::invh);
    const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
    const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
    // const R hei = Th[0].lenEdge(0);
    const R invh = 1./hei;

    // We define fh on the cutSpace
    Fun_h fv(Wh, fun_force);
    Fun_h fq(CPh, fun_div);
    // Fun_h u0(Wh, fun_dirichlet);
    Fun_h p0(Ph, fun_neumann);
    Fun_h phat(Wh, fun_interfacePr);

    Normal n;
    Tangent t;
    FunTest p(Wh,1,2), q(Wh,1,2), u(Wh,2), v(Wh,2);

    // u=1e-2 p=1e1 REALLY GOOD for jump problem!
    double uPenParam = 1e0;//1e-2; //cont 1e-1`
    double pPenParam = 1e0;//1e1; // cont 1e2
    double jumpParam = 1e0; // [anything<1e0 doesn't give full u convergence]

    double theta = 0; // [extension param] 1e-2 good

    // [ASSEMBLY]

    //:: a(u,v)_Omega
    darcyCutMix.addBilinear(
      innerProduct(mu*u, v)
      -innerProduct(p, div(v))
      +innerProduct(div(u), q)
    );

    darcyCutMix.addBilinear(
      innerProduct((1-theta)*mu_G*average(u.t()*n), average(v.t()*n))
      +innerProduct((1-theta)*xi0*mu_G*jump(u.t()*n), jump(v.t()*n)) // b(p,v)-b(q,u) bdry terms
      ,interface
    );

    // l(v)_Omega
    darcyCutMix.addLinear(
      innerProduct(fq.expression(), q)
    );

    // darcyCutMix.addBilinearFormExtDomain(
    //   innerProduct(div(u), q)
    //   ,macro,1
    // );
    // darcyCutMix.addLinearFormExtDomain(
    //   innerProduct(fq.expression(), q)
    //   ,macro,1
    // );


    darcyCutMix.addLinear(
      -innerProduct(p0.expression(), v*n) // Only on Gamma_N (pressure)
      // + innerProduct(u0, invh*(v.t()*n)*lambdaB) // Only on Gamma_D (vel normal comp)
      // - innerProduct(u0, q) // Only on Gamma_D (vel normal comp)
      , boundary
    );

    // ExpressionFunFEM<Mesh2> pphat(phat, 2, op_id,0,0);
    FunTest ppphat(phat, 2, 1);
    // FunTest ppphat(CPh,pphat);
    darcyCutMix.addLinear(
      -innerProduct((1-theta)*ppphat, jump(v*n))
      // -innerProduct(phat.expression(), jump(v.t()*n))
      ,interface
    );

    int N = Wh.NbDoF();



    {  // Diagonal Preconditionning

      std::map<std::pair<int,int>,R> P;
      darcyCutMix.setMatrixTo(P);
      darcyCutMix.addBilinear(
        innerProduct(mu*u, v)
        + innerProduct(p, q)
      );
      darcyCutMix.addBilinear(
        innerProduct((1-theta)*mu_G*average(u*n), average(v*n))
        +innerProduct((1-theta)*xi0*mu_G*jump(u*n), jump(v*n)) // b(p,v)-b(q,u) bdry terms
        ,interface
      );

      Scotty_diagonal_preconditioner(N, P);

      // darcyCutMix.preconditionning(P);

      darcyCutMix.pmat = &darcyCutMix.DF;

      // matlab::Export(darcyCutMix.DF, "matP"+to_string(i)+".dat");
      // matlab::Export(P, "matA"+to_string(i)+".dat");
      // // // darcyCutMix.solve();
      // return 0;
    }

  // matlab::Export(darcyCutMix.mat, "matA"+to_string(i)+".dat");

  // Full/macro Stabilization
  FunTest grad2un = grad(grad(u)*n)*n;
  darcyCutMix.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
  // innerProduct(uPenParam*hei*jump(u*t), jump(v*t)) // [Method 1: Remove jump in vel]
  // +innerProduct(uPenParam*pow(hei,3)*jump(grad(u)*n), jump(grad(v)*n))
  // // +innerProduct(uPenParam*pow(hei,5)*jump(grad2un), jump(grad2un))
  // +innerProduct(pPenParam*hei*jump(p), jump(q))
  // +innerProduct(pPenParam*pow(hei,3)*jump(grad(p)), jump(grad(q)))

   innerProduct(uPenParam*hei*jump(u*t), jump(v*t)) // [Method 1: Remove jump in vel]
  +innerProduct(uPenParam*pow(hei,3)*jump(grad(u)*n), jump(grad(v)*n))
  +innerProduct(uPenParam*pow(hei,5)*jump(grad2un), jump(grad2un))
  +innerProduct(pPenParam*hei*jump(p), jump(div(v)))
  +innerProduct(pPenParam*hei*jump(div(u)), jump(q))
  +innerProduct(pPenParam*pow(hei,3)*jump(grad(p)), jump(grad(div(v))))
  +innerProduct(pPenParam*pow(hei,3)*jump(grad(div(v))) , jump(grad(q)))
  , macro
);

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


  // matlab::Export(darcyCutMix.mat, "matA.dat");
  // matlab::Export(darcyCutMix.mat, "matA"+to_string(i)+".dat");

  // int N = Wh.NbDoF();
  // std::map<std::pair<int,int>,double> Res;
  // multiply(N, macro.St, darcyCutMix.mat, Res);
  // matlab::Export(Res, "matStA.dat");
  // extension.do_extension();


  // matlab::Export(darcyCutMix.mat, "matB"+to_string(i)+".dat");
  // nx += 100;
  // ny += 100;
  // continue;
  // matlab::Export(darcyCutMix.mat, "matA"+to_string(i)+".dat");
  // extension.erase_rows_to_fix_RT0();
  // if(MPIcf::IamMaster()){
  //   std::cout << "dof = \t" << Wh.NbDoF()<< std::endl;
  // matlab::Export(darcyCutMix.mat, "matB"+to_string(i)+".dat");
  // nx *= 2;
  // ny *= 2;
  // continue;

  // return 0;

  // extension.make_S();
  // Rn b(darcyCutMix.rhs.size());
  // extension.precond(b);
  // matlab::Export(darcyCutMix.mat, "matB"+to_string(i)+".dat");
  // continue;
  // extension.do_extension();
  // extension.solve_weak_extension("mumps");
  // extension.make_S();
  // extension.solve_weak_extension("mumps");

  // matlab::Export(darcyCutMix.mat, "matB"+to_string(i)+".dat");
  // { // EXPORT SAVE AND LOAD WITH MATLAB
  //   matlab::Export(darcyCutMix.mat, "Amat.dat");
  //   matlab::Export(darcyCutMix.rhs, "rhs.dat");
  //
  //   std::cout << " Export matrix and rhs " << darcyCutMix.rhs.size() << std::endl;
  //   getchar();
  //   matlab::loadVector(darcyCutMix.rhs, darcyCutMix.rhs.size(),  "rhs.txt");
  //   std::cout << " Vector loaded "<< darcyCutMix.rhs.size() << std::endl;
  // }
  // darcyCutMix.solve("umfpack");
  // darcyCutMix.solve();



  // Rn tmp(N);
  // Rn sigma(N);
  // int M = N-macro.dof2rm.size(); // = darcyCutMix.rhs.size()
  // multiply(N,M, macro.Pt, darcyCutMix.rhs, tmp);
  // multiply(N,N, macro.S, tmp, sigma);
  // Fun_h femSolh(Wh, sigma);
  Fun_h femSolh(Wh, darcyCutMix.rhs);

  // macro.imposeDirichlet(darcyCutMix.mat,darcyCutMix.rhs);

  // matlab::Export(macro.P, "matP.dat");

  // Fun_h femSolh(Wh, tmp);
  // darcyCutMix.solve();
  //
  // matlab::Export(darcyCutMix.rhs, "myRHS.dat");


  // L2 norm vel
  R errU = L2normCut(femSolh,fun_exact,0,2);

  // R errP = L2normCut(femSolh,fun_exact,0,2, &macro);//L2normCut(femSolh,fun_exact,2,1);
  R errP = L2normCut(femSolh,fun_exact,2,1);

  // L2 norm div
  ExpressionFunFEM<Mesh2> femSol_0dx(femSolh, 0, op_dx);
  ExpressionFunFEM<Mesh2> femSol_1dy(femSolh, 1, op_dy);
  R errDiv = L2normCut(femSol_0dx+femSol_1dy,fun_div,Wh);

  // std::cout << errDiv << std::endl;

  // R errDivLoc = L2normCut(femSol_0dx+femSol_1dy,fun_div,Wh,&macro);//L2normCutLoc(femSol_0dx+femSol_1dy,fun_div,Wh,macro);

  // Max error div calculation
  R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy,fun_div,Wh);

  // [PLOTTING]
  {
    Paraview2 writerS(Wh, levelSet, "darcyPaperPlotExample1_"+to_string(i)+".vtk");
    writerS.add(femSolh, "velocity" , 0, 2);
    // writerS.add(femSolh, "velocityExt", 0, 2, macro);
    writerS.add(femSolh, "pressure" , 2, 1);
    writerS.add(femSol_0dx+femSol_1dy, "divergence");
    // writerS.add(femSol_0dx+femSol_1dy, "divergenceExt", macro);
    Fun_h solh(Wh, fun_exact);

    // [For looking at the error]
    solh.v -= femSolh.v;
    solh.v.map(fabs);
    writerS.add(solh, "velocityError" , 0, 2);

    Fun_h divSolh(Wh, fun_div);
    ExpressionFunFEM<Mesh2> femDiv(divSolh, 0, op_id);

      // Paraview2 writerS2(Wh, levelSet, "darcyCutMixedEx"+to_string(i)+".vtk");
    // writerS2.add(solh, "velocity", 0, 2);
    // writerS2.add(solh, "pressure", 2, 1);
    // writerS2.add(divSolh, "divergence",0,1);
    // writerS.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceErrorExt", macro);
    writerS.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceError");

  }

  pPrint.push_back(errP);
  uPrint.push_back(errU);
  // continue;
  divPrint.push_back(errDiv);
  // divPrintLoc.push_back(errDivLoc);

  maxDivPrint.push_back(maxErrDiv);
  h.push_back(1./nx);

  if(i==0) {convpPr.push_back(0);convuPr.push_back(0);convdivPr.push_back(0);convdivPrLoc.push_back(0);convmaxdivPr.push_back(0);}
  else {
    convpPr.push_back( log(pPrint[i]/pPrint[i-1])/log(h[i]/h[i-1]));
    convuPr.push_back( log(uPrint[i]/uPrint[i-1])/log(h[i]/h[i-1]));
    convdivPr.push_back( log(divPrint[i]/divPrint[i-1])/log(h[i]/h[i-1]));
    // convdivPrLoc.push_back( log(divPrintLoc[i]/divPrintLoc[i-1])/log(h[i]/h[i-1]));

    convmaxdivPr.push_back( log(maxDivPrint[i]/maxDivPrint[i-1])/log(h[i]/h[i-1]));
  }
  // nx += 100;
  // ny += 100;
  // nx = 2*nx - 1;
  // ny = 2*ny - 1;
  // int nn = iters/5;
  nx *= 2;
  ny *= 2;
  // nx = (int)round( (1+0.2*i)*nx/2 )*2; // Makes a nonuniform refinement to an EVEN integer
  // ny = (int)round( (1+0.2*i)*ny/2 )*2;
  // std::cout << nx << std::endl;
  // shift = 0.5+(i+1)*hei/iters; // moves one grid cell over entire span

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

// FunTest ff(femSolh, 0, 2);
// R ee = integralCut(
//   innerProduct(div(u), div(v))
//   , femSolh, femSolh
//   , Wh
// );
// std::cout << ee << std::endl;
}

#endif

#ifdef PROBLEM_CUT_MIXED_DARCY_BLOCK

namespace Data_CutMixedDarcy {
  bool solHasJump = true;
  R shift = 0.5;
  R interfaceRad = 0.2501; // not exactly 1/4 to avoid interface cutting exaclty a vertex
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
  R fun_exactU(const R2 P, int compInd, int dom) {
    Diff<R,2> X(P.x,0), Y(P.y,1);
    Diff<R,2> r2 = (X-shift)*(X-shift) + (Y-shift)*(Y-shift);
    R radius2 = interfaceRad*interfaceRad;
    R cst = (dom==0)*3./2;
    R mul = (dom==0)*2 + (dom==1)*1;
    Diff<R, 2> val = r2/(mul*radius2) + cst;
    return -val.d[compInd];
  }
  R fun_exactP(const R2 P, int compInd, int dom) {
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

static void exactRHSintegration(const CutFESpace2& Vh, Rn& rhs, R (*f)(const R2)){
  typedef typename Mesh2::BorderElement BorderElement;
  typedef typename FESpace2::FElement FElement;
  typedef typename Mesh2::Element Element;
  typedef typename FElement::QFB QFB;
  typedef typename QFB::QuadraturePoint QuadraturePoint;


  KNMK<double> fv(Vh[0].NbDoF(),3,1); //  the value for basic fonction
  What_d Fop = Fwhatd(1);
  const QFB &qfb(*QF_Simplex<R1>(5));

  for( int ifac = Vh.first_boundary_element(); ifac < Vh.last_boundary_element(); ifac+=Vh.next_boundary_element()) {

    const BorderElement & BE(Vh.Th.be(ifac)); // The face
    int ifaceK; // index of face of triangle corresp to edge (0,1,2)
    const int kb = Vh.Th.BoundaryElement(ifac, ifaceK); // index of element (triangle), ifaceK gets modified inside
    const int k = Vh.idxElementFromBackMesh(kb, 0);

    const FElement& FK((Vh)[k]);
    R2 normal = FK.T.N(ifaceK);
    const R meas = BE.mesure();

    int nb_face_onB = 0;
    for(int i=0;i<Element::nea;++i){
      int ib = i;
      if(Vh.Th.ElementAdj(kb,ib) == -1) nb_face_onB += 1;
    }
    assert(nb_face_onB > 0);
    double measOnB = nb_face_onB*meas;
    const int kv = Vh.idxElementFromBackMesh(kb, 0);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const R2 mip = BE((R1)ip); // mip is in the global edge
      const R Cint = meas * ip.getWeight();

      FK.BF(Fop,FK.T.toKref(mip), fv);
      double val_fh = f(mip);

      for(int d=0;d<2;++d){
        for(int i = FK.dfcbegin(d); i < FK.dfcend(d); ++i) {
          rhs(FK.loc2glb(i)) +=  -Cint * val_fh * fv(i,d,op_id)*normal[d];
        }
      }
    }
  }
}

int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx =8; // 6
  int ny =8; // 6

  vector<double> uPrint,pPrint,divPrint,divPrintLoc,maxDivPrint,h,convuPr,convpPr,convdivPr,convdivPrLoc,convmaxdivPr;

  int iters = 10;

  for(int i=0;i<iters;++i) {
    Mesh2 Th(nx, ny, 0., 0., 1., 1.);
    // Mesh2 Th("../mesh/RTmesh20.msh");
    FESpace2 Lh(Th, DataFE<Mesh2>::P1);
    Fun_h levelSet(Lh, fun_levelSet);

    Interface2 interface(Th, levelSet.v);

    KN<const GTypeOfFE<Mesh2>* > arrayFE(2);
    arrayFE(0) = &DataFE<Mesh2>::RT0;
    arrayFE(1) = &DataFE<Mesh2>::P0;
    GTypeOfFESum<Mesh2> mixedFE(arrayFE);

    FESpace2 Vh(Th, DataFE<Mesh2>::RT0);
    FESpace2 Ph(Th, DataFE<Mesh2>::P0) ; // for the RHS


    FESpace2 P2h (Th, DataFE<Mesh2>::P2); // for the RHS

    // CutFEM stuff
    CutFESpace2 Wh(Vh , interface, {1,-1});
    CutFESpace2 WPh(Ph , interface, {1,-1});
    CutFESpace2 CPh(P2h, interface, {1,-1});

    // macro.tag_extension_edges();
    // // macro.tag_exhaust_edges();

    // gnuplot::save(Th);
    // gnuplot::save(interface);
    // gnuplot::save(macro);
    // gnuplot::save(Wh,0,"Vh1.dat");
    // gnuplot::save(Wh,1,"Vh2.dat");

    // getchar();
    CutFEM<Mesh2> darcyCutMix(Wh);
    darcyCutMix.add(WPh);

    MacroElement macro(Wh,  1e-2);

    Extension extension(darcyCutMix);
    extension.tag_extension_edges(macro, Wh);
    // extension.tag_extension_edges(macro, WPh);

    // gnuplot::save(macroU, extension);
    // gnuplot::save(interface);


    R xi = 3./4; // > 1/2
    R xi0 = (xi-0.5)/2.;
    CutFEM_Parameter mu("mu",1.,1.);
    // const CutFEM_Parameter& invh(Parameter::invh);
    const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
    const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
    const R hei = Th[0].lenEdge(0);
    const R invh = 1./hei;

    // We define fh on the cutSpace
    Fun_h fv(Wh, fun_force);
    Fun_h fq(CPh, fun_div);
    Fun_h u0(Wh, fun_dirichlet);
    Fun_h p0(P2h, fun_neumann);
    Fun_h phat(WPh, fun_interfacePr);


    int idx_p0 = Wh.NbDoF();
    std::cout << idx_p0 << std::endl;
    KN_<double> data_uh(darcyCutMix.rhs(SubArray(Wh.NbDoF(),0)));
    KN_<double> data_ph(darcyCutMix.rhs(SubArray(WPh.NbDoF(),idx_p0)));
    Fun_h uh(Wh, data_uh);
    Fun_h ph(WPh, data_ph);


    Normal n;
    Tangent t;
    FunTest p(WPh,1,0), q(WPh,1,0), u(Wh,2), v(Wh,2);

    // u=1e-2 p=1e1 REALLY GOOD for jump problem!
    double uPenParam = 1e-1;//1e-2; //cont 1e-1`
    double pPenParam = 1e1;//1e1; // cont 1e2
    double jumpParam = 1e0; // [anything<1e0 doesn't give full u convergence]

    double theta = 0; // [extension param] 1e-2 good

    // [ASSEMBLY]

    //:: a(u,v)_Omega
    darcyCutMix.addBilinearFormExtDomain(
      innerProduct(mu*u, v)
      -innerProduct(p, div(v))
      +innerProduct(div(u), q)
      ,theta);

      darcyCutMix.addBilinear(
         innerProduct((1-theta)*mu_G*average(u.t()*n), average(v.t()*n))
        +innerProduct((1-theta)*xi0*mu_G*jump(u.t()*n), jump(v.t()*n)) // b(p,v)-b(q,u) bdry terms
        ,interface
      );

      // l(v)_Omega
      darcyCutMix.addLinearFormExtDomain(
        innerProduct(fq.expression(), q)
        ,theta);

        darcyCutMix.addLinear(
          -innerProduct(p0.expression(), v*n) // Only on Gamma_N (pressure)
          // + innerProduct(u0, invh*(v.t()*n)*lambdaB) // Only on Gamma_D (vel normal comp)
          // - innerProduct(u0, q) // Only on Gamma_D (vel normal comp)
          , boundary
        );

        // ExpressionFunFEM<Mesh2> pphat(phat, 2, op_id,0,0);
        FunTest ppphat(phat, 0, 1);
        // FunTest ppphat(CPh,pphat);
        darcyCutMix.addLinear(
          -innerProduct((1-theta)*ppphat, jump(v.t()*n))
          // -innerProduct(phat.expression(), jump(v.t()*n))
          ,interface
        );

        int N = Wh.NbDoF();
        /*
        {
        darcyCutMix.saveMatrix();
        darcyCutMix.addBilinear(
        innerProduct(mu*u, v)
        + innerProduct(p, q)
      );
      darcyCutMix.addBilinear(
      innerProduct((1-theta)*mu_G*average(u.t()*n), average(v.t()*n))
      +innerProduct((1-theta)*xi0*mu_G*jump(u.t()*n), jump(v.t()*n)) // b(p,v)-b(q,u) bdry terms
      ,interface
    );

    macro.precondDiag(N, darcyCutMix.NL, darcyCutMix.DF, darcyCutMix.rhs);

    darcyCutMix.pmat = &darcyCutMix.DF;
    matlab::Export(darcyCutMix.mat, "matA"+to_string(i)+".dat");
    darcyCutMix.solve();
  }
  */

  // matlab::Export(darcyCutMix.mat, "matA"+to_string(i)+".dat");

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



  // extension.do_extension();
  // matlab::Export(darcyCutMix.mat, "matB"+to_string(i)+".dat");
  // matlab::Export(darcyCutMix.mat, "matBlock.dat");
  // int N = Wh.NbDoF();
  // std::map<std::pair<int,int>,double> Res;
  // multiply(N, macro.St, darcyCutMix.mat, Res);
  // matlab::Export(Res, "matStA.dat");
  // macro.make_S();
  // macro.precond(darcyCutMix.mat,darcyCutMix.rhs);  // changing indices
  // matlab::Export(darcyCutMix.mat, "matB"+to_string(i)+".dat");


  extension.solve();
  // darcyCutMix.solve();

  // Rn tmp(N);
  // Rn sigma(N);
  // int M = N-macro.dof2rm.size(); // = darcyCutMix.rhs.size()
  // multiply(N,M, macro.Pt, darcyCutMix.rhs, tmp);
  // multiply(N,N, macro.S, tmp, sigma);
  // darcyCutMix.rhs = sigma;
  // Fun_h femSolh(Wh, sigma);
  // Fun_h femSolh(Wh, darcyCutMix.rhs);

  // macro.imposeDirichlet(darcyCutMix.mat,darcyCutMix.rhs);

  // matlab::Export(macro.P, "matP.dat");

  // Fun_h femSolh(Wh, tmp);
  // darcyCutMix.solve();
  //
  // matlab::Export(darcyCutMix.rhs, "myRHS.dat");


  // L2 norm vel
  R errU = L2normCut(uh,fun_exactU,0,2);
  // R errP = L2normCut(femSolh,fun_exact,0,2, &macro);//L2normCut(femSolh,fun_exact,2,1);
  R errP = L2normCut(ph,fun_exactP,0,1);

  // L2 norm div
  ExpressionFunFEM<Mesh2> femSol_0dx(uh, 0, op_dx);
  ExpressionFunFEM<Mesh2> femSol_1dy(uh, 1, op_dy);
  R errDiv = L2normCut(femSol_0dx+femSol_1dy,fun_div,Wh);
  // R errDivLoc = L2normCut(femSol_0dx+femSol_1dy,fun_div,Wh,&macro);//L2normCutLoc(femSol_0dx+femSol_1dy,fun_div,Wh,macro);

  // Max error div calculation
  R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy,fun_div,Wh);

  // [PLOTTING]
  // {
  //   Paraview2 writerS(Wh, levelSet, "darcyCutMixedBDM1_"+to_string(i)+".vtk");
  //   writerS.add(femSolh, "velocity" , 0, 2);
  //   writerS.add(femSolh, "velocityExt", 0, 2, macro);
  //   writerS.add(femSolh, "pressure" , 2, 1);
  //   writerS.add(femSol_0dx+femSol_1dy, "divergence");
  //   writerS.add(femSol_0dx+femSol_1dy, "divergenceExt", macro);
  //   Fun_h solh(Wh, fun_exact);
  //   // Fun_h divSolh(Wh, fun_div); // [WARNING PERHAPS WORKED WHEN DIV ONLY IN COMP 2]
  //   // ExpressionFunFEM<Mesh2> divSolEFF(divSolh, 0, op_id);
  //
  //   // [For looking at the error]
  //   solh.v -= femSolh.v;
  //   solh.v.map(fabs);
  //   writerS.add(solh, "velocityError" , 0, 2);
  //
  //   Fun_h divSolh(Wh, fun_div);
  //   ExpressionFunFEM<Mesh2> femDiv(divSolh, 0, op_id);
  //
  //     // Paraview2 writerS2(Wh, levelSet, "darcyCutMixedEx"+to_string(i)+".vtk");
  //   // writerS2.add(solh, "velocity", 0, 2);
  //   // writerS2.add(solh, "pressure", 2, 1);
  //   // writerS2.add(divSolh, "divergence",0,1);
  //   writerS.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceErrorExt", macro);
  //   writerS.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceError");
  //
  // }

  pPrint.push_back(errP);
  uPrint.push_back(errU);
  // continue;
  divPrint.push_back(errDiv);
  // divPrintLoc.push_back(errDivLoc);

  maxDivPrint.push_back(maxErrDiv);
  h.push_back(1./nx);

  if(i==0) {convpPr.push_back(0);convuPr.push_back(0);convdivPr.push_back(0);convdivPrLoc.push_back(0);convmaxdivPr.push_back(0);}
  else {
    convpPr.push_back( log(pPrint[i]/pPrint[i-1])/log(h[i]/h[i-1]));
    convuPr.push_back( log(uPrint[i]/uPrint[i-1])/log(h[i]/h[i-1]));
    convdivPr.push_back( log(divPrint[i]/divPrint[i-1])/log(h[i]/h[i-1]));
    // convdivPrLoc.push_back( log(divPrintLoc[i]/divPrintLoc[i-1])/log(h[i]/h[i-1]));

    convmaxdivPr.push_back( log(maxDivPrint[i]/maxDivPrint[i-1])/log(h[i]/h[i-1]));
  }

  // nx *= 2;
  // ny *= 2;
  nx += 1;
  ny += 1;
  // nx = (int)round( (1+0.2*i)*nx/2 )*2; // Makes a nonuniform refinement to an EVEN integer
  // ny = (int)round( (1+0.2*i)*ny/2 )*2;
  // std::cout << nx << std::endl;
  // shift = 0.5+(i+1)*hei/iters; // moves one grid cell over entire span

}

std::cout << "\n" << std::left
<< std::setw(10) << std::setfill(' ') << "h"
<< std::setw(15) << std::setfill(' ') << "err_new u"
<< std::setw(15) << std::setfill(' ') << "conv p"
<< std::setw(15) << std::setfill(' ') << "err u"
<< std::setw(15) << std::setfill(' ') << "conv u"
<< std::setw(15) << std::setfill(' ') << "err divu"
<< std::setw(15) << std::setfill(' ') << "conv divu"
<< std::setw(15) << std::setfill(' ') << "err_new divu"
<< std::setw(15) << std::setfill(' ') << "convLoc divu"
<< std::setw(15) << std::setfill(' ') << "err maxdivu"
<< std::setw(15) << std::setfill(' ') << "conv maxdivu"
<< "\n" << std::endl;
for(int i=0;i<uPrint.size();++i) {
  std::cout << std::left
  << std::setw(10) << std::setfill(' ') << h[i]
  << std::setw(15) << std::setfill(' ') << pPrint[i]
  << std::setw(15) << std::setfill(' ') << convpPr[i]
  << std::setw(15) << std::setfill(' ') << uPrint[i]
  << std::setw(15) << std::setfill(' ') << convuPr[i]
  << std::setw(15) << std::setfill(' ') << divPrint[i]
  << std::setw(15) << std::setfill(' ') << convdivPr[i]
  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
  << std::setw(15) << std::setfill(' ') << maxDivPrint[i]
  << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
  << std::endl;
}
}

#endif

#ifdef DARCY_EXAMPLE_PUPPI
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

int main(int argc, char** argv ) {
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx =10; // 6
  int ny =10; // 6

  vector<double> uPrint,pPrint,divPrint,divPrintLoc,maxDivPrint,h,convuPr,convpPr,convdivPr,convdivPrLoc,convmaxdivPr;

  int iters = 6;

  for(int i=0;i<iters;++i) {
    Mesh2 Th(nx, ny, 0., 0., d_x, d_y);
    // Mesh2 Th("../mesh/RTmesh20.msh");
    FESpace2 Lh(Th, DataFE<Mesh2>::P1);
    Fun_h levelSet(Lh, fun_levelSet);
    Interface2 interface(Th, levelSet.v);

    KN<const GTypeOfFE<Mesh2>* > arrayFE(2);
    arrayFE(0) = &DataFE<Mesh2>::RT0;
    arrayFE(1) = &DataFE<Mesh2>::P0;
    GTypeOfFESum<Mesh2> mixedFE(arrayFE);
    FESpace2 Vh(Th, mixedFE);
    FESpace2 P0h(Th, DataFE<Mesh2>::P0); // for the RHS
    FESpace2 Ph(Th, DataFE<Mesh2>::P2); // for the RHS

    // // [Artifical interface]
    // FESpace2 Wh(Th, interface, mixedFE);
    // FESpace2 CPh(Th, DataFE<Mesh2>::P2);
    // CutFEM<Mesh2> darcyCutMix(Wh);

    // [CutFEM]
    CutFESpace2 Wh(Vh , interface, {1});
    CutFESpace2 CPh(Ph, interface, {1});
    //
    // MacroElement macro(Wh, 1e-2);
    // macro.tag_extension_edges();
    // // // macro.tag_exhaust_edges();
    //
    // // gnuplot::save(Th);
    // // gnuplot::save(interface);
    // // gnuplot::save(macro);
    // // gnuplot::save(Wh,0,"Vh1.dat");
    // // gnuplot::save(Wh,1,"Vh2.dat");
    //
    CutFEM<Mesh2> darcyCutMix(Wh);
    const R hei = Th[0].lenEdge(0);
    // const R invh = 1./hei;


    // MacroElement macro(Wh, 0.1);
    // Extension extension(darcyCutMix);
    // // extension.tag_extension_edges(macro);
    // // extension.tag_extension_edges(macro, innerEdge);
    // // extension.tag_extension_edges(macro, innerEdge, boundary);
    // extension.tag_exhaust_edges(macro);



    R xi = 3./4; // > 1/2
    R xi0 = (xi-0.5)/2.;
    CutFEM_Parameter mu("mu",1.,1.);
    // const CutFEM_Parameter& h(Parameter::h);
    const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
    const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
    const CutFEM_Parameter& invh(Parameter::invmeas);

    // We define fh on the cutSpace
    Fun_h fv(Wh, fun_force);
    Fun_h fq(CPh, fun_div);


    Fun_h p0(CPh, fun_natural);
    // Fun_h u0W(CPh, fun_dirichlet_west);
    Fun_h u0(Wh, fun_enforced);
    // Fun_h p0(CPh, fun_neumann);
    // Fun_h u0(CPh, fun_dirichlet);

    Fun_h phat(P0h, fun_interfacePr);

    Normal n;
    Tangent t;
    FunTest p(Wh,1,2), q(Wh,1,2), u(Wh,2), v(Wh,2);

    // u=1e-2 p=1e1 REALLY GOOD for jump problem!
    double uPenParam = 1e0;//1e-2; //cont 1e-1`
    double pPenParam = 1e0;//1e1; // cont 1e2
    double jumpParam = 1e0; // [anything<1e0 doesn't give full u convergence]
    double boundaryParam = 1e-3;
    // [ASSEMBLY]
    //:: a(u,v)_Omega
    darcyCutMix.addBilinear(
      innerProduct(mu*u, v)
      -innerProduct(p, div(v))
      +innerProduct(div(u), q)
    );

    // darcyCutMix.addBilinear(
    //   innerProduct(mu_G*average(u*n), average(v*n))
    //   +innerProduct(xi0*mu_G*jump(u*n), jump(v*n)) // b(p,v)-b(q,u) bdry terms
    // ,interface);

    // l(v)_Omega
    darcyCutMix.addLinear(
      innerProduct(fq.expression(), q)
    );

    // darcyCutMix.addBilinearFormExtDomain(
    //   innerProduct(div(u), q)
    //   ,macro,1
    // );
    // darcyCutMix.addLinearFormExtDomain(
    //   innerProduct(fq.expression(), q)
    //   ,macro,1
    // );


    darcyCutMix.addBilinear(
      innerProduct(boundaryParam*u, invh*(v))
      // - innerProduct(u*n,q)
      + innerProduct(p, v*n)
      ,boundary
      ,{1,4});


    darcyCutMix.addBilinear(
      +innerProduct(boundaryParam*u, invh*(v))
      // - innerProduct(u*n,q)
      + innerProduct(p, v*n)
    ,interface);



    // darcyCutMix.addBilinear(
    //   +innerProduct(p, invh*(v*n)) //*lambdaB Only on Gamma_N (vel normal comp)
    //   // -innerProduct(p0E.expression(), q) // Only on Gamma_N (vel normal comp)
    // ,boundary,{3,2});

    darcyCutMix.addLinear(
      -innerProduct(p0.expression(), v*n)
      ,boundary
      ,{2,3}
    );
    darcyCutMix.addLinear(
      +innerProduct(u0.expression(2), invh*(v)*boundaryParam) //*lambdaB Only on Gamma_N (vel normal comp)
      // -innerProduct(u0*n, q)
      ,boundary
      ,{1, 4}
    );
    darcyCutMix.addLinear(
      +innerProduct(u0.expression(2), invh*(v)*boundaryParam) //*lambdaB Only on Gamma_N (vel normal comp)
      // -innerProduct(u0*n, q)
      ,interface
    );

    // Full/macro Stabilization
    // darcyCutMix.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
    //     innerProduct(uPenParam*pow(hei,3)*jump(grad(u)), jump(grad(v)))
    //    +innerProduct(uPenParam*hei*jump(u*t), jump(v*t)) // [Method 1: Remove jump in vel]
    //    +innerProduct(pPenParam*pow(hei,1)*jump(grad(p)), jump(grad(q)))
    //    +innerProduct(pPenParam*invh*jump(p), jump(q))
    //    // , macro
    //  );



    // darcyCutMix.addLinear(
    //   -innerProduct(p0.expression(), v*n) // Gamma_N
    //   +innerProduct(u0.expression(), invh*(v*n)) //*lambdaB Only on Gamma_D (vel normal comp)
    //   -innerProduct(u0.expression(), q) // Only on Gamma_D (vel normal comp)
    // ,boundary);

    // darcyCutMix.addLinear(
    //   -innerProduct(phat.expression(), jump(v*n))
    // ,interface);

    // matlab::Export(darcyCutMix.mat, "matA"+to_string(i)+".dat");
    int N = Wh.NbDoF();

    // macro.precondDiag(N, darcyCutMix.NL, darcyCutMix.DF, darcyCutMix.rhs);
    // darcyCutMix.pmat = &darcyCutMix.DF;
    // matlab::Export(darcyCutMix.mat, "matA"+to_string(i)+".dat");

    // std::map<std::pair<int,int>,double> Res;
    // multiply(N, macro.St, darcyCutMix.mat, Res);
    // matlab::Export(Res, "matStA.dat");

    // macro.make_S();
    // macro.precond(darcyCutMix.mat,darcyCutMix.rhs);  // changing indices
    // extension.erase_rows_to_fix_RT0();

    // matlab::Export(darcyCutMix.mat, "matB"+to_string(i)+".dat");
    // // darcyCutMix.solve();
    //
    // nx*=2;
    // ny*=2;
    // continue;


    darcyCutMix.solve();
    // darcyCutMix.solve("umfpack");

    // Rn tmp(N);
    // Rn sigma(N);
    // int M = N-macro.dof2rm.size(); // = darcyCutMix.rhs.size()
    // multiply(N,M, macro.Pt, darcyCutMix.rhs, tmp);
    // multiply(N,N, macro.S, tmp, sigma);
    // Fun_h femSolh(Wh, sigma);

    Fun_h femSolh(Wh, darcyCutMix.rhs);

    // macro.imposeDirichlet(darcyCutMix.mat,darcyCutMix.rhs);
    // matlab::Export(macro.P, "matP.dat");

    // Fun_h femSolh(Wh, tmp);
    // darcyCutMix.solve();
    // matlab::Export(darcyCutMix.rhs, "myRHS.dat");


    // L2 norm vel
    R errU = L2normCut(femSolh,fun_exact,0,2);

    // FunTest ff(femSolh, 0, 2);
    // R ee = integralCut(
    //   innerProduct(div(u), div(v))
    //   , femSolh, femSolh
    //   , Wh
    // );
    // std::cout << ee << std::endl;

    // R errP = L2normCut(femSolh,fun_exact,0,2, &macro);//L2normCut(femSolh,fun_exact,2,1);
    R errP = L2normCut(femSolh,fun_exact,2,1);

    // L2 norm div
    ExpressionFunFEM<Mesh2> femSol_0dx(femSolh, 0, op_dx);
    ExpressionFunFEM<Mesh2> femSol_1dy(femSolh, 1, op_dy);
    R errDiv = L2normCut(femSol_0dx+femSol_1dy,fun_div,Wh);
    // R errDivLoc = L2normCut(femSol_0dx+femSol_1dy,fun_div,Wh,&macro);//L2normCutLoc(femSol_0dx+femSol_1dy,fun_div,Wh,macro);

    // Max error div calculation
    R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy,fun_div,Wh);

    // [PLOTTING]
    if(MPIcf::IamMaster()){
      Paraview2 writerS(Wh, levelSet, "darcyPuppiCut_"+to_string(i)+".vtk");
      // Paraview2 writerS(Wh, "darcyPuppi_"+to_string(i)+".vtk");

      writerS.add(femSolh, "velocity" , 0, 2);
      // writerS.add(femSolh, "velocityExt", 0, 2, macro);
      writerS.add(femSolh, "pressure" , 2, 1);
      writerS.add(femSol_0dx+femSol_1dy, "divergence");
      // writerS.add(femSol_0dx+femSol_1dy, "divergenceExt", macro);

      Fun_h solh(Wh, fun_exact);
      Fun_h divSolh(Wh, fun_div);

      // [Exact solution]
      // writerS.add(solh, "velocityEx", 0, 2);
      // writerS.add(solh, "pressureEx", 2, 1);
      // writerS.add(divSolh, "divergenceEx",0,1);

      // [For looking at the error]
      solh.v -= femSolh.v;
      solh.v.map(fabs);
      writerS.add(solh, "velocityError" , 0, 2);

      ExpressionFunFEM<Mesh2> femDiv(divSolh, 0, op_id);
      writerS.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceError");
      // writerS.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceErrorExt", macro);

        // Paraview2 writerS2(Wh, levelSet, "darcyCutMixedEx"+to_string(i)+".vtk");
      // writerS2.add(solh, "velocity", 0, 2);
      // writerS2.add(solh, "pressure", 2, 1);
      // writerS2.add(divSolh, "divergence",0,1);
    }

    pPrint.push_back(errP);
    uPrint.push_back(errU);
    // continue;
    divPrint.push_back(errDiv);
    // divPrintLoc.push_back(errDivLoc);

    maxDivPrint.push_back(maxErrDiv);
    h.push_back(1./nx);

    if(i==0) {convpPr.push_back(0);convuPr.push_back(0);convdivPr.push_back(0);convdivPrLoc.push_back(0);convmaxdivPr.push_back(0);}
    else {
      convpPr.push_back( log(pPrint[i]/pPrint[i-1])/log(h[i]/h[i-1]));
      convuPr.push_back( log(uPrint[i]/uPrint[i-1])/log(h[i]/h[i-1]));
      convdivPr.push_back( log(divPrint[i]/divPrint[i-1])/log(h[i]/h[i-1]));
      // convdivPrLoc.push_back( log(divPrintLoc[i]/divPrintLoc[i-1])/log(h[i]/h[i-1]));

      convmaxdivPr.push_back( log(maxDivPrint[i]/maxDivPrint[i-1])/log(h[i]/h[i-1]));
    }

    nx *= 2;
    ny *= 2;
    // nx += 1;
    // ny += 1;
    // nx = (int)round( (1+0.2*i)*nx/2 )*2; // Makes a nonuniform refinement to an EVEN integer
    // ny = (int)round( (1+0.2*i)*ny/2 )*2;
    // std::cout << nx << std::endl;
    // shift = 0.5+(i+1)*hei/iters; // moves one grid cell over entire span
  }
  std::cout << "\n" << std::left
  << std::setw(10) << std::setfill(' ') << "h"
  << std::setw(15) << std::setfill(' ') << "err_new u"
  << std::setw(15) << std::setfill(' ') << "conv p"
  << std::setw(15) << std::setfill(' ') << "err u"
  << std::setw(15) << std::setfill(' ') << "conv u"
  << std::setw(15) << std::setfill(' ') << "err divu"
  << std::setw(15) << std::setfill(' ') << "conv divu"
  << std::setw(15) << std::setfill(' ') << "err_new divu"
  << std::setw(15) << std::setfill(' ') << "convLoc divu"
  << std::setw(15) << std::setfill(' ') << "err maxdivu"
  << std::setw(15) << std::setfill(' ') << "conv maxdivu"
  << "\n" << std::endl;
  for(int i=0;i<uPrint.size();++i) {
    std::cout << std::left
    << std::setw(10) << std::setfill(' ') << h[i]
    << std::setw(15) << std::setfill(' ') << pPrint[i]
    << std::setw(15) << std::setfill(' ') << convpPr[i]
    << std::setw(15) << std::setfill(' ') << uPrint[i]
    << std::setw(15) << std::setfill(' ') << convuPr[i]
    << std::setw(15) << std::setfill(' ') << divPrint[i]
    << std::setw(15) << std::setfill(' ') << convdivPr[i]
    // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
    // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
    << std::setw(15) << std::setfill(' ') << maxDivPrint[i]
    << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
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

#ifdef PROBLEM_RT_PROJECTION

/*
COMPLETE+
INFO-
TODO*
_________
+ RT0 projects ok for stabilisation with derivative
- Divergence not zero with stabilisation, and no convergence
* Is it correctly coded?

*/

namespace Data_RTProjection {
  R interfaceRad = 0.25;
  R shift = 0.5;

  R fun_radius2(const R2 P){
    return (P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift);
  }
  R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x-shift)*(P.x-shift) + (P.y-shift)*(P.y-shift)) - interfaceRad;
  }
  R fun_rhs(const R2 P, int compInd, int dom) { // also the exact solution, since we project
//    if (dom == 0) {
      if (compInd == 0) {
        return 4*P.x+4*cos(P.y)+P.y*P.y;
      } else {
        return -4*P.y-4*cos(P.x);
      }
//    } else return 0;
  }
  R fun_div(const R2 P, int compInd, int dom) {
    return 0;
  }
}

using namespace Data_RTProjection;

int main(int argc, char** argv )
{
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 10;
  int ny = 10;
  int iters = 5;

  vector<double> uPrint,divPrint,maxDivPrint,h,convuPr,convmaxdivPr,convdivPr;

  for(int i=0;i<iters;++i) {
    Mesh2 Th(nx, ny, 0., 0., 1., 1.);

    FESpace2 Lh(Th, DataFE<Mesh2>::P1);
    Fun_h levelSet(Lh, fun_levelSet);
    Interface2 interface(Th, levelSet.v);

    FESpace2 Vh(Th, DataFE<Mesh2>::BDM1);

    // CutFEM stuff
    SubDomain2 Vh1(Vh, levelSet, 1);
    SubDomain2 Vh2(Vh, levelSet,-1);

    KN<SubDomain2*> subDomainsVh(2);
    subDomainsVh(0) = &Vh1;  // domain index 0
    subDomainsVh(1) = &Vh2;  // domain index 1
    CutFESpace2 Wh(subDomainsVh, interface);

    CutFEM<Mesh2,Interface2> RTproj(Wh, interface);

    // const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
    // const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
    const R hei = Th[0].lenEdge(0);
    // const R invh = 1./hei;

    // We define fh on the cutSpace
    Fun_h fv(Wh, fun_rhs);

    Normal n;
    Tangent t;
    FunTest u(Wh,2), v(Wh,2);
    // FunTest v1(Wh, 1, 0), v2(Wh, 1, 1);
    // a(u,v)_Omega
    RTproj.addBilinearFormDomain(
      innerProduct(u, v)
    );

    // int penParam = 1; //1e-2;
    // RTproj.addFaceStabilization( // h^(2k+1)
    //   // innerProduct(penParam*pow(hei,1)*jump(u.t()*t), jump(v.t()*t)) // [Method 1: Remove jump in vel]
    //   innerProduct(penParam*pow(hei,3)*jump(grad(u).t()*n), jump(grad(v).t()*n))
    //   + innerProduct(penParam*pow(hei,3)*jump(grad(u).t()*t), jump(grad(v).t()*t))
    //
    //   // + innerProduct(penParam*pow(hei,5)*jump(DDU(u)), jump(DDU(v))) //[SECOND DERIV]
    //   // + innerProduct(penParam*pow(hei,3)*jump(div(u)), jump(div(v)))
    // );
    // RTproj.addBilinearFormInterface(
    //   innerProduct(penParam*jump(u.t()*t),jump(v.t()*t))
    //   + innerProduct(penParam*jump(u.t()*n),jump(v.t()*n))
    //   + innerProduct(penParam*pow(hei,2)*jump(grad(u).t()*n), jump(grad(v).t()*n))
    //   + innerProduct(penParam*pow(hei,2)*jump(grad(u).t()*t), jump(grad(v).t()*t))
    //   // + innerProduct(penParam*pow(hei,2)*jump(div(u)), jump(div(v)))
    // );

    // ExpressionFunFEM<Mesh2> fvdx(fv, 0, op_dx);
    // ExpressionFunFEM<Mesh2> fvdy(fv, 1, op_dy);
    // l(v)_Omega
    RTproj.addLinearFormDomain(
      innerProduct(fv,v)
      // innerProduct(fvdx, v1) + innerProduct(fvdy, v2)
    );

    std::cout << "Check: Assembly is passed." << std::endl;
    matlab::Export(RTproj.mat, "mat"+std::to_string(i)+"Cut.dat");



    RTproj.solve();

    // L2 Error calculations
    // Rn solVec(Wh.NbDoF());
    // interpolate(Wh, solVec, fun_rhs);
    // R errU = RTproj.L2norm(solVec,0,2); // L2 norm of velocity

    Fun_h femSolh(Wh, RTproj.rhs);

    // L2 norm vel
    R errU = L2normCut(femSolh,fun_rhs);

    // L2 norm div
    ExpressionFunFEM<Mesh2> femSol_0dx(femSolh, 0, op_dx);
    ExpressionFunFEM<Mesh2> femSol_1dy(femSolh, 1, op_dy);
    R errDiv = L2normCut(femSol_0dx+femSol_1dy,fun_div,Wh);

    // Max error div calculation
    R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy,fun_div, Wh);

    // Plotting
    Rn solVec(Wh.NbDoF());
    interpolate(Wh, solVec, fun_rhs); // R errU = RTproj.L2norm(solVec,0,2); //remnants of prev L2norm

    Fun2 femSol2(Wh, RTproj.rhs);
    VTKcutWriter2 writer(femSol2, levelSet, "RTproj_"+to_string(i)+".vtk");
    writer.add("velocity", 0, 2);
    Fun2 sol2(Wh, solVec);
    VTKcutWriter2 writerS2(sol2, levelSet, "RTprojEx"+to_string(i)+".vtk");
    writerS2.add("velocity", 0, 2);

    divPrint.push_back(errDiv); // uPrint,divPrint,maxDivPrint,h,convuPr,convmaxdivPr,convdivPr;
    maxDivPrint.push_back(maxErrDiv);
    uPrint.push_back(errU);
    h.push_back(1./nx);

    if(i==0) {convdivPr.push_back(0);convmaxdivPr.push_back(0);convuPr.push_back(0);}
    else {
     convdivPr.push_back( log(divPrint[i]/divPrint[i-1])/log(h[i]/h[i-1]));
     convmaxdivPr.push_back( log(maxDivPrint[i]/maxDivPrint[i-1])/log(h[i]/h[i-1]));
     convuPr.push_back( log(uPrint[i]/uPrint[i-1])/log(h[i]/h[i-1]));
    }

    // shift = 0.5+(i+1)*hei/iters; // moves one grid cell over entire span
    nx *= 2;
    ny *= 2;
   }

   std::cout <<
    "\n\n h \t err div \t\t conv div \t\t err maxDiv \t\t conv maxDiv \t\t err u \t\t conv u \n\n"
    << std::endl;
   for(int i=0;i<uPrint.size();++i) {
     std::cout << std::left
     << std::setw(10) << std::setfill(' ') << h[i]
     << std::setw(24) << std::setfill(' ') << divPrint[i]
     << std::setw(24) << std::setfill(' ') << convdivPr[i]
     << std::setw(24) << std::setfill(' ') << maxDivPrint[i]
     << std::setw(20) << std::setfill(' ') << convmaxdivPr[i]
     << std::setw(20) << std::setfill(' ') << uPrint[i]
     << std::setw(20) << std::setfill(' ') << convuPr[i] << std::endl;
   }

}
#endif

#ifdef PROBLEM_MIXED_DARCY

/*
COMPLETE+
INFO-
TODO*
_________

+ Konvergensordning 1 korrekt för RT0/P0
+ Konvergensordning 2 korrekt för RT1/P1dc
+ Konvergensordning 1 korrekt för RT0/P0 med artificiell rand
+ Konvergensordning 2 korrekt för RT1/P1dc med artificiell rand

*/

namespace Data_MixedDarcy {
  bool artificialInterface = false; // [Change to true/false depending on artificial interface used/not]

  R interfaceRad = 0.25;
  R shift = 0.5;

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
     return r2/(2*radius2)+3./2.;
  }
  R fun_force(const R2 P, int compInd) {
    return 0;
  }
  R fun_div(const R2 P, int compInd, int dom) {
    R radius2 = interfaceRad*interfaceRad;
    return -2./radius2;
  }
  R fun_exact(const R2 P, int compInd, int dom) {
    Diff<R,2> X(P.x,0), Y(P.y,1);
    Diff<R,2> r2 = (X-0.5)*(X-0.5) + (Y-0.5)*(Y-0.5);
    R radius2 = interfaceRad*interfaceRad;
    Diff<R, 2> val = r2/(2*radius2) + 3./2;
    if (compInd==2) {
      return val.val;
    }
    else {
      return -val.d[compInd];
    }
  }
}

using namespace Data_MixedDarcy;

int main(int argc, char** argv )
{
  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh2> Fun_h;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 10;
  int ny = 10;
  int nIter = 5;

  vector<double> uPrint,pPrint,divPrint,maxDivPrint,h,convuPr,convpPr,convdivPr,convmaxdivPr;

  for(int i=0;i<nIter;++i) {
    Mesh2 Th(nx, ny, 0., 0., 1., 1.);

    KN<const GTypeOfFE<Mesh2>* > arrayFE(2);
    arrayFE(0) = &DataFE<Mesh2>::BDM1;
    arrayFE(1) = &DataFE<Mesh2>::P0;
    GTypeOfFESum<Mesh2> mixedFE(arrayFE);
    FESpace2 Vh(Th, mixedFE);
    FESpace2 Ph(Th, DataFE<Mesh2>::P2); // for the RHS

    FEM<Mesh2> darcyMixed(Vh);

    Paraview2 writerS(Vh, "darcyMixed_"+to_string(i)+".vtk"); // for plots
    Paraview2 writerS2(Vh, "darcyMixedEx"+to_string(i)+".vtk");

    if (artificialInterface) {
      FESpace2 Lh(Th, DataFE<Mesh2>::P1);
      Fun_h levelSet(Lh, fun_levelSet);
      Interface2 interface(Th, levelSet.v);
      CutFEM<Mesh2,Interface2> darcyMixed(Vh, interface); // [Overwrites darcyMixed]

      Paraview2 writerS(Vh, levelSet, "darcyMixed_"+to_string(i)+".vtk"); // for plots
      Paraview2 writerS2(Vh, levelSet, "darcyMixedEx"+to_string(i)+".vtk");
    }

    R xi = 3./4; // > 1/2
    R xi0 = (xi-0.5)/2.;
    R mu = 1;// CutFEM_Parameter mu("mu",100.,1.);

    // const CutFEM_Parameter& invh(Parameter::invh);
    // const CutFEM_Parameter& lambdaG(Parameter::lambdaG);
    // const CutFEM_Parameter& lambdaB(Parameter::lambdaB);
    const R hei = Th[0].lenEdge(0);
    const R invh = 1./hei;

    Fun_h fv(Vh, fun_force);
    Fun_h fq(Vh, fun_div);
    Fun_h u0(Vh, fun_dirichlet);
    Fun_h p0(Ph, fun_neumann);

    int d = 2;
    Normal n;
    Tangent t;
    FunTest p(Vh,1,d), q(Vh,1,d), u(Vh,d), v(Vh,d);

    // a(u,v)_Omega
    darcyMixed.addBilinearFormDomain(
      innerProduct(mu*u, v)
      - innerProduct(p, div(v)) + innerProduct(div(u), q) // b(p,v)-b(q,u)
    );
    // Only on Gamma_D (vel normal comp)
    // darcyMixed.addBilinearFormBorder(
    //   innerProduct(lambdaB*u.t()*n, invh*v.t()*n)
    //   + innerProduct(p, v.t()*n) - innerProduct(u.t()*n, q)
    // );

    ExpressionFunFEM<Mesh2> ffq(fq, 0, op_id);

    // l(v)_Omega
    darcyMixed.addLinearFormDomain(
      innerProduct(fv, v)
      + innerProduct(ffq, q)
    );

    // ExpressionFunFEM<Mesh2> u0(gh, d, op_id);
    darcyMixed.addLinearFormBorder(
      - innerProduct(p0, v.t()*n) // Only on Gamma_N (pressure)
      // + innerProduct(u0, invh*(v.t()*n)*lambdaB) // Only on Gamma_D (vel normal comp)
      // - innerProduct(u0, q) // Only on Gamma_D (vel normal comp)
    );

    std::cout << "Check: Assembly is passed." << std::endl;
    matlab::Export(darcyMixed.mat, "mat"+std::to_string(i)+"Cut.dat");
    darcyMixed.solve();

    // [L2 ERRORS]
    Fun_h femSolh(Vh, darcyMixed.rhs);

    // L2 norm vel
    R errU = L2norm(femSolh,fun_exact,0,2);
    R errP = L2norm(femSolh,fun_exact,2,1);

    // L2 norm div
    ExpressionFunFEM<Mesh2> femSol_0dx(femSolh, 0, op_dx);
    ExpressionFunFEM<Mesh2> femSol_1dy(femSolh, 1, op_dy);
    R errDiv = L2norm(femSol_0dx+femSol_1dy,fun_div,Vh);

    // Max error div calculation
    R maxErrDiv = maxNorm(femSol_0dx+femSol_1dy,fun_div,Vh);

    // Rn solExVec(Vh.NbDoF());
    // interpolate(Vh, solExVec, fun_exact);
    // R errU = darcyMixed.L2norm(solExVec,0,2);
    // R errP = darcyMixed.L2norm(solExVec,2,1);

    // [PLOTTING]
    // Defined writers in the beginning if statement (for scope reasons)
    if (i<5) {
      writerS.add(femSolh, "velocity", 0, 2);
      writerS.add(femSolh, "pressure", 2, 1);
      writerS.add(femSol_0dx+femSol_1dy, "divergence");

      Fun_h solh(Vh, fun_exact);
      Fun_h divSolh(Vh, fun_div);

      writerS2.add(solh, "velocity", 0, 2);
      writerS2.add(solh, "pressure", 2, 1);
      writerS2.add(divSolh, "divergence", 0, 1);
      // ExpressionFunFEM<Mesh2> divSolEFF(divSolh, 0, op_id);
      // writerS2.add((femSol_0dx+femSol_1dy-divSolEFF)*(femSol_0dx+femSol_1dy-divSolEFF), "divergence"); // Plots error
    }

    pPrint.push_back(errP);
    uPrint.push_back(errU);
    divPrint.push_back(errDiv);
    maxDivPrint.push_back(maxErrDiv);
    h.push_back(1./nx);

    if(i==0) {convpPr.push_back(0);convuPr.push_back(0);convdivPr.push_back(0);convmaxdivPr.push_back(0);}
    else {
     convpPr.push_back( log(pPrint[i]/pPrint[i-1])/log(h[i]/h[i-1]));
     convuPr.push_back( log(uPrint[i]/uPrint[i-1])/log(h[i]/h[i-1]));
     convdivPr.push_back( log(divPrint[i]/divPrint[i-1])/log(h[i]/h[i-1]));
     convmaxdivPr.push_back( log(maxDivPrint[i]/maxDivPrint[i-1])/log(h[i]/h[i-1]));
    }

    // nx *= 2;
    // ny *= 2;
    nx += 1;
    ny += 1;
   }

   std::cout << "\n" << std::left
   << std::setw(10) << std::setfill(' ') << "h"
   << std::setw(15) << std::setfill(' ') << "err p"
   << std::setw(15) << std::setfill(' ') << "conv p"
   << std::setw(15) << std::setfill(' ') << "err u"
   << std::setw(15) << std::setfill(' ') << "conv u"
   << std::setw(15) << std::setfill(' ') << "err divu"
   << std::setw(15) << std::setfill(' ') << "conv divu"
   << std::setw(15) << std::setfill(' ') << "err maxdivu"
   << std::setw(15) << std::setfill(' ') << "conv maxdivu"
   << "\n" << std::endl;
   for(int i=0;i<uPrint.size();++i) {
     std::cout << std::left
     << std::setw(10) << std::setfill(' ') << h[i]
     << std::setw(15) << std::setfill(' ') << pPrint[i]
     << std::setw(15) << std::setfill(' ') << convpPr[i]
     << std::setw(15) << std::setfill(' ') << uPrint[i]
     << std::setw(15) << std::setfill(' ') << convuPr[i]
     << std::setw(15) << std::setfill(' ') << divPrint[i]
     << std::setw(15) << std::setfill(' ') << convdivPr[i]
     << std::setw(15) << std::setfill(' ') << maxDivPrint[i]
     << std::setw(15) << std::setfill(' ') << convmaxdivPr[i] << std::endl;
   }
 }
#endif
