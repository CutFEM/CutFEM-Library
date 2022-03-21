#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "../util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif

// #include "baseCutProblem.hpp"
// #include "levelSet.hpp"
// #include "extension.cpp"
// #include "../num/gnuplot.hpp"
#include "interface_levelSet.hpp"
#include "cut_mesh.hpp"
#include "baseCutProblem.hpp"

// #include "gnuplot.hpp"
#define TEST_3D



#ifdef TEST_2D
R fun_expr(const R2 P,const int i) {return 10;}
R fun_div0(const R2 P,const int i, int d) {return (i==0)?4*sin(P[0]):4*cos(P[1]);}
R fdiv(const R2 P,const int i, int d) {return 4*cos(P[0]) - 4*sin(P[1]);}


R fun_levelSet(const R2 P, const int i) {
  double x=P.x, y=P.y;
  return 0.5*(sqrt(x*x + y*y) -0.1*cos(5*atan2(y,x)+2) - 0.25);
}

R fun_levelSet2_1(const R2 P, const int i) {
  R shiftX = 0.5;
  R shiftY = 0.5;
  R r = 0.250001;//0.3;
  double x=P.x, y=P.y;
  // return 0.5*(sqrt(x*x + y*y) -0.075*cos(5*atan2(y,x)+2) - 0.25);
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY)) - r;
  // return P.x-P.y + 0.1197;
  // return P.y - 0.550766;

}

R fun_levelSet2_2(const R2 P, const int i) {
  R shiftX = 1.5;
  R shiftY = 0.5;
  R r = 0.3;
  // return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY)) - r;
  return P.x - 2*P.y - .2197;
  // return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY)) - r;
  // return P.x+P.y - 0.40753;

  // return P.x+P.y - 0.40753;
  // return P.x - 0.50737;
}
R fun_levelSet2_3(const R2 P, const int i) {
  R shiftX = 1.5;
  R shiftY = 1.5;
  R r = 0.3;
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY)) - r;
  // return P.x-P.y + 0.07064;
}

R fun_levelSet2_4(const R2 P, const int i) {
  R shiftX = 0.5;
  R shiftY = 1.5;
  R r = 0.3;
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY)) - r;
  // return P.x+P.y - 0.40753;
}

R fun_id(const R2 P,const int i) {return i;}//P[i];}
R fun_time(const R2 P,const int i, const R t) {return 1+t;}
// R fun_test(const R2 P,const int i) {return P.x + P.y;}
R fun_test(const R2 P,const int i, int d) {return d;}

//
// template<((double)(*pfun)(double)) f>
// double myfun(double x) { return f(x);}
//
// double myx2(double x) {return x*x;}

int main(int argc, char** argv ){
  // typedef MeshQuad2 Mesh;
  // typedef Cut_MeshQ2 CutMesh;
  // typedef FESpaceQ2       Space;
  // typedef CutFESpaceQ2 CutSpace;
  typedef Mesh2 Mesh;
  typedef Cut_MeshT2 CutMesh;
  typedef FESpace2       Space;
  typedef CutFESpaceT2 CutSpace;


  typedef typename Space::FElement FElement;
  typedef FunFEM<Mesh> Fun_h;
  //
  // const QuadratureFormular1d& QFB = QF_GaussLegendre2;

  const double cpubegin = CPUtime();

  MPIcf cfMPI(argc,argv);


  int nx = 20;
  int ny = 20;
  int nz = 3;
  // MeshHexa Th(nx, ny, nz, 0., 0., 0., 1., 1., 1.);
  // Mesh Th(nx, ny, -0.5, -0.5, 1., 1.);  // Flower shape
  // Mesh Th(nx, ny, 0., 0., 2., 2.);
  Mesh Th(nx, ny, 0., 0., 1., 1.);
  Th.info();

  // FESpaceQ3 Lh(Th, DataFE<Mesh>::P1);
  // Fun_h levelSet(Lh, fun_levelSet3);
  Space Lh(Th, DataFE<Mesh>::P1);
  Fun_h levelSet1(Lh, fun_levelSet2_1);
  // Fun_h levelSet2(Lh, fun_levelSet2_2);
  // Fun_h levelSet3(Lh, fun_levelSet2_3);
  // Fun_h levelSet4(Lh, fun_levelSet2_4);

  // KN<Fun_h> ls(2);
  // ls(0).init(Lh, fun_levelSet);
  // ls(1).init(Lh, fun_levelSet2_2);

  Interface_LevelSet<Mesh> interface1(Th, levelSet1);
  // Interface_LevelSet<Mesh> interface2(Th, levelSet2);
  // Interface_LevelSet<Mesh> interface3(Th, levelSet3);
  // Interface_LevelSet<Mesh> interface4(Th, levelSet4);
  // KN<Interface_LevelSet<Mesh>*> ls(2);
  // ls(0) = &interface1;
  // ls(1) = &interface2;


  // Time_Interface<Mesh> interface(2);
  // interface.init(0,Th, ls[0]);
  // interface.init(1,Th, ls[1]);

  // interface.init(Th, &levelSet1);


  Cut_Mesh<Mesh> cutTh(Th,interface1);
  // cutTh.truncate(interface1, -1);
  // cutTh.create_surface_mesh(interface1);
  // cutTh.add(interface1);
  // cutTh.add(interface2);
  // cutTh.truncate(interface3, -1);
  // cutTh.truncate(interface4, -1);
  // cutTh.add(interface1);


  // cutTh.add(interface4, -1);
  // CutMesh cutTh(Th, interface);

  // CutSpace Vh(cutTh, Lh);
  // Fun_h ftest(Vh, fun_test);
  MacroElementCL<Mesh> macro(cutTh, 2e-1);
  // MacroElementSurfaceCL<Mesh> macro_surface(interface1, 1e-1);


  // const MeshHexa::Element& T(Th[0]);
  // std::vector<int> list_cut;

  // Rn ls(8);
  // for(int i=0;i<8;++i) ls[i] = fun_levelSet3(T[i],0);
  // std::cout << ls << std::endl;
  // const SignPattern<MeshHexa::Element> signPattern(ls);

  // signPattern.get_ordered_list_node(list_cut);

  // for(int i=0;i<list_cut.size();++i) std::cout << list_cut[i] << std::endl;

  // Interface_LevelSet<MeshHexa> interface(Th, levelSet.v);

  // KN<const GTypeOfFE<MeshQuad2>* > arrayFE(2);
  // arrayFE(0) = &DataFE<MeshQuad2>::P1;
  // arrayFE(1) = &DataFE<MeshQuad2>::P1;
  // GTypeOfFESum<MeshQuad2> FE(arrayFE);
  // FESpaceQ2 Ph(Th, FE);
  // Fun_h test(Ph, fun_id);

  // const FElement& FK(Lh[0]);
  // KNMK<double> bf(FK.NbDoF(),FK.N,3);
  // R2 x(0.0, 0.0);
  // FK.BF(Fop_Dall, x, bf);
  // std::cout << bf<< std::endl;

  // ParaviewCut<Mesh> writer1(cutTh, "macro_active_mesh.vtk");
  // writer1.write_active_mesh();
  // writer1.add(levelSet1, "levelSet1", 0, 1);
  // ParaviewCut<Mesh> writer1(cutTh, "macro_backMesh.vtk");
  // writer1.write_active_mesh();
  // writer1.add(levelSet1, "levelSet1", 0, 1);

  // ParaviewCut<Mesh> writer (cutTh, "macro_surface_element.vtk");
  // writer.write_macro_element(macro_surface);
  // ParaviewCut<Mesh> writer3 (cutTh, "macro_surface_outter_edge.vtk");
  // writer3.writeMacroOutterEdge(macro_surface);
  // ParaviewCut<Mesh> writer4 (cutTh, "macro_surface_inner_edge.vtk");
  // writer4.writeMacroInnerEdge(macro_surface);

  ParaviewCut<Mesh> writer (cutTh, "macro_element.vtk");
  writer.write_macro_element(macro, 0);
  ParaviewCut<Mesh> writer2 (cutTh, "macro_inner_edge.vtk");
  writer2.writeMacroInnerEdge(macro, 0);
  ParaviewCut<Mesh> writer3 (cutTh, "macro_outter_edge.vtk");
  writer3.writeMacroOutterEdge(macro, 0);



  // writer.add(levelSet1, "levelSet1", 0, 1);
  // writer.add(levelSet2, "levelSet2", 0, 1);
  // writer.add(levelSet3, "levelSet3", 0, 1);
  // writer.add(levelSet4, "levelSet4", 0, 1);


  // ParaviewCut<Mesh> writer(cutTh,  "background_mesh_quad.vtk");
  // writer.add(ftest, "test_fun", 0, 1);
  // writer.add(levelSet1, "levelSet1", 0, 1);
  // writer.add(levelSet2, "levelSet2", 0, 1);
  // Paraview<Mesh> writer(Vh, levelSet, "quad_Test.vtk");
  // Paraview<Mesh> writers(Lh, nullptr, "hexa_Test2.vtk");
  // writer.add(levelSet, "levelSet", 0, 1);
  // writer.add(ftest, "test_fun", 0, 1);

  // getchar();

  std::cout << " need to check case 3 and 4 " << std::endl;


  return 0;

}

#endif
#ifdef TEST_3D

R fun_levelSet(const R3 P, const int i) {
  R shiftX = 0.;
  R shiftY = 0.;
  R shiftZ = 0.;
  R r = 1.0001;
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY) + (P.z-shiftZ)*(P.z-shiftZ)) - r;
  // return -0.25*P.x - 0.25*P.y + 0.25*P.z - 0.125;  // 3 cuts
  // return P.x - 0.5;  // 4 cuts
  // return 0.5*P.x +0.5*P.y - P.z - 0.2;   // 5 cuts
  // return 0.5*P.x +0.5*P.y - P.z - 0.002;   // 5 cuts
  // return 0.5*P.x +0.3*P.y - 0.5*P.z - 0.2;   // 6 cuts
  // return 0.5*P.x -0.5*P.y + 0.5*P.z - 0.2;   // 6 cuts
  // return -0.5*P.x +0.5*P.y + 0.5*P.z - 0.2;   // 6 cuts
}
R fun_levelSet1(const R3 P, const int i) {
  R shiftX = 0.75;
  R shiftY = 0.75;
  R shiftZ = 0.75;
  R r = 0.6;
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY) + (P.z-shiftZ)*(P.z-shiftZ)) - r;
}
R fun_levelSet2(const R3 P, const int i) {
  R shiftX = -0.75;
  R shiftY = -0.75;
  R shiftZ = -0.75;
  R r = 0.6;
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY) + (P.z-shiftZ)*(P.z-shiftZ)) - r;
}
R fun_test(const R3 P,const int i, int d) {return d;}
R fun_velocity(const R3 P, const int i){
  if(i == 0) return P.z/10;
  else if(i == 1) return 0;
  else       return -P.x/10;
}
double fun_kappa_E1(int i, double hh, double meas, double measK, double meas_Cut) {
    return meas_Cut/measK;
}
double fun_kappa_E2(int i, double hh, double meas, double measK, double meas_Cut) {
    return 1-meas_Cut/measK;
}

int main(int argc, char** argv )
{
  typedef Mesh3       Mesh;
  typedef CutMeshT3   CutMesh;
  typedef FESpace3    Space;
  typedef CutFESpaceT3 CutSpace;
  typedef TestFunction<3> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  const double cpubegin = CPUtime();

  MPIcf cfMPI(argc,argv);

  // MESH PARAMETERS
  // ---------------------------------------------------------------------------
  int    nx = 10, ny = 10, nz = 10;
  double lx = 3., ly = 3., lz = 3.;
  Mesh Th(nx, ny, nz, -1.5, -1.5, -1.5, lx, ly, lz);
  Th.info();
  double h = lx/(nx-1);

  // PROBLEM PARAMETERS
  // ---------------------------------------------------------------------------
  double A0 = 1.0, A1 = 1.0, A2 = 0.5;
  double kappa1 = 2.0, kappa2 = 0.5, kappa01 = 1, kappa02 = 2;
  double kappaTilde1 = kappa1/kappa01, kappaTilde2 = kappa2/kappa02;
  CutFEM_Parameter kappaTilde("kappaTilde", kappaTilde1, kappaTilde2);
  CutFEM_Parameter kappaTildeA("kappaTildeA", kappaTilde1*A1, kappaTilde2*A2);

  // Local weights for average function across bulk faces
  CutFEM_Parameter kappa_E1("kappa_E1", fun_kappa_E1);
  CutFEM_Parameter kappa_E2("kappa_E2", fun_kappa_E2);

  // FULLSTAB PARAMETERS
  double tau_a0 = 1e0, tau_b0 = 0;
  double tau_a1 = 1e1, tau_b1 = 0;
  double tau_a2 = 1e1, tau_b2 = 0;
  double lambdaB = A1*1e1/h;
  double lambdaA0 = A0*5e3/h;


  // LEVEL-SET & INTERFACE CONSTRUCTION
  // ---------------------------------------------------------------------------
  Space Lh(Th, DataFE<Mesh>::P1);
  Fun_h levelSet(Lh, fun_levelSet);
  Interface_LevelSet<Mesh> interface(Th, levelSet);
  // Fun_h levelSet1(Lh, fun_levelSet1);
  // Interface_LevelSet<Mesh> interface1(Th, levelSet1);
  // Fun_h levelSet2(Lh, fun_levelSet2);
  // Interface_LevelSet<Mesh> interface2(Th, levelSet2);

  // CUTMESH && SURFACE MESH CONSTRUCTION
  // ---------------------------------------------------------------------------
  Cut_Mesh<Mesh> Khi(Th, interface);
  // Khi.add(interface1);
  // Khi.add(interface2);
  Khi.info();
  Cut_Mesh<Mesh> Kh0(Th);
  Kh0.create_surface_mesh(interface);

  // FINITE ELEMENT SPACE
  // ---------------------------------------------------------------------------
  Space Vh(Th, DataFE<Mesh>::P1dc);
  CutSpace Wh(Khi, Vh);
  Wh.info();
  CutSpace Sh(Kh0, Vh);
  Sh.info();


  // INTERPOLATION OF THE VELOCITY FIELD
  // ---------------------------------------------------------------------------
  Lagrange3 FEvelocity(1);
  Space velVh(Th, FEvelocity);
  Fun_h vel(velVh, fun_velocity);

  // BUILDING THE SYSTEM
  // ---------------------------------------------------------------------------
  CutFEM<Mesh> problem({&Wh, &Sh});
  FunTest u(Wh,1), v(Wh,1);         // Bulk
  FunTest u0(Sh,1), v0(Sh,1);       // Surface
  Normal n;
  // BEGIN ASSEMBLY
  // ---------------------------------------------------------------------------
  std::cout << " -------------------  ASSEMBLY  -------------------" << std::endl;
  double t0 = MPIcf::Wtime();
  // problem.addBilinear(
  //           innerProduct(kappaTildeA*grad(u),grad(v))                  // (3.12)
  //         + innerProduct(kappaTilde*(vel.expression(),grad(u)), v)*0.5 // (3.14)
  //         - innerProduct(kappaTilde*u,(vel.expression(),grad(v)))*0.5  // (3.14)
  //         , Khi
  // );

  problem.addBilinear(
          - innerProduct(kappaTildeA*average(grad(u)*n),jump(v))      // (3.12)
          - innerProduct(kappaTildeA*jump(u), average(grad(v)*n))         // (3.13)
          + innerProduct(average((vel*n)*u), kappaTilde*jump(v))*0.5         // (3.15)
          - innerProduct(jump((vel*n)*u),    kappaTilde*average(v))*0.5      // (3.15)
          , Khi
          , innerFace
  );



  // problem.addLinear(
  //         innerProduct(1,v)
  //       , Khi
  // );
  // matlab::Export(problem.get_matrix(), "A.dat");
  // matlab::Export(problem.get_rhs(), "rhs.dat");

  double t1 = MPIcf::Wtime();
  std::cout << " Time assembly \t" << t1-t0 << std::endl;
  //
  // Fun_h uh(Wh, problem.get_rhs());
  Rn_ data_uh = problem.rhs_(SubArray(Wh.get_nb_dof(),0));
  Fun_h uh(Wh, problem.get_rhs());
  Fun_h ftest(Wh, fun_test);
  // ftest.print();

  ParaviewCut<Mesh> writer(Khi,"hexa_Test.vtk");
  // writer.add(levelSet2, "levelSet", 0, 1);
  writer.add(uh, "sol", 0, 1);
  writer.add(ftest, "domain", 0, 1);

}


// SEEMS THAT THE 3 CUTS CASE IS PROBLEMATIC

// {
  // Cut_Mesh2 cutTh(Th, levelSet.v);




  // double x =2.;
  // std::cout << myfun<myx2>(x) << std::endl;

  // int nx = 2;
  // int ny = 2;
  // Mesh Th(nx, ny, 0., 0., 1, 1.);
  // FESpace2 Vh(Th, DataFE<Mesh2>::RT1);
  // const Element& K(Th[0]);   // The reference element
  // const FElement& FK(Vh[0]);
  //
  // KNMK<double> bf(FK.NbDoF(),FK.N,1);
  //
  // const QuadratureFormular1d& QF = QF_GaussLegendre4;
  // const QuadratureFormular2d& QFK = QuadratureFormular_T_5;
  //
  // // wanna compute int_Fi phi_BDM1_i . n * phi_P1_j(Fi)
  // for(int e=0;e<3;++e){
  //   std::cout << " edge " << e << std::endl;
  //
  //   for(int df=0;df<8;++df){
  //
  //     // int df = 3;
  //     double val = 0, val2=0;
  //     double meas = K.lenEdge(e);
  //     int eOrientation = K.EdgeOrientation(e);
  //     R2 normal = K.EdgeOrientation(e)*K.N(e);
  //     // R2 normal = -K.Edge(e).perp();
  //
  //     for(int iq=0;iq<QF.getNbrOfQuads();++iq) {
  //       QuadraturePoint1d ip_1d(QF[iq]);
  //       R2 ip_hat = K.toKref(ip_1d, e);
  //       FK.BF(Fop_D0, ip_hat, bf);
  //
  //       // R l0 = QF[iq].x, l1 = 1 - QF[iq].x;
  //       // R p0 = (2 * l0 - l1) * 2;     // poly othogonaux to \lambda_1
  //       // R p1 = (2 * l1 - l0) * 2;     // poly othogonaux to \lambda_0
  //       // R lambda1 = eOrientation * p0 * QF[iq].a;    // [some quadrature function?]
  //       // R lambda0 = eOrientation * p1 * QF[iq].a;    //
  //       // if(eOrientation < 0) {
  //       //   Exchange(lambda1, lambda0);
  //       // }
  //       //
  //       // val  += lambda0*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y) ;
  //       // val2 += lambda1*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y) ;
  //
  //       double lambda1 = 6*ip_1d.x-2;
  //       double lambda0 = 4-6*ip_1d.x;
  //       if(eOrientation < 0) {
  //         Exchange(lambda1, lambda0);
  //       }
  //       val += meas*ip_1d.a*lambda0*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y) ;
  //       val2 += meas*ip_1d.a*lambda1*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y) ;
  //
  //       // val  += (-6*ip_1d.x+3)  *meas*ip_1d.getWeight()*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y);
  //       // val2 += 1*meas*ip_1d.getWeight()*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y);
  //     }
  //     std::cout << " bf " << df << "\t" << val << "\t" << val2 << std::endl;
  //   }
  //   std::cout << "----------------"<< std::endl;
  // }
  // {
  //   R2 B[2] = {K.Edge(1), K.Edge(2)};
  //   B[0] = B[0].perp();
  //   B[1] = B[1].perp();
  //   std::cout << "bubbles "<< std::endl;
  //
  //   for(int df=0;df<8;++df){
  //
  //     double val = 0, val2=0;
  //
  //     for(int iq=0;iq<QFK.getNbrOfQuads();++iq) {
  //
  //       QuadraturePoint2d ip_2d(QFK[iq]);
  //       R2 mip_Ks   = FK.map(ip_2d);
  //       R2 ip_hat = K.toKref(mip_Ks);
  //
  //       FK.BF(Fop_D0, ip_hat, bf);
  //       double w = QFK[iq].a * 0.5;
  //       val  += w*(bf(df,0,0)*B[0].x + bf(df,1,0)*B[0].y) ;
  //       val2 += w*(bf(df,0,0)*B[1].x + bf(df,1,0)*B[1].y) ;
  //     }
  //
  //   std::cout << " bf " << df << "\t" << val << "\t" << val2 << std::endl;
  // }
  //
  // std::cout << "----------------"<< std::endl;
  // }

  //
  // int nx = 10;
  // int ny = 10;
  // Mesh Th(nx, ny, 0., 0., 1, 1.);
  // FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  // Fun_h levelSet(Lh, fun_levelSet);
  // Interface2 interface(Th, levelSet.v);
  //
  // Lagrange2 FEtest(1);
  // FESpace2 Vh(Th, FEtest);
  //
  // CutFESpace2 Wh(Lh, interface, {1,-1});
  //
  // CutFEM<Mesh> poisson(Wh);
  // FunTest u(Wh,1), v(Wh,1);
  // const double t0= MPIcf::Wtime();
  // Normal n;
  // Fun_h fq(Vh, fun_id);
  // CutFEM_Parameter kappaTilde("kappaTilde", 1,2);
  //
  // poisson.addBilinear(
  //   // innerProduct(ln(u),v)
  //   // innerProduct(average(u, kappaTilde), jump(v)),
  //   // innerProduct(kappaTilde*jump(u), jump(v)),
  //   // innerProduct(fabs(fq*n)*jump(u), jump(v)),
  //   // innerEdge
  //   interface
  // );
  //
  //
  // getchar();
  //
  // std::cout << "time \t" << MPIcf::Wtime() - t0 << std::endl;


  //
  // Fun_h levelSet(Lh, fun_levelSet);
  // Interface2 interface(Th, levelSet.v);
  //
  // CutFESpace2 Wh(Lh , interface, {1,-1});
  // MacroElement macro(Wh, 0.01);
  //
  //
  // KN<const GTypeOfFE<Mesh2>* > myFE(2);
  // myFE(0) = &DataFE<Mesh2>::P2;
  // myFE(1) = &DataFE<Mesh2>::P2;
  // GTypeOfFESum<Mesh2> FEvh(myFE);
  // FESpace2 Vh(Th, FEvh);
  // FunTest du(Vh, 2);
  // std::cout << gradS(average2(du)) << std::endl;

  // std::cout << gradS(average2(du)) << std::endl;

  // gnuplot::save(Th);
  // gnuplot::save(interface);
  // gnuplot::save(macro);
  // gnuplot::save(Wh,0,"Vh1.dat");
  // gnuplot::save(Wh,1,"Vh2.dat");
//   std::map<std::pair<int,int>,double> A, S, St, C;
//   Rn b(3);
//   Rn x(3);
//   b(0) = 1;
//   b(1) = 2;
//   b(2) = 3;
//   A[make_pair(0,0)] = 1;
//   A[make_pair(0,1)] = 2;
//   A[make_pair(0,2)] = 3;
//   A[make_pair(1,0)] = 3;
//   A[make_pair(1,1)] = 1;
//   A[make_pair(1,2)] = 5;
//   A[make_pair(2,0)] = 3;
//   A[make_pair(2,1)] = 5;
//   A[make_pair(2,2)] = 2;
//
//   S[make_pair(0,0)] = 1;
//   S[make_pair(1,2)] = 1;
//   S[make_pair(2,2)] = 1;
//
//   St[make_pair(0,0)] = 1;
//   St[make_pair(2,1)] = 1;
//   St[make_pair(2,2)] = 1;
//
//   multiply(3, St, A, C);
//   multiply(3, C,  S, A);
//   for(auto it = A.begin(); it!= A.end();++it){
//     std::cout << "( " << it->first.first << " , "
//               << it->first.second << " )  = \t"
//               << it->second << std::endl;
//             }
//   multiply(3,3, St, b, x);
//   std::cout << x << std::endl;
//
//   removeDF(3, A, x, {1});
//
//   for(auto it = A.begin(); it!= A.end();++it){
//     std::cout << "( " << it->first.first << " , "
//               << it->first.second << " )  = \t"
//               << it->second << std::endl;
//             }
//   std::cout << x << std::endl;
//
// set<int> tt;
// tt.insert(2);
// tt.insert(1);
// tt.insert(5);
// tt.insert(1);
// show("test" , tt);
// }


#endif
