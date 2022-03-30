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
#include "paraview.hpp"

#include "../num/matlab.hpp"
#define TEST_2D



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
  typedef ActiveMeshT2 CutMesh;
  typedef FESpace2       Space;
  typedef CutFESpaceT2 CutSpace;


  typedef typename Space::FElement FElement;
  typedef FunFEM<Mesh> Fun_h;
  typedef TestFunction<2> FunTest;

  // const QuadratureFormular1d& QFB = QF_GaussLegendre2

  const double cpubegin = CPUtime();

  MPIcf cfMPI(argc,argv);

  Mesh Th(10, 10, 0., 0., 1., 1.);
  Th.info();
  Space Vh(Th, DataFE<Mesh2>::P1);
  FEM<Mesh2> problem(Vh);

  FunTest u(Vh,1), v(Vh,1);
  //----------------------------------------------
  // problem.addBilinear((u,v) + (grad(u), grad(v)) , Th);
  // problem.addLinear((1,v), Th);
  // matlab::Export(problem.mat_, "matK_new.dat");
  // matlab::Export(problem.rhs_, "rhsK_new.dat");
  // problem.cleanMatrix();
  // problem.rhs_ = 0.;
  //
  // //----------------------------------------------
  // problem.addBilinear((u,v) + (grad(u), grad(v)) , Th, boundary);
  // problem.addLinear((1,v), Th, boundary);
  // matlab::Export(problem.mat_, "matdw_new.dat");
  // matlab::Export(problem.rhs_, "rhsdw_new.dat");
  // problem.cleanMatrix();
  // problem.rhs_ = 0.;

  //----------------------------------------------
  problem.addBilinear((jump(u),jump(v)), Th, innerEdge);
  problem.addLinear((1,jump(v)), Th, innerEdge);
  matlab::Export(problem.mat_, "matdK_new.dat");
  matlab::Export(problem.rhs_, "rhsdK_new.dat");
  problem.cleanMatrix();
  problem.rhs_ = 0.;


  return 0;


  // const double cpubegin = CPUtime();
  //
  // MPIcf cfMPI(argc,argv);
  //
  //
  // int nx = 20;
  // int ny = 20;
  // int nz = 3;
  // // MeshHexa Th(nx, ny, nz, 0., 0., 0., 1., 1., 1.);
  // // Mesh Th(nx, ny, -0.5, -0.5, 1., 1.);  // Flower shape
  // // Mesh Th(nx, ny, 0., 0., 2., 2.);
  // Mesh Th(nx, ny, 0., 0., 1., 1.);
  // Th.info();

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

  InterfaceLevelSet<Mesh> interface1(Th, levelSet1);
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


  ActiveMesh<Mesh> cutTh(Th,interface1);
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
  MacroElement<Mesh> macro(cutTh, 2e-1);
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

  // Paraview<Mesh> writer;// (cutTh, "macro_element.vtk");
  // writer.writeMacroElement(macro, 0, "macro_element.vtk");
  // // ParaviewCut<Mesh> writer2 (cutTh, "macro_inner_edge.vtk");
  // writer.writeMacroInnerEdge(macro, 0, "macro_inner_edge.vtk");
  // // ParaviewCut<Mesh> writer3 (cutTh, "macro_outter_edge.vtk");
  // writer.writeMacroOutterEdge(macro, 0, "macro_outter_edge.vtk");
  //


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
  // return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY) + (P.z-shiftZ)*(P.z-shiftZ)) - r;
  // return -0.25*P.x - 0.25*P.y + 0.25*P.z - 0.125;  // 3 cuts
  // return P.x - 0.5;  // 4 cuts
  // return 0.5*P.x +0.5*P.y - P.z - 0.2;   // 5 cuts
  // return 0.5*P.x +0.5*P.y - P.z - 0.002;   // 5 cuts
  // return 0.5*P.x +0.3*P.y - 0.5*P.z - 0.2;   // 6 cuts
  // return 0.5*P.x -0.5*P.y + 0.5*P.z - 0.2;   // 6 cuts
  // return -0.5*P.x +0.5*P.y + 0.5*P.z - 0.2;   // 6 cuts

  return -P.x - P.y - 0.03876;

}
R fun_levelSet1(const R3 P, const int i) {
  R shiftX = 0.6;
  R shiftY = 0.6;
  R shiftZ = 0.6;
  R r = 0.4;
  return -sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY) + (P.z-shiftZ)*(P.z-shiftZ)) + r;
}
R fun_levelSet2(const R3 P, const int i) {
  R shiftX = -0.6;
  R shiftY = -0.6;
  R shiftZ = -0.6;
  R r = 0.4;
  return -sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY) + (P.z-shiftZ)*(P.z-shiftZ)) + r;
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

double fun_parameter_test(int domain, double h, double meas, double measK, double meas_cut) {
  double mu = CutFEM_ParameterList::listParameter["mu"]->evaluate(domain,h,meas,measK,meas_cut);
  return 2*mu;
}

class  TestParameter : public Virtual_Parameter {
public:
  const ParameterCutFEM& mu;
  TestParameter(const ParameterCutFEM& m) : mu(m){}
  double evaluate(int domain, double h, double meas, double measK, double meas_Cut) const {
    return 2*mu.evaluate(domain,h,meas,measK,meas_Cut);
  };
};

int main(int argc, char** argv )
{
  typedef Mesh3       Mesh;
  typedef ActiveMeshT3  ActiveMesh;
  typedef FESpace3    Space;
  typedef CutFESpaceT3 CutSpace;
  typedef TestFunction<3> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  const double cpubegin = CPUtime();

  MPIcf cfMPI(argc,argv);

  // MESH PARAMETERS
  // ---------------------------------------------------------------------------
  int    nx = 20, ny = 20, nz = 20;
  double lx = 3., ly = 3., lz = 3.;
  Mesh Th(nx, ny, nz, -1.5, -1.5, -1.5, lx, ly, lz);
  Th.info();
  double h = lx/(nx-1);

  // PROBLEM PARAMETERS
  // ---------------------------------------------------------------------------
  double A0 = 1.0, A1 = 1.0, A2 = 0.5;
  double kappa1 = 2.0, kappa2 = 0.5, kappa01 = 1, kappa02 = 2;
  double kappaTilde1 = kappa1/kappa01, kappaTilde2 = kappa2/kappa02;


  std::vector<double> data_mu = {1,2,3,4};
  ParameterCutFEM mu(data_mu);
  ParameterCutFEM rho(fun_parameter_test);

  std::vector<R3> data_beta = {R3(1,1,1), R3(2,2,2), R3(3,3,3)};
  CutFEM_R3 beta(data_beta);

  TestParameter mu2(mu);

  // ParameterCutFEM kappaTilde("kappaTilde", kappaTilde1, kappaTilde2);
  // ParameterCutFEM kappaTildeA("kappaTildeA", kappaTilde1*A1, kappaTilde2*A2);
  //
  // // Local weights for average function across bulk faces
  // ParameterCutFEM kappa_E1("kappa_E1", fun_kappa_E1);
  // ParameterCutFEM kappa_E2("kappa_E2", fun_kappa_E2);

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
  InterfaceLevelSet<Mesh> interface(Th, levelSet);

  Fun_h levelSet1(Lh, fun_levelSet1);
  InterfaceLevelSet<Mesh> interface1(Th, levelSet1);
  Fun_h levelSet2(Lh, fun_levelSet2);
  InterfaceLevelSet<Mesh> interface2(Th, levelSet2);

  // CUTMESH && SURFACE MESH CONSTRUCTION
  // ---------------------------------------------------------------------------
  ActiveMesh Khi(Th, interface);
  Khi.add(interface1);
  Khi.add(interface2);
  Khi.info();
  ActiveMesh Kh0(Th);
  Kh0.create_surface_mesh(interface1);

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
  // FunTest u(Wh,1,0,2), v(Wh,1,0,2);         // Bulk
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

  // problem.addBilinear(
  //         - innerProduct(kappaTildeA*average(grad(u)*n),jump(v))      // (3.12)
  //         - innerProduct(kappaTildeA*jump(u), average(grad(v)*n))         // (3.13)
  //         + innerProduct(average((vel*n)*u), kappaTilde*jump(v))*0.5         // (3.15)
  //         - innerProduct(jump((vel*n)*u),    kappaTilde*average(v))*0.5      // (3.15)
  //         , Khi
  //         , innerFace
  // );

  // problem.addBilinear(
  //         innerProduct(A0*gradS(u0),gradS(v0))             // (3.12)
  //         // + innerProduct((vel.expression(),gradS(u0)),v0)*0.5     // (3.14)
  //         // - innerProduct(u0,(vel.expression(),gradS(v0)))*0.5     // (3.14)
  //         , interface1
  // );

  // problem.addLinear(
  //         innerProduct(1,mu*u)
  //         , Khi
  //         , boundary
  //     );


  // problem.addLinear(
  //         innerProduct(1,mu*u0)
  //       , interface1
  //     );

  problem.addLinear(
          innerProduct(1,(1./mu2)*v)
        , Khi
  );
  // problem.addLinear(
  //         innerProduct(1,jump(v0))
  //       , Kh0
  //       , innerFace
  // );
  // matlab::Export(problem.get_matrix(), "A.dat");
  // matlab::Export(problem.get_rhs(), "rhs.dat");

  double t1 = MPIcf::Wtime();
  std::cout << " Time assembly \t" << t1-t0 << std::endl;
  //
  // Fun_h uh(Wh, problem.get_rhs());
  int idx0_s = Wh.get_nb_dof();
  Rn_ data_uh = problem.rhs_(SubArray(Wh.get_nb_dof(),0));
  Rn_ data_sh = problem.rhs_(SubArray(Sh.get_nb_dof(),idx0_s));



  // Fun_h uh(Wh, problem.get_rhs());
  Fun_h uh(Wh, data_uh);
  Fun_h us(Sh, data_sh);
  Fun_h ftest(Wh, fun_test);
  // ftest.print();

  ParaviewCut<Mesh> writer(Khi,"hexa_Test.vtk");
  // writer.add(levelSet2, "levelSet", 0, 1);
  writer.add(uh, "sol", 0, 1);
  // writer.add(us, "sol_surf", 0, 1);
  writer.add(ftest, "domain", 0, 1);

  // ParaviewCut<Mesh> writerS(Kh0,"surface_Test.vtk");
  // writerS.add(us, "sol", 0, 1);
  // writerS.add(levelSet1, "levelSet", 0, 1);


  // ParaviewCut<Mesh> writer(Khi,"hexa_Test.vtk");
  // // writer.add(levelSet2, "levelSet", 0, 1);
  // writer.add(uh, "sol", 0, 1);
}





#endif
