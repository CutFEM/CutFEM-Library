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
  R shiftX = 0.45;
  R shiftY = 0.45;
  R r = 0.250001;//0.3;
  double x=P.x, y=P.y;
  // return 0.5*(sqrt(x*x + y*y) -0.075*cos(5*atan2(y,x)+2) - 0.25);
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY)) - r;
  // return P.x-P.y + 0.1197;
  // return P.y - 0.550766;

}

R fun_levelSet2_2(const R2 P, const int i) {
  R shiftX = 0.55;
  R shiftY = 0.55;
  R r = 0.250001;//0.3;
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY)) - r;
  // return P.x - 2*P.y - .2197;
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


  int nx = 18;
  int ny = 18;
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
  Fun_h levelSet2(Lh, fun_levelSet2_2);
  // Fun_h levelSet3(Lh, fun_levelSet2_3);
  // Fun_h levelSet4(Lh, fun_levelSet2_4);

  // KN<Fun_h> ls(2);
  // ls(0).init(Lh, fun_levelSet);
  // ls(1).init(Lh, fun_levelSet2_2);

  // Interface_LevelSet<Mesh> interface1(Th, levelSet1);
  // Interface_LevelSet<Mesh> interface2(Th, levelSet2);
  // Interface_LevelSet<Mesh> interface3(Th, levelSet3);
  // Interface_LevelSet<Mesh> interface4(Th, levelSet4);
  // KN<Interface_LevelSet<Mesh>*> ls(2);
  // ls(0) = &interface1;
  // ls(1) = &interface2;


  Time_Interface<Mesh> interface1(2);
  interface1.init(0,Th, levelSet1);
  interface1.init(1,Th, levelSet2);

  // interface.init(Th, &levelSet1);


  Cut_Mesh<Mesh> cutTh(Th);//,interface1);
  cutTh.truncate(interface1, 1);
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
  // MacroElementCL<Mesh> macro(cutTh, 1);
  // MacroElementSurfaceCL<Mesh> macro(interface1, 1e-1);


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


  Paraview<Mesh> writer0 (Lh, nullptr, "backgroundMesh.vtk");
  writer0.add(levelSet1, "levelSet1", 0, 1);
  writer0.add(levelSet2, "levelSet2", 0, 1);

  ParaviewCut<Mesh> writer1,writer2,writer3,writer4;
  writer1.writeActiveMesh(cutTh, "activeMesh.vtk");
  writer2.writeFaceStab(cutTh, 0, "faceStab.vtk");
  // writer2.writeMacroElement(macro, 0, "macroElement.vtk");
  // writer3.writeMacroInnerEdge(macro, 0, "macroInnerEdge.vtk");
  // writer4.writeMacroOutterEdge(macro, 0, "macroOutterEdge.vtk");

  // ParaviewCut<Mesh> writer (cutTh, "macro_element.vtk");
  // writer.write_macro_element(macro, 0);
  // ParaviewCut<Mesh> writer2 (cutTh, "macro_inner_edge.vtk");
  // writer2.writeMacroInnerEdge(macro, 0);
  // ParaviewCut<Mesh> writer3 (cutTh, "macro_outter_edge.vtk");
  // writer3.writeMacroOutterEdge(macro, 0);



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
R fun_id(const R3 P,const int i) {return P[i];}
R fun_1(const R3 P,const int i) {return 1;}


int main(int argc, char** argv )
{
  typedef MeshHexa     Mesh;
  typedef Cut_MeshQ3   CutMesh;
  typedef FESpaceQ3    Space;
  typedef CutFESpaceQ3 CutSpace;

  typedef typename Space::FElement FElement;
  typedef FunFEM<Mesh> Fun_h;
  const double cpubegin = CPUtime();

  MPIcf cfMPI(argc,argv);


  int    nx = 10, ny = 10, nz = 10;
  double lx = 3., ly = 3., lz = 3.;
  MeshHexa Th(nx, ny, nz, -1.5, -1.5, -1.5, lx, ly, lz);
  Th.info();


  Space Lh(Th, DataFE<Mesh>::P1);
  Fun_h levelSet(Lh, fun_levelSet);
  Interface_LevelSet<Mesh> interface(Th, levelSet);


  Cut_Mesh<Mesh> Kh(Th, interface);
  Cut_Mesh<Mesh> Sh(Th);
  Sh.create_surface_mesh(interface);




  ParaviewCut<Mesh> writer(Th,"hexa_Test.vtk");
  writer.add(levelSet, "levelSet", 0, 1);

  // writer.add(ftest, "test_fun", 0, 1);

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
