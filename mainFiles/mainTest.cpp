#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "../util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif
#include "interface_levelSet.hpp"
#include "cut_mesh.hpp"
#include "baseProblem.hpp"
#include "paraview.hpp"
#include <unordered_map>
#include "omp.h"

#include "../num/matlab.hpp"
#define TEST_PERFORMANCE_ASSEMBLY


#ifdef TEST_PERFORMANCE_ASSEMBLY


struct hash_pair {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

typedef Mesh2 Mesh;
typedef ActiveMeshT2 CutMesh;
typedef FESpace2       Space;
typedef CutFESpaceT2 CutSpace;
typedef typename Space::FElement FElement;
typedef FunFEM<Mesh> Fun_h;
typedef TestFunction<2> FunTest;

int main(int argc, char** argv ){

  MPIcf cfMPI(argc,argv);


  Mesh Th(200, 200, 0., 0., 1., 1.);
  // FESpace2 Vh(Th, DataFE<Mesh2>::P1);
  Lagrange2 FEu(2);
  Space Vh(Th, FEu);
  Mesh1 Qh(2, 0, 1);
  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  const TimeSlab& In(Ih[0]);

  FunTest u(Vh,2), v(Vh,2);
  ListItemVF<2> kernel = contractProduct(Eps(u),Eps(v));
  // std::cout << kernel << std::endl;


  // std::map<std::pair<int,int>, double> A;
  // std::map<std::pair<int,int>, double> B0;
  // std::map<std::pair<int,int>, double> B1;
  // std::map<std::pair<int,int>, double> B2;
  // std::map<std::pair<int,int>, double> B3;
  // std::map<std::pair<int,int>, double> *A;

  int thread_count = 2;
  cout << "Threads: ";
  cin >> thread_count;
  MPIcf::Bcast(thread_count, MPIcf::Master(), 1);
  std::unordered_map<std::pair<int,int>, double,hash_pair> A[thread_count];

  // std::map<std::pair<int,int>, double> A[thread_count];


  // progress bar0("Build sparsity", Th.last_element());
  //
  // for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {
  //   bar0 += Th.next_element();
  //
  //   const FElement& FK(Vh[k]);
  //
  //   for(int l=0;l<Vh.N;++l){
  //
  //     // if(it==0 && jt == 0) A = &B0;
  //     // else if(it==0 && jt == 1) A = &B1;
  //     // else if(it==1 && jt == 0) A = &B2;
  //     // else A = &B3;
  //     for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
  //       for(int jt=In.dfcbegin(0); jt<In.dfcend(0); ++jt) {
  //         for(int i = FK.dfcbegin(l); i < FK.dfcend(l); ++i) {
  //           for(int j = FK.dfcbegin(l); j < FK.dfcend(l); ++j) {
  //             // (*A)[std::make_pair(FK.loc2glb(i),FK.loc2glb(j))] = 0;
  //             for(int ii=1;ii<thread_count;++ii){
  //
  //               A[ii][std::make_pair(FK.loc2glb(i, it),FK.loc2glb(j, jt))] = 0;
  //
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }
  // bar0.end();



  omp_set_num_threads(thread_count);
  // #pragma omp parallel for
// #pragma omp parallel
double t0 = MPIcf::Wtime();

    int iam = 0, np = 1;
    // #pragma omp parallel
    // omp_set_num_threads(2);
#pragma omp parallel default(shared) private(iam, np)
{
  np = omp_get_num_threads();
  iam = omp_get_thread_num();
  // printf("Hello from thread %d out of %d from process %d out of %d \n",
  // iam, np, MPIcf::my_rank(), MPIcf::size());
// #pragma omp  for private(index_i0, index_j0)
//default(shared)
// std::map<std::pair<int,int>, double> A;
progress bar("Add Bilinear Mesh", Th.last_element()/thread_count);

  // private(index_i0, index_j0), firstprivate(localContributionMatrix)
  #pragma omp for
  for(int k=Th.first_element(); k<Th.last_element(); k+= Th.next_element()) {
    bar += Th.next_element();

    // printf("Hello from thread %d out of %d from process %d out of %d \n",
    // iam, np, MPIcf::my_rank(), MPIcf::size());

    // printf("thread %d do element %d \n",iam, k);

    const FElement& FK(Vh[k]);

    for(int l=0;l<kernel.size();++l){

      // for(int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
      //   for(int j = FK.dfcbegin(0); j < FK.dfcend(0); ++j) {
      //     A[std::make_pair(FK.loc2glb(i),FK.loc2glb(j))] +=  i*j;
      //   }
      // }

      for(int it=In.dfcbegin(0); it<In.dfcend(0); ++it) {
        for(int jt=In.dfcbegin(0); jt<In.dfcend(0); ++jt) {
          for(int i = FK.dfcbegin(kernel[l].cv); i < FK.dfcend(kernel[l].cv); ++i) {
            for(int j = FK.dfcbegin(kernel[l].cu); j < FK.dfcend(kernel[l].cu); ++j) {
// #pragma omp critical
              (A[iam])[std::make_pair(FK.loc2glb(i, it),FK.loc2glb(j, jt))] += i*j;
              // (*A)[std::make_pair(FK.loc2glb(i),FK.loc2glb(j))] += i*j;
              // B[it][jt][std::make_pair(FK.loc2glb(i),FK.loc2glb(j))] += i*j;

            }
          }
        }
      }
    }
  }
bar.end();
}
std::cout << " real time " << MPIcf::Wtime() - t0 << std::endl;
// std::cout << A[0].size() << std::endl;

// matlab::Export(A[0], "mat0.dat");
// std::map<std::pair<int,int>, double> B = A[0];
std::map<std::pair<int,int>, double> B;
double tt0 = MPIcf::Wtime();
for(int i=0;i<thread_count;++i){
  for(const auto& it : A[i]) {
    B[it.first] += it.second;
  }
}
std::cout << " time finalizing " << MPIcf::Wtime() - tt0 << std::endl;

// matlab::Export(A[0], "mat1.dat");
  return 0;
}

#endif


#ifdef TEST_2D
R fun_expr(const R2 P,const int i) {return 10;}
R fun_div0(const R2 P,const int i, int d) {return (i==0)?4*sin(P[0]):4*cos(P[1]);}
R fdiv(const R2 P,const int i, int d) {return 4*cos(P[0]) - 4*sin(P[1]);}


// R fun_levelSet(const R2 P, const int i) {
//   double x=P.x, y=P.y;
//   // return 0.5*(sqrt(x*x + y*y) -0.1*cos(5*atan2(y,x)+2) - 0.25);
//   return x - 0.25;
// }
// R fun_levelSet(const R2 P, const int i, const double t) {
//   double x=P.x, y=P.y;
//
//   return x - 0.25 - t;
// }
double fun_levelSet(const R2 P, const int i, const R t) {
    R xc = 0.5 + 0.28*sin(M_PI*t), yc = 0.5 - 0.28*cos(M_PI*t);
    return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.22)*(P.y-0.22)) - 0.17;

    // return sqrt((P.x - xc)*(P.x - xc) + (P.y - yc)*(P.y - yc)) - 0.17;
}

// Level-set function initial
double fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x-0.5)*(P.x-0.5) + (P.y-0.22)*(P.y-0.22)) - 0.17;
}
R fun_velocity(const R2 P, const int i) {
    if (i == 0) return M_PI*(0.5-P.y);
    else return M_PI*(P.x-0.5);
}


R fun_levelSet2_1(const R2 P, const int i) {
  R shiftX = 0.5;
  R shiftY = 0.5;
  R r = 0.250001;//0.3;
  double x=P.x, y=P.y;
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY)) - r;
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

R fun_levelSet_t0(const R2 P, const int i) {
  R shiftX = 0.45;
  R shiftY = 0.45;
  R r = 0.250001;//0.3;
  double x=P.x, y=P.y;
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY)) - r;
}
R fun_levelSet_t1(const R2 P, const int i) {
  R shiftX = 0.55;
  R shiftY = 0.55;
  R r = 0.250001;//0.3;
  double x=P.x, y=P.y;
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY)) - r;
}

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

  // Mesh Th(14, 14, 0., 0., 1., 1.);
  Mesh Th(10, 10, 0., 0., 1., 1.);
  Th.info();
  double h = 1. / 10;
  double dt = h / 3;
  int Nt = (int) (0.25 / dt);
  Mesh1 Qh(Nt,0,0.25);  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  dt = 0.25 / Nt;
  const QuadratureFormular1d& qTime(*Lobatto(3));
  TimeInterface<Mesh> interface(qTime);
  std::cout << dt << std::endl;

  Lagrange2 FEvelocity(2);
  FESpace2 VelVh(Th, FEvelocity);
  Fun_h vel(VelVh, fun_velocity);
  Normal N;


  FESpace2 Lh(Th, DataFE<Mesh2>::P1);
  const Uint nbTime = qTime.n;
  const TimeSlab& In(Ih[0]);
  vector<Fun_h> ls(nbTime);
  for (int i=0; i<nbTime; i++) {
    R tt = In.Pt(R1(qTime(i).x));
    ls[i].init(Lh, fun_levelSet, tt);
    interface.init(i, Th, ls[i]);
  }


//
// return 0;

  CutFEM<Mesh> convdiff(qTime);

  ActiveMesh<Mesh> Kh2(Th);
  Kh2.truncate(interface, 1);


  // Paraview<Mesh> writer(Kh2, "test_mesh.vtk");
  // writer.add(ls[0], "ls0", 0, 1);
  // writer.add(ls[1], "ls1", 0, 1);
  // writer.add(ls[2], "ls2", 0, 1);

  FESpace2 Vh(Th, DataFE<Mesh>::P1dc);
  CutSpace Wh(Kh2, Vh);
  Wh.info();



  convdiff.initSpace(Wh, In);
  double kappaTilde2 = 1;
  FunTest u(Wh, 1), v(Wh, 1);
  convdiff.addBilinear(
    innerProduct((vel*N)*u, kappaTilde2*v)
      , interface
      , In
  );

  matlab::Export(convdiff.mat_, "mat0.dat");


  // Space Vh(Th, DataFE<Mesh2>::P1);
  // FEM<Mesh2> problem(Vh);
  // FunTest u(Vh,1), v(Vh,1);


  // Space Lh(Th, DataFE<Mesh>::P1);
  // Fun_h levelSet1(Lh, fun_levelSet2_1);
  //
  // Fun_h levelSet_t0(Lh, fun_levelSet_t0);
  // Fun_h levelSet_t1(Lh, fun_levelSet_t1);

  // Fun_h levelSet4(Lh, fun_levelSet2_4);

  // KN<Fun_h> ls(2);
  // ls(0).init(Lh, fun_levelSet);
  // ls(1).init(Lh, fun_levelSet2_2);

  // InterfaceLevelSet<Mesh> interface1(Th, levelSet1);
  // Interface_LevelSet<Mesh> interface2(Th, levelSet2);
  // Interface_LevelSet<Mesh> interface3(Th, levelSet3);
  // Interface_LevelSet<Mesh> interface4(Th, levelSet4);
  // KN<Interface_LevelSet<Mesh>*> ls(2);
  // ls(0) = &interface1;
  // ls(1) = &interface2;


  // TimeInterface<Mesh> interface(2);
  // interface.init(0,Th, levelSet_t0);
  // interface.init(1,Th, levelSet_t1);



  // ActiveMesh<Mesh> Kh0(Th);
  // Kh0.createSurfaceMesh(interface);
  // ActiveMesh<Mesh> Kh1(Th);
  // Kh1.truncate(interface, -1);
  // ActiveMesh<Mesh> Kh2(Th);
  // Kh2.truncate(interface, 1);


  // MacroElement<Mesh> macro1(Kh1, 0.2);
  // MacroElement<Mesh> macro1full(Kh1, 1);
  // MacroElement<Mesh> macro2(Kh2, 0.2);
  // MacroElement<Mesh> macro2full(Kh2, 1);
  // MacroElementSurface<Mesh> macro_surface(interface1, 0.2);


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


  // Paraview<Mesh> writer(Th, "backMesh.vtk");
  // writer.writeActiveMesh(Kh0, "active_mesh0.vtk");
  // writer.writeActiveMesh(Kh1, "active_mesh1.vtk");
  // writer.writeActiveMesh(Kh2, "active_mesh2.vtk");

  //
  // writer.writeFaceStab(Kh0, 0, "fullstab_face0.vtk");
  // writer.writeFaceStab(Kh1, 0, "fullstab_face1.vtk");
  // writer.writeFaceStab(Kh2, 0, "fullstab_face2.vtk");
  //
  //
  // writer.writeNonStabMesh(Kh1, "Omega1_nonstabElement.vtk");
  // writer.writeNonStabMesh(Kh2, "Omega2_nonstabElement.vtk");
  // writer.writeNonStabMeshEdge(Kh1, "Omega1_nonstabEdge.vtk");
  // writer.writeNonStabMeshEdge(Kh2, "Omega2_nonstabEdge.vtk");

  // writer.writeMacroInnerEdge(macro1, 0, "macro_inner_edge1.vtk");
  // writer.writeMacroOutterEdge(macro1, 0, "macro_outter_edge1.vtk");
  // writer.writeNonStabMesh(macro1, 0, "nonStab_macro_element1.vtk");
  // writer.writeMacroInnerEdge(macro1full, 0, "macrofull_inner_edge1.vtk");
  // writer.writeMacroOutterEdge(macro1full, 0, "macrofull_outter_edge1.vtk");
  // writer.writeNonStabMesh(macro1full, 0, "nonStab_macroFull_element1.vtk");
  //
  // writer.writeMacroInnerEdge(macro2, 0, "macro_inner_edge2.vtk");
  // writer.writeMacroOutterEdge(macro2, 0, "macro_outter_edge2.vtk");
  // writer.writeNonStabMesh(macro2, 0, "nonStab_macro_element2.vtk");
  // writer.writeMacroInnerEdge(macro2full, 0, "macrofull_inner_edge2.vtk");
  // writer.writeMacroOutterEdge(macro2full, 0, "macrofull_outter_edge2.vtk");
  // writer.writeNonStabMesh(macro2full, 0, "nonStab_macroFull_element2.vtk");
  //
  // writer.writeMacroElement(macro, 0, "macro_element.vtk");
  // writer.writeMacroElement(macro, 0, "macro_element.vtk");
  // writer.writeMacroElement(macro_surface, "macro_element0.vtk");
  // writer.writeMacroOutterEdge(macro_surface, "macro_outter_edge0.vtk");
  // writer.writeMacroInnerEdge(macro_surface, "macro_inner_edge0.vtk");
  // writer.writeNonStabMesh(macro_surface, "nonStab_maxro_element0.vtk");
  //

  // Paraview<Mesh> writer(Th, "backMesh_smooth.vtk");
  // writer.add(levelSet_t0, "levelSet1", 0, 1);
  // writer.add(levelSet_t1, "levelSet2", 0, 1);
  // writer.add(levelSet1  , "levelSet", 0, 1);
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


class  TestParameter : public VirtualParameter {
public:
  const CutFEMParameter& mu;
  TestParameter(const CutFEMParameter& m) : mu(m){}
  double evaluate(int domain, double h, double meas, double measK, double meas_Cut) const {
    return 2*mu.evaluate(domain,h,meas,measK,meas_Cut);
  }
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
  CutFEMParameter mu(data_mu);

  std::vector<R3> data_beta = {R3(1,1,1), R3(2,2,2), R3(3,3,3)};
  CutFEM_R3 beta(data_beta);

  TestParameter mu2(mu);

  // CutFEMParameter kappaTilde("kappaTilde", kappaTilde1, kappaTilde2);
  // CutFEMParameter kappaTildeA("kappaTildeA", kappaTilde1*A1, kappaTilde2*A2);
  //
  // // Local weights for average function across bulk faces
  // CutFEMParameter kappa_E1("kappa_E1", fun_kappa_E1);
  // CutFEMParameter kappa_E2("kappa_E2", fun_kappa_E2);

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
  Kh0.createSurfaceMesh(interface1);

  // FINITE ELEMENT SPACE
  // ---------------------------------------------------------------------------
  Space Vh(Th, DataFE<Mesh>::P1dc);
  CutSpace Wh(Khi, Vh);
  Wh.info();
  CutSpace Sh(Kh0, Vh);
  Sh.info();


  // INTERPOLATION OF THE VELOCITY FIELD
  // ---------------------------------------------------------------------------
  Lagrange3<Mesh> FEvelocity(1);
  Space velVh(Th, FEvelocity);
  Fun_h vel(velVh, fun_velocity);

  // BUILDING THE SYSTEM
  // ---------------------------------------------------------------------------
  CutFEM<Mesh> problem(Wh); problem.add(Sh);
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

  Paraview<Mesh> writer(Khi,"hexa_Test.vtk");
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
