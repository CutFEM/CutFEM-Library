#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "../util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif

#include "baseCutProblem.hpp"
#include "levelSet.hpp"
#include "../num/gnuplot.hpp"
// #include "gnuplot.hpp"
#define TEST_2D



#ifdef TEST_2D
R fun_expr(const R2 P,const int i) {return 10;}
R fun_div0(const R2 P,const int i, int d) {return (i==0)?4*sin(P[0]):4*cos(P[1]);}
R fdiv(const R2 P,const int i, int d) {return 4*cos(P[0]) - 4*sin(P[1]);}

R fun_levelSet(const R2 P, const int i) {
  R shiftX = 0.5;
  R shiftY = 0.5;
  R r = 0.3;
  return sqrt((P.x-shiftX)*(P.x-shiftX) + (P.y-shiftY)*(P.y-shiftY)) - r;
}
R fun_id(const R2 P,const int i) {return (i==0);}//P[i];}
R fun_time(const R2 P,const int i, const R t) {return 1+t;}
R fun_test(const R2 P,const int i) {return P.x + P.y;}

//
// template<((double)(*pfun)(double)) f>
// double myfun(double x) { return f(x);}
//
// double myx2(double x) {return x*x;}

int main(int argc, char** argv )
{
  typedef Mesh2 Mesh;
  typedef FESpace2 FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename Mesh2::Element Element;

  typedef TestFunction<2> FunTest;
  typedef FunFEM<Mesh> Fun_h;

  const QuadratureFormular1d& QFB = QF_GaussLegendre2;

  const double cpubegin = CPUtime();

  MPIcf cfMPI(argc,argv);

  // double x =2.;
  // std::cout << myfun<myx2>(x) << std::endl;

  int nx = 2;
  int ny = 2;
  Mesh Th(nx, ny, 0., 0., 1, 1.);
  FESpace2 Vh(Th, DataFE<Mesh2>::RT1);
  const Element& K(Th[0]);   // The reference element
  const FElement& FK(Vh[0]);

  KNMK<double> bf(FK.NbDoF(),FK.N,1);

  const QuadratureFormular1d& QF = QF_GaussLegendre4;
  const QuadratureFormular2d& QFK = QuadratureFormular_T_5;

  // wanna compute int_Fi phi_BDM1_i . n * phi_P1_j(Fi)
  for(int e=0;e<3;++e){
    std::cout << " edge " << e << std::endl;

    for(int df=0;df<8;++df){

      // int df = 3;
      double val = 0, val2=0;
      double meas = K.lenEdge(e);
      int eOrientation = K.EdgeOrientation(e);
      R2 normal = K.EdgeOrientation(e)*K.N(e);
      // R2 normal = -K.Edge(e).perp();

      for(int iq=0;iq<QF.getNbrOfQuads();++iq) {
        QuadraturePoint1d ip_1d(QF[iq]);
        R2 ip_hat = K.toKref(ip_1d, e);
        FK.BF(Fop_D0, ip_hat, bf);

        // R l0 = QF[iq].x, l1 = 1 - QF[iq].x;
        // R p0 = (2 * l0 - l1) * 2;     // poly othogonaux to \lambda_1
        // R p1 = (2 * l1 - l0) * 2;     // poly othogonaux to \lambda_0
        // R lambda1 = eOrientation * p0 * QF[iq].a;    // [some quadrature function?]
        // R lambda0 = eOrientation * p1 * QF[iq].a;    //
        // if(eOrientation < 0) {
        //   Exchange(lambda1, lambda0);
        // }
        //
        // val  += lambda0*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y) ;
        // val2 += lambda1*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y) ;

        double lambda1 = 6*ip_1d.x-2;
        double lambda0 = 4-6*ip_1d.x;
        if(eOrientation < 0) {
          Exchange(lambda1, lambda0);
        }
        val += meas*ip_1d.a*lambda0*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y) ;
        val2 += meas*ip_1d.a*lambda1*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y) ;

        // val  += (-6*ip_1d.x+3)  *meas*ip_1d.getWeight()*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y);
        // val2 += 1*meas*ip_1d.getWeight()*(bf(df,0,0)*normal.x + bf(df,1,0)*normal.y);
      }
      std::cout << " bf " << df << "\t" << val << "\t" << val2 << std::endl;
    }
    std::cout << "----------------"<< std::endl;
  }
  {
    R2 B[2] = {K.Edge(1), K.Edge(2)};
    B[0] = B[0].perp();
    B[1] = B[1].perp();
    std::cout << "bubbles "<< std::endl;

    for(int df=0;df<8;++df){

      double val = 0, val2=0;

      for(int iq=0;iq<QFK.getNbrOfQuads();++iq) {

        QuadraturePoint2d ip_2d(QFK[iq]);
        R2 mip_Ks   = FK.map(ip_2d);
        R2 ip_hat = K.toKref(mip_Ks);

        FK.BF(Fop_D0, ip_hat, bf);
        double w = QFK[iq].a * 0.5;
        val  += w*(bf(df,0,0)*B[0].x + bf(df,1,0)*B[0].y) ;
        val2 += w*(bf(df,0,0)*B[1].x + bf(df,1,0)*B[1].y) ;
      }

    std::cout << " bf " << df << "\t" << val << "\t" << val2 << std::endl;
  }

  std::cout << "----------------"<< std::endl;
  }

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

  return 0;

}

#endif
#ifdef TEST_3D

// R fun_expr(const R2 P,const int i) {return 10;}
// R fun_div0(const R2 P,const int i, int d) {return (i==0)?4*sin(P[0]):4*cos(P[1]);}
// R fdiv(const R2 P,const int i, int d) {return 4*cos(P[0]) - 4*sin(P[1]);}
R fun_levelSet(const R3 P, const int i) {
  return sqrt(P.x*P.x + P.y*P.y + P.z*P.z) - 0.25;
}
R fun_id(const R3 P,const int i) {return P[i];}
R fun_1(const R3 P,const int i) {return 1;}

// R fun_time(const R2 P,const int i, const R t) {return 1+t;}


int main(int argc, char** argv )
{
  const int d = 3;
  typedef Mesh3 Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef Interface3 Interface;
  typedef TestFunction<d> FunTest;
  typedef FunFEM<Mesh> Fun_h;
  typedef GLevelSet<Mesh> LevelSet;
  // typedef GCurvature<Mesh> Curvature;

  MPIcf cfMPI(argc,argv);
  const double cpubegin = CPUtime();

  int nx = 40;
  int ny = 40;
  int nz = 40;

  Mesh3 Th(nx, ny, nz, -1, -1, -1, 2., 2., 2.);
  double meshSize = (4./nx);
  double dT = meshSize/4;//1./100;
  GTime::time_step = dT;
  GTime::total_number_iteration = 1;//int(0.25/dT);
  Mesh1 Qh(GTime::total_number_iteration+1, GTime::t0, GTime::final_time());


  FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
  const QuadratureFormular1d& qTime(*Lobatto(3));
  const Uint nbTime = qTime.n;
  const Uint ndfTime = Ih[0].NbDoF();

  std::cout << " nx        \t" << nx << std::endl;
  std::cout << " nb nodes  \t" << Th.nv << std::endl;
  std::cout << " nb elements\t" << Th.nt << std::endl;
  std::cout << " Mesh size \t" << meshSize << std::endl;
  std::cout << " Time Step \t" << dT << std::endl;
  std::cout << " Number of iteration \t" << GTime::total_number_iteration << std::endl;
  std::cout << " We are using the Simpson rule (3 quadrature points)" << std::endl;
  std::cout << " number of quadrature points in time : \t" << nbTime << std::endl;
  std::cout << " number of dof in time per time slab : \t" << ndfTime << std::endl;


  // levelSet stuff
  // ----------------------------------------------
  std::cout << "\n --------------------------------------- \n " << std::endl;
  std::cout << " We use a P2 levelSet function and project it to a piecewise linear space \n in order to create the interface " << std::endl;
  FESpace Lh  (Th, DataFE<Mesh>::P1);
  Fun_h levelSet(Lh, fun_levelSet);

  Interface3 interface(Th, levelSet.v);


  SubDomain3 Lh1(Lh, levelSet, 1);
  SubDomain3 Lh2(Lh, levelSet,-1);
  KN<SubDomain3*> subDomainsLh(2);
  subDomainsLh(0) = &Lh1;  // domain index 0
  subDomainsLh(1) = &Lh2;  // domain index 1
  CutFESpace3 Zh(subDomainsLh, interface);

  Fun_h solCut1(Zh, fun_1);
  std::cout << integral(solCut1, 0, 1) << "\t" << 4./3*M_PI*pow(0.25,3)<< std::endl;
  std::cout << integralSurf(solCut1, 0, 0) << "\t" << 4*M_PI*pow(0.25,2) << std::endl;


  // Paraview3 writerS(Zh, levelSet, "testPlotCut3D.vtk");
  // writerS.add(levelSet, "levelSet", 0, 1);
}



#endif
