#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
//#include "../util/cputime.h"
#ifdef USE_MPI
#  include "cfmpi.hpp"
#endif

#include "finiteElement.hpp"
#include "baseProblem.hpp"
#include "paraview.hpp"
#include "../num/matlab.hpp"

// #include "../num/gnuplot.hpp"

/*
NOTES:
- compared to FANYI YANG AND XIAOPING XIE formulation -crlcrlu-k²u-gradp = f, we use -crlcrlu-k²u+gradp = f similar to Brezzi

*/




// #define PROBLEM_FITTED_MAXWELL

#define PROBLEM_UNFITTED_MAXWELL

//#define PROBLEM_UNFITTED_MAXWELL_PRIMAL


// {lacking mesh for circle fitted example, haven't found square example}
#ifdef PROBLEM_FITTED_MAXWELL

  namespace Erik_Data_FITTED_MAXWELL {
    R mu = 1;
    R eps_r = 1;
    R k = 1;
    
    R fun_rhs(double *P, int i) {
      R x = P[0];
      R y = P[1];
      if(i==0) return 0;
      else return 0;
    }
    R fun_exact_u(double *P, int i) {
      R x = P[0];
      R y = P[1];
      if(i==0)      return  20*x*pow(y,3);
      else if(i==1) return 5*pow(x,4)-5*pow(y,4);
    }
    R fun_exact_p(double *P, int i) {
      R x = P[0];
      R y = P[1];
      return 20*(3*pow(x,2)*y-pow(y,3));//-2.21664/sqrt(M_PI*interfaceRad*interfaceRad);
    }
    R fun_div(double *P, int i) {
      R x = P[0];
      R y = P[1];
      return 0;
    }
  }
  using namespace Erik_Data_FITTED_MAXWELL;

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

    std::vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 5;
    for(int i=0;i<iters;++i) { 

      std::cout << "\n ------------------------------------- " << std::endl;
      Mesh Kh(nx, ny, 0., 0., 1., 1.);
      const R hi = 1./(nx-1);
      const R penaltyParam = 1e3/hi; // interesting: 1e-1/sqrt(hi) gets conv 1/1 down to 0.00625, otherwise 1e-1/hi p variable goes down to 0.5


      Lagrange2 FEvelocity(2);
      Space VELh(Kh, FEvelocity);
      Space SCAh(Kh, DataFE<Mesh>::P2);

      Space Uh(Kh, DataFE<Mesh>::P1);
      Space Vh(Kh, DataFE<Mesh2>::RT0); 
      Space Wh(Kh, DataFE<Mesh2>::P0);

      Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
      Fun_h u0(VELh, fun_exact_u);
      Fun_h p0(SCAh, fun_exact_p); ExpressionFunFEM<Mesh> exactp(p0,0,op_id); 

      CutFEM<Mesh2> maxwell(Uh); maxwell.add(Vh); maxwell.add(Wh);

      Normal n;
      Tangent t;
      /* Syntax:
      FunTest (fem space, #components, place in space)
      */
      FunTest w(Uh,1,0), tau(Uh,1,0), u(Vh,2,0), v(Vh,2,0), p(Wh,1,0), q(Wh,1,0);

      {
      // [Bulk]
      maxwell.addBilinear( // w = curl u 
        innerProduct(1./mu*w, tau)
        - innerProduct(u, rotgrad(tau))       // + -> w = -curl u
        , Kh
      );
      maxwell.addBilinear( // -mu Delta u + grad p
        innerProduct(rotgrad(w), v)           // + -> -
        - innerProduct(k*k*eps_r*u, v)
        - innerProduct(p, div(v))
        , Kh
      );
      maxwell.addLinear(
        innerProduct(fh.expression(2), v)
        , Kh
      );
      maxwell.addBilinear(
        innerProduct(div(u), q)
        , Kh
      );
      // [Vorticity BC 2 (natural)]
      maxwell.addLinear(
        - innerProduct(u0*t,tau)              // + -> +
        - innerProduct(p0.expression(), v*n)
        , Kh, INTEGRAL_BOUNDARY
      );
      }

      // std::cout << integral(Khi,exactp,0) << std::endl;
      matlab::Export(maxwell.mat_, "mat"+std::to_string(i)+"Cut.dat");
      maxwell.solve();

      // EXTRACT SOLUTION
      int nb_vort_dof = Uh.get_nb_dof();
      int nb_flux_dof = Vh.get_nb_dof();
      Rn_ data_wh = maxwell.rhs_(SubArray(nb_vort_dof,0));
      Rn_ data_uh = maxwell.rhs_(SubArray(nb_flux_dof,nb_vort_dof));// Rn_ data_uh = maxwell.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));
      Rn_ data_ph = maxwell.rhs_(SubArray(Wh.get_nb_dof(),nb_vort_dof+nb_flux_dof));// Rn_ data_ph = maxwell.rhs_(SubArray(maxwell_.get_nb_dof(),nb_vort_dof+nb_flux_dof));
      Fun_h wh(Uh, data_wh);
      Fun_h uh(Vh, data_uh);
      Fun_h ph(Wh, data_ph);

      // [Post process pressure]
      R meanP = integral(Kh,exactp,0);
      ExpressionFunFEM<Mesh> fem_p(ph,0,op_id);
      R meanPfem = integral(Kh,fem_p,0);
      // std::cout << meanP << std::endl;
      CutFEM<Mesh2> post(Wh);
      post.addLinear(
        innerProduct(1,q)
        , Kh
      ); 
      R area = post.rhs_.sum();
      ph.v -= meanPfem/area;
      ph.v += meanP/area;


      ExpressionFunFEM<Mesh> dx_uh0(uh, 0, op_dx);
      ExpressionFunFEM<Mesh> dy_uh1(uh, 1, op_dy);

      // [Errors]
      {
        // Fun_h solw(Uh, fun_exact_w);
        Fun_h solu(Vh, fun_exact_u); Fun_h soluErr(Vh, fun_exact_u);
        Fun_h solp(Wh, fun_exact_p);
        soluErr.v -= uh.v;
        soluErr.v.map(fabs);
        // Fun_h divSolh(Wh, fun_div);
        // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

        Paraview<Mesh> writer(Kh, "maxwell_"+std::to_string(i)+".vtk");
        writer.add(wh, "vorticity" , 0, 1);
        writer.add(uh, "velocity" , 0, 2);
        writer.add(ph, "pressure" , 0, 1);
        writer.add(dx_uh0+dy_uh1, "divergence");
        // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
        writer.add(solu, "velocityExact" , 0, 2);
        writer.add(soluErr, "velocityError" , 0, 2);
        // writer.add(solh, "velocityError" , 0, 2);

        // writer.add(fabs(femDiv, "divergenceError");
      }

      // R errW      = L2normCut(wh,fun_exact_w,0,1);
      R errU      = L2norm(uh,fun_exact_u,0,2);
      R errP      = L2norm(ph,fun_exact_p,0,1);
      R errDiv    = L2norm(dx_uh0+dy_uh1,Kh);
      R maxErrDiv = maxNorm(dx_uh0+dy_uh1,Kh);
      // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
      // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);


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
    << std::setw(15) << std::setfill(' ') << "err p"
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



#ifdef PROBLEM_UNFITTED_MAXWELL

  namespace Erik_Data_UNFITTED_MAXWELL {

    R shift = 0;
    R rad = 0.7-1e-12;
    R fun_levelSet(double *P, const int i) {
      return sqrt((P[0]-shift)*(P[0]-shift) + (P[1]-shift)*(P[1]-shift)) - rad;
    }

    R mu_r = 1;
    R eps_r = 1;
    R k = 1;

    R fun_div(double *P, int i, int dom) {
      R x = P[0];
      R y = P[1];
      return 0;
    }
    R fun_rhs(double *P, int i, int dom) {
      R x = P[0];
      R y = P[1];
      if(i==0) return 2*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y) - eps_r*k*k*cos(M_PI*x)*sin(M_PI*y) + 2*x;
      else if(i==1) return -2*M_PI*M_PI*cos(M_PI*y)*sin(M_PI*x) + eps_r*k*k*sin(M_PI*x)*cos(M_PI*y) + 2*y;
    }
    R fun_exact_u(double *P, int i, int dom) {
      R x = P[0];
      R y = P[1];
      if(i==0) return   cos(M_PI*x)*sin(M_PI*y);
      else if(i==1) return -sin(M_PI*x)*cos(M_PI*y);
    }
    R fun_exact_p(double *P, int i, int dom ) {
      R x = P[0];
      R y = P[1];
      return x*x+y*y-rad*rad;// - 670.171/(pi*interfaceRad*interfaceRad);
    }
  }
  using namespace Erik_Data_UNFITTED_MAXWELL;

  int main(int argc, char** argv ) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2   Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc,argv);

    int nx = 11;
    int ny = 11;
    // int d = 2;

    std::vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 4;
    for(int i=0;i<iters;++i) { 

      std::cout << "\n ------------------------------------- " << std::endl;
      Mesh Kh(nx, ny, -1.0, -1.0, 2.0, 2.0); // [start x, start y, length x, length y]
      const R hi = 1./(nx-1); // 1./(nx-1)
      const R penaltyParam = 1e1; // 4e3, 1e2, 4e2, 4e0

      Space Lh(Kh, DataFE<Mesh2>::P1);
      Fun_h levelSet(Lh, fun_levelSet);
      InterfaceLevelSet<Mesh> interface(Kh, levelSet);


      Lagrange2 FEvelocity(4);
      Space VELh_(Kh, FEvelocity);
      Space SCAh_(Kh, DataFE<Mesh>::P4);

      Space Uh_(Kh, DataFE<Mesh>::P1);    // vorticity
      Space Vh_(Kh, DataFE<Mesh2>::RT0);  // velocity
      Space Wh_(Kh, DataFE<Mesh2>::P0);   // pressure

      // [Remove exterior]
      ActiveMesh<Mesh> Khi(Kh);
      Khi.truncate(interface, 1);

      MacroElement<Mesh> macro(Khi, 1);


      CutSpace VELh(Khi, VELh_);
      CutSpace SCAh(Khi, SCAh_);

      CutSpace Uh(Khi, Uh_);
      CutSpace Vh(Khi, Vh_);
      CutSpace Wh(Khi, Wh_);

      Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
      Fun_h u0(VELh, fun_exact_u);
      Fun_h p0(SCAh, fun_exact_p); ExpressionFunFEM<Mesh> exactp(p0,0,op_id); 

      CutFEM<Mesh2> maxwell(Uh); maxwell.add(Vh); maxwell.add(Wh);

      Normal n;
      Tangent t;
      /* Syntax:
      FunTest (fem space, #components, place in space)
      */
      FunTest w(Uh,1,0), tau(Uh,1,0), u(Vh,2,0), v(Vh,2,0), p(Wh,1,0), q(Wh,1,0);
      {
      // [Bulk]
      maxwell.addBilinear( // w = curl u 
        innerProduct(1./mu_r*w, tau)
        - innerProduct(u, rotgrad(tau))       // + -> w = -curl u
        , Khi
      );
      maxwell.addBilinear( // -mu Delta u + grad p
        innerProduct(rotgrad(w), v)           // + -> -
        - innerProduct(k*k*eps_r*u, v)
        - innerProduct(p, div(v))
        , Khi
      );
      maxwell.addLinear(
        innerProduct(fh.exprList(2), v)
        , Khi
      );
      maxwell.addBilinear(
        innerProduct(div(u), q)
        , Khi
      );
      // [Vorticity BC 2 (natural)]
      maxwell.addLinear(
        + innerProduct(u0*t,tau)              // + -> +
        - innerProduct(p0.expr(), v*n)
        , interface
      );
      double wPenParam = 1e0; 
      double uPenParam = 1e0;
      double pPenParam = 1e0;
      // // FunTest grad2un = grad(grad(u)*n)*n;
      // // // FunTest grad2pn = grad(grad(p)*n)*n;
      // // // FunTest grad2divun = grad(grad(div(u))*n)*n;
      // maxwell.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
      //   // [More primal stab]
      //   // // innerProduct(uPenParam*pow(hi,0)*jump(w), jump(tau)) // [w in P1, continuous]
      //   // +innerProduct(pPenParam*pow(hi,3)*jump(rotgrad(w)), jump(rotgrad(tau)))
      //   // +innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v))
      //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
      //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*t), jump(grad(v)*t))
      //   // // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
      //   // -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
      //   // +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
      //   // // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
      //   // // -innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))) , jump(grad(q)))

      //   // [Mixed stab (works great here!)]
      maxwell.addFaceStabilization(
        +innerProduct(pPenParam*pow(hi,1)*jump(rotgrad(w)), jump(v))
        -innerProduct(uPenParam*pow(hi,1)*jump(u), jump(rotgrad(tau)))
        +innerProduct(uPenParam*pow(hi,3)*jump(grad(rotgrad(w))*n), jump(grad(v)*n))
        -innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(rotgrad(tau))*n))
        -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
        +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
        -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
        +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))), jump(grad(q)))
        , Khi, macro
      );
      }

      // std::cout << integral(Khi,exactp,0) << std::endl;
      matlab::Export(maxwell.mat_, "mat"+std::to_string(i)+"Cut.dat");
      maxwell.solve();

      // EXTRACT SOLUTION
      int nb_vort_dof = Uh.get_nb_dof();
      int nb_flux_dof = Vh.get_nb_dof();
      Rn_ data_wh = maxwell.rhs_(SubArray(nb_vort_dof,0));
      Rn_ data_uh = maxwell.rhs_(SubArray(nb_flux_dof,nb_vort_dof));// Rn_ data_uh = maxwell.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));
      Rn_ data_ph = maxwell.rhs_(SubArray(Wh.get_nb_dof(),nb_vort_dof+nb_flux_dof));// Rn_ data_ph = maxwell.rhs_(SubArray(maxwell_.get_nb_dof(),nb_vort_dof+nb_flux_dof));
      Fun_h wh(Uh, data_wh);
      Fun_h uh(Vh, data_uh);
      Fun_h ph(Wh, data_ph);

      // [Post process pressure]
      // R meanP = integral(Khi,exactp,0);
      // ExpressionFunFEM<Mesh> fem_p(ph,0,op_id);
      // R meanPfem = integral(Khi,fem_p,0);
      // // std::cout << meanP << std::endl;
      // CutFEM<Mesh2> post(Wh);
      // post.addLinear(
      //   innerProduct(1,q)
      //   , Khi
      // ); 
      // R area = post.rhs_.sum();
      // ph.v -= meanPfem/area;
      // ph.v += meanP/area;


      ExpressionFunFEM<Mesh> dx_uh0(uh, 0, op_dx);
      ExpressionFunFEM<Mesh> dy_uh1(uh, 1, op_dy);
      auto uh_0dx = dx(uh.expr(0));

      auto uh_1dy = dy(uh.expr(1));

  
      // [Errors]
      {
        // Fun_h solw(Uh, fun_exact_w);
        Fun_h solu(Vh, fun_exact_u); Fun_h soluErr(Vh, fun_exact_u);
        Fun_h solp(Wh, fun_exact_p);
        soluErr.v -= uh.v;
        soluErr.v.map(fabs);
        // Fun_h divSolh(Wh, fun_div);
        // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

        Paraview<Mesh> writer(Khi, "maxwell_"+std::to_string(i)+".vtk");
        writer.add(wh, "vorticity" , 0, 1);
        writer.add(uh, "velocity" , 0, 2);
        writer.add(ph, "pressure" , 0, 1);
        writer.add(uh_0dx + uh_1dy, "divergence"); 

        // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
        writer.add(solp, "pressureExact" , 0, 1);
        writer.add(solu, "velocityExact" , 0, 2);
        writer.add(soluErr, "velocityError" , 0, 2);
        // writer.add(solh, "velocityError" , 0, 2);

        // writer.add(fabs(femDiv, "divergenceError");
      }

      // R errW      = L2normCut(wh,fun_exact_w,0,1);
      R errU      = L2normCut(uh,fun_exact_u,0,2);
      R errP      = L2normCut(ph,fun_exact_p,0,1);
      
      double errDiv = L2normCut(uh_0dx + uh_1dy, Khi);

      double maxErrDiv = maxNormCut(uh_0dx + uh_1dy, Khi);


      ul2.push_back(errU);
      pl2.push_back(errP);
      divl2.push_back(errDiv);
      divmax.push_back(maxErrDiv);
      h.push_back(hi);
      if(i==0) {convu.push_back(0); convp.push_back(0);}
      else {
        convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
        convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
      }

      if(i==4) {
        nx = nx+5;
        ny = ny+5;
      } else {
        nx = 2*nx-1;
        ny = 2*ny-1;
      }
    }
    std::cout << "\n" << std::left
    << std::setw(10) << std::setfill(' ') << "h"
    << std::setw(15) << std::setfill(' ') << "err p"
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


#ifdef PROBLEM_UNFITTED_MAXWELL_PRIMAL

  namespace Erik_Data_UNFITTED_MAXWELL {

    R shift = 0;
    R rad = 0.7-1e-12;
    R fun_levelSet(double *P, const int i) {
      return sqrt((P[0]-shift)*(P[0]-shift) + (P[1]-shift)*(P[1]-shift)) - rad;
    }

    R mu_r = 1;
    R eps_r = 1;
    R k = 1;

    R fun_div(double *P, int i, int dom) {
      R x = P[0];
      R y = P[1];
      return 0;
    }
    R fun_rhs(double *P, int i, int dom) {
      R x = P[0];
      R y = P[1];
      if(i==0) return   2*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y) - eps_r*k*k*cos(M_PI*x)*sin(M_PI*y) + 2*x;
      else if(i==1) return -2*M_PI*M_PI*cos(M_PI*y)*sin(M_PI*x) + eps_r*k*k*sin(M_PI*x)*cos(M_PI*y) + 2*y;
    }
    R fun_exact_u(double *P, int i, int dom) {
      R x = P[0];
      R y = P[1];
      if(i==0) return   cos(M_PI*x)*sin(M_PI*y);
      else if(i==1) return -sin(M_PI*x)*cos(M_PI*y);
    }
    R fun_exact_p(double *P, int i, int dom ) {
      R x = P[0];
      R y = P[1];
      return x*x+y*y-rad*rad;// - 670.171/(pi*interfaceRad*interfaceRad);
    }
  }
  using namespace Erik_Data_UNFITTED_MAXWELL;

  int main(int argc, char** argv ) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2   Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc,argv);

    int nx = 11;
    int ny = 11;
    // int d = 2;

    std::vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 4;
    for(int i=0;i<iters;++i) { 

      std::cout << "\n ------------------------------------- " << std::endl;
      Mesh Kh(nx, ny, -1.0, -1.0, 2.0, 2.0); // [start x, start y, length x, length y]
      const R hi = 1./(nx-1); // 1./(nx-1)
      const R penaltyParam = 1e1; // 4e3, 1e2, 4e2, 4e0

      Space Lh(Kh, DataFE<Mesh2>::P1);
      Fun_h levelSet(Lh, fun_levelSet);
      InterfaceLevelSet<Mesh> interface(Kh, levelSet);


      Lagrange2 FEvelocity(4);
      Space VELh_(Kh, FEvelocity);
      Space SCAh_(Kh, DataFE<Mesh>::P4);

      Space Uh_(Kh, DataFE<Mesh>::P1);
      Space Vh_(Kh, DataFE<Mesh2>::RT0); 
      Space Wh_(Kh, DataFE<Mesh2>::P0);

      // [Remove exterior]
      ActiveMesh<Mesh> Khi(Kh);
      Khi.truncate(interface, 1);

      MacroElement<Mesh> macro(Khi, 1);


      CutSpace VELh(Khi, VELh_);
      CutSpace SCAh(Khi, SCAh_);

      CutSpace Uh(Khi, Uh_);
      CutSpace Vh(Khi, Vh_);
      CutSpace Wh(Khi, Wh_);

      Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
      Fun_h u0(VELh, fun_exact_u);
      Fun_h p0(SCAh, fun_exact_p); ExpressionFunFEM<Mesh> exactp(p0,0,op_id); 

      CutFEM<Mesh2> maxwell(Uh); maxwell.add(Vh); maxwell.add(Wh);

      Normal n;
      Tangent t;
      /* Syntax:
      FunTest (fem space, #components, place in space)
      */
      FunTest w(Uh,1,0), tau(Uh,1,0), u(Vh,2,0), v(Vh,2,0), p(Wh,1,0), q(Wh,1,0);
      {
      // [Bulk]
      maxwell.addBilinear( // w = curl u 
        innerProduct(1./mu_r*w, tau)
        - innerProduct(u, rotgrad(tau))       // + -> w = -curl u
        , Khi
      );
      maxwell.addBilinear( // -mu Delta u + grad p
        innerProduct(rotgrad(w), v)           // + -> -
        - innerProduct(k*k*eps_r*u, v)
        - innerProduct(p, div(v))
        , Khi
      );
      maxwell.addLinear(
        innerProduct(fh.expression(2), v)
        , Khi
      );
      maxwell.addBilinear(
        innerProduct(div(u), q)
        , Khi
      );
      // [Vorticity BC 2 (natural)]
      maxwell.addLinear(
        + innerProduct(u0*t,tau)              // + -> +
        - innerProduct(p0.expression(), v*n)
        , interface
      );
      double wPenParam = 1e0; 
      double uPenParam = 1e0;
      double pPenParam = 1e0;
      // // FunTest grad2un = grad(grad(u)*n)*n;
      // // // FunTest grad2pn = grad(grad(p)*n)*n;
      // // // FunTest grad2divun = grad(grad(div(u))*n)*n;
      // maxwell.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
      //   // [More primal stab]
      //   // // innerProduct(uPenParam*pow(hi,0)*jump(w), jump(tau)) // [w in P1, continuous]
      //   // +innerProduct(pPenParam*pow(hi,3)*jump(rotgrad(w)), jump(rotgrad(tau)))
      //   // +innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v))
      //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
      //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*t), jump(grad(v)*t))
      //   // // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
      //   // -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
      //   // +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
      //   // // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
      //   // // -innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))) , jump(grad(q)))

      //   // [Mixed stab (works great here!)]
      maxwell.addFaceStabilization(
        +innerProduct(pPenParam*pow(hi,1)*jump(rotgrad(w)), jump(v))
        -innerProduct(uPenParam*pow(hi,1)*jump(u), jump(rotgrad(tau)))
        +innerProduct(uPenParam*pow(hi,3)*jump(grad(rotgrad(w))*n), jump(grad(v)*n))
        -innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(rotgrad(tau))*n))
        -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
        +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
        -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
        +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))), jump(grad(q)))

        , Khi, macro
      );
      }

      // std::cout << integral(Khi,exactp,0) << std::endl;
      matlab::Export(maxwell.mat_, "mat"+std::to_string(i)+"Cut.dat");
      maxwell.solve();

      // EXTRACT SOLUTION
      int nb_vort_dof = Uh.get_nb_dof();
      int nb_flux_dof = Vh.get_nb_dof();
      Rn_ data_wh = maxwell.rhs_(SubArray(nb_vort_dof,0));
      Rn_ data_uh = maxwell.rhs_(SubArray(nb_flux_dof,nb_vort_dof));// Rn_ data_uh = maxwell.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));
      Rn_ data_ph = maxwell.rhs_(SubArray(Wh.get_nb_dof(),nb_vort_dof+nb_flux_dof));// Rn_ data_ph = maxwell.rhs_(SubArray(maxwell_.get_nb_dof(),nb_vort_dof+nb_flux_dof));
      Fun_h wh(Uh, data_wh);
      Fun_h uh(Vh, data_uh);
      Fun_h ph(Wh, data_ph);

      // [Post process pressure]
      // R meanP = integral(Khi,exactp,0);
      // ExpressionFunFEM<Mesh> fem_p(ph,0,op_id);
      // R meanPfem = integral(Khi,fem_p,0);
      // // std::cout << meanP << std::endl;
      // CutFEM<Mesh2> post(Wh);
      // post.addLinear(
      //   innerProduct(1,q)
      //   , Khi
      // ); 
      // R area = post.rhs_.sum();
      // ph.v -= meanPfem/area;
      // ph.v += meanP/area;


      ExpressionFunFEM<Mesh> dx_uh0(uh, 0, op_dx);
      ExpressionFunFEM<Mesh> dy_uh1(uh, 1, op_dy);

      // [Errors]
      {
        // Fun_h solw(Uh, fun_exact_w);
        Fun_h solu(Vh, fun_exact_u); Fun_h soluErr(Vh, fun_exact_u);
        Fun_h solp(Wh, fun_exact_p);
        soluErr.v -= uh.v;
        soluErr.v.map(fabs);
        // Fun_h divSolh(Wh, fun_div);
        // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

        Paraview<Mesh> writer(Khi, "maxwell_"+std::to_string(i)+".vtk");
        writer.add(wh, "vorticity" , 0, 1);
        writer.add(uh, "velocity" , 0, 2);
        writer.add(ph, "pressure" , 0, 1);
        writer.add(dx_uh0+dy_uh1, "divergence");
        // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
        writer.add(solp, "pressureExact" , 0, 1);
        writer.add(solu, "velocityExact" , 0, 2);
        writer.add(soluErr, "velocityError" , 0, 2);
        // writer.add(solh, "velocityError" , 0, 2);

        // writer.add(fabs(femDiv, "divergenceError");
      }

      // R errW      = L2normCut(wh,fun_exact_w,0,1);
      R errU      = L2normCut(uh,fun_exact_u,0,2);
      R errP      = L2normCut(ph,fun_exact_p,0,1);
      R errDiv    = L2normCut(dx_uh0+dy_uh1,fun_div,Khi);
      R maxErrDiv = maxNormCut(dx_uh0+dy_uh1,fun_div,Khi);
      // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
      // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);


      ul2.push_back(errU);
      pl2.push_back(errP);
      divl2.push_back(errDiv);
      divmax.push_back(maxErrDiv);
      h.push_back(hi);
      if(i==0) {convu.push_back(0); convp.push_back(0);}
      else {
        convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
        convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
      }

      if(i==4) {
        nx = nx+5;
        ny = ny+5;
      } else {
        nx = 2*nx-1;
        ny = 2*ny-1;
      }
    }
    std::cout << "\n" << std::left
    << std::setw(10) << std::setfill(' ') << "h"
    << std::setw(15) << std::setfill(' ') << "err p"
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