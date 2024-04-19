#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "../util/cputime.h"
#include "cfmpi.hpp"

#include "finiteElement.hpp"
#include "baseProblem.hpp"
#include "paraview.hpp"
#include "../num/matlab.hpp"

// #include "../num/gnuplot.hpp"

// #define PROBLEM_FITTED_STOKESRT

// #define PROBLEM_FITTED_STOKES_VORTICITY

// #define POISSEUILLE_EXAMPLE

// #define DYNAMIC_DROP_EXAMPLE

// #define DYNAMIC_DROP_EXAMPLE_VORTICITY

// #define PROBLEM_UNFITTED_STOKESRT

// #define PROBLEM_UNFITTED_STOKESRT_SYMMETRIC_ALT // (2023 autumn)

// #define PROBLEM_UNFITTED_STOKES_VORTICITY

// #define PROBLEM_FITTED_CORIOLIS_STOKESRT

// need to add
#define PROBLEM_UNFITTED_CORIOLIS_STOKESRT // [not pressure robust with "Puppi" stab 2k-1/2k+1 or 2k+1/2k+1! But
// pressure robust with ours!]
// need to add
// #define PROBLEM_UNFITTED_CORIOLIS_STOKES_VORTICITY

// need to add
// #define PROBLEM_UNFITTED_PRESROB2_STOKES

// need to add
// #define PROBLEM_UNFITTED_PRESROB2_STOKES_VORTICITY

// if I have motivation
// #define PROBLEM_UNFITTED_STOKES3D

#ifdef PROBLEM_FITTED_STOKESRT

namespace Erik_Data_STOKESRT {

// R fun_rhs(const R2 P, int i) {
//   R x = P.x;
//   R y = P.y;
//    if(i==0) return -2*y;
//    else return 2*x;
// }
// R fun_exact_u(const R2 P, int i) {
//   R x = P.x;
//   R y = P.y;
//   // if(i==0) return      2*x*y*(x-1)*(y-1)*x*(1-x)*(2*y-1);
//   // else if(i==1) return 2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x-1);
//   if(i==0)      return  x*x*y;
//   else if(i==1) return -x*y*y;
// }
// R fun_exact_p(const R2 P, int i) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }

// [This example needs ability to compute mean of exact p, its not 0]
R fun_rhs(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 0;
    else
        return 0;
}
R fun_exact_u(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 20 * x * pow(y, 3);
    else if (i == 1)
        return 5 * pow(x, 4) - 5 * pow(y, 4);
}
R fun_exact_p(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    return 60 * pow(x, 2) * y - 20 * pow(y, 3);
}
R fun_div(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    return 0;
}
} // namespace Erik_Data_STOKESRT
using namespace Erik_Data_STOKESRT;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 10;
    int ny = 10;
    // int d = 2;
    // ProblemOption optionStokes;
    // optionStokes.order_space_element_quadrature_ = 9;
    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    for (int i = 0; i < 3; ++i) { // i<3 umfpack..

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Th(nx, ny, 0., 0., 1., 1.);

        Lagrange2 FEvelocity(4);
        Space Uh(Th, FEvelocity);

        // Space Vh(Th, FEvelocity);
        Space Vh(Th, DataFE<Mesh2>::BDM1);
        Space Qh(Th, DataFE<Mesh2>::P0);

        // const R hi = Th[0].lenEdge(0);
        const R hi           = 1. / (nx - 1);
        const R penaltyParam = 4e3 / hi; // 1e1/pow(hi,3);
        const R sigma        = 1e1;      // 1e1 or 1e-1// HIGHER 1e2,1e3,etc WORSENS THE SOL [??]

        Fun_h fh(Uh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h gh(Vh, fun_exact_u);
        // Fun_h gh(Qh, fun_div); // THIS IS IF div u != 0

        FEM<Mesh2> stokes(Vh);
        stokes.add(Qh);
        // FEM<Mesh2> stokes(Vh, optionStokes); stokes.add(Qh);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest u(Vh, 2, 0), p(Qh, 1, 0), v(Vh, 2, 0), q(Qh, 1, 0);
        R mu = 1;
        // const MeshParameter& h(Parameter::h);
        // const MeshParameter& invlEdge(Parameter::invmeas);

        {
            stokes.addBilinear(
                contractProduct(mu * grad(u), grad(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q), Th);
            stokes.addLinear(innerProduct(fh.expression(2), v), Th);
            stokes.addBilinear(-innerProduct(average(grad(u * t) * n, 0.5, 0.5), jump(v * t)) -
                                   innerProduct(jump(u * t), average(grad(v * t) * n, 0.5, 0.5)) +
                                   innerProduct(1. / hi * (sigma * jump(u * t)), jump(v * t)),
                               Th, INTEGRAL_INNER_EDGE_2D);
            stokes.addBilinear(-innerProduct(grad(u) * n, v) // natural
                                   - innerProduct(u, grad(v) * n) + innerProduct(penaltyParam * u, v)

                                   + innerProduct(p, v * n) // natural
                               ,
                               Th, INTEGRAL_BOUNDARY);
            stokes.addLinear(-innerProduct(gh.expression(2), grad(v) * n) +
                                 innerProduct(gh.expression(2), penaltyParam * v),
                             Th, INTEGRAL_BOUNDARY);
            // stokes.addBilinear(
            //   - innerProduct(grad(u)*n, v)  //natural
            //   - innerProduct(u, grad(v)*n)
            //   + innerProduct(penaltyParam*u*t, v*t)

            //   // + innerProduct(p, v*n)  // natural
            //   , Th
            //   , INTEGRAL_BOUNDARY
            // );
            // stokes.addLinear(
            //   - innerProduct(gh.expression(2), grad(v)*n)
            //   + innerProduct(gh*t, penaltyParam*v*t)
            //   , Th
            //   , INTEGRAL_BOUNDARY
            // );
            // stokes.setDirichlet(gh,Th);
            // R meanP = integral(Th,exactp,0);
            stokes.addLagrangeMultiplier(innerProduct(1., p), 0, Th);
        }

        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        // nx = 2*nx;
        // ny = 2*ny;
        // continue;
        stokes.solve("mumps");

        // EXTRACT SOLUTION
        int idx0_s  = Vh.get_nb_dof();
        Rn_ data_uh = stokes.rhs_(SubArray(Vh.get_nb_dof(), 0));
        Rn_ data_ph = stokes.rhs_(SubArray(Qh.get_nb_dof(), idx0_s));
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Qh, data_ph);
        ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

        Space Rh(Th, DataFE<Mesh>::P2);
        Fun_h pph(Rh, fun_exact_p);
        ExpressionFunFEM<Mesh> tru_p(pph, 0, op_id);
        R meanP = integral(Th, tru_p);
        ExpressionFunFEM<Mesh> fem_p(ph, 0, op_id);
        R meanPfem = integral(Th, fem_p);
        // std::cout << meanP << std::endl;

        ph.v -= meanPfem;
        ph.v += meanP;

        {
            Fun_h soluErr(Vh, fun_exact_u);
            Fun_h soluh(Vh, fun_exact_u);
            soluErr.v -= uh.v;
            soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);
            Paraview<Mesh> writer(Th, "stokes_" + to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(femSol_0dx + femSol_1dy, "divergence");
            writer.add(soluh, "velocityExact", 0, 2);
            writer.add(soluErr, "velocityError", 0, 2);
            // writer.add(solh, "velocityError" , 0, 2);
            // writer.add(fabs(femDiv, "divergenceError");
        }

        // Rn solExVec(mixedSpace.NbDoF());
        // interpolate(mixedSpace, solExVec, fun_exact);
        // R errU = L2norm(sol,fun_exact,0,2);
        // R errP = L2norm(sol,fun_exact,2,1);

        R errU      = L2norm(uh, fun_exact_u, 0, 2);
        R errP      = L2norm(ph, fun_exact_p, 0, 1);
        R errDiv    = L2norm(femSol_0dx + femSol_1dy, Th);
        R maxErrDiv = maxNorm(femSol_0dx + femSol_1dy, Th);

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
        h.push_back(1. / nx);
        if (i == 0) {
            convu.push_back(0);
            convp.push_back(0);
        } else {
            convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
        }

        nx = 2 * nx;
        ny = 2 * ny;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ') << "err_p"
              << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ') << "err u"
              << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
                  << pl2[i] << std::setw(15) << std::setfill(' ') << convp[i] << std::setw(15) << std::setfill(' ')
                  << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
#endif

#ifdef PROBLEM_FITTED_STOKES_VORTICITY

namespace Erik_Data_FITTED_STOKES_VORTICITY {
R fun_rhs(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 0;
    else
        return 0;
}
R fun_exact_u(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 20 * x * pow(y, 3);
    else if (i == 1)
        return 5 * pow(x, 4) - 5 * pow(y, 4);
}
R fun_exact_p(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    return 20 * (3 * pow(x, 2) * y - pow(y, 3)); //-2.21664/sqrt(M_PI*interfaceRad*interfaceRad);
}
R fun_div(const R2 P, int i) {
    R x = P.x;
    R y = P.y;
    return 0;
}
} // namespace Erik_Data_FITTED_STOKES_VORTICITY
using namespace Erik_Data_FITTED_STOKES_VORTICITY;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 10;
    int ny = 10;
    // int d = 2;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 5;
    for (int i = 0; i < iters; ++i) {

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, 0., 0., 1., 1.);
        const R hi           = 1. / (nx - 1);
        const R penaltyParam = 1e3 / hi; // interesting: 1e-1/sqrt(hi) gets conv 1/1 down to 0.00625, otherwise 1e-1/hi
                                         // p variable goes down to 0.5

        Lagrange2 FEvelocity(2);
        Space VELh(Kh, FEvelocity);
        Space SCAh(Kh, DataFE<Mesh>::P2);

        Space Uh(Kh, DataFE<Mesh>::P1);
        Space Vh(Kh, DataFE<Mesh2>::RT0);
        Space Wh(Kh, DataFE<Mesh2>::P0);

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h u0(VELh, fun_exact_u);
        Fun_h p0(SCAh, fun_exact_p);
        ExpressionFunFEM<Mesh> exactp(p0, 0, op_id);

        CutFEM<Mesh2> stokes(Uh);
        stokes.add(Vh);
        stokes.add(Wh);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest w(Uh, 1, 0), tau(Uh, 1, 0), u(Vh, 2, 0), v(Vh, 2, 0), p(Wh, 1, 0), q(Wh, 1, 0);

        R mu = 1;
        {
            // [Bulk]
            stokes.addBilinear(                                                // w = curl u
                innerProduct(1. / mu * w, tau) - innerProduct(u, rotgrad(tau)) // + -> w = -curl u
                ,
                Kh);
            stokes.addBilinear(             // -mu Delta u + grad p
                innerProduct(rotgrad(w), v) // + -> -
                    - innerProduct(p, div(v)),
                Kh);
            stokes.addLinear(innerProduct(fh.expression(2), v), Kh);
            stokes.addBilinear(innerProduct(div(u), q), Kh);
            // [Velocity Dirichlet BC]
            {
                stokes.addBilinear(
                    // + innerProduct(p, v*n)
                    // + innerProduct(u*t, tau)
                    +innerProduct(u * n, penaltyParam * v * n), Kh, INTEGRAL_BOUNDARY);
                stokes.addLinear(              // [TANGENT/NORMAL POINTS IN OTHER DIRECTION FOR FITTED VS UNFITTED?!?!]
                                               // + innerProduct(u0*t, tau)
                    -innerProduct(u0 * t, tau) // [moved over to RHS like a Neuman condition..?]
                        + innerProduct(u0 * n, penaltyParam * v * n),
                    Kh, INTEGRAL_BOUNDARY);
                // [Sets uniqueness of the pressure]
                // R meanP = integral(Kh,exactp,0);
                // stokes.addLagrangeMultiplier(
                //   innerProduct(1.,p), 0
                //   , Kh
                // );
                // [Sets uniqueness of the pressure in another way such that divu = 0]
                FEM<Mesh2> lagr(Uh);
                lagr.add(Vh);
                lagr.add(Wh);
                Rn zero_vec = lagr.rhs_;
                lagr.addLinear(innerProduct(1., p), Kh);
                Rn lag_row(lagr.rhs_);
                lagr.rhs_ = zero_vec;
                lagr.addLinear(innerProduct(1, v * n), Kh, INTEGRAL_BOUNDARY);
                stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);
            }
            // [Vorticity BC 2 (natural)]
            {
                // stokes.addLinear(
                //   - innerProduct(u0*t,tau)              // + -> +
                //   - innerProduct(p0.expression(), v*n)
                //   , Kh, INTEGRAL_BOUNDARY
                // );
            }
        }

        // [grad -> eps]
        {
            // stokes.addBilinear(
            //   contractProduct(2*mu*Eps(u),Eps(v))
            //   - innerProduct(p,div(v))
            //   + innerProduct(div(u),q)
            //   , Khi
            // );
            // stokes.addLinear(
            //   innerProduct(fh.expression(2),v)
            //   , Khi
            // );
            // stokes.addBilinear(
            //   - innerProduct(average(2*mu*(Eps(u)*t)*n,0.5,0.5), jump(v*t))
            //   + innerProduct(jump(u*t), average(2*mu*(Eps(v)*t)*n,0.5,0.5))
            //   + innerProduct(1./hi*sigma*(jump(u*t)), jump(v*t))
            //   , Khi
            //   , INTEGRAL_INNER_EDGE_2D
            // );
            // stokes.addBilinear(
            //   - innerProduct(2*mu*Eps(u)*n, v)  //natural
            //   + innerProduct(u, 2*mu*Eps(v)*n)
            //   + innerProduct(penaltyParam*u, v)

            //   // + innerProduct(p, v*n)  // natural
            //   // - innerProduct(u*n, q)  // essential

            //   // + innerProduct(penaltyParam*u*n, v*n)
            //   , interface
            // );
            // stokes.addLinear(
            //   + innerProduct(gh.expression(2), 2*mu*Eps(v)*n)
            //   + innerProduct(gh.expression(2), penaltyParam*v)

            //   // - innerProduct(gh.expression(2), q*n)

            //   // + innerProduct(eps()*n, 2*mu*v) // HOW TO EPS
            //   // + innerProduct(u*n, penaltyParam*v*n)
            //   , interface
            // );

            // FunTest Eps2un = (grad(Eps(u)*n)+grad(Eps(u)*n).t())*n;
            // FunTest grad2pn = grad(grad(p)*n)*n;
            // // FunTest grad2divun = grad(grad(div(u))*n)*n;
            // stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
            //   //  innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v)) // [Method 1: Remove jump in vel]
            //   // +innerProduct(uPenParam*pow(hi,3)*jump(Eps(u)*n), jump(Eps(v)*n))
            //   // +innerProduct(uPenParam*pow(hi,5)*jump(Eps2un), jump(Eps2un))
            //   // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
            //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

            //   +innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v)) // [Method 1: Remove jump in vel]
            //   +innerProduct(uPenParam*pow(hi,3)*jump(Eps(u)*n), jump(Eps(v)*n))
            //   // +innerProduct(uPenParam*pow(hi,5)*jump(Eps2un), jump(Eps2un))

            //   -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
            //   +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
            //   -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
            //   +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(v))) , jump(grad(q)))
            //   // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
            //   // +innerProduct(pPenParam*pow(hi,5)*jump(grad2divun) , jump(grad2pn))
            //   , Khi
            // );

            // // R meanP = integral(Khi,exactp,0); std::cout << meanP << std::endl;
            // // stokes.addLagrangeMultiplier(
            // //   innerProduct(1.,p), 0
            // //   // innerProduct(1.,p)+innerProduct(1.,div(u)), meanP
            // //   , Khi, INTEGRAL_EXTENSION, 1
            // // );
            // // // std::cout << Vh.get_nb_dof() << " " << Ph.get_nb_dof() << std::endl;
            // // // stokes.mat_[std::make_pair(50,Vh.get_nb_dof()+1)] = 1e30;
        }

        // std::cout << integral(Khi,exactp,0) << std::endl;
        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_vort_dof = Uh.get_nb_dof();
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_wh     = stokes.rhs_(SubArray(nb_vort_dof, 0));
        Rn_ data_uh     = stokes.rhs_(SubArray(
            nb_flux_dof, nb_vort_dof)); // Rn_ data_uh = stokes.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));
        Rn_ data_ph     = stokes.rhs_(SubArray(
            Wh.get_nb_dof(),
            nb_vort_dof +
                nb_flux_dof)); // Rn_ data_ph = stokes.rhs_(SubArray(stokes_.get_nb_dof(),nb_vort_dof+nb_flux_dof));
        Fun_h wh(Uh, data_wh);
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Wh, data_ph);

        // [Post process pressure]
        R meanP = integral(Kh, exactp, 0);
        ExpressionFunFEM<Mesh> fem_p(ph, 0, op_id);
        R meanPfem = integral(Kh, fem_p, 0);
        // std::cout << meanP << std::endl;
        CutFEM<Mesh2> post(Wh);
        post.addLinear(innerProduct(1, q), Kh);
        R area = post.rhs_.sum();
        ph.v -= meanPfem / area;
        ph.v += meanP / area;

        ExpressionFunFEM<Mesh> dx_uh0(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> dy_uh1(uh, 1, op_dy);

        // [Errors]
        {
            // Fun_h solw(Uh, fun_exact_w);
            Fun_h solu(Vh, fun_exact_u);
            Fun_h soluErr(Vh, fun_exact_u);
            Fun_h solp(Wh, fun_exact_p);
            soluErr.v -= uh.v;
            soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Kh, "stokes_" + to_string(i) + ".vtk");
            writer.add(wh, "vorticity", 0, 1);
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(dx_uh0 + dy_uh1, "divergence");
            // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
            writer.add(solu, "velocityExact", 0, 2);
            writer.add(soluErr, "velocityError", 0, 2);
            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        // R errW      = L2normCut(wh,fun_exact_w,0,1);
        R errU      = L2norm(uh, fun_exact_u, 0, 2);
        R errP      = L2norm(ph, fun_exact_p, 0, 1);
        R errDiv    = L2norm(dx_uh0 + dy_uh1, Kh);
        R maxErrDiv = maxNorm(dx_uh0 + dy_uh1, Kh);
        // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

        ul2.push_back(errU);
        pl2.push_back(errP);
        divl2.push_back(errDiv);
        divmax.push_back(maxErrDiv);
        h.push_back(1. / nx);
        if (i == 0) {
            convu.push_back(0);
            convp.push_back(0);
        } else {
            convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
        }

        nx *= 2;
        ny *= 2;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ') << "err p"
              << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ') << "err u"
              << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
                  << pl2[i] << std::setw(15) << std::setfill(' ') << convp[i] << std::setw(15) << std::setfill(' ')
                  << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
#endif

#ifdef POISSEUILLE_EXAMPLE

namespace Data_POISSEUILLE {
double mu1 = 1;
double mu2 = 1;
R2 shift(0., 0.);
double sigma = 1;
R fun_levelSet(const R2 P, int i) {
    return sqrt((P.x - shift.x) * (P.x - shift.x) + (P.y - shift.y) * (P.y - shift.y)) - 1.00001;
}

R fun_exact_p(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 0.;
}
} // namespace Data_POISSEUILLE
using namespace Data_POISSEUILLE;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 10;
    int ny = 10;
    // int d = 2;

    ProblemOption optionStokes;
    optionStokes.order_space_element_quadrature_ = 5;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    for (int i = 0; i < 4; ++i) { // i<3

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, -3., -2., 6., 4.);
        const R hi           = 2. / (nx - 1);
        const R penaltyParam = 1e1 / hi;
        // const R sigma = 1;
        double uPenParam     = 1e0;
        double pPenParam     = 1e0;
        double jumpParam     = 1e0;

        // CutFEMParameter mu(1e-1,1e-3);
        // CutFEMParameter invmu(1e1,1e3);
        CutFEMParameter mu(mu1, mu2);
        // double sigma = 0;//24.5;//700;
        double kappa1 = 0.5; // mu(1)/(mu(0)+mu(1));
        double kappa2 = 0.5; // mu(0)/(mu(0)+mu(1));

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(2);
        Space UUh(Kh, FEvelocity);

        // Space Qh(Kh, DataFE<Mesh>::P1);

        Space Wh(Kh, DataFE<Mesh2>::RT1); // REMOVE KhESE TWO LATER
        // Space Wh(Kh, FEvelocity);
        Space Qh(Kh, DataFE<Mesh2>::P1dc); // FOR MIXEDSPACE

        ActiveMesh<Mesh> Khi(Kh, interface);

        ActiveMesh<Mesh> Kh0(Kh);
        Kh0.createSurfaceMesh(interface);
        ActiveMesh<Mesh> Kh1(Kh);
        Kh1.truncate(interface, -1);
        ActiveMesh<Mesh> Kh2(Kh);
        Kh2.truncate(interface, 1);

        // Khi.truncate(interface, 1);
        CutSpace Vh(Khi, Wh);
        CutSpace Ph(Khi, Qh);
        CutSpace Uh(Khi, UUh);

        // Fun_h fh(Uh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        // Fun_h gh(Uh, fun_exact_u);
        // Fun_h gh(Qh, fun_div); // THIS IS IF div u != 0

        CutFEM<Mesh2> stokes(Vh, optionStokes);
        stokes.add(Ph);

        Normal n;
        Tangent t;
        FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0), p1(Ph, 1, 0, 0);

        stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                           Khi);

        stokes.addBilinear(-innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2), jump(v * t)) +
                               innerProduct(jump(u * t), average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
                               innerProduct(1. / hi * (jump(u * t)), jump(v * t)),
                           Khi, INTEGRAL_INNER_EDGE_2D);
        stokes.addBilinear(innerProduct(jump(u), -2 * mu * average(Eps(v) * n, kappa1, kappa2)) +
                               innerProduct(-2 * mu * average(Eps(u) * n, kappa1, kappa2), jump(v)) +
                               innerProduct(1. / hi / hi * jump(u), jump(v)) +
                               innerProduct(average(p, kappa1, kappa2), jump(v * n))
                           // - innerProduct(jump(u*n), average(q,0.5,0.5))
                           ,
                           interface);
        stokes.addLinear(innerProduct(3. / 2, average(v * n, kappa2, kappa1)) * sigma, interface);

        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v)                               // natural
                               + innerProduct(u, 2 * mu * Eps(v) * n) + innerProduct(p, v * n) // natural
                               + innerProduct(penaltyParam * u, v)
                           // - innerProduct(u*n, q)  // essential
                           ,
                           Khi, INTEGRAL_BOUNDARY, {1, 3});
        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) // natural
                                                                 // + innerProduct(u, 2*mu*Eps(v)*n)
                                                                 // + innerProduct(p, v*n)  // natural
                                                                 // + innerProduct(penaltyParam*u, v)
                                                                 // - innerProduct(u*n, q)  // essential
                           ,
                           Khi, INTEGRAL_BOUNDARY, {2, 4});
        stokes.addLinear(innerProduct(1, v * n) // natural
                         ,
                         Khi, INTEGRAL_BOUNDARY, {4});
        stokes.addLinear(innerProduct(-1, v * n) // natural
                         ,
                         Khi, INTEGRAL_BOUNDARY, {2});

        // l(v)_Omega
        // stokes.addLinear(
        //   innerProduct(fh.expression(2),v)
        //   , Khi
        // );

        FunTest grad2un = grad(grad(u) * n) * n;
        // stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
        //  innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v)) // [Method 1: Remove jump in vel]
        // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
        // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
        // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
        // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

        // // +innerProduct(uPenParam*pow(hi,1)*jump(u), mu*jump(v)) // [Method 1: Remove jump in vel]
        // // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), mu*jump(grad(v)*n))
        // // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))

        // // -innerProduct(pPenParam*hi*jump(p), (1./mu)*jump(div(v)))
        // // +innerProduct(pPenParam*hi*jump(div(u)), (1./mu)*jump(q))
        // // -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), (1./mu)*jump(grad(div(v))))
        // // +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))) , (1./mu)*jump(grad(q)))

        // // +innerProduct(uPenParam*pow(hi,3)*jump(div(u)), mu*jump(div(v))) // [Method 1: Remove jump in vel]
        // // +innerProduct(uPenParam*pow(hi,5)*jump(grad(div(u))), mu*jump(grad(div(v))))
        //   , Khi
        // // , macro
        // );

        // Sets uniqueness of the pressure
        // double tt0 = MPIcf::Wtime();
        // int N = stokes.get_nb_dof();
        // std::map<int, double> df2rm;
        // R2 P = Qh[0].Pt(0);
        // double val = fun_exact_p(P, 0, 0);
        // df2rm.insert({Vh.get_nb_dof(), val});
        //
        // eraseRow(N, stokes.mat_, stokes.rhs_, df2rm);
        // std::cout << "Time setting p val " << MPIcf::Wtime() - tt0 << std::endl;
        // stokes.addLagrangeMultiplier(
        //   innerProduct(1.,p), 0., Khi
        // );

        stokes.solve();

        // EXTRACT SOLUTION
        int idx0_s  = Vh.get_nb_dof();
        Rn_ data_uh = stokes.rhs_(SubArray(Vh.get_nb_dof(), 0));
        Rn_ data_ph = stokes.rhs_(SubArray(Ph.get_nb_dof(), idx0_s));
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Ph, data_ph);
        ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

        {
            // Fun_h soluErr(Vh, fun_exact_u);
            // Fun_h soluh(Vh, fun_exact_u);
            // soluErr.v -= uh.v;
            // soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Khi, "poisseuilleSp_" + to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(fabs(femSol_0dx + femSol_1dy), "divergence");
            // writer.add(soluh, "velocityExact" , 0, 2);
            // writer.add(soluErr, "velocityError" , 0, 2);

            // writer.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceError");

            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        R errU      = 0.; // L2normCut(uh,fun_exact_u,0,2);
        R errP      = 0.; // L2normCut(ph,fun_exact_p,0,1);
        R errDiv1   = L2normCut(femSol_0dx + femSol_1dy, Khi, 1);
        R errDiv0   = L2normCut(femSol_0dx + femSol_1dy, Khi, 0);
        R maxErrDiv = maxNormCut(femSol_0dx + femSol_1dy, Khi);

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
        h.push_back(1. / nx);
        if (i == 0) {
            convu.push_back(0);
            convp.push_back(0);
        } else {
            convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
        }

        nx *= 2;
        ny *= 2;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ')
              << "err_p u" << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ')
              << "err u" << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
              << "err divu0" << std::setw(15) << std::setfill(' ')
              << "err divu1"

              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              // << std::setw(15) << std::setfill(' ') << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
                  << pl2[i] << std::setw(15) << std::setfill(' ') << convp[i] << std::setw(15) << std::setfill(' ')
                  << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
#endif

#ifdef DYNAMIC_DROP_EXAMPLE
namespace DynamicDropData {
R2 shift(0, 0);
double sigma = 0;
double mu1   = 1;
double mu2   = 1;

R fun_levelSet(const R2 P, int i) {
    return sqrt((P.x - shift.x) * (P.x - shift.x) + (P.y - shift.y) * (P.y - shift.y)) - 2. / 3;
}

static R falpha1(const R r) {
    // const R MU1 = 1e-1;
    // const R MU2 = 1e-3;
    const R MU1 = mu1;
    const R MU2 = mu2;
    const R r0  = 2. / 3;
    return 1. / MU1 + (1. / MU2 - 1. / MU1) * exp(r * r - r0 * r0);
}
static R falpha2(const R r) {
    R MU2 = mu2; // 1e-3;
    return 1. / MU2;
}

static R falpha(const R r) {
    const R r0 = 2. / 3;
    return (r < r0) ? falpha2(r) : falpha1(r);
}

static R fun_rhs(const R2 P, int i) {
    const R r2 = Norme2_2(P - shift);
    const R s  = exp(-r2);
    R2 R(4 * s * (r2 - 2) * P.y + 3 * P.x * P.x, -4 * s * (r2 - 2) * P.x);
    return (i < 2) ? R[i] : 0;
}
static R fun_boundary(const R2 P, const int i) {
    const R r = Norme2(P - shift);
    R2 R(-P.y, P.x);
    R = falpha(r) * exp(-r * r) * R;
    return (i < 2) ? R[i] : 0;
}

static R2 fun_velocity1(const R2 P) {
    R r = Norme2(P - shift);
    R2 R(-P.y, P.x);
    R = falpha1(r) * exp(-r * r) * R;
    return R;
}
static R2 fun_velocity2(const R2 P) {
    R r = Norme2(P - shift);
    R2 R(-P.y, P.x);
    R = falpha2(r) * exp(-r * r) * R;
    return R;
}
static R fun_pressure1(const R2 P) { return pow(P.x, 3); }
static R fun_pressure2(const R2 P) {
    // R sigma = 1;//24.5;//700;
    return pow(P.x, 3) + sigma * 3. / 2.;
}

// static R fun_solution(const R2 P, int ci, int domain) {
//   if(domain == 0) return (ci < 2)? fun_velocity1(P)[ci] : fun_pressure1(P);
//   else return (ci < 2)? fun_velocity2(P)[ci] : fun_pressure2(P);
// }
static R fun_exact_u(const R2 P, int ci, int domain) {
    if (domain == 0)
        return fun_velocity1(P)[ci];
    else
        return fun_velocity2(P)[ci];
}
static R fun_exact_p(const R2 P, int ci, int domain) {
    if (domain == 0)
        return fun_pressure1(P);
    else
        return fun_pressure2(P);
}

R2 fparam(double t) { return R2(2. / 3 * cos(t + 1. / 3), 2. / 3 * sin(t + 1. / 3)); }

} // namespace DynamicDropData
namespace Erik_Data_StokesDiv {
double mu1   = 1;
double mu2   = 1;
double sigma = 0;

// R fun_rhs(const R2 P, int i) {
//   return 0;
// }
// R fun_exact(const R2 P, int i) {
//   if(i==0) return 0;
//   else if(i==1) return P.x;
//   else return 0;
// }
R shift        = 0.5;
R interfaceRad = 0.250001; // 2./3; // not exactly 1/4 to avoid interface cutting exaclty a vertex
R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x - shift) * (P.x - shift) + (P.y - shift) * (P.y - shift)) - interfaceRad;
}

R fun_rhs(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return -2 * y;
    else
        return 2 * x;
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
    if (i == 0)
        return x * x * y;
    else if (i == 1)
        return -x * y * y;
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
//   if(i==0) return
//   -2*mu*(pow(x,4)*(6*y-3)+pow(x,3)*(6-12*y)+3*pow(x,2)*(4*pow(y,3)-6*pow(y,2)+4*y-1)-6*x*y*(2*pow(y,2)-3*y+1)+y*(2*pow(y,2)-3*y+1))
//   + 10*(3*pow(x-0.5,2)*pow(y,2)-3*pow(1-x,2)*pow(y-0.5,3)); else if(i==1) return
//   2*mu*(2*pow(x,3)*(6*pow(y,2)-6*y+1)-3*pow(x,2)*(6*pow(y,2)-6*y+1)+x*(6*pow(y,4)-12*pow(y,3)+12*pow(y,2)-6*y+1)-3*pow(y-1,2)*pow(y,2))
//   + 10*(2*pow(x-0.5,3)*y+3*pow(1-x,3)*pow(y-0.5,2)); else return 0;
// }
// R fun_exact(const R2 P, int i) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return pow(x,2)*pow(1-x,2)*y*(1-y)*(1-2*y);
//   else if(i==1) return (-x)*(1-x)*(1-2*x)*pow(y,2)*pow(1-y,2);
//   else return 10*(pow(x-0.5,3)*pow(y,2)+pow(1-x,3)*pow(y-0.5,3));
// }
} // namespace Erik_Data_StokesDiv
// using namespace Erik_Data_StokesDiv;
using namespace DynamicDropData;

// class LambdaGamma : public VirtualParameter { public :
// const CutFEMParameter& mu_;
//
// LambdaBoundary(const CutFEM_Parameter& mu, double G, double H) : mu_(mu) , G_(G) , H_(H) {} double evaluate(int
// domain, double h, double meas, double measK, double measCut) const { double gamma = meas / h ; double alpha = measK /
// h / h; double val = mu_.evaluate(domain,h,meas,measK,measCut) return val/h(G_+H_*gamma/alpha) ;}
// };

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 10;
    int ny = 10;
    // int d = 2;

    ProblemOption optionStokes;
    optionStokes.order_space_element_quadrature_ = 9;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    for (int i = 0; i < 3; ++i) { // i<3

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, -1., -1., 2., 2.);
        const R hi           = 2. / (nx - 1);
        const R penaltyParam = 1e-1 / hi;
        // const R sigma = 1;
        double uPenParam     = 1e-1;
        double pPenParam     = 1e-1;
        double jumpParam     = 1e-1;

        // CutFEMParameter mu(1e-1,1e-3);
        // CutFEMParameter invmu(1e1,1e3);
        CutFEMParameter mu(mu1, mu2);
        // double sigma = 0;//24.5;//700;
        double kappa1 = mu(1) / (mu(0) + mu(1));
        double kappa2 = mu(0) / (mu(0) + mu(1));

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(2);
        Space UUh(Kh, FEvelocity);

        Space Wh(Kh, DataFE<Mesh2>::BDM1); // REMOVE KhESE TWO LATER
        // Space Wh(Kh, FEvelocity);
        Space Qh(Kh, DataFE<Mesh2>::P0); // FOR MIXEDSPACE

        ActiveMesh<Mesh> Khi(Kh, interface);

        MacroElement<Mesh> macro(Khi, 1);

        CutSpace Vh(Khi, Wh);
        CutSpace Ph(Khi, Qh);
        CutSpace Uh(Khi, UUh);

        Fun_h fh(Uh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h gh(Vh, fun_exact_u);

        CutFEM<Mesh2> stokes(Vh, optionStokes);
        stokes.add(Ph);
        CutSpace Ph2(Khi, Qh); // stokes.add(Ph2); FunTest phi(Ph2,1,0), psi(Ph2,1,0);

        Normal n;
        Tangent t;
        FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0), p1(Ph, 1, 0, 0);

        stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                           Khi);
        stokes.addLinear(+innerProduct(fh.expression(2), v), Khi);
        stokes.addBilinear(-innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2), jump(v * t)) -
                               innerProduct(jump(u * t), average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
                               innerProduct(jumpParam * 1. / hi * (jump(u * t)), jump(v * t)),
                           Khi, INTEGRAL_INNER_EDGE_2D);
        stokes.addBilinear(-innerProduct(2 * average(mu * Eps(u) * n, kappa1, kappa2), jump(v)) // [(2muEps(u)-pI)n]=0
                               - innerProduct(jump(u), 2 * average(mu * Eps(v) * n, kappa1,
                                                                   kappa2)) // sym term = 0 but in Dirichlet fashion
                               + innerProduct(average(p, kappa1, kappa2), jump(v * n)) // [(2muEps(u)-pI)n]=0
                               // - innerProduct(jump(u*n), average(q,0.5,0.5))                  // sym term, ignored
                               + innerProduct(jumpParam * 1. / hi * jump(u), jump(v)) // [u]=0
                           ,
                           interface);
        stokes.addLinear(+innerProduct(3. / 2, average(v * n, kappa2, kappa1)) * sigma, interface);

        // stokes.addBilinear(
        //   - innerProduct(2*mu*Eps(u)*n, v)  //natural
        //   - innerProduct(u, 2*mu*Eps(v)*n)
        //   + innerProduct(p, v*n)  // natural
        //   + innerProduct(penaltyParam*u, v)
        //   , Khi
        //   , INTEGRAL_BOUNDARY
        // );
        // stokes.addLinear(
        //   - innerProduct(gh.expression(2), 2*mu*Eps(v)*n)
        //   + innerProduct(gh.expression(2), penaltyParam*v)
        //   , Khi
        //   , INTEGRAL_BOUNDARY
        // );
        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) // natural
                               - innerProduct(u, 2 * mu * Eps(v) * n)
                               // + innerProduct(p, v*n)  // natural
                               + innerProduct(penaltyParam * u * t, v * t),
                           Khi, INTEGRAL_BOUNDARY);
        stokes.addLinear(-innerProduct(gh.expression(2), 2 * mu * Eps(v) * n) +
                             innerProduct(gh * t, penaltyParam * v * t),
                         Khi, INTEGRAL_BOUNDARY);
        stokes.setDirichlet(gh, Khi);

        FunTest grad2un = grad(grad(u) * n) * n;
        stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
                                     //  innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v)) // [Method 1: Remove jump in
                                     //  vel]
                                     // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
                                     // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
                                     // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
                                     // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

            +innerProduct(uPenParam * pow(hi, -1) * jump(u), jump(v)) // [Method 1: Remove jump in vel]
                + innerProduct(uPenParam * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
                innerProduct(uPenParam * pow(hi, 3) * jump(grad2un), jump(grad2un))

                - innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
                innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
                innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(v))), jump(grad(q)))

            // +innerProduct(uPenParam*pow(hi,3)*jump(div(u)), mu*jump(div(v))) // [Method 1: Remove jump in vel]
            // +innerProduct(uPenParam*pow(hi,5)*jump(grad(div(u))), mu*jump(grad(div(v))))

            // -innerProduct(pPenParam*hi*jump(phi), (1./mu)*jump(div(v)))
            // +innerProduct(pPenParam*hi*jump(div(u)), (1./mu)*jump(psi))
            // +innerProduct(pPenParam*pow(hi,1)*jump(phi), jump(q))
            // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(psi))

            ,
            Khi, macro);

        // [Set uniqueness of pressure]
        stokes.addLagrangeMultiplier(innerProduct(1., p1), 0., Khi);
        // CutFEM<Mesh2> lagr(Vh); lagr.add(Ph);
        // Rn zero_vec = lagr.rhs_;
        // lagr.addLinear(
        //   innerProduct(1.,p)
        //   , Khi
        // );
        // Rn lag_row(lagr.rhs_); lagr.rhs_ = zero_vec;
        // lagr.addLinear(
        //   innerProduct(1, v*n)
        //   , Khi, INTEGRAL_BOUNDARY
        // );
        // stokes.addLagrangeVecToRowAndCol(lag_row,lagr.rhs_,0);

        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int idx0_s  = Vh.get_nb_dof();
        Rn_ data_uh = stokes.rhs_(SubArray(Vh.get_nb_dof(), 0));
        Rn_ data_ph = stokes.rhs_(SubArray(Ph.get_nb_dof(), idx0_s));
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Ph, data_ph);
        ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

        {
            Fun_h soluErr(Vh, fun_exact_u);
            Fun_h soluh(Vh, fun_exact_u);
            soluErr.v -= uh.v;
            soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Khi, "stokes_" + to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(fabs(femSol_0dx + femSol_1dy), "divergence");
            writer.add(soluh, "velocityExact", 0, 2);
            writer.add(soluErr, "velocityError", 0, 2);
            // writer.add(fabs((femSol_0dx+femSol_1dy)-femDiv), "divergenceError");

            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        R errU      = L2normCut(uh, fun_exact_u, 0, 2);
        R errP      = L2normCut(ph, fun_exact_p, 0, 1);
        R errDiv1   = L2normCut(femSol_0dx + femSol_1dy, Khi, 1);
        R errDiv0   = L2normCut(femSol_0dx + femSol_1dy, Khi, 0);
        R maxErrDiv = maxNormCut(femSol_0dx + femSol_1dy, Khi);

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
        // divmax.push_back(errDiv1);
        divmax.push_back(maxErrDiv);

        h.push_back(1. / nx);
        if (i == 0) {
            convu.push_back(0);
            convp.push_back(0);
        } else {
            convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
        }

        nx *= 2;
        ny *= 2;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ') << "err p"
              << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ') << "err u"
              << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
              << "err divu0"
              // << std::setw(15) << std::setfill(' ') << "err divu1"

              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
                  << pl2[i] << std::setw(15) << std::setfill(' ') << convp[i] << std::setw(15) << std::setfill(' ')
                  << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
#endif

#ifdef DYNAMIC_DROP_EXAMPLE_VORTICITY
namespace DynamicDropData {
R2 shift(0, 0);
double sigma = 0;
double mu1   = 1;
double mu2   = 1;

R fun_levelSet(const R2 P, int i) {
    return sqrt((P.x - shift.x) * (P.x - shift.x) + (P.y - shift.y) * (P.y - shift.y)) - 2. / 3;
}

static R falpha1(const R r) {
    // const R MU1 = 1e-1;
    // const R MU2 = 1e-3;
    const R MU1 = mu1;
    const R MU2 = mu2;
    const R r0  = 2. / 3;
    return 1. / MU1 + (1. / MU2 - 1. / MU1) * exp(r * r - r0 * r0);
}
static R falpha2(const R r) {
    R MU2 = mu2; // 1e-3;
    return 1. / MU2;
}

static R falpha(const R r) {
    const R r0 = 2. / 3;
    return (r < r0) ? falpha2(r) : falpha1(r);
}

static R fun_rhs(const R2 P, int i) {
    const R r2 = Norme2_2(P - shift);
    const R s  = exp(-r2);
    R2 R(4 * s * (r2 - 2) * P.y + 3 * P.x * P.x, -4 * s * (r2 - 2) * P.x);
    return (i < 2) ? R[i] : 0;
}
static R fun_boundary(const R2 P, const int i) {
    const R r = Norme2(P - shift);
    R2 R(-P.y, P.x);
    R = falpha(r) * exp(-r * r) * R;
    return (i < 2) ? R[i] : 0;
}

static R2 fun_velocity1(const R2 P) {
    R r = Norme2(P - shift);
    R2 R(-P.y, P.x);
    R = falpha1(r) * exp(-r * r) * R;
    return R;
}
static R2 fun_velocity2(const R2 P) {
    R r = Norme2(P - shift);
    R2 R(-P.y, P.x);
    R = falpha2(r) * exp(-r * r) * R;
    return R;
}
static R fun_pressure1(const R2 P) { return pow(P.x, 3); }
static R fun_pressure2(const R2 P) {
    // R sigma = 1;//24.5;//700;
    return pow(P.x, 3) + sigma * 3. / 2.;
}

// static R fun_solution(const R2 P, int ci, int domain) {
//   if(domain == 0) return (ci < 2)? fun_velocity1(P)[ci] : fun_pressure1(P);
//   else return (ci < 2)? fun_velocity2(P)[ci] : fun_pressure2(P);
// }
static R fun_exact_u(const R2 P, int ci, int domain) {
    if (domain == 0)
        return fun_velocity1(P)[ci];
    else
        return fun_velocity2(P)[ci];
}
static R fun_exact_p(const R2 P, int ci, int domain) {
    if (domain == 0)
        return fun_pressure1(P);
    else
        return fun_pressure2(P);
}

R2 fparam(double t) { return R2(2. / 3 * cos(t + 1. / 3), 2. / 3 * sin(t + 1. / 3)); }

} // namespace DynamicDropData
namespace Erik_Data_StokesDiv {
double mu1   = 1;
double mu2   = 1;
double sigma = 0;

// R fun_rhs(const R2 P, int i) {
//   return 0;
// }
// R fun_exact(const R2 P, int i) {
//   if(i==0) return 0;
//   else if(i==1) return P.x;
//   else return 0;
// }
R shift        = 0.5;
R interfaceRad = 0.250001; // 2./3; // not exactly 1/4 to avoid interface cutting exaclty a vertex
R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x - shift) * (P.x - shift) + (P.y - shift) * (P.y - shift)) - interfaceRad;
}

R fun_rhs(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return -2 * y;
    else
        return 2 * x;
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
    if (i == 0)
        return x * x * y;
    else if (i == 1)
        return -x * y * y;
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
//   if(i==0) return
//   -2*mu*(pow(x,4)*(6*y-3)+pow(x,3)*(6-12*y)+3*pow(x,2)*(4*pow(y,3)-6*pow(y,2)+4*y-1)-6*x*y*(2*pow(y,2)-3*y+1)+y*(2*pow(y,2)-3*y+1))
//   + 10*(3*pow(x-0.5,2)*pow(y,2)-3*pow(1-x,2)*pow(y-0.5,3)); else if(i==1) return
//   2*mu*(2*pow(x,3)*(6*pow(y,2)-6*y+1)-3*pow(x,2)*(6*pow(y,2)-6*y+1)+x*(6*pow(y,4)-12*pow(y,3)+12*pow(y,2)-6*y+1)-3*pow(y-1,2)*pow(y,2))
//   + 10*(2*pow(x-0.5,3)*y+3*pow(1-x,3)*pow(y-0.5,2)); else return 0;
// }
// R fun_exact(const R2 P, int i) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return pow(x,2)*pow(1-x,2)*y*(1-y)*(1-2*y);
//   else if(i==1) return (-x)*(1-x)*(1-2*x)*pow(y,2)*pow(1-y,2);
//   else return 10*(pow(x-0.5,3)*pow(y,2)+pow(1-x,3)*pow(y-0.5,3));
// }
} // namespace Erik_Data_StokesDiv
// using namespace Erik_Data_StokesDiv;
using namespace DynamicDropData;

// class LambdaGamma : public VirtualParameter { public :
// const CutFEMParameter& mu_;
//
// LambdaBoundary(const CutFEM_Parameter& mu, double G, double H) : mu_(mu) , G_(G) , H_(H) {} double evaluate(int
// domain, double h, double meas, double measK, double measCut) const { double gamma = meas / h ; double alpha = measK /
// h / h; double val = mu_.evaluate(domain,h,meas,measK,measCut) return val/h(G_+H_*gamma/alpha) ;}
// };

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 11;
    int ny = 11;
    // int d = 2;

    ProblemOption optionStokes;
    optionStokes.order_space_element_quadrature_ = 9;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    for (int i = 0; i < 3; ++i) { // i<3

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, -1., -1., 2., 2.);
        const R hi           = 2. / (nx - 1);
        const R penaltyParam = 1e1 / hi;
        // const R sigma = 1;
        double uPenParam     = 1e0;
        double pPenParam     = 1e0;
        double cc            = 1e0;

        // CutFEMParameter mu(1e-1,1e-3);
        // CutFEMParameter invmu(1e1,1e3);
        CutFEMParameter mu(mu1, mu2);
        // double sigma = 0;//24.5;//700;
        double kappa1 = mu(1) / (mu(0) + mu(1));
        double kappa2 = mu(0) / (mu(0) + mu(1));

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(2);
        Space VELh_(Kh, FEvelocity);
        Space SCAh_(Kh, DataFE<Mesh>::P2);

        Space Uh_(Kh, DataFE<Mesh>::P2);
        Space Vh_(Kh, DataFE<Mesh2>::RT1);
        Space Wh_(Kh, DataFE<Mesh2>::P1dc);

        ActiveMesh<Mesh> Khi(Kh, interface);

        CutSpace VELh(Khi, VELh_);
        CutSpace SCAh(Khi, SCAh_);

        CutSpace Uh(Khi, Uh_);
        CutSpace Vh(Khi, Vh_);
        CutSpace Wh(Khi, Wh_);

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h u0(VELh, fun_exact_u);
        Fun_h p0(SCAh, fun_exact_p);
        ExpressionFunFEM<Mesh> exactp(p0, 0, op_id);

        CutFEM<Mesh2> stokes(Uh);
        stokes.add(Vh);
        stokes.add(Wh);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest w(Uh, 1, 0), tau(Uh, 1, 0), u(Vh, 2, 0), v(Vh, 2, 0), p(Wh, 1, 0), q(Wh, 1, 0);

        // [Bulk]
        stokes.addBilinear( // w = curl u
            +innerProduct(1. / mu * w, tau) - innerProduct(u, rotgrad(tau)), Khi);
        stokes.addBilinear( // mu Delta u + grad p
            +innerProduct(rotgrad(w), v) - innerProduct(p, div(v)), Khi);
        stokes.addLinear(+innerProduct(fh.expression(2), v), Khi);
        stokes.addBilinear(+innerProduct(div(u), q), Khi);
        // [Interface]
        stokes.addBilinear(
            // + innerProduct(jump(u*t), average(tau,kappa1,kappa2)) // natural = 0 Dirichlet, so leave alone
            -innerProduct(average(u * t, kappa1, kappa2), jump(tau)) // natural

                // - innerProduct(average(w,kappa1,kappa2), jump(v*t)) // antisym 1
                + innerProduct(jump(w * t), cc * average(v, kappa1, kappa2)) // antisym 2 (baked into itf cond)

                // + innerProduct(average(u,kappa1,kappa2), 2*jump(mu*grad(v)*n)) // made up 1
                + innerProduct(2 * jump(mu * grad(u) * n),
                               cc * average(v, kappa1, kappa2)) // made up 2 (baked into itf cond)

                + innerProduct(average(p, kappa1, kappa2), jump(v * n)) // natural
                // + innerProduct(jump(p*n), cc*average(v,kappa1,kappa2)) // natural (baked into itf cond)

                + innerProduct(penaltyParam * (jump(u)), jump(v)),
            interface);
        stokes.addLinear(+innerProduct(3. / 2, cc * average(v * n, kappa2, kappa1)) * sigma, interface);
        // [Boundary natural]
        stokes.addLinear(-innerProduct(u0 * t, -1 * tau) // tangent points in wrong direction
                             - innerProduct(p0.expression(), v * n),
                         Khi, INTEGRAL_BOUNDARY);
        FunTest grad2un = grad(grad(u) * n) * n;
        // // FunTest grad2pn = grad(grad(p)*n)*n;
        // // FunTest grad2divun = grad(grad(div(u))*n)*n;
        stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
            +innerProduct(pPenParam * pow(hi, -1) * jump(rotgrad(w)), jump(v)) -
                innerProduct(pPenParam * pow(hi, -1) * jump(u), jump(rotgrad(tau))) +
                innerProduct(uPenParam * pow(hi, 1) * jump(grad(rotgrad(w)) * n), jump(grad(v) * n)) -
                innerProduct(uPenParam * pow(hi, 1) * jump(grad(u) * n), jump(grad(rotgrad(tau)) * n))

                - innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
                innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q))
            // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
            // -innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))) , jump(grad(q)))

            ,
            Khi
            // , macro
        );
        // [Boundary enforced]
        // stokes.addBilinear( // int_Omg grad(p)*v = int_itf p v*t - int_Omg p div(v)
        //   + innerProduct(p, v*n) // natural
        //   + innerProduct(penaltyParam*u, v) // stability
        //   , Khi, INTEGRAL_BOUNDARY
        // );
        // stokes.addLinear(
        //   - innerProduct(u0*t,tau) // [wtf why is + now correct..?]
        //   + innerProduct(u0.expression(2), penaltyParam*v)
        //   , Khi, INTEGRAL_BOUNDARY
        // );
        // R meanP = integral(Khi,exactp,0);
        // stokes.addLagrangeMultiplier(
        //   innerProduct(1.,p), meanP
        //   , Khi
        // );
        // CutFEM<Mesh2> lagr(Uh); lagr.add(Vh); lagr.add(Wh);
        // Rn zero_vec = lagr.rhs_;
        // lagr.addLinear(
        //   innerProduct(1.,p)
        //   , Khi
        // );
        // Rn lag_row(lagr.rhs_); lagr.rhs_ = zero_vec;
        // // lagr.addLinear(
        // //   -innerProduct(1, div(v))
        // //   , Kh_i, INTEGRAL_EXTENSION, 1
        // // );
        // // lagr.addLinear(
        // //   innerProduct(1, div(v))
        // //   , Kh_i
        // // );
        // lagr.addLinear(
        //   innerProduct(1, v*n)
        //   , interface
        // );
        // stokes.addLagrangeVecToRowAndCol(lag_row,lagr.rhs_,0);
        // FunTest grad2un = grad(grad(u)*n)*n;
        // // // FunTest grad2pn = grad(grad(p)*n)*n;
        // // // FunTest grad2divun = grad(grad(div(u))*n)*n;
        // stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
        //   // innerProduct(uPenParam*pow(hi,0)*jump(w), jump(tau)) // [w in P1, continuous]
        //   // +innerProduct(pPenParam*pow(hi,3)*jump(rotgrad(w)), jump(rotgrad(tau)))
        //   // +innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(v)) // [maybe should be 2k-1 if can scale pressure
        //   also]
        //   // +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(v)*n))
        //   // +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*t), jump(grad(v)*t))
        //   // // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
        //   // -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
        //   // +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
        //   // // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
        //   // // -innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))) , jump(grad(q)))

        //   +innerProduct(pPenParam*pow(hi,-1)*jump(rotgrad(w)), jump(v))
        //   -innerProduct(pPenParam*pow(hi,-1)*jump(u), jump(rotgrad(tau)))
        //   +innerProduct(uPenParam*pow(hi,1)*jump(grad(rotgrad(w))*n), jump(grad(v)*n))
        //   -innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(rotgrad(tau))*n))

        //   -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
        //   +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
        //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
        //   // -innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))) , jump(grad(q)))

        //   , Khi
        //   // , macro
        // );

        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_vort_dof = Uh.get_nb_dof();
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_wh     = stokes.rhs_(SubArray(nb_vort_dof, 0));
        Rn_ data_uh     = stokes.rhs_(SubArray(
            nb_flux_dof, nb_vort_dof)); // Rn_ data_uh = stokes.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));
        Rn_ data_ph     = stokes.rhs_(SubArray(
            Wh.get_nb_dof(),
            nb_vort_dof +
                nb_flux_dof)); // Rn_ data_ph = stokes.rhs_(SubArray(stokes_.get_nb_dof(),nb_vort_dof+nb_flux_dof));
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
            Fun_h solu(Vh, fun_exact_u);
            Fun_h soluErr(Vh, fun_exact_u);
            Fun_h solp(Wh, fun_exact_p);
            soluErr.v -= uh.v;
            soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Khi, "stokes_" + to_string(i) + ".vtk");
            writer.add(wh, "vorticity", 0, 1);
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(dx_uh0 + dy_uh1, "divergence");
            // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
            writer.add(solu, "velocityExact", 0, 2);
            writer.add(soluErr, "velocityError", 0, 2);
            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        // R errW      = L2normCut(wh,fun_exact_w,0,1);
        R errU      = L2normCut(uh, fun_exact_u, 0, 2);
        R errP      = L2normCut(ph, fun_exact_p, 0, 1);
        R errDiv    = L2normCut(dx_uh0 + dy_uh1, Khi);
        R maxErrDiv = maxNormCut(dx_uh0 + dy_uh1, Khi);
        // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

        ul2.push_back(errU);
        pl2.push_back(errP);
        divl2.push_back(errDiv);
        divmax.push_back(maxErrDiv);
        h.push_back(hi);
        if (i == 0) {
            convu.push_back(0);
            convp.push_back(0);
        } else {
            convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
        }

        nx = 2 * nx - 1;
        ny = 2 * ny - 1;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ') << "err p"
              << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ') << "err u"
              << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
                  << pl2[i] << std::setw(15) << std::setfill(' ') << convp[i] << std::setw(15) << std::setfill(' ')
                  << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
#endif

#ifdef PROBLEM_UNFITTED_STOKESRT

namespace Erik_Data_UNFITTED_STOKESRT {

R shift        = 0.5;
// R interfaceRad = 0.25+1e-10;//2./3; // not exactly 1/4 to avoid interface cutting exaclty a vertex
R interfaceRad = sqrt(0.25) - 1e-12; // [<-- Olshanskii example sqrt(0.25)=0.5 ]
R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x - shift) * (P.x - shift) + (P.y - shift) * (P.y - shift)) - interfaceRad;
}

// [Olshanskii example]
R fun_div(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 0;
}
R fun_rhs(const R2 P, int i, int dom) {
    // R mu=1;
    R x = P.x;
    R y = P.y;
    // if(i==0)      return  40*x*(x*x - y*y) - (16*y - 8);
    // else if(i==1) return (16*x - 8) - 40*y*(x*x - y*y);
    if (i == 0)
        return 40 * x * x * x - 40 * x * y * y - 32 * y + 16;
    else if (i == 1)
        return -40 * x * x * y + 32 * x + 40 * y * y * y - 16;
}
R fun_exact_u(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 2 * (x * x - x + 0.25 + y * y - y) * (2 * y - 1);
    else
        return -2 * (x * x - x + 0.25 + y * y - y) * (2 * x - 1);
}
R fun_exact_p(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 10 * (x * x - y * y) * (x * x - y * y); // - 670.171/(pi*interfaceRad*interfaceRad);
}

// [??]
// R fun_rhs(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return 0;
//   else return 0;
// }
// R fun_exact_u(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0)      return  20*x*pow(y,3);
//   else if(i==1) return 5*pow(x,4)-5*pow(y,4);
// }
// R fun_exact_p(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 20*(3*pow(x,2)*y-pow(y,3));//-2.21664/sqrt(M_PI*interfaceRad*interfaceRad);
// }
// R fun_div(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }

// [Ex: divu = 0, inhomog. BC]
// R fun_rhs(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//    if(i==0) return -2*y;
//    else return 2*x;
//   // if(i==0) return -4*(6*x*x - 6*x + 1)*y*(2*y*y - 3*y + 1) - 12*(x-1)*(x-1)*x*x*(2*y - 1);
//   // else if(i==1) return 4*x*(2*x*x - 3 *x + 1)*(6*y*y - 6*y + 1) + 12*(2*x - 1)*(y - 1)*(y-1)*y*y;
//   // else return 0;
// }
// R fun_exact_u(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   // if(i==0) return      2*x*y*(x-1)*(y-1)*x*(1-x)*(2*y-1);
//   // else if(i==1) return 2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x-1);
//   if(i==0)      return  x*x*y;
//   else if(i==1) return -x*y*y;
// }
// R fun_exact_p(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0.;
// }
// R fun_div(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }

// R fun_rhs(const R2 P, int i, int dom) {
//   R mu=1;
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return
//   -2*mu*(pow(x,4)*(6*y-3)+pow(x,3)*(6-12*y)+3*pow(x,2)*(4*pow(y,3)-6*pow(y,2)+4*y-1)-6*x*y*(2*pow(y,2)-3*y+1)+y*(2*pow(y,2)-3*y+1))
//   + 10*(3*pow(x-0.5,2)*pow(y,2)-3*pow(1-x,2)*pow(y-0.5,3)); else if(i==1) return
//   2*mu*(2*pow(x,3)*(6*pow(y,2)-6*y+1)-3*pow(x,2)*(6*pow(y,2)-6*y+1)+x*(6*pow(y,4)-12*pow(y,3)+12*pow(y,2)-6*y+1)-3*pow(y-1,2)*pow(y,2))
//   + 10*(2*pow(x-0.5,3)*y+3*pow(1-x,3)*pow(y-0.5,2)); else return 0;
// }
// R fun_exact_u(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return pow(x,2)*pow(1-x,2)*y*(1-y)*(1-2*y);
//   else return (-x)*(1-x)*(1-2*x)*pow(y,2)*pow(1-y,2);
// }
// R fun_exact_p(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 10*(pow(x-0.5,3)*pow(y,2)+pow(1-x,3)*pow(y-0.5,3));
// }
// R fun_div(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }
} // namespace Erik_Data_UNFITTED_STOKESRT
using namespace Erik_Data_UNFITTED_STOKESRT;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 11;
    int ny = 11;
    // int d = 2;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 4;
    for (int i = 0; i < iters; ++i) { // i<3

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, 0., 0., 1., 1.);
        const R hi       = 1. / (nx - 1);
        const R sigma    = 1e-2; // 1e-2; // 1e-2
        const R lamun    = 400;  // sigma
        const R lamut    = 4e3;  // 4e3, 8e2
        double uPenParam = 1e0;  // sigma, 6e2, 6e0
        double pPenParam = 1e0;  // sigma, 4e2, 4e0

        // const R vnPen = 4e3; // 4e3, 8e2

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(4);
        Space VELh_(Kh, FEvelocity);
        Space SCAh_(Kh, DataFE<Mesh>::P2);

        Space Wh(Kh, DataFE<Mesh2>::BDM1);
        Space Qh(Kh, DataFE<Mesh2>::P0); // FOR MIXEDSPACE

        ActiveMesh<Mesh> Khi(Kh);
        Khi.truncate(interface, 1);

        MacroElement<Mesh> macro(Khi, 1);

        CutSpace VELh(Khi, VELh_);
        CutSpace SCAh(Khi, SCAh_);

        CutSpace Vh(Khi, Wh);
        CutSpace Ph(Khi, Qh);

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h gh(VELh, fun_exact_u);
        Fun_h exactph(SCAh, fun_exact_p);
        ExpressionFunFEM<Mesh> exactp(exactph, 0, op_id);

        CutFEM<Mesh2> stokes(Vh);
        stokes.add(Ph);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);

        R mu = 1;
        {
            stokes.addBilinear(
                contractProduct(mu * grad(u), grad(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q), Khi);
            stokes.addLinear(innerProduct(fh.expression(2), v), Khi);
            // [NECESSARY FOR RT / BDM ELEMENTS]
            stokes.addBilinear(-innerProduct(mu * average(grad(u * t) * n, 0.5, 0.5), jump(v * t)) +
                                   innerProduct(jump(u * t), mu * average(grad(v * t) * n, 0.5, 0.5)) +
                                   innerProduct(1. / hi * sigma * jump(u * t), jump(v * t)),
                               Khi, INTEGRAL_INNER_EDGE_2D);
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v)      // natural
                                   + innerProduct(u, mu * grad(v) * n) // symmetry
                                   // + innerProduct(1./hi*lamut*u*t, v*t) // stability
                                   // + innerProduct(1./hi*lamun*u*n, v*n) // stability
                                   + innerProduct(1. / hi * lamut * u, v) // stability

                                   + innerProduct(p, v * n) // natural
                               // - innerProduct(u*n, q)
                               ,
                               interface);
            stokes.addLinear(+innerProduct(gh.expression(2), mu * grad(v) * n)
                                 // + innerProduct(gh*t, 1./hi*lamut*v*t)
                                 // + innerProduct(gh*n, 1./hi*lamun*v*n)
                                 + innerProduct(gh.expression(2), 1. / hi * lamut * v)

                             // - innerProduct(gh*n, q)
                             ,
                             interface);
            // [Sets uniqueness of the pressure]
            // R meanP = integral(Khi,exactp,0);
            // stokes.addLagrangeMultiplier(
            //   innerProduct(1, p), 0
            //   , Khi
            // );
            CutFEM<Mesh2> lagr(Vh);
            lagr.add(Ph);
            Rn zero_vec = lagr.rhs_;
            lagr.addLinear(innerProduct(1., p), Khi);
            Rn lag_row(lagr.rhs_);
            lagr.rhs_ = zero_vec;
            lagr.addLinear(innerProduct(1, v * n), interface);
            stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);

            FunTest grad2un = grad(grad(u) * n) * n;
            // // FunTest grad2pn = grad(grad(p)*n)*n;
            // // FunTest grad2divun = grad(grad(div(u))*n)*n;
            stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
                                         //  innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v)) // [Method 1: Remove
                                         //  jump in vel]
                                         // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
                                         // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
                                         // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
                                         // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

                +innerProduct(uPenParam * pow(hi, -1) * jump(u), jump(v)) // [Method 1: Remove jump in vel]
                    + innerProduct(uPenParam * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
                    innerProduct(uPenParam * pow(hi, 3) * jump(grad2un), jump(grad2un)) -
                    innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
                    innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(v))), jump(grad(q)))
                // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
                // +innerProduct(pPenParam*pow(hi,5)*jump(grad2divun) , jump(grad2pn))
                ,
                Khi, macro);
        }

        // [grad -> eps]
        {
            // stokes.addBilinear(
            //   contractProduct(2*mu*Eps(u),Eps(v))
            //   - innerProduct(p,div(v))
            //   + innerProduct(div(u),q)
            //   , Khi
            // );
            // stokes.addLinear(
            //   innerProduct(fh.expression(2),v)
            //   , Khi
            // );
            // stokes.addBilinear(
            //   - innerProduct(average(2*mu*(Eps(u)*t)*n,0.5,0.5), jump(v*t))
            //   + innerProduct(jump(u*t), average(2*mu*(Eps(v)*t)*n,0.5,0.5))
            //   + innerProduct(1./hi*sigma*(jump(u*t)), jump(v*t))
            //   , Khi
            //   , INTEGRAL_INNER_EDGE_2D
            // );
            // stokes.addBilinear(
            //   - innerProduct(2*mu*Eps(u)*n, v)  //natural
            //   + innerProduct(u, 2*mu*Eps(v)*n)
            //   + innerProduct(1./hi*penaltyParam*u, v)

            //   + innerProduct(p, v*n)  // natural
            //   , interface
            // );
            // stokes.addLinear(
            //   + innerProduct(gh.expression(2), 2*mu*Eps(v)*n)
            //   + innerProduct(gh.expression(2), 1./hi*penaltyParam*v)

            //   , interface
            // );

            // // [Sets uniqueness of the pressure]
            // // R meanP = integral(Khi,exactp,0);
            // // stokes.addLagrangeMultiplier(
            // //   innerProduct(1, p), 0
            // //   , Khi
            // // );
            // CutFEM<Mesh2> lagr(Vh); lagr.add(Ph);
            // Rn zero_vec = lagr.rhs_;
            // lagr.addLinear(
            //   innerProduct(1, p)
            //   , Khi
            // );
            // Rn lag_row(lagr.rhs_); lagr.rhs_ = zero_vec;
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , interface
            // );
            // stokes.addLagrangeVecToRowAndCol(lag_row,lagr.rhs_,0);

            // FunTest grad2un = grad(grad(u)*n)*n;
            // // FunTest grad2pn = grad(grad(p)*n)*n;
            // // FunTest grad2divun = grad(grad(div(u))*n)*n;
            // stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
            //   //  innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v))
            //   // +innerProduct(uPenParam*pow(hi,3)*jump(Eps(u)*n), jump(Eps(v)*n))
            //   // +innerProduct(uPenParam*pow(hi,5)*jump(Eps2un), jump(Eps2un))
            //   // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
            //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

            //   // +innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v))
            //   // +innerProduct(uPenParam*pow(hi,3)*jump(Eps(u)*n), jump(Eps(v)*n))
            //   // // +innerProduct(uPenParam*pow(hi,5)*jump(Eps2un), jump(Eps2un))
            //   // -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
            //   // +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
            //   // -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
            //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(v))) , jump(grad(q)))
            //   // // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
            //   // // +innerProduct(pPenParam*pow(hi,5)*jump(grad2divun) , jump(grad2pn))

            //   +innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v))
            //   +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
            //   +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
            //   -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
            //   +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
            //   -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
            //   +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(v))) , jump(grad(q)))
            //   // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
            //   // +innerProduct(pPenParam*pow(hi,5)*jump(grad2divun) , jump(grad2pn))
            //   , Khi
            //   , macro
            // );
        }

        // [divu extension alternative]
        // stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
        //   +innerProduct(uPenParam*hi*jump(u), jump(v)) // [Method 1: Remove jump in vel]
        //   +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
        //   +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
        //   -innerProduct(pPenParam*hi*jump(p), jump(div(v)))
        //   -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
        //   // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
        //   , Khi
        //   // , macro
        // );
        // stokes.addBilinearOtherSide(
        //   innerProduct(div(u), q)
        // , Khi, 1
        // );
        // FunTest p1(Ph,1,0,0);
        // stokes.addLagrangeMultiplierBothSides(
        //   innerProduct(1.,p1), 0.
        //   , Khi, 1
        // );

        // std::cout << integral(Khi,exactp,0) << std::endl;
        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_uh     = stokes.rhs_(SubArray(nb_flux_dof, 0));
        Rn_ data_ph     = stokes.rhs_(SubArray(Ph.get_nb_dof(), nb_flux_dof));
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Ph, data_ph);

        // [Post process pressure]
        R meanP = integral(Khi, exactp, 0);
        ExpressionFunFEM<Mesh> fem_p(ph, 0, op_id);
        R meanPfem = integral(Khi, fem_p, 0);
        // std::cout << meanP << std::endl;
        CutFEM<Mesh2> post(Qh);
        post.addLinear(innerProduct(1, q), Khi);
        R area = post.rhs_.sum();
        ph.v -= meanPfem / area;
        ph.v += meanP / area;

        // [Other post process fix w/o Lag mult]
        // R meanP = integral(Khi,exactp,0); std::cout << meanP << std::endl; std::cout
        // << 2.21664/sqrt(M_PI*interfaceRad*interfaceRad) << std::endl; Rn meanp_vec(data_ph); for (int m=0;
        // m<meanp_vec.size(); m++) {
        //   meanp_vec[m] = 2.21664/sqrt(M_PI*interfaceRad*interfaceRad);
        // }
        // ph.v += meanp_vec;

        // [divu post process "fix"]
        // R data_lagrange = stokes.rhs_[1+nb_flux_dof+Ph.get_nb_dof()]; // [stokes.rhs_.size() =
        // 1+idx0_s+Ph.get_nb_dof()] std::cout << data_lagrange << std::endl; Rn lagrange_vec(data_ph); for (int m=0;
        // m<lagrange_vec.size(); m++) {
        //   // lagrange_vec[m] = 0;
        //   lagrange_vec[m] = data_lagrange;
        // }
        // Fun_h lambdah(Ph, lagrange_vec);
        // ExpressionFunFEM<Mesh> fflambdah(lambdah, 0, op_id);

        // [Post processing u (instead of divu)]
        // Space SSh(Kh, DataFE<Mesh2>::RT0);
        // CutSpace Sh(Khi, SSh);

        // [L2 projection]
        // CutFEM<Mesh2> post_proj(Sh);
        // FunTest rho(Sh,2,0), gamma(Sh,2,0);
        // post_proj.addBilinear(
        //   innerProduct(rho, gamma)
        //   +innerProduct(div(rho), div(gamma))
        //   , Khi
        // );
        // post_proj.addLinear(
        //   innerProduct(uh.expression(2), gamma)
        //   -innerProduct(lambdah.expression(), div(gamma))
        //   , Khi
        // );
        // post_proj.solve();
        // Fun_h rhoh(Sh, post_proj.rhs_);

        // Rn rhoh_vec = rhoh.v;
        // std::cout << rhoh_vec.size() << std::endl;
        // Rn rhoh_rt1vec(nb_flux_dof); // uh.v.size()
        // for(int idx=0; idx<nb_flux_dof; idx++) rhoh_rt1vec[idx]=0;
        // // std::cout << rhoh_vec << std::endl;

        // typedef GFESpace<Mesh> FESpace;
        // typedef typename FESpace::FElement FElement;
        // for(int k=0;k<Khi.get_nb_element();++k) {
        //   const FElement& FK(Vh[k]);
        //   const FElement& FK_p(Sh[k]);

        //   for (int idx_loc=0; idx_loc<3; idx_loc++) {
        //     int idx = FK.loc2glb(idx_loc*2);
        //     int idx_p = FK_p.loc2glb(idx_loc);
        //     if (abs(rhoh_rt1vec[idx]) != 0) continue;

        //     rhoh_rt1vec[idx] = rhoh_vec[idx_p];
        //   }
        // }
        // uh.v -= rhoh_rt1vec;

        // ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
        // ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

        // [Hdiv interpolation]
        // Fun_h rhoh(Sh, data_uh); // [does not get correct size 81?? gets size of Vh]

        // CutFEM<Mesh2> post_interp(Sh);
        // FunTest rho(Sh,2,0), gamma(Sh,2,0);
        // post_interp.addBilinear(
        //   innerProduct(jump(rho*n), average(gamma*n))
        //   +innerProduct(average(rho*n), jump(gamma*n))
        //   , Khi, innerFacet
        // );
        // post_interp.addLinear(
        //   innerProduct(uh*n, jump(gamma*n))
        //   // +innerProduct(jump(uh*n), average(gamma*n))
        //   , Khi, innerFacet
        // );
        // post_interp.addBilinear(
        //   innerProduct(rho*n, gamma*n)
        //   , interface
        // );
        // post_interp.addLinear(
        //   innerProduct(uh*n, gamma*n)
        //   , interface
        // );
        // post_interp.solve();
        // Fun_h rhoh(Sh, post_interp.rhs_);

        // typedef GFESpace<Mesh> FESpace;
        // typedef typename FESpace::FElement FElement;
        // for(int k=0;k<Khi.get_nb_element();++k) {
        //   const FElement& FK(Vh[k]);

        //   // [BASICALLY WANT TO DO CHAIN ALGORITHM]
        //   for (int idx_loc=0; idx_loc<3; idx_loc++) {
        //     // int idx = FK.loc2glb(idx_loc*2);
        //     // int idx_p = FK_p.loc2glb(idx_loc);
        //     // if (abs(rt0_vec[FK.loc2glb(idx_loc*2)]-uh.v[FK.loc2glb(idx_loc*2)]) > 1e-15) {
        //     //   assert(idx_loc != 2);
        //     //   continue;
        //     // }
        //     // else {
        //     //   rt0_vec[FK.loc2glb(idx_loc*2)] = uh.v[FK.loc2glb(idx_loc*2)]+data_lagrange*FK.T.measure();
        //     //   break;
        //     // }

        //     rt0_vec[FK.loc2glb(idx_loc*2+1)] = 0;
        //     // rhoh_vec[FK.loc2glb(idx_loc*2+2)] = 0;
        //     // rhoh_vec[FK.loc2glb(idx_loc*2+3)] = 0;
        //   }
        // }

        // // uh.v = rt0_vec;
        // uh.v -= rt0_vec;

        ExpressionFunFEM<Mesh> uh_0dx(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> uh_1dy(uh, 1, op_dy);

        // Fun_h test(Vh, rt0_vec);
        // ExpressionFunFEM<Mesh> femSol_0dx(test, 0, op_dx);
        // ExpressionFunFEM<Mesh> femSol_1dy(test, 1, op_dy);

        // [Errors]
        {
            Fun_h soluErr(Vh, fun_exact_u);
            Fun_h soluh(Vh, fun_exact_u);
            soluErr.v -= uh.v;
            soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Khi, "stokes_" + to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(uh_0dx + uh_1dy, "divergence");
            writer.add(soluh, "velocityExact", 0, 2);
            writer.add(soluErr, "velocityError", 0, 2);
            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        R errU      = L2normCut(uh, fun_exact_u, 0, 2);
        R errP      = L2normCut(ph, fun_exact_p, 0, 1);
        R errDiv    = L2normCut(uh_0dx + uh_1dy, fun_div, Khi);
        R maxErrDiv = maxNormCut(uh_0dx + uh_1dy, fun_div, Khi);
        // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

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
        h.push_back(hi);
        if (i == 0) {
            convu.push_back(0);
            convp.push_back(0);
        } else {
            convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
        }

        nx = 2 * nx - 1;
        ny = 2 * ny - 1;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ') << "err_p"
              << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ') << "err u"
              << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
                  << pl2[i] << std::setw(15) << std::setfill(' ') << convp[i] << std::setw(15) << std::setfill(' ')
                  << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
#endif

#ifdef PROBLEM_UNFITTED_STOKESRT_SYMMETRIC_ALT

namespace Erik_Data_UNFITTED_STOKESRT {

R shift        = 0.5;
// R interfaceRad = 0.25+1e-10;//2./3; // not exactly 1/4 to avoid interface cutting exaclty a vertex
R interfaceRad = sqrt(0.25) - 1e-12; // [<-- Olshanskii example sqrt(0.25)=0.5 ]
R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x - shift) * (P.x - shift) + (P.y - shift) * (P.y - shift)) - interfaceRad;
}

// [Olshanskii example]
R fun_div(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 0;
}
R fun_rhs(const R2 P, int i, int dom) {
    // R mu=1;
    R x = P.x;
    R y = P.y;
    // if(i==0)      return  40*x*(x*x - y*y) - (16*y - 8);
    // else if(i==1) return (16*x - 8) - 40*y*(x*x - y*y);
    if (i == 0)
        return 40 * x * x * x - 40 * x * y * y - 32 * y + 16;
    else if (i == 1)
        return -40 * x * x * y + 32 * x + 40 * y * y * y - 16;
}
R fun_exact_u(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 2 * (x * x - x + 0.25 + y * y - y) * (2 * y - 1);
    else
        return -2 * (x * x - x + 0.25 + y * y - y) * (2 * x - 1);
}
R fun_exact_p(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 10 * (x * x - y * y) * (x * x - y * y); // - 670.171/(pi*interfaceRad*interfaceRad);
}

// [??]
// R fun_rhs(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return 0;
//   else return 0;
// }
// R fun_exact_u(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0)      return  20*x*pow(y,3);
//   else if(i==1) return 5*pow(x,4)-5*pow(y,4);
// }
// R fun_exact_p(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 20*(3*pow(x,2)*y-pow(y,3));//-2.21664/sqrt(M_PI*interfaceRad*interfaceRad);
// }
// R fun_div(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }

// [Ex: divu = 0, inhomog. BC]
// R fun_rhs(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//    if(i==0) return -2*y;
//    else return 2*x;
//   // if(i==0) return -4*(6*x*x - 6*x + 1)*y*(2*y*y - 3*y + 1) - 12*(x-1)*(x-1)*x*x*(2*y - 1);
//   // else if(i==1) return 4*x*(2*x*x - 3 *x + 1)*(6*y*y - 6*y + 1) + 12*(2*x - 1)*(y - 1)*(y-1)*y*y;
//   // else return 0;
// }
// R fun_exact_u(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   // if(i==0) return      2*x*y*(x-1)*(y-1)*x*(1-x)*(2*y-1);
//   // else if(i==1) return 2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x-1);
//   if(i==0)      return  x*x*y;
//   else if(i==1) return -x*y*y;
// }
// R fun_exact_p(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0.;
// }
// R fun_div(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }

// R fun_rhs(const R2 P, int i, int dom) {
//   R mu=1;
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return
//   -2*mu*(pow(x,4)*(6*y-3)+pow(x,3)*(6-12*y)+3*pow(x,2)*(4*pow(y,3)-6*pow(y,2)+4*y-1)-6*x*y*(2*pow(y,2)-3*y+1)+y*(2*pow(y,2)-3*y+1))
//   + 10*(3*pow(x-0.5,2)*pow(y,2)-3*pow(1-x,2)*pow(y-0.5,3)); else if(i==1) return
//   2*mu*(2*pow(x,3)*(6*pow(y,2)-6*y+1)-3*pow(x,2)*(6*pow(y,2)-6*y+1)+x*(6*pow(y,4)-12*pow(y,3)+12*pow(y,2)-6*y+1)-3*pow(y-1,2)*pow(y,2))
//   + 10*(2*pow(x-0.5,3)*y+3*pow(1-x,3)*pow(y-0.5,2)); else return 0;
// }
// R fun_exact_u(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return pow(x,2)*pow(1-x,2)*y*(1-y)*(1-2*y);
//   else return (-x)*(1-x)*(1-2*x)*pow(y,2)*pow(1-y,2);
// }
// R fun_exact_p(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 10*(pow(x-0.5,3)*pow(y,2)+pow(1-x,3)*pow(y-0.5,3));
// }
// R fun_div(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }
} // namespace Erik_Data_UNFITTED_STOKESRT
using namespace Erik_Data_UNFITTED_STOKESRT;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 11;
    int ny = 11;
    // int d = 2;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 4;
    for (int i = 0; i < iters; ++i) { // i<3

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, 0., 0., 1., 1.);
        const R hi          = 1. / (nx - 1);
        const R sigma       = 1e0;  // 1e0
        const R itfPenParam = 4e2;  // 4e2
        double uStabParam   = 1e0;  // 1e0
        double pStabParam   = 1e0;  // 1e0
        double itfStabParam = 1e-2; // 1e-2

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(4);
        Space VELh_(Kh, FEvelocity);
        Space SCAh_(Kh, DataFE<Mesh>::P2);

        Space Wh(Kh, DataFE<Mesh2>::RT1);
        Space Qh(Kh, DataFE<Mesh2>::P1dc); // FOR MIXEDSPACE

        // ACTIVE MESH
        ActiveMesh<Mesh> Khi(Kh);
        Khi.truncate(interface, 1);
        MacroElement<Mesh> macro(Khi, 1);

        CutSpace VELh(Khi, VELh_);
        CutSpace SCAh(Khi, SCAh_);

        CutSpace Vh(Khi, Wh);
        CutSpace Ph(Khi, Qh);

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h gh(VELh, fun_exact_u);
        Fun_h exactph(SCAh, fun_exact_p);
        ExpressionFunFEM<Mesh> exactp(exactph, 0, op_id);

        // SURFACE MESH
        ActiveMesh<Mesh> Kh_itf(Kh);
        Kh_itf.createSurfaceMesh(interface);
        CutSpace Ph_itf(Kh_itf, Qh);
        MacroElementSurface<Mesh> macro_itf(Kh_itf, interface, 0.6); // 0.3

        // PROBLEM SETUP
        CutFEM<Mesh2> stokes(Vh);
        stokes.add(Ph);
        stokes.add(Ph_itf);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);
        FunTest p_itf(Ph_itf, 1, 0), q_itf(Ph_itf, 1, 0);

        R mu = 1;
        {
            stokes.addBilinear(
                contractProduct(mu * grad(u), grad(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q), Khi);
            stokes.addLinear(innerProduct(fh.expression(2), v), Khi);
            // [NECESSARY FOR RT / BDM ELEMENTS]
            stokes.addBilinear(-innerProduct(mu * average(grad(u * t) * n, 0.5, 0.5), jump(v * t)) +
                                   innerProduct(jump(u * t), mu * average(grad(v * t) * n, 0.5, 0.5)) +
                                   innerProduct(1. / hi * sigma * jump(u * t), jump(v * t)),
                               Khi, INTEGRAL_INNER_EDGE_2D);
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v)                       // natural
                                   + innerProduct(u, mu * grad(v) * n)                  // symmetry
                                   + innerProduct(u * t, 1. / hi * itfPenParam * v * t) // stability

                               ,
                               interface);
            stokes.addLinear(+innerProduct(gh.expression(2), mu * grad(v) * n) +
                                 innerProduct(gh * t, 1. / hi * itfPenParam * v * t)

                                 ,
                             interface);
            stokes.addBilinear(+innerProduct(p_itf, v * n) + innerProduct(u * n, q_itf), interface);
            stokes.addLinear(+innerProduct(gh * n, q_itf), interface);
            // [Sets uniqueness of the pressure]
            // R meanP = integral(Khi,exactp,0);
            stokes.addLagrangeMultiplier(+innerProduct(1, p), 0, Khi);

            FunTest grad2un = grad(grad(u) * n) * n;
            // // FunTest grad2pn = grad(grad(p)*n)*n;
            // // FunTest grad2divun = grad(grad(div(u))*n)*n;
            stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]

                +innerProduct(uStabParam * pow(hi, -1) * jump(u), jump(v)) // [Method 1: Remove jump in vel]
                    + innerProduct(uStabParam * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
                    innerProduct(uStabParam * pow(hi, 3) * jump(grad2un), jump(grad2un)) -
                    innerProduct(pStabParam * pow(hi, 1) * jump(p), jump(div(v))) +
                    innerProduct(pStabParam * pow(hi, 1) * jump(div(u)), jump(q)) -
                    innerProduct(pStabParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                    innerProduct(pStabParam * pow(hi, 3) * jump(grad(div(u))), jump(grad(q)))
                // -innerProduct(pStabParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
                // +innerProduct(pStabParam*pow(hi,5)*jump(grad2divun) , jump(grad2pn))
                ,
                Khi, macro);
            stokes.addFaceStabilization(
                -innerProduct(itfStabParam * pow(hi, -1) * jump(p_itf), jump(q_itf))
                // -innerProduct(itfStabParam*pow(hi,1)*jump(grad(p_itf)*n), jump(grad(q_itf)*n))
                ,
                Kh_itf
                // , macro_itf // somehow fails at last iteration when not using macro (due to umfpack maybe?)
            );
            stokes.addBilinear(-innerProduct(itfStabParam * pow(hi, 1) * grad(p_itf) * n, grad(q_itf) * n), Kh_itf);
        }

        // [divu extension alternative]
        // stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
        //   +innerProduct(uPenParam*hi*jump(u), jump(v)) // [Method 1: Remove jump in vel]
        //   +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
        //   +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
        //   -innerProduct(pPenParam*hi*jump(p), jump(div(v)))
        //   -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
        //   // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
        //   , Khi
        //   // , macro
        // );
        // stokes.addBilinearOtherSide(
        //   innerProduct(div(u), q)
        // , Khi, 1
        // );
        // FunTest p1(Ph,1,0,0);
        // stokes.addLagrangeMultiplierBothSides(
        //   innerProduct(1.,p1), 0.
        //   , Khi, 1
        // );

        // std::cout << integral(Khi,exactp,0) << std::endl;
        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_uh     = stokes.rhs_(SubArray(nb_flux_dof, 0));
        Rn_ data_ph     = stokes.rhs_(SubArray(Ph.get_nb_dof(), nb_flux_dof));
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Ph, data_ph);

        // [Post process pressure]
        R meanP = integral(Khi, exactp, 0);
        ExpressionFunFEM<Mesh> fem_p(ph, 0, op_id);
        R meanPfem = integral(Khi, fem_p, 0);
        // std::cout << meanP << std::endl;
        CutFEM<Mesh2> post(Qh);
        post.addLinear(innerProduct(1, q), Khi);
        R area = post.rhs_.sum();
        ph.v -= meanPfem / area;
        ph.v += meanP / area;

        // [Other post process fix w/o Lag mult]
        // R meanP = integral(Khi,exactp,0); std::cout << meanP << std::endl; std::cout
        // << 2.21664/sqrt(M_PI*interfaceRad*interfaceRad) << std::endl; Rn meanp_vec(data_ph); for (int m=0;
        // m<meanp_vec.size(); m++) {
        //   meanp_vec[m] = 2.21664/sqrt(M_PI*interfaceRad*interfaceRad);
        // }
        // ph.v += meanp_vec;

        // [divu post process "fix"]
        // R data_lagrange = stokes.rhs_[1+nb_flux_dof+Ph.get_nb_dof()]; // [stokes.rhs_.size() =
        // 1+idx0_s+Ph.get_nb_dof()] std::cout << data_lagrange << std::endl; Rn lagrange_vec(data_ph); for (int m=0;
        // m<lagrange_vec.size(); m++) {
        //   // lagrange_vec[m] = 0;
        //   lagrange_vec[m] = data_lagrange;
        // }
        // Fun_h lambdah(Ph, lagrange_vec);
        // ExpressionFunFEM<Mesh> fflambdah(lambdah, 0, op_id);

        // [Post processing u (instead of divu)]
        // Space SSh(Kh, DataFE<Mesh2>::RT0);
        // CutSpace Sh(Khi, SSh);

        // [L2 projection]
        // CutFEM<Mesh2> post_proj(Sh);
        // FunTest rho(Sh,2,0), gamma(Sh,2,0);
        // post_proj.addBilinear(
        //   innerProduct(rho, gamma)
        //   +innerProduct(div(rho), div(gamma))
        //   , Khi
        // );
        // post_proj.addLinear(
        //   innerProduct(uh.expression(2), gamma)
        //   -innerProduct(lambdah.expression(), div(gamma))
        //   , Khi
        // );
        // post_proj.solve();
        // Fun_h rhoh(Sh, post_proj.rhs_);

        // Rn rhoh_vec = rhoh.v;
        // std::cout << rhoh_vec.size() << std::endl;
        // Rn rhoh_rt1vec(nb_flux_dof); // uh.v.size()
        // for(int idx=0; idx<nb_flux_dof; idx++) rhoh_rt1vec[idx]=0;
        // // std::cout << rhoh_vec << std::endl;

        // typedef GFESpace<Mesh> FESpace;
        // typedef typename FESpace::FElement FElement;
        // for(int k=0;k<Khi.get_nb_element();++k) {
        //   const FElement& FK(Vh[k]);
        //   const FElement& FK_p(Sh[k]);

        //   for (int idx_loc=0; idx_loc<3; idx_loc++) {
        //     int idx = FK.loc2glb(idx_loc*2);
        //     int idx_p = FK_p.loc2glb(idx_loc);
        //     if (abs(rhoh_rt1vec[idx]) != 0) continue;

        //     rhoh_rt1vec[idx] = rhoh_vec[idx_p];
        //   }
        // }
        // uh.v -= rhoh_rt1vec;

        // ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
        // ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

        // [Hdiv interpolation]
        // Fun_h rhoh(Sh, data_uh); // [does not get correct size 81?? gets size of Vh]

        // CutFEM<Mesh2> post_interp(Sh);
        // FunTest rho(Sh,2,0), gamma(Sh,2,0);
        // post_interp.addBilinear(
        //   innerProduct(jump(rho*n), average(gamma*n))
        //   +innerProduct(average(rho*n), jump(gamma*n))
        //   , Khi, innerFacet
        // );
        // post_interp.addLinear(
        //   innerProduct(uh*n, jump(gamma*n))
        //   // +innerProduct(jump(uh*n), average(gamma*n))
        //   , Khi, innerFacet
        // );
        // post_interp.addBilinear(
        //   innerProduct(rho*n, gamma*n)
        //   , interface
        // );
        // post_interp.addLinear(
        //   innerProduct(uh*n, gamma*n)
        //   , interface
        // );
        // post_interp.solve();
        // Fun_h rhoh(Sh, post_interp.rhs_);

        // typedef GFESpace<Mesh> FESpace;
        // typedef typename FESpace::FElement FElement;
        // for(int k=0;k<Khi.get_nb_element();++k) {
        //   const FElement& FK(Vh[k]);

        //   // [BASICALLY WANT TO DO CHAIN ALGORITHM]
        //   for (int idx_loc=0; idx_loc<3; idx_loc++) {
        //     // int idx = FK.loc2glb(idx_loc*2);
        //     // int idx_p = FK_p.loc2glb(idx_loc);
        //     // if (abs(rt0_vec[FK.loc2glb(idx_loc*2)]-uh.v[FK.loc2glb(idx_loc*2)]) > 1e-15) {
        //     //   assert(idx_loc != 2);
        //     //   continue;
        //     // }
        //     // else {
        //     //   rt0_vec[FK.loc2glb(idx_loc*2)] = uh.v[FK.loc2glb(idx_loc*2)]+data_lagrange*FK.T.measure();
        //     //   break;
        //     // }

        //     rt0_vec[FK.loc2glb(idx_loc*2+1)] = 0;
        //     // rhoh_vec[FK.loc2glb(idx_loc*2+2)] = 0;
        //     // rhoh_vec[FK.loc2glb(idx_loc*2+3)] = 0;
        //   }
        // }

        // // uh.v = rt0_vec;
        // uh.v -= rt0_vec;

        ExpressionFunFEM<Mesh> uh_0dx(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> uh_1dy(uh, 1, op_dy);

        // Fun_h test(Vh, rt0_vec);
        // ExpressionFunFEM<Mesh> femSol_0dx(test, 0, op_dx);
        // ExpressionFunFEM<Mesh> femSol_1dy(test, 1, op_dy);

        // [Errors]
        {
            Fun_h soluErr(Vh, fun_exact_u);
            Fun_h soluh(Vh, fun_exact_u);
            soluErr.v -= uh.v;
            soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Khi, "stokes_" + to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(uh_0dx + uh_1dy, "divergence");
            writer.add(soluh, "velocityExact", 0, 2);
            writer.add(soluErr, "velocityError", 0, 2);
            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        R errU      = L2normCut(uh, fun_exact_u, 0, 2);
        R errP      = L2normCut(ph, fun_exact_p, 0, 1);
        R errDiv    = L2normCut(uh_0dx + uh_1dy, fun_div, Khi);
        R maxErrDiv = maxNormCut(uh_0dx + uh_1dy, fun_div, Khi);
        // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

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
        h.push_back(hi);
        if (i == 0) {
            convu.push_back(0);
            convp.push_back(0);
        } else {
            convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
        }

        nx = 2 * nx - 1;
        ny = 2 * ny - 1;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ') << "err_p"
              << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ') << "err u"
              << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
                  << pl2[i] << std::setw(15) << std::setfill(' ') << convp[i] << std::setw(15) << std::setfill(' ')
                  << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
#endif

#ifdef PROBLEM_UNFITTED_STOKES_VORTICITY

namespace Erik_Data_UNFITTED_STOKES_VORTICITY {

R shift = 0.5;

// [Olshanskii]
R interfaceRad = sqrt(0.25) - 1e-12; // [<-- Olshanskii example+1e-12]
R fun_div(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 0;
}
R fun_rhs(const R2 P, int i, int dom) {
    // R mu=1;
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 40 * x * x * x - 40 * x * y * y - 32 * y + 16;
    else if (i == 1)
        return -40 * x * x * y + 32 * x + 40 * y * y * y - 16;
}
R fun_exact_u(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 2 * (x * x - x + 0.25 + y * y - y) * (2 * y - 1);
    else
        return -2 * (x * x - x + 0.25 + y * y - y) * (2 * x - 1);
}
R fun_exact_p(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 10 * (x * x - y * y) * (x * x - y * y); // - 670.171/(pi*interfaceRad*interfaceRad);
}

// [??]
// R interfaceRad = 0.25+1e-8;//2./3; // not exactly 1/4 to avoid interface cutting exaclty a vertex
// R fun_rhs(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return 0;
//   else return 0;
// }
// R fun_exact_u(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0)      return  20*x*pow(y,3);
//   else if(i==1) return 5*pow(x,4)-5*pow(y,4);
// }
// R fun_exact_p(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 20*(3*pow(x,2)*y-pow(y,3));//-2.21664/sqrt(M_PI*interfaceRad*interfaceRad);
// }
// R fun_div(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }

// [Ex: divu = 0, inhomog. BC]
// R fun_rhs(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//    if(i==0) return -2*y;
//    else return 2*x;
//   // if(i==0) return -4*(6*x*x - 6*x + 1)*y*(2*y*y - 3*y + 1) - 12*(x-1)*(x-1)*x*x*(2*y - 1);
//   // else if(i==1) return 4*x*(2*x*x - 3 *x + 1)*(6*y*y - 6*y + 1) + 12*(2*x - 1)*(y - 1)*(y-1)*y*y;
//   // else return 0;
// }
// R fun_exact_u(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   // if(i==0) return      2*x*y*(x-1)*(y-1)*x*(1-x)*(2*y-1);
//   // else if(i==1) return 2*x*y*(x-1)*(y-1)*y*(y-1)*(2*x-1);
//   if(i==0)      return  x*x*y;
//   else if(i==1) return -x*y*y;
// }
// R fun_exact_p(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0.;
// }
// R fun_div(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }

// R fun_rhs(const R2 P, int i, int dom) {
//   R mu=1;
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return
//   -2*mu*(pow(x,4)*(6*y-3)+pow(x,3)*(6-12*y)+3*pow(x,2)*(4*pow(y,3)-6*pow(y,2)+4*y-1)-6*x*y*(2*pow(y,2)-3*y+1)+y*(2*pow(y,2)-3*y+1))
//   + 10*(3*pow(x-0.5,2)*pow(y,2)-3*pow(1-x,2)*pow(y-0.5,3)); else if(i==1) return
//   2*mu*(2*pow(x,3)*(6*pow(y,2)-6*y+1)-3*pow(x,2)*(6*pow(y,2)-6*y+1)+x*(6*pow(y,4)-12*pow(y,3)+12*pow(y,2)-6*y+1)-3*pow(y-1,2)*pow(y,2))
//   + 10*(2*pow(x-0.5,3)*y+3*pow(1-x,3)*pow(y-0.5,2)); else return 0;
// }
// R fun_exact_u(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   if(i==0) return pow(x,2)*pow(1-x,2)*y*(1-y)*(1-2*y);
//   else return (-x)*(1-x)*(1-2*x)*pow(y,2)*pow(1-y,2);
// }
// R fun_exact_p(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 10*(pow(x-0.5,3)*pow(y,2)+pow(1-x,3)*pow(y-0.5,3));
// }
// R fun_div(const R2 P, int i, int dom) {
//   R x = P.x;
//   R y = P.y;
//   return 0;
// }

R fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x - shift) * (P.x - shift) + (P.y - shift) * (P.y - shift)) - interfaceRad;
}
} // namespace Erik_Data_UNFITTED_STOKES_VORTICITY
using namespace Erik_Data_UNFITTED_STOKES_VORTICITY;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 11;
    int ny = 11;
    // int d = 2;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 4;
    for (int i = 0; i < iters; ++i) {

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, 0.0, 0.0, 1.0, 1.0);
        const R hi           = 1. / (nx - 1); // 1./(nx-1)
        const R penaltyParam = 8e2;           // 4e3, 8e2

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(4);
        Space VELh_(Kh, FEvelocity);
        Space SCAh_(Kh, DataFE<Mesh>::P4);

        Space Uh_(Kh, DataFE<Mesh>::P1); // Nedelec order 0 type 1
        Space Vh_(Kh, DataFE<Mesh2>::RT0);
        Space Wh_(Kh, DataFE<Mesh2>::P0);

        // [Remove exterior]
        ActiveMesh<Mesh> Khi(Kh);
        Khi.truncate(interface, 1);

        MacroElement<Mesh> macro(Khi, 1); // we use 0.25 for vorticity BC2

        CutSpace VELh(Khi, VELh_);
        CutSpace SCAh(Khi, SCAh_);

        CutSpace Uh(Khi, Uh_);
        CutSpace Vh(Khi, Vh_);
        CutSpace Wh(Khi, Wh_);

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h u0(VELh, fun_exact_u);
        Fun_h p0(SCAh, fun_exact_p);
        ExpressionFunFEM<Mesh> exactp(p0, 0, op_id);

        CutFEM<Mesh2> stokes(Uh);
        stokes.add(Vh);
        stokes.add(Wh);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest w(Uh, 1, 0), tau(Uh, 1, 0), u(Vh, 2, 0), v(Vh, 2, 0), p(Wh, 1, 0), q(Wh, 1, 0);

        R mu = 1;
        {
            // [Bulk]
            stokes.addBilinear( // w = curl u
                innerProduct(1. / mu * w, tau) - innerProduct(u, rotgrad(tau)), Khi);
            stokes.addBilinear( // mu Delta u + grad p
                innerProduct(rotgrad(w), v) - innerProduct(p, div(v)), Khi);
            stokes.addLinear(+innerProduct(fh.expression(2), v), Khi);
            stokes.addBilinear(+innerProduct(div(u), q), Khi);
            // [Dirichlet Velocity BC]
            {
                const MeshParameter &itf_h(Parameter::measureIntegral);
                stokes.addBilinear( // int_Omg grad(p)*v = int_itf p v*t - int_Omg p div(v)
                                    // + innerProduct(p, v*n)
                    +innerProduct(1. / hi * penaltyParam * u * n, v * n)
                    // - innerProduct(u*t, tau)
                    // + innerProduct(w, v*t)
                    // + innerProduct(1./hi*penaltyParam*u, v)
                    ,
                    interface);
                stokes.addLinear(+innerProduct(u0 * t, tau) // [wtf why is + now correct..?]
                                     + innerProduct(u0 * n, 1. / hi * penaltyParam * v * n)
                                 // - innerProduct(u0*t,tau)
                                 // + innerProduct(u0.expression(2), 1./hi*penaltyParam*v)
                                 ,
                                 interface);

                // [Sets uniqueness of the pressure]
                // R meanP = integral(Khi,exactp,0);
                // stokes.addLagrangeMultiplier(
                //   innerProduct(1, p), 0
                //   , Khi
                // );
                // [Sets uniqueness of the pressure in another way such that divu = 0]
                CutFEM<Mesh2> lagr(Uh);
                lagr.add(Vh);
                lagr.add(Wh);
                Rn zero_vec = lagr.rhs_;
                lagr.addLinear(innerProduct(1, p), Khi);
                Rn lag_row(lagr.rhs_);
                lagr.rhs_ = zero_vec;
                lagr.addLinear(innerProduct(1, v * n), interface);
                stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);
                // [Stabilization]
                double wPenParam = 1e1; // 1e1
                double uPenParam = 1e1; // 1e-1 ~ 1/penParam (2e0 for (0,lamm,0))
                double pPenParam = 1e1; // 1e0 (2e0 for (0,lamm,0))
                FunTest grad2un  = grad(grad(u) * n) * n;
                FunTest grad2wn  = grad(grad(w) * n) * n;
                // // FunTest grad2pn = grad(grad(p)*n)*n;
                // // FunTest grad2divun = grad(grad(div(u))*n)*n;
                stokes.addFaceStabilization(
                    /* "Primal" stab: (lw,0,la) */
                    // innerProduct(uPenParam*pow(hi,1)*jump(w), jump(tau)) // [w in P1, continuous]
                    +innerProduct(wPenParam * pow(hi, 3) * jump(grad(w) * n), jump(grad(tau) * n)) +
                        innerProduct(uPenParam * pow(hi, 5) * jump(grad2wn), jump(grad2wn)) +
                        innerProduct(uPenParam * pow(hi, 1) * jump(u),
                                     jump(v)) // [maybe should be 2k-1 if can scale pressure also]
                        + innerProduct(uPenParam * pow(hi, 3) * jump(grad(u) * n), jump(grad(v) * n)) +
                        innerProduct(uPenParam * pow(hi, 5) * jump(grad2un), jump(grad2un)) -
                        innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
                        innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
                        innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                        innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(u))), jump(grad(q)))

                    /* Mixed stab: (0,lm,0) */
                    // innerProduct(uPenParam*pow(hi,1)*jump(w), jump(tau)) // [w in P1, continuous]
                    // +innerProduct(wPenParam*pow(hi,3)*jump(grad(w)*n), jump(grad(tau)*n))
                    // +innerProduct(wPenParam*pow(hi,5)*jump(grad2wn), jump(grad2wn))
                    // +innerProduct(uPenParam*pow(hi,1)*jump(rotgrad(w)), jump(v))
                    // -innerProduct(uPenParam*pow(hi,1)*jump(u), jump(rotgrad(tau)))
                    // +innerProduct(uPenParam*pow(hi,3)*jump(grad(rotgrad(w))), jump(grad(v)))
                    // -innerProduct(uPenParam*pow(hi,3)*jump(grad(u)), jump(grad(rotgrad(tau))))
                    // -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
                    // +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
                    // -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
                    // +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))), jump(grad(q)))

                    //   /* 21/12/22 The last two terms below DESTROY the divergence despite not entering the q-row!
                    //   WHY?! -innerProduct(pPenParam*pow(hi,-1)*jump(u*t), jump(tau))
                    //   +innerProduct(pPenParam*pow(hi,-1)*jump(w), jump(v*t))
                    //   -innerProduct(pPenParam*pow(hi,1)*jump(grad(u*t)), jump(rotgrad(tau)))
                    //   +innerProduct(pPenParam*pow(hi,1)*jump(rotgrad(w)), jump(grad(v*t)))
                    //   */

                    ,
                    Khi, macro);
            }
            // [Vorticity BC 2 (natural)] (n x u = t . u = u0, p = p0)
            {
                // stokes.addLinear( // [normal points in "wrong" direction, t = -t & n = -n]
                //   + innerProduct(u0*t, tau)
                //   + innerProduct(p0.expression(), v*n)
                //   , interface
                // );
                // double wPenParam = 1e0;
                // double uPenParam = 1e0;
                // double pPenParam = 1e0;
                // FunTest grad2un = grad(grad(u)*n)*n;
                // // // FunTest grad2pn = grad(grad(p)*n)*n;
                // // // FunTest grad2divun = grad(grad(div(u))*n)*n;
                // stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
                //   // [More primal stab]
                //   // // // innerProduct(uPenParam*pow(hi,0)*jump(w), jump(tau)) // [w in P1, continuous]
                //   // +innerProduct(pPenParam*pow(hi,3)*jump(rotgrad(w)), jump(rotgrad(tau)))
                //   // +innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v))
                //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
                //   // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
                //   // -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
                //   // +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
                //   // -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
                //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))) , jump(grad(q)))

                //   // [Mixed stab (works great here!)]
                //   +innerProduct(pPenParam*pow(hi,1)*jump(rotgrad(w)), jump(v))
                //   -innerProduct(uPenParam*pow(hi,1)*jump(u), jump(rotgrad(tau)))
                //   +innerProduct(uPenParam*pow(hi,3)*jump(grad(rotgrad(w))*n), jump(grad(v)*n))
                //   -innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(rotgrad(tau))*n))
                //   -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
                //   +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
                //   -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
                //   +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))), jump(grad(q)))

                //   , Khi, macro
                // );
            }
        }

        // std::cout << integral(Khi,exactp,0) << std::endl;
        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_vort_dof = Uh.get_nb_dof();
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_wh     = stokes.rhs_(SubArray(nb_vort_dof, 0));
        Rn_ data_uh     = stokes.rhs_(SubArray(
            nb_flux_dof, nb_vort_dof)); // Rn_ data_uh = stokes.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));
        Rn_ data_ph     = stokes.rhs_(SubArray(
            Wh.get_nb_dof(),
            nb_vort_dof +
                nb_flux_dof)); // Rn_ data_ph = stokes.rhs_(SubArray(stokes_.get_nb_dof(),nb_vort_dof+nb_flux_dof));
        Fun_h wh(Uh, data_wh);
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Wh, data_ph);

        // [Post process pressure]
        R meanP = integral(Khi, exactp, 0);
        ExpressionFunFEM<Mesh> fem_p(ph, 0, op_id);
        R meanPfem = integral(Khi, fem_p, 0);
        // std::cout << meanP << std::endl;
        CutFEM<Mesh2> post(Wh);
        post.addLinear(innerProduct(1, q), Khi);
        R area = post.rhs_.sum();
        ph.v -= meanPfem / area;
        ph.v += meanP / area;

        ExpressionFunFEM<Mesh> dx_uh0(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> dy_uh1(uh, 1, op_dy);

        // [Errors]
        {
            // Fun_h solw(Uh, fun_exact_w);
            Fun_h solu(Vh, fun_exact_u);
            Fun_h soluErr(Vh, fun_exact_u);
            Fun_h solp(Wh, fun_exact_p);
            soluErr.v -= uh.v;
            soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Khi, "stokes_" + to_string(i) + ".vtk");
            writer.add(wh, "vorticity", 0, 1);
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(dx_uh0 + dy_uh1, "divergence");
            // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
            writer.add(solp, "pressureExact", 0, 1);
            writer.add(solu, "velocityExact", 0, 2);
            writer.add(soluErr, "velocityError", 0, 2);
            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        // R errW      = L2normCut(wh,fun_exact_w,0,1);
        R errU      = L2normCut(uh, fun_exact_u, 0, 2);
        R errP      = L2normCut(ph, fun_exact_p, 0, 1);
        R errDiv    = L2normCut(dx_uh0 + dy_uh1, fun_div, Khi);
        R maxErrDiv = maxNormCut(dx_uh0 + dy_uh1, fun_div, Khi);
        // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

        ul2.push_back(errU);
        pl2.push_back(errP);
        divl2.push_back(errDiv);
        divmax.push_back(maxErrDiv);
        h.push_back(hi);
        if (i == 0) {
            convu.push_back(0);
            convp.push_back(0);
        } else {
            convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
        }

        // if(i==4) {
        //   nx = nx+5;
        //   ny = ny+5;
        // } else {
        nx = 2 * nx - 1;
        ny = 2 * ny - 1;
        // }
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ') << "err p"
              << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ') << "err u"
              << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
                  << pl2[i] << std::setw(15) << std::setfill(' ') << convp[i] << std::setw(15) << std::setfill(' ')
                  << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
#endif

#ifdef PROBLEM_FITTED_CORIOLIS_STOKESRT

namespace Erik_Data_CORIOLIS_STOKESRT {
R sdBox(const R2 p, const R2 b) {
    R2 absp(abs(p.x), abs(p.y));
    R2 d = absp - b;
    R2 maxd0(max(d.x, 0.0), max(d.y, 0.0));
    return maxd0.norm() + min(max(d.x, d.y), 0.0);
    // return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}
R fun_levelSet(const R2 P, const int i) {
    // const R2 b(2.,1.); // b.x=width, b.y=height
    // const R2 shift(2.,0.);
    const R2 b(2.0, 1.0); // b.x=width, b.y=height
    const R2 shift(4., 0.);
    return sdBox(P - shift, b) - 0.01;
}

// [Coriolis example]
R fun_y(const R2 P, int i, int dom) { return P.y; }
R fun_div(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 0;
}
R fun_rhs(const R2 P, int i, int dom) {
    // R mu=1;
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 0; // 100*(2-x)*(2-x);
    else if (i == 1)
        return 0;
}
R fun_uin(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return y * (2 - y) / 2; // original paper 4*(2-y)*(y-1)
    else
        return 0;
}
R fun_uout(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 1. / 4 * (2 - y) * (y - 1); // original paper 4*(2-y)*(y-1)
    // if(i==0)    return 1./4*(2-y)*y; // [for a rectangle mesh]
    else
        return 0;
}
R fun_0(const R2 P, int i) { return 0; }
} // namespace Erik_Data_CORIOLIS_STOKESRT

using namespace Erik_Data_CORIOLIS_STOKESRT;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 11;
    int ny = 11;
    // int d = 2;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 1;
    for (int i = 0; i < iters; ++i) { // i<3

        // std::cout << "\n ------------------------------------- " << std::endl;
        // Mesh Kh(nx, ny, 0., 0., 4., 2.);
        // Mesh Kh("../Lmesh"+to_string(i)+".msh"); // [doesnt work.....]
        // Mesh Kh("../Rectmesh.msh");
        Mesh Kh("../Lmesh3.msh");

        // const R hi = 1./(nx-1);
        R hi = 0;
        for (int k = 0; k < Kh.nt; k++) {
            const auto &K(Kh[k]);
            hi = std::max(hi, K.hMax());
        }
        std::cout << hi << std::endl;

        Lagrange2 FEvelocity(4);
        Space VELh(Kh, FEvelocity);
        // Space VELh(Kh, DataFE<Mesh2>::BDM1);
        Space SCAh(Kh, DataFE<Mesh>::P2);

        Space Vh(Kh, DataFE<Mesh2>::BDM1);
        // Space Wh(Kh, FEvelocity);
        Space Ph(Kh, DataFE<Mesh2>::P0); // FOR MIXEDSPACE

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h uhin(VELh, fun_uin);
        Fun_h uhout(VELh, fun_uout);
        Fun_h zerofunh(VELh, fun_0);
        // Fun_h yfunh(VELh, fun_y);
        // Fun_h exactph(SCAh, fun_exact_p); ExpressionFunFEM<Mesh> exactp(exactph,0,op_id);

        FEM<Mesh2> stokes(Vh);
        stokes.add(Ph);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);
        FunTest u1(Vh, 1, 0), u2(Vh, 1, 1), v1(Vh, 1, 0), v2(Vh, 1, 1);

        R mu    = 0.01;
        R omega = 100; // 50 // original paper 100;

        {
            stokes.addBilinear(contractProduct(mu * grad(u), grad(v)) - innerProduct(p, div(v)) +
                                   innerProduct(div(u), q)

                                   - innerProduct(2 * omega * u2, v1) + innerProduct(2 * omega * u1, v2),
                               Kh);
            stokes.addLinear(innerProduct(fh.expression(2), v), Kh);
            const R sigma      = 4e0;   // 8e-1
            const R pp_tangent = sigma; // 1e2, 4e2
            const R pp_normal  = 8e-1;  // 1e2, 4e3
            // [NECESSARY FOR RT / BDM ELEMENTS]
            stokes.addBilinear(-innerProduct(mu * average(grad(u * t) * n, 0.5, 0.5), jump(v * t)) +
                                   innerProduct(jump(u * t), mu * average(grad(v * t) * n, 0.5, 0.5)) +
                                   innerProduct(1. / hi * (sigma * jump(u * t)), jump(v * t)),
                               Kh, INTEGRAL_INNER_EDGE_2D);
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v) // u=0
                                   + innerProduct(u, mu * grad(v) * n) +
                                   innerProduct(1. / hi * pp_tangent * u * t, v * t) +
                                   innerProduct(1. / hi * pp_normal * u * n, v * n) + innerProduct(p, v * n),
                               Kh, INTEGRAL_BOUNDARY, {1, 2, 3, 5} // , {1,3}
            );
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v) // natural
                                   + innerProduct(u, mu * grad(v) * n) +
                                   innerProduct(1. / hi * pp_tangent * u * t, v * t) +
                                   innerProduct(1. / hi * pp_normal * u * n, v * n) + innerProduct(p, v * n),
                               Kh, INTEGRAL_BOUNDARY, {4, 6} // , {2,4}
            );
            stokes.addLinear(+innerProduct(uhin.expression(2), mu * grad(v) * n) +
                                 innerProduct(uhin * t, 1. / hi * pp_tangent * v * t) +
                                 innerProduct(uhin * n, 1. / hi * pp_normal * v * n),
                             Kh, INTEGRAL_BOUNDARY, {6} // , {4}
            );
            stokes.addLinear(+innerProduct(uhout.expression(2), mu * grad(v) * n) +
                                 innerProduct(uhout * t, 1. / hi * pp_tangent * v * t) +
                                 innerProduct(uhout * n, 1. / hi * pp_normal * v * n),
                             Kh, INTEGRAL_BOUNDARY, {4} // , {2}
            );
            // stokes.setDirichlet(uhin,Kh,{6});
            // stokes.setDirichlet(uhout,Kh,{4});
            // stokes.setDirichlet(zerofunh,Kh,{1,2,3,5});
            // [Sets uniqueness of the pressure]
            // R meanP = integral(Kh,p0,0);
            // stokes.addLagrangeMultiplier(
            //   innerProduct(1.,p), 0
            //   , Kh
            // );
            FEM<Mesh2> lagr(Vh);
            lagr.add(Ph);
            Rn zero_vec = lagr.rhs_;
            lagr.addLinear(innerProduct(1, p), Kh);
            Rn lag_row(lagr.rhs_);
            lagr.rhs_ = zero_vec;
            lagr.addLinear(innerProduct(1, v * n), Kh, INTEGRAL_BOUNDARY);
            stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);
        }

        // [grad -> eps]
        {
            // stokes.addBilinear(
            //   contractProduct(2*mu*Eps(u),Eps(v))
            //   - innerProduct(p,div(v))
            //   + innerProduct(div(u),q)

            //   - innerProduct(2*omega*u2,v1)
            //   + innerProduct(2*omega*u1,v2)
            //   , Khi
            // );
            // stokes.addLinear(
            //   innerProduct(fh.expression(2),v)
            //   , Khi
            // );
            // stokes.addBilinear(
            //   - innerProduct(average(2*mu*(Eps(u)*t)*n,0.5,0.5), jump(v*t))
            //   + innerProduct(jump(u*t), average(2*mu*(Eps(v)*t)*n,0.5,0.5))
            //   + innerProduct(1./hi*sigma*(jump(u*t)), jump(v*t))
            //   , Khi
            //   , INTEGRAL_INNER_EDGE_2D
            // );
            // stokes.addBilinear(
            //   // - innerProduct(2*mu*Eps(u)*n, v)  //natural
            //   // + innerProduct(u, 2*mu*Eps(v)*n)
            //   + innerProduct(1./hi*penaltyParam*u, v)
            //   , Khi, INTEGRAL_BOUNDARY, {1,3}
            // );
            // stokes.addBilinear(
            //   - innerProduct(2*mu*Eps(u)*n, v)  //natural
            //   + innerProduct(u, 2*mu*Eps(v)*n)
            //   + innerProduct(1./hi*penaltyParam*u, v)
            //   , Khi, INTEGRAL_BOUNDARY, {2,4}
            // );
            // stokes.addLinear(
            //   + innerProduct(uhout.expression(2), 2*mu*Eps(v)*n)
            //   + innerProduct(uhout.expression(2), 1./hi*penaltyParam*v)

            //   , Khi, INTEGRAL_BOUNDARY, {2}
            // );
            // stokes.addLinear(
            //   + innerProduct(uhin.expression(2), 2*mu*Eps(v)*n)
            //   + innerProduct(uhin.expression(2), 1./hi*penaltyParam*v)

            //   , Khi, INTEGRAL_BOUNDARY, {4}
            // );

            // // [Sets uniqueness of the pressure]
            // // R meanP = integral(Khi,exactp,0);
            // // stokes.addLagrangeMultiplier(
            // //   innerProduct(1.,p), 0
            // //   , Khi
            // // );
            // CutFEM<Mesh2> lagr(Vh); lagr.add(Ph);
            // Rn zero_vec = lagr.rhs_;
            // lagr.addLinear(
            //   innerProduct(1.,p)
            //   , Khi
            // );
            // Rn lag_row(lagr.rhs_); lagr.rhs_ = zero_vec;
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , interface
            // );
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , Khi, INTEGRAL_BOUNDARY
            // );
            // stokes.addLagrangeVecToRowAndCol(lag_row,lagr.rhs_,0);

            // // FunTest grad2un = grad(grad(u)*n)*n;
            // // FunTest grad2pn = grad(grad(p)*n)*n;
            // // // FunTest grad2divun = grad(grad(div(u))*n)*n;
            // // stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
            // //   // +innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(v))
            // //   // +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(v)*n))
            // //   // // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
            // //   // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
            // //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

            // //   +innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(v))
            // //   +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(v)*n))
            // //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad2un), jump(grad2un))
            // //   -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
            // //   +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
            // //   -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
            // //   +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(v))) , jump(grad(q)))
            // //   // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
            // //   // +innerProduct(pPenParam*pow(hi,5)*jump(grad2divun) , jump(grad2pn))
            // //   , Khi
            // //   , macro
            // // );
        }

        // std::cout << integral(Khi,exactp,0) << std::endl;
        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_uh     = stokes.rhs_(SubArray(nb_flux_dof, 0));
        Rn_ data_ph     = stokes.rhs_(SubArray(Ph.get_nb_dof(), nb_flux_dof));
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Ph, data_ph);

        // [Post process pressure]
        // R meanP = integral(Khi,p0,0);
        // ExpressionFunFEM<Mesh> fem_p(ph,0,op_id);
        // R meanPfem = integral(Khi,fem_p,0);
        // // std::cout << meanP << std::endl;
        // CutFEM<Mesh2> post(Qh);
        // post.addLinear(
        //   innerProduct(1,q)
        //   , Khi
        // );
        // R area = post.rhs_.sum();
        // ph.v -= meanPfem/area;
        // ph.v += meanP/area;

        ExpressionFunFEM<Mesh> uh_0dx(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> uh_1dy(uh, 1, op_dy);

        // Fun_h test(Vh, rt0_vec);
        // ExpressionFunFEM<Mesh> femSol_0dx(test, 0, op_dx);
        // ExpressionFunFEM<Mesh> femSol_1dy(test, 1, op_dy);

        // [Errors]
        {
            // Fun_h soluErr(Vh, fun_exact_u);
            // Fun_h soluh(Vh, fun_exact_u);
            // soluErr.v -= uh.v;
            // soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Kh, "stokes_" + to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(uh_0dx + uh_1dy, "divergence");
            // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
            // writer.add(soluh, "velocityExact" , 0, 2);
            // writer.add(soluErr, "velocityError" , 0, 2);
            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        // R errU      = L2normCut(uh,fun_exact_u,0,2);
        // R errP      = L2normCut(ph,fun_exact_p,0,1);
        R errDiv    = L2norm(uh_0dx + uh_1dy, Kh);
        R maxErrDiv = maxNorm(uh_0dx + uh_1dy, Kh);
        // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

        // // solExVec -= stokesDiv.rhs;
        // for(int i=0;i<solExVec.size();++i){
        //   solExVec(i) = fabs(solExVec(i)-stokesDiv.rhs(i));
        // }
        //
        // Fun_h solEx(mixedSpace, solExVec);
        // writerS.add(solEx,"uh_err", 0, 2);

        // writerS.add(solEx,"uh_ex", 0, 2);
        // writerS.add(solEx,"ph_ex", 2, 1);

        // ul2.push_back(errU);
        // pl2.push_back(errP);
        divl2.push_back(errDiv);
        divmax.push_back(maxErrDiv);
        h.push_back(hi);
        // if(i==0) {convu.push_back(0); convp.push_back(0);}
        // else {
        //   convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
        //   convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
        // }

        nx = 2 * nx - 1;
        ny = 2 * ny - 1;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ')
              << "h"
              // << std::setw(15) << std::setfill(' ') << "err_p"
              // << std::setw(15) << std::setfill(' ') << "conv p"
              // << std::setw(15) << std::setfill(' ') << "err u"
              // << std::setw(15) << std::setfill(' ') << "conv u"
              << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ')
                  << h[i]
                  // << std::setw(15) << std::setfill(' ') << pl2[i]
                  // << std::setw(15) << std::setfill(' ') << convp[i]
                  // << std::setw(15) << std::setfill(' ') << ul2[i]
                  // << std::setw(15) << std::setfill(' ') << convu[i]
                  << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
#endif

#ifdef PROBLEM_UNFITTED_CORIOLIS_STOKESRT

namespace Erik_Data_CORIOLIS_STOKESRT {
R sdBox(const R2 p, const R2 b) {
    R2 absp(abs(p.x), abs(p.y));
    R2 d = absp - b;
    R2 maxd0(max(d.x, 0.0), max(d.y, 0.0));
    return maxd0.norm() + min(max(d.x, d.y), 0.0);
    // return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}
R fun_levelSet(const R2 P, const int i) {
    // const R2 b(2.,1.); // b.x=width, b.y=height
    // const R2 shift(2.,0.);
    const R2 b(2.0, 1.0); // b.x=width, b.y=height
    const R2 shift(4., 0.);
    return sdBox(P - shift, b) - 0.01;
}

// [Coriolis example]
R fun_y(const R2 P, int i, int dom) { return P.y; }
R fun_div(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 0;
}
R fun_rhs(const R2 P, int i, int dom) {
    // R mu=1;
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 0; // 100*(2-x)*(2-x);
    else if (i == 1)
        return 0;
}
R fun_uin(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return y * (2 - y) / 2; // original paper 4*(2-y)*(y-1)
    else
        return 0;
}
R fun_uout(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 1. / 4 * (2 - y) * (y - 1); // original paper 4*(2-y)*(y-1)
    else
        return 0;
}
// R fun_p(const R2 P, int i, int dom ) {
//   R x = P.x;
//   R y = P.y;
//   return 10*(x*x-y*y)*(x*x-y*y);// - 670.171/(pi*interfaceRad*interfaceRad);
// }
} // namespace Erik_Data_CORIOLIS_STOKESRT

using namespace Erik_Data_CORIOLIS_STOKESRT;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 11;
    int ny = 11;
    // int d = 2;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 4;
    for (int i = 0; i < iters; ++i) { // i<3

        // std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, 0., 0., 4., 2.);
        // Mesh Kh("../Lmesh"+to_string(i)+".msh");
        // if(i==1) Mesh Kh("../Lmesh1.msh");
        // else if(i==2) Mesh Kh("../Lmesh2.msh");

        const R hi = 1. / (nx - 1);

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(4);
        Space VELh_(Kh, FEvelocity);
        Space SCAh_(Kh, DataFE<Mesh>::P2);

        Space Wh(Kh, DataFE<Mesh2>::BDM1);
        // Space Wh(Kh, FEvelocity);
        Space Qh(Kh, DataFE<Mesh2>::P0); // FOR MIXEDSPACE

        ActiveMesh<Mesh> Khi(Kh);
        Khi.truncate(interface, -1);

        MacroElement<Mesh> macro(Khi, 1.);

        CutSpace VELh(Khi, VELh_);
        CutSpace SCAh(Khi, SCAh_);

        CutSpace Vh(Khi, Wh);
        CutSpace Ph(Khi, Qh);

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h uhin(VELh, fun_uin);
        Fun_h uhout(VELh, fun_uout);
        Fun_h yfunh(VELh, fun_y);
        // Fun_h exactph(SCAh, fun_exact_p); ExpressionFunFEM<Mesh> exactp(exactph,0,op_id);

        CutFEM<Mesh2> stokes(Vh);
        stokes.add(Ph);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);
        FunTest u1(Vh, 1, 0), u2(Vh, 1, 1), v1(Vh, 1, 0), v2(Vh, 1, 1);

        R mu    = 0.01;
        R omega = 100; // 50 // original paper 100;

        // [grad] DOESNT WORK AT ALL!!!?
        {
            stokes.addBilinear(contractProduct(mu * grad(u), grad(v)) - innerProduct(p, div(v)) +
                                   innerProduct(div(u), q)

                                   - innerProduct(2 * omega * u2, v1) + innerProduct(2 * omega * u1, v2),
                               Khi);
            stokes.addLinear(innerProduct(fh.expression(2), v), Khi);
            const R sigma      = 1e0; // 4
            const R pp_tangent = 1e0; // 4=sigma, 400, 400
            const R pp_normal  = 4e2; // 0.8, 4, 40
            double uPenParam   = 1e0; // 4, 1e-1, 1e1
            double pPenParam   = 1e0; // 4, 1e-1, 1e1
            // [NECESSARY FOR RT / BDM ELEMENTS]
            stokes.addBilinear(-innerProduct(average(grad(u * t) * n, 0.5, 0.5), jump(v * t)) +
                                   innerProduct(jump(u * t), average(grad(v * t) * n, 0.5, 0.5)) +
                                   innerProduct(1. / hi * (sigma * jump(u * t)), jump(v * t)),
                               Khi, INTEGRAL_INNER_EDGE_2D);
            // const R pp_13 = 1e2; // 1e2, 4e2
            // const R pp_24 = 1e2; // 1e2, 4e3
            // const R pp_itf = 1e2; // 1e2, 1e2
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v) // natural
                                   + innerProduct(u, mu * grad(v) * n) +
                                   innerProduct(1. / hi * pp_tangent * u * t, v * t) +
                                   innerProduct(1. / hi * pp_normal * u * n, v * n) + innerProduct(p, v * n),
                               interface);
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v) // natural
                                   + innerProduct(u, mu * grad(v) * n) +
                                   innerProduct(1. / hi * pp_tangent * u * t, v * t) +
                                   innerProduct(1. / hi * pp_normal * u * n, v * n) + innerProduct(p, v * n),
                               Khi, INTEGRAL_BOUNDARY, {2, 4});
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v) // natural
                                   + innerProduct(u, mu * grad(v) * n) +
                                   innerProduct(1. / hi * pp_tangent * u * t, v * t) +
                                   innerProduct(1. / hi * pp_normal * u * n, v * n) + innerProduct(p, v * n),
                               Khi, INTEGRAL_BOUNDARY, {1, 3});
            stokes.addLinear(+innerProduct(uhin.expression(2), mu * grad(v) * n) +
                                 innerProduct(uhin * t, 1. / hi * pp_tangent * v * t) +
                                 innerProduct(uhin * n, 1. / hi * pp_normal * v * n),
                             Khi, INTEGRAL_BOUNDARY, {4});
            stokes.addLinear(+innerProduct(uhout.expression(2), mu * grad(v) * n) +
                                 innerProduct(uhout * t, 1. / hi * pp_tangent * v * t) +
                                 innerProduct(uhout * n, 1. / hi * pp_normal * v * n),
                             Khi, INTEGRAL_BOUNDARY, {2});
            // [Sets uniqueness of the pressure]
            // stokes.addLagrangeMultiplier(
            //   innerProduct(1.,p), 0
            //   , Khi
            // );
            CutFEM<Mesh2> lagr(Vh);
            lagr.add(Ph);
            Rn zero_vec = lagr.rhs_;
            lagr.addLinear(innerProduct(1, p), Khi);
            Rn lag_row(lagr.rhs_);
            lagr.rhs_ = zero_vec;
            lagr.addLinear(innerProduct(1, v * n), interface);
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , Khi, INTEGRAL_BOUNDARY
            // );
            stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);

            FunTest grad2un = grad(grad(u) * n) * n;
            // FunTest grad2pn = grad(grad(p)*n)*n;
            // FunTest grad2divun = grad(grad(div(u))*n)*n;
            stokes.addFaceStabilization(
                // [2k-1/2k+1 and changing to 2k+1/2k+1 gives the same bad effect!]
                // +innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(v))
                // +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(v)*n))
                // +innerProduct(uPenParam*pow(hi,3)*jump(grad2un), jump(grad2un))
                // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
                // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

                // [Pressure robust!]
                +innerProduct(uPenParam * pow(hi, -1) * jump(u), jump(v)) +
                    innerProduct(uPenParam * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
                    innerProduct(uPenParam * pow(hi, 3) * jump(grad2un), jump(grad2un)) -
                    innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
                    innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(v))), jump(grad(q)))
                // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
                // +innerProduct(pPenParam*pow(hi,5)*jump(grad2divun) , jump(grad2pn))
                ,
                Khi, macro);
        }

        // [grad -> eps]
        {
            // stokes.addBilinear(
            //   contractProduct(2*mu*Eps(u),Eps(v))
            //   - innerProduct(p,div(v))
            //   + innerProduct(div(u),q)

            //   - innerProduct(2*omega*u2,v1)
            //   + innerProduct(2*omega*u1,v2)
            //   , Khi
            // );
            // stokes.addLinear(
            //   innerProduct(fh.expression(2),v)
            //   , Khi
            // );
            // stokes.addBilinear(
            //   - innerProduct(average(2*mu*(Eps(u)*t)*n,0.5,0.5), jump(v*t))
            //   + innerProduct(jump(u*t), average(2*mu*(Eps(v)*t)*n,0.5,0.5))
            //   + innerProduct(1./hi*sigma*(jump(u*t)), jump(v*t))
            //   , Khi
            //   , INTEGRAL_INNER_EDGE_2D
            // );
            // stokes.addBilinear(
            //   // - innerProduct(2*mu*Eps(u)*n, v)  //natural
            //   // + innerProduct(u, 2*mu*Eps(v)*n)
            //   + innerProduct(10*1./hi*penaltyParam*u, v)
            //   , interface
            // );
            // stokes.addBilinear(
            //   // - innerProduct(2*mu*Eps(u)*n, v)  //natural
            //   // + innerProduct(u, 2*mu*Eps(v)*n)
            //   + innerProduct(1./hi*penaltyParam*u, v)
            //   , Khi, INTEGRAL_BOUNDARY, {1,3}
            // );
            // stokes.addBilinear(
            //   - innerProduct(2*mu*Eps(u)*n, v)  //natural
            //   + innerProduct(u, 2*mu*Eps(v)*n)
            //   + innerProduct(1./hi*penaltyParam*u, v)
            //   , Khi, INTEGRAL_BOUNDARY, {2,4}
            // );
            // stokes.addLinear(
            //   + innerProduct(uhout.expression(2), 2*mu*Eps(v)*n)
            //   + innerProduct(uhout.expression(2), 1./hi*penaltyParam*v)

            //   , Khi, INTEGRAL_BOUNDARY, {2}
            // );
            // stokes.addLinear(
            //   + innerProduct(uhin.expression(2), 2*mu*Eps(v)*n)
            //   + innerProduct(uhin.expression(2), 1./hi*penaltyParam*v)

            //   , Khi, INTEGRAL_BOUNDARY, {4}
            // );

            // // [Sets uniqueness of the pressure]
            // // R meanP = integral(Khi,exactp,0);
            // // stokes.addLagrangeMultiplier(
            // //   innerProduct(1.,p), 0
            // //   , Khi
            // // );
            // CutFEM<Mesh2> lagr(Vh); lagr.add(Ph);
            // Rn zero_vec = lagr.rhs_;
            // lagr.addLinear(
            //   innerProduct(1.,p)
            //   , Khi
            // );
            // Rn lag_row(lagr.rhs_); lagr.rhs_ = zero_vec;
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , interface
            // );
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , Khi, INTEGRAL_BOUNDARY
            // );
            // stokes.addLagrangeVecToRowAndCol(lag_row,lagr.rhs_,0);

            // FunTest grad2un = grad(grad(u)*n)*n;
            // FunTest grad2pn = grad(grad(p)*n)*n;
            // // FunTest grad2divun = grad(grad(div(u))*n)*n;
            // stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
            //   // +innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(v))
            //   // +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(v)*n))
            //   // // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
            //   // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
            //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

            //   +innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(v))
            //   +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(v)*n))
            //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad2un), jump(grad2un))
            //   -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
            //   +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
            //   -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
            //   +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(v))) , jump(grad(q)))
            //   // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
            //   // +innerProduct(pPenParam*pow(hi,5)*jump(grad2divun) , jump(grad2pn))
            //   , Khi
            //   , macro
            // );
        }

        // std::cout << integral(Khi,exactp,0) << std::endl;
        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_uh     = stokes.rhs_(SubArray(nb_flux_dof, 0));
        Rn_ data_ph     = stokes.rhs_(SubArray(Ph.get_nb_dof(), nb_flux_dof));
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Ph, data_ph);

        // [Post process pressure]
        // R meanP = integral(Khi,p0,0);
        // ExpressionFunFEM<Mesh> fem_p(ph,0,op_id);
        // R meanPfem = integral(Khi,fem_p,0);
        // // std::cout << meanP << std::endl;
        // CutFEM<Mesh2> post(Qh);
        // post.addLinear(
        //   innerProduct(1,q)
        //   , Khi
        // );
        // R area = post.rhs_.sum();
        // ph.v -= meanPfem/area;
        // ph.v += meanP/area;

        ExpressionFunFEM<Mesh> uh_0dx(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> uh_1dy(uh, 1, op_dy);

        // Fun_h test(Vh, rt0_vec);
        // ExpressionFunFEM<Mesh> femSol_0dx(test, 0, op_dx);
        // ExpressionFunFEM<Mesh> femSol_1dy(test, 1, op_dy);

        // [Errors]
        {
            // Fun_h soluErr(Vh, fun_exact_u);
            // Fun_h soluh(Vh, fun_exact_u);
            // soluErr.v -= uh.v;
            // soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Khi, "stokes_" + to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(uh_0dx + uh_1dy, "divergence");
            // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
            // writer.add(soluh, "velocityExact" , 0, 2);
            // writer.add(soluErr, "velocityError" , 0, 2);
            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        // R errU      = L2normCut(uh,fun_exact_u,0,2);
        // R errP      = L2normCut(ph,fun_exact_p,0,1);
        R errDiv    = L2normCut(uh_0dx + uh_1dy, Khi);
        R maxErrDiv = maxNormCut(uh_0dx + uh_1dy, Khi);
        // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

        // // solExVec -= stokesDiv.rhs;
        // for(int i=0;i<solExVec.size();++i){
        //   solExVec(i) = fabs(solExVec(i)-stokesDiv.rhs(i));
        // }
        //
        // Fun_h solEx(mixedSpace, solExVec);
        // writerS.add(solEx,"uh_err", 0, 2);

        // writerS.add(solEx,"uh_ex", 0, 2);
        // writerS.add(solEx,"ph_ex", 2, 1);

        // ul2.push_back(errU);
        // pl2.push_back(errP);
        divl2.push_back(errDiv);
        divmax.push_back(maxErrDiv);
        h.push_back(hi);
        // if(i==0) {convu.push_back(0); convp.push_back(0);}
        // else {
        //   convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
        //   convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
        // }

        nx = 2 * nx - 1;
        ny = 2 * ny - 1;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ')
              << "h"
              // << std::setw(15) << std::setfill(' ') << "err_p"
              // << std::setw(15) << std::setfill(' ') << "conv p"
              // << std::setw(15) << std::setfill(' ') << "err u"
              // << std::setw(15) << std::setfill(' ') << "conv u"
              << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ')
                  << h[i]
                  // << std::setw(15) << std::setfill(' ') << pl2[i]
                  // << std::setw(15) << std::setfill(' ') << convp[i]
                  // << std::setw(15) << std::setfill(' ') << ul2[i]
                  // << std::setw(15) << std::setfill(' ') << convu[i]
                  << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
#endif

#ifdef PROBLEM_UNFITTED_CORIOLIS_STOKES_VORTICITY

namespace Erik_Data_CORIOLIS_STOKES {

R sdBox(const R2 p, const R2 b) {
    R2 absp(abs(p.x), abs(p.y));
    R2 d = absp - b;
    R2 maxd0(max(d.x, 0.0), max(d.y, 0.0));
    return maxd0.norm() + min(max(d.x, d.y), 0.0);
}
R fun_levelSet(const R2 P, const int i) {
    const R2 b(2.0, 1.0); // b.x=width, b.y=height
    const R2 shift(4., 0.);
    return sdBox(P - shift, b) - 0.01;
}

// [Coriolis example]
R fun_y(const R2 P, int i, int dom) { return P.y; }
R fun_div(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 0;
}
R fun_rhs(const R2 P, int i, int dom) {
    // R mu=1;
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 0; // 100*(2-x)*(2-x);
    else if (i == 1)
        return 0;
}
R fun_uin(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return y * (2 - y) / 2; // original paper 4*(2-y)*(y-1)
    else
        return 0;
}
R fun_uout(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 1. / 4 * (2 - y) * (y - 1); // original paper 4*(2-y)*(y-1)
    else
        return 0;
}
} // namespace Erik_Data_CORIOLIS_STOKES

using namespace Erik_Data_CORIOLIS_STOKES;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 11;
    int ny = 11;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 4;
    for (int i = 0; i < iters; ++i) { // i<3

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, 0., 0., 4., 2.);
        const R hi           = 1. / (nx - 1);
        const R penaltyParam = 400; // 400
        double uPenParam     = 1e0; // 1e-1, 1e1
        double pPenParam     = 1e0; // 1e-1, 1e1

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(4);
        Space VELh_(Kh, FEvelocity);
        Space SCAh_(Kh, DataFE<Mesh>::P2);

        Space Uh_(Kh, DataFE<Mesh>::P2);
        Space Vh_(Kh, DataFE<Mesh2>::RT1);
        Space Wh_(Kh, DataFE<Mesh2>::P1dc);

        ActiveMesh<Mesh> Khi(Kh);
        Khi.truncate(interface, -1);

        MacroElement<Mesh> macro(Khi, 1.);

        CutSpace VELh(Khi, VELh_);
        CutSpace SCAh(Khi, SCAh_);

        CutSpace Uh(Khi, Uh_);
        CutSpace Vh(Khi, Vh_);
        CutSpace Wh(Khi, Wh_);

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h uhin(VELh, fun_uin);
        Fun_h uhout(VELh, fun_uout);
        Fun_h yfunh(VELh, fun_y);

        CutFEM<Mesh2> stokes(Uh);
        stokes.add(Vh);
        stokes.add(Wh);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest w(Uh, 1, 0), tau(Uh, 1, 0), u(Vh, 2, 0), v(Vh, 2, 0), p(Wh, 1, 0), q(Wh, 1, 0);
        FunTest u1(Vh, 1, 0), u2(Vh, 1, 1), v1(Vh, 1, 0), v2(Vh, 1, 1);

        R mu    = 0.01;
        R omega = 100; // 50 // original paper 100;
        {
            // [Bulk]
            stokes.addBilinear( // coriolis
                -innerProduct(2 * omega * u2, v1) + innerProduct(2 * omega * u1, v2), Khi);
            stokes.addBilinear( // w = curl u
                innerProduct(1. / mu * w, tau) - innerProduct(u, rotgrad(tau)), Khi);
            stokes.addBilinear( // mu Delta u + grad p
                innerProduct(rotgrad(w), v) - innerProduct(p, div(v)), Khi);
            stokes.addLinear(innerProduct(fh.expression(2), v), Khi);
            stokes.addBilinear(+innerProduct(div(u), q), Khi);
            stokes.addBilinear(+innerProduct(p, v * n) +
                                   innerProduct(1. / hi * penaltyParam * u * n, v * n) // stability
                               ,
                               interface);
            stokes.addBilinear(+innerProduct(p, v * n) +
                                   innerProduct(1. / hi * penaltyParam * u * n, v * n) // stability
                               ,
                               Khi, INTEGRAL_BOUNDARY);
            stokes.addLinear(+innerProduct(uhout * t, tau) + innerProduct(uhout * n, 1. / hi * penaltyParam * v * n),
                             Khi, INTEGRAL_BOUNDARY, {2});
            stokes.addLinear(+innerProduct(uhin * t, tau) + innerProduct(uhin * n, 1. / hi * penaltyParam * v * n), Khi,
                             INTEGRAL_BOUNDARY, {4});
            // [Sets uniqueness of the pressure]
            // stokes.addLagrangeMultiplier(
            //   innerProduct(1.,p), 0
            //   , Khi
            // );
            // [Sets uniqueness of the pressure in another way such that divu = 0]
            CutFEM<Mesh2> lagr(Uh);
            lagr.add(Vh);
            lagr.add(Wh);
            Rn zero_vec = lagr.rhs_;
            lagr.addLinear(innerProduct(1., p), Khi);
            Rn lag_row(lagr.rhs_);
            lagr.rhs_ = zero_vec;
            lagr.addLinear(innerProduct(1, v * n), interface);
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , Khi, INTEGRAL_BOUNDARY
            // );
            stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);
            // [Stabilization]
            FunTest grad2un = grad(grad(u) * n) * n;
            FunTest grad2wn = grad(grad(w) * n) * n;
            // // FunTest grad2pn = grad(grad(p)*n)*n;
            // // FunTest grad2divun = grad(grad(div(u))*n)*n;
            stokes.addFaceStabilization(
                +innerProduct(uPenParam * pow(hi, 3) * jump(grad(w) * n), jump(grad(tau) * n)) +
                    innerProduct(uPenParam * pow(hi, 5) * jump(grad2wn), jump(grad2wn)) +
                    innerProduct(uPenParam * pow(hi, 1) * jump(u), jump(v)) +
                    innerProduct(uPenParam * pow(hi, 3) * jump(grad(u) * n), jump(grad(v) * n)) +
                    innerProduct(uPenParam * pow(hi, 5) * jump(grad2un), jump(grad2un))

                    // +innerProduct(pPenParam*pow(hi,-1)*jump(rotgrad(w)), jump(v))
                    // -innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(rotgrad(tau)))
                    // +innerProduct(uPenParam*pow(hi,1)*jump(grad(rotgrad(w))*n), jump(grad(v)*n))
                    // -innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(rotgrad(tau))*n))

                    - innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
                    innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(u))), jump(grad(q))) +
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v))))

                // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
                // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))
                ,
                Khi, macro);
        }

        // std::cout << integral(Khi,exactp,0) << std::endl;
        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_vort_dof = Uh.get_nb_dof();
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_wh     = stokes.rhs_(SubArray(nb_vort_dof, 0));
        Rn_ data_uh     = stokes.rhs_(SubArray(
            nb_flux_dof, nb_vort_dof)); // Rn_ data_uh = stokes.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));
        Rn_ data_ph     = stokes.rhs_(SubArray(
            Wh.get_nb_dof(),
            nb_vort_dof +
                nb_flux_dof)); // Rn_ data_ph = stokes.rhs_(SubArray(stokes_.get_nb_dof(),nb_vort_dof+nb_flux_dof));
        Fun_h wh(Uh, data_wh);
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Wh, data_ph);

        // [Post process pressure]
        // R meanP = integral(Khi,exactp,0);
        // ExpressionFunFEM<Mesh> fem_p(ph,0,op_id);
        // R meanPfem = integral(Khi,fem_p,0);
        // // std::cout << meanP << std::endl;
        // CutFEM<Mesh2> post(Qh);
        // post.addLinear(
        //   innerProduct(1,q)
        //   , Khi
        // );
        // R area = post.rhs_.sum();
        // ph.v -= meanPfem/area;
        // ph.v += meanP/area;

        ExpressionFunFEM<Mesh> uh_0dx(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> uh_1dy(uh, 1, op_dy);

        // Fun_h test(Vh, rt0_vec);
        // ExpressionFunFEM<Mesh> femSol_0dx(test, 0, op_dx);
        // ExpressionFunFEM<Mesh> femSol_1dy(test, 1, op_dy);

        // [Errors]
        {
            // Fun_h soluErr(Vh, fun_exact_u);
            // Fun_h soluh(Vh, fun_exact_u);
            // soluErr.v -= uh.v;
            // soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Khi, "stokes_" + to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(uh_0dx + uh_1dy, "divergence");
            // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
            // writer.add(soluh, "velocityExact" , 0, 2);
            // writer.add(soluErr, "velocityError" , 0, 2);
            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        // R errU      = L2normCut(uh,fun_exact_u,0,2);
        // R errP      = L2normCut(ph,fun_exact_p,0,1);
        R errDiv    = L2normCut(uh_0dx + uh_1dy, Khi);
        R maxErrDiv = maxNormCut(uh_0dx + uh_1dy, Khi);
        // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

        // ul2.push_back(errU);
        // pl2.push_back(errP);
        divl2.push_back(errDiv);
        divmax.push_back(maxErrDiv);
        h.push_back(hi);
        // if(i==0) {convu.push_back(0); convp.push_back(0);}
        // else {
        //   convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
        //   convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
        // }

        nx = 2 * nx - 1;
        ny = 2 * ny - 1;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ')
              << "h"
              // << std::setw(15) << std::setfill(' ') << "err_p"
              // << std::setw(15) << std::setfill(' ') << "conv p"
              // << std::setw(15) << std::setfill(' ') << "err u"
              // << std::setw(15) << std::setfill(' ') << "conv u"
              << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ')
                  << h[i]
                  // << std::setw(15) << std::setfill(' ') << pl2[i]
                  // << std::setw(15) << std::setfill(' ') << convp[i]
                  // << std::setw(15) << std::setfill(' ') << ul2[i]
                  // << std::setw(15) << std::setfill(' ') << convu[i]
                  << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
#endif

#ifdef PROBLEM_UNFITTED_PRESROB2_STOKES

namespace Erik_Data_UNFITTED_STOKESRT {

R Ra = 1e6;

R fun_levelSet(const R2 P, const int i) { return P.y - 1.; }

// [Example 1 from Neilan pressure robust paper]
R fun_div(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 0;
}
R fun_rhs(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 0;
    else if (i == 1)
        return Ra * (1 - y + 3 * y * y);
}
R fun_exact_u(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 0;
    else
        return 0;
}
R fun_exact_p(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return Ra * (y * y * y - y * y / 2 + y - 7 / 12);
}
} // namespace Erik_Data_UNFITTED_STOKESRT
using namespace Erik_Data_UNFITTED_STOKESRT;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    ProblemOption optionStokes;
    optionStokes.order_space_element_quadrature_ = 9;

    int nx = 11;
    int ny = 11;
    // int d = 2;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp, gradul2;

    int iters = 4;
    for (int i = 0; i < iters; ++i) { // i<3

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, 0., 0., 1., 1.251);
        const R hi       = 1. / (nx - 1);
        const R sigma    = 1e0;      // 1e-2
        const R lambdau  = Ra * 1e5; // Ra*1e8; // 4e3, 8e2
        double uPenParam = 1e0;      // sigma, 6e2, 6e0
        double pPenParam = 1e0;      // sigma, 4e2, 4e0

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(4);
        Space VELh_(Kh, FEvelocity);
        Space SCAh_(Kh, DataFE<Mesh>::P2);

        Space Wh(Kh, DataFE<Mesh2>::BDM1);
        // Space Wh(Kh, FEvelocity);
        Space Qh(Kh, DataFE<Mesh2>::P0); // FOR MIXEDSPACE

        ActiveMesh<Mesh> Khi(Kh);
        Khi.truncate(interface, 1);

        MacroElement<Mesh> macro(Khi, 1);

        CutSpace VELh(Khi, VELh_);
        CutSpace SCAh(Khi, SCAh_);

        CutSpace Vh(Khi, Wh);
        CutSpace Ph(Khi, Qh);

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h gh(VELh, fun_exact_u);
        Fun_h exactph(SCAh, fun_exact_p);
        ExpressionFunFEM<Mesh> exactp(exactph, 0, op_id);

        CutFEM<Mesh2> stokes(Vh);
        stokes.add(Ph);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);

        R mu = 1;
        {
            stokes.addBilinear(
                contractProduct(mu * grad(u), grad(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q), Khi);
            stokes.addLinear(innerProduct(fh.expression(2), v), Khi);
            // [NECESSARY FOR RT / BDM ELEMENTS]
            stokes.addBilinear(-innerProduct(mu * average(grad(u * t) * n, 0.5, 0.5), jump(v * t)) +
                                   innerProduct(jump(u * t), mu * average(grad(v * t) * n, 0.5, 0.5)) +
                                   innerProduct(1. / hi * sigma * jump(u * t), jump(v * t)),
                               Khi, INTEGRAL_INNER_EDGE_2D);
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v)           // natural
                                   + innerProduct(u, mu * grad(v) * n)      // symmetry
                                   + innerProduct(1. / hi * lambdau * u, v) // stability

                                   + innerProduct(p, v * n) // natural
                               ,
                               interface);
            // stokes.addLinear(
            //   + innerProduct(gh.expression(2), mu*grad(v)*n)
            //   + innerProduct(gh.expression(2), 1./hi*lambdau*v)
            //   , interface
            // );
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v)           // natural
                                   + innerProduct(u, mu * grad(v) * n)      // symmetry
                                   + innerProduct(1. / hi * lambdau * u, v) // stability

                                   + innerProduct(p, v * n) // natural
                               ,
                               Khi, INTEGRAL_BOUNDARY);
            // stokes.addLinear(
            //   + innerProduct(gh.expression(2), mu*grad(v)*n)
            //   + innerProduct(gh.expression(2), 1./hi*lambdau*v)
            //   , Khi, INTEGRAL_BOUNDARY
            // );
            // [Sets uniqueness of the pressure]
            R meanP = integral(Khi, exactp, 0);
            // stokes.addLagrangeMultiplier(
            //   innerProduct(1, p), meanP
            //   , Khi
            // );
            CutFEM<Mesh2> lagr(Vh);
            lagr.add(Ph);
            lagr.addLinear(innerProduct(1, p), Khi);
            Rn lag_row(lagr.rhs_);
            lagr.rhs_ = 0.;
            lagr.addLinear(innerProduct(1, v * n), interface);
            lagr.addLinear(innerProduct(1, v * n), Khi, INTEGRAL_BOUNDARY);
            stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, meanP);

            FunTest grad2un = grad(grad(u) * n) * n;
            // // FunTest grad2pn = grad(grad(p)*n)*n;
            // // FunTest grad2divun = grad(grad(div(u))*n)*n;
            stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
                                         //  innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(v))
                                         // +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(v)*n))
                                         // +innerProduct(uPenParam*pow(hi,3)*jump(grad2un), jump(grad2un))
                                         // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
                                         // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

                +innerProduct(uPenParam * pow(hi, -1) * jump(u), jump(v)) +
                    innerProduct(uPenParam * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
                    innerProduct(uPenParam * pow(hi, 3) * jump(grad2un), jump(grad2un)) -
                    innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
                    innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(v))), jump(grad(q))),
                Khi, macro);
        }

        // std::cout << integral(Khi,exactp,0) << std::endl;
        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_uh     = stokes.rhs_(SubArray(nb_flux_dof, 0));
        Rn_ data_ph     = stokes.rhs_(SubArray(Ph.get_nb_dof(), nb_flux_dof));
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Ph, data_ph);

        ExpressionFunFEM<Mesh> uh_0dx(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> uh_1dy(uh, 1, op_dy);
        ExpressionFunFEM<Mesh> uh_0dy(uh, 0, op_dy);
        ExpressionFunFEM<Mesh> uh_1dx(uh, 1, op_dx);

        // Fun_h test(Vh, rt0_vec);
        // ExpressionFunFEM<Mesh> femSol_0dx(test, 0, op_dx);
        // ExpressionFunFEM<Mesh> femSol_1dy(test, 1, op_dy);

        // [Errors]
        Fun_h soluErr(Vh, fun_exact_u);
        Fun_h soluh(Vh, fun_exact_u);
        soluErr.v -= uh.v;
        soluErr.v.map(fabs);
        // Fun_h divSolh(Wh, fun_div);
        // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

        Paraview<Mesh> writer(Khi, "stokes_" + to_string(i) + ".vtk");
        writer.add(uh, "velocity", 0, 2);
        writer.add(ph, "pressure", 0, 1);
        writer.add(uh_0dx + uh_1dy, "divergence");
        writer.add(soluh, "velocityExact", 0, 2);
        writer.add(soluErr, "velocityError", 0, 2);
        // writer.add(solh, "velocityError" , 0, 2);

        // writer.add(fabs(femDiv, "divergenceError");

        R errU      = L2normCut(uh, fun_exact_u, 0, 2);
        R errGradU  = integral(Khi, uh_0dx * uh_0dx + uh_0dy * uh_0dy + uh_1dx * uh_1dx + uh_1dy * uh_1dy, 0);
        R errP      = L2normCut(ph, fun_exact_p, 0, 1);
        R errDiv    = L2normCut(uh_0dx + uh_1dy, fun_div, Khi);
        R maxErrDiv = maxNormCut(uh_0dx + uh_1dy, fun_div, Khi);
        // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

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
        gradul2.push_back(errGradU);
        pl2.push_back(errP);
        divl2.push_back(errDiv);
        divmax.push_back(maxErrDiv);
        h.push_back(hi);
        if (i == 0) {
            convu.push_back(0);
            convp.push_back(0);
        } else {
            convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
        }

        nx = 2 * nx - 1;
        ny = 2 * ny - 1;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ') << "err_p"
              << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ') << "err u"
              << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ') << "err maxdivu" << std::setw(15) << std::setfill(' ')
              << "err gradu" << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
                  << pl2[i] << std::setw(15) << std::setfill(' ') << convp[i] << std::setw(15) << std::setfill(' ')
                  << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ') << divmax[i] << std::setw(15) << std::setfill(' ') << gradul2[i]
                  << std::endl;
    }
}
#endif

#ifdef PROBLEM_UNFITTED_PRESROB2_STOKES_VORTICITY

namespace Erik_Data_UNFITTED_STOKES_VORTICITY {

R Ra = 1e16;

R fun_levelSet(const R2 P, const int i) { return P.y - 1.; }

// [Example 1 from Neilan pressure robust paper]
R fun_div(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 0;
}
R fun_rhs(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 0;
    else if (i == 1)
        return Ra * (1 - y + 3 * y * y);
}
R fun_exact_u(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 0;
    else
        return 0;
}
R fun_exact_p(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return Ra * (y * y * y - y * y / 2 + y - 7 / 12);
}
} // namespace Erik_Data_UNFITTED_STOKES_VORTICITY
using namespace Erik_Data_UNFITTED_STOKES_VORTICITY;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 11;
    int ny = 11;
    // int d = 2;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp, gradul2;

    int iters = 4;
    for (int i = 0; i < iters; ++i) {

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, 0.0, 0.0, 1.0, 1.25001);
        const R hi           = 1. / (nx - 1); // 1./(nx-1)
        const R penaltyParam = Ra * 1e5;      // 4e3, 8e2

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(4);
        Space VELh_(Kh, FEvelocity);
        Space SCAh_(Kh, DataFE<Mesh>::P4);

        Space Uh_(Kh, DataFE<Mesh>::P1); // Nedelec order 0 type 1
        Space Vh_(Kh, DataFE<Mesh2>::RT0);
        Space Wh_(Kh, DataFE<Mesh2>::P0);

        // [Remove exterior]
        ActiveMesh<Mesh> Khi(Kh);
        Khi.truncate(interface, 1);

        MacroElement<Mesh> macro(Khi, 1); // we use 0.25 for vorticity BC2

        CutSpace VELh(Khi, VELh_);
        CutSpace SCAh(Khi, SCAh_);

        CutSpace Uh(Khi, Uh_);
        CutSpace Vh(Khi, Vh_);
        CutSpace Wh(Khi, Wh_);

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h u0(VELh, fun_exact_u);
        Fun_h p0(SCAh, fun_exact_p);
        ExpressionFunFEM<Mesh> exactp(p0, 0, op_id);

        CutFEM<Mesh2> stokes(Uh);
        stokes.add(Vh);
        stokes.add(Wh);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest w(Uh, 1, 0), tau(Uh, 1, 0), u(Vh, 2, 0), v(Vh, 2, 0), p(Wh, 1, 0), q(Wh, 1, 0);

        R mu = 1;
        // [Bulk]
        stokes.addBilinear( // w = curl u
            innerProduct(1. / mu * w, tau) - innerProduct(u, rotgrad(tau)), Khi);
        stokes.addBilinear( // mu Delta u + grad p
            innerProduct(rotgrad(w), v) - innerProduct(p, div(v)), Khi);
        stokes.addLinear(+innerProduct(fh.expression(2), v), Khi);
        stokes.addBilinear(+innerProduct(div(u), q), Khi);
        // [Dirichlet Velocity BC u0 = 0]
        {
            stokes.addBilinear(+innerProduct(p, v * n) +
                                   innerProduct(1. / hi * penaltyParam * u * n, v * n) // stability
                               ,
                               Khi, INTEGRAL_BOUNDARY);
            stokes.addBilinear(+innerProduct(p, v * n) +
                                   innerProduct(1. / hi * penaltyParam * u * n, v * n) // stability
                               ,
                               interface);

            // // [Sets uniqueness of the pressure]
            R meanP = integral(Khi, exactp, 0);
            // std::cout << meanP << std::endl;
            // stokes.addLagrangeMultiplier(
            //   innerProduct(1, p), meanP
            //   , Khi
            // );
            // [Sets uniqueness of the pressure in another way such that divu = 0]
            CutFEM<Mesh2> lagr(Uh);
            lagr.add(Vh);
            lagr.add(Wh);
            lagr.addLinear(innerProduct(1, p), Khi);
            Rn lag_row(lagr.rhs_);
            lagr.rhs_ = 0.;
            lagr.addLinear(innerProduct(1, v * n), interface);
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , Khi, INTEGRAL_BOUNDARY
            // );
            stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, meanP);
            // [Stabilization]
            double wPenParam = 1e0; // 1e1
            double uPenParam = 1e0; // 1e-1 ~ 1/penParam (2e0 for (0,lamm,0))
            double pPenParam = 1e0; // 1e0 (2e0 for (0,lamm,0))
            FunTest grad2un  = grad(grad(u) * n) * n;
            FunTest grad2wn  = grad(grad(w) * n) * n;
            // // FunTest grad2pn = grad(grad(p)*n)*n;
            // // FunTest grad2divun = grad(grad(div(u))*n)*n;
            stokes.addFaceStabilization(
                /* "Primal" stab: (lw,0,la) */
                // innerProduct(uPenParam*pow(hi,1)*jump(w), jump(tau)) // [w in P1, continuous]
                +innerProduct(wPenParam * pow(hi, 3) * jump(grad(w) * n), jump(grad(tau) * n)) +
                    innerProduct(uPenParam * pow(hi, 5) * jump(grad2wn), jump(grad2wn)) +
                    innerProduct(uPenParam * pow(hi, 1) * jump(u),
                                 jump(v)) // [maybe should be 2k-1 if can scale pressure also]
                    + innerProduct(uPenParam * pow(hi, 3) * jump(grad(u) * n), jump(grad(v) * n)) +
                    innerProduct(uPenParam * pow(hi, 5) * jump(grad2un), jump(grad2un))

                    - innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
                    innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(u))), jump(grad(q)))

                // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
                // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

                ,
                Khi, macro);
        }

        // std::cout << integral(Khi,exactp,0) << std::endl;
        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_vort_dof = Uh.get_nb_dof();
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_wh     = stokes.rhs_(SubArray(nb_vort_dof, 0));
        Rn_ data_uh     = stokes.rhs_(SubArray(
            nb_flux_dof, nb_vort_dof)); // Rn_ data_uh = stokes.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));
        Rn_ data_ph     = stokes.rhs_(SubArray(
            Wh.get_nb_dof(),
            nb_vort_dof +
                nb_flux_dof)); // Rn_ data_ph = stokes.rhs_(SubArray(stokes_.get_nb_dof(),nb_vort_dof+nb_flux_dof));
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

        ExpressionFunFEM<Mesh> uh0_dx(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> uh0_dy(uh, 0, op_dy);
        ExpressionFunFEM<Mesh> uh1_dx(uh, 1, op_dx);
        ExpressionFunFEM<Mesh> uh1_dy(uh, 1, op_dy);
        // [Errors]
        {
            // Fun_h solw(Uh, fun_exact_w);
            Fun_h solu(Vh, fun_exact_u);
            Fun_h soluErr(Vh, fun_exact_u);
            Fun_h solp(Wh, fun_exact_p);
            soluErr.v -= uh.v;
            soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Khi, "stokes_" + to_string(i) + ".vtk");
            writer.add(wh, "vorticity", 0, 1);
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(uh0_dx + uh1_dy, "divergence");
            // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
            writer.add(solp, "pressureExact", 0, 1);
            writer.add(solu, "velocityExact", 0, 2);
            writer.add(soluErr, "velocityError", 0, 2);
            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        // R errW      = L2normCut(wh,fun_exact_w,0,1);
        R errU      = L2normCut(uh, fun_exact_u, 0, 2);
        R errGradU  = integral(Khi, uh0_dx * uh0_dx + uh0_dy * uh0_dy + uh1_dx * uh1_dx + uh1_dy * uh1_dy, 0);
        R errP      = L2normCut(ph, fun_exact_p, 0, 1);
        R errDiv    = L2normCut(uh0_dx + uh1_dy, fun_div, Khi);
        R maxErrDiv = maxNormCut(uh0_dx + uh1_dy, fun_div, Khi);
        // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

        ul2.push_back(errU);
        pl2.push_back(errP);
        divl2.push_back(errDiv);
        divmax.push_back(maxErrDiv);
        gradul2.push_back(errGradU);
        h.push_back(hi);
        if (i == 0) {
            convu.push_back(0);
            convp.push_back(0);
        } else {
            convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
        }
        nx = 2 * nx - 1;
        ny = 2 * ny - 1;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ') << "err p"
              << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ') << "err u"
              << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ') << "err maxdivu" << std::setw(15) << std::setfill(' ')
              << "conv err gradu" << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
                  << pl2[i] << std::setw(15) << std::setfill(' ') << convp[i] << std::setw(15) << std::setfill(' ')
                  << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ') << divmax[i] << std::setw(15) << std::setfill(' ') << gradul2[i]
                  << std::endl;
    }
}
#endif

#ifdef PROBLEM_UNFITTED_STOKES3D

namespace Erik_Data_UNFITTED_STOKES3D {

R3 shift(0., 0., 0.);
R fun_levelSet(const R3 P, int i) {
    return sqrt((P.x - shift.x) * (P.x - shift.x) + (P.y - shift.y) * (P.y - shift.y) +
                (P.z - shift.z) * (P.z - shift.z)) -
           2. / 3;
}
R fun_rhs(const R3 P, int i) { return 0; }
R fun_boundary(const R3 P, int i) { return (i == 0) ? 0.5 * P.z : 0; }

R fun_kkk(const R3 P, int i) { return 0.5 * P.z; }
} // namespace Erik_Data_UNFITTED_STOKES3D
using namespace Erik_Data_UNFITTED_STOKES3D;

int main(int argc, char **argv) {
    typedef TestFunction<3> FunTest;
    typedef FunFEM<Mesh3> Fun_h;
    typedef Mesh3 Mesh;
    typedef ActiveMeshT3 CutMesh;
    typedef FESpace3 Space;
    typedef CutFESpaceT3 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    const int d = 3;
    int nx      = 10;
    int ny      = 10;
    int nz      = 10;
    Mesh3 Kh(nx, ny, nz, -1., -1., -1., 2., 2., 2.);
    const R hi           = 1. / (nx - 1); // 1./(nx-1)
    const R penaltyParam = 4e3;           // 4e3, 8e2

    Space Uh_(Kh, DataFE<Mesh>::Ned0); // Nedelec order 0 type 1
    Space Vh_(Kh, DataFE<Mesh>::RT0);
    Space Wh_(Kh, DataFE<Mesh>::P0);

    Fun_h fh0(Uh_, fun_kkk);

    // std::cout << fh0.v << std::endl;

    //   Paraview<Mesh> writer(Kh, "stokes3D_"+to_string(0)+".vtk");
    //   writer.add(fh0, "kkk" , 0, 3);

    // return 0;

    // FEM<Mesh3> stokes3D_({&Uh_, &Vh_, &Wh_}); std::getchar();

    Space Lh(Kh, DataFE<Mesh>::P1);
    Fun_h levelSet(Lh, fun_levelSet);
    InterfaceLevelSet<Mesh> interface(Kh, levelSet);

    // [Remove exterior]
    ActiveMesh<Mesh> Khi(Kh);
    Khi.truncate(interface, 1);

    CutSpace Uh(Khi, Uh_);
    CutSpace Vh(Khi, Vh_);
    CutSpace Wh(Khi, Wh_);

    // Interpolate data
    Fun_h fh(Vh, fun_rhs);
    Fun_h u0(Vh, fun_boundary);

    // Init system matrix & assembly
    CutFEM<Mesh> stokes3D(Uh);
    stokes3D.add(Vh);
    stokes3D.add(Wh);
    // CutFEM<Mesh> stokes3D(Vh); stokes3D.add(Wh);

    Normal n;
    /* Syntax:
    FunTest (fem space, #components, place in space)
    */
    FunTest w(Uh, 3, 0), tau(Uh, 3, 0);
    FunTest u(Vh, 3, 0), v(Vh, 3, 0), p(Wh, 1, 0), q(Wh, 1, 0);
    R mu = 1;

    // // [Bulk]
    stokes3D.addBilinear( // w = curl u
        innerProduct(1. / mu * w, tau) - innerProduct(u, curl(tau)), Khi);
    stokes3D.addBilinear( // mu Delta u + grad p
        innerProduct(curl(w), v) - innerProduct(p, div(v)), Khi);
    stokes3D.addLinear(+innerProduct(fh.expression(3), v), Khi);
    stokes3D.addBilinear(+innerProduct(div(u), q), Khi);
    // [Dirichlet Velocity BC]
    // const MeshParameter &itf_h(Parameter::measureIntegral);
    stokes3D.addBilinear( // int_Omg grad(p)*v = int_itf p v*t - int_Omg p div(v)
        +innerProduct(p, v * n) + innerProduct(1. / hi * penaltyParam * u * n, v * n)
        // - innerProduct(u*t, tau)
        // + innerProduct(w, v*t)
        // + innerProduct(1./hi*penaltyParam*u, v)
        ,
        interface);
    stokes3D.addLinear(
        // + innerProduct(cross(n,u0), tau) // [wtf why is + now correct..?]
        // + innerProduct(u0*n, 1./hi*penaltyParam*v*n)
        // - innerProduct(u0*t,tau)
        +innerProduct(u0 * n, 1. / hi * penaltyParam * v), interface);

    // [Sets uniqueness of the pressure]
    // R meanP = integral(Khi,p0,0);
    stokes3D.addLagrangeMultiplier(innerProduct(1, p), 0, Khi);
    // // [Sets uniqueness of the pressure in another way such that divu = 0]
    // CutFEM<Mesh> lagr(Uh); lagr.add(Vh); lagr.add(Wh);
    // Rn zero_vec = lagr.rhs_;
    // lagr.addLinear(
    //   innerProduct(1, p)
    //   , Khi
    // );
    // Rn lag_row(lagr.rhs_); lagr.rhs_ = zero_vec;
    // lagr.addLinear(
    //   innerProduct(1, v*n)
    //   , interface
    // );
    // stokes3D.addLagrangeVecToRowAndCol(lag_row,lagr.rhs_,0);
    // // // [Stabilization]
    // // double wPenParam = 1e1; // 1e1
    // // double uPenParam = 1e1; // 1e-1 ~ 1/penParam (2e0 for (0,lamm,0))
    // // double pPenParam = 1e1; // 1e0 (2e0 for (0,lamm,0))
    // // FunTest grad2un = grad(grad(u)*n)*n;
    // // FunTest grad2wn = grad(grad(w)*n)*n;
    // // // // FunTest grad2pn = grad(grad(p)*n)*n;
    // // // // FunTest grad2divun = grad(grad(div(u))*n)*n;
    // // stokes3D.addFaceStabilization(
    // //   /* "Primal" stab: (lw,0,la) */
    // //   // innerProduct(uPenParam*pow(hi,1)*jump(w), jump(tau)) // [w in P1, continuous]
    // //   +innerProduct(wPenParam*pow(hi,3)*jump(grad(w)*n), jump(grad(tau)*n))
    // //   +innerProduct(uPenParam*pow(hi,5)*jump(grad2wn), jump(grad2wn))
    // //   +innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v)) // [maybe should be 2k-1 if can scale pressure also]
    // //   +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
    // //   +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
    // //   -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
    // //   +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
    // //   -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
    // //   +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))) , jump(grad(q)))

    // //   /* Mixed stab: (0,lm,0) */
    // //   // innerProduct(uPenParam*pow(hi,1)*jump(w), jump(tau)) // [w in P1, continuous]
    // //   // +innerProduct(wPenParam*pow(hi,3)*jump(grad(w)*n), jump(grad(tau)*n))
    // //   // +innerProduct(wPenParam*pow(hi,5)*jump(grad2wn), jump(grad2wn))
    // //   // +innerProduct(uPenParam*pow(hi,1)*jump(curl(w)), jump(v))
    // //   // -innerProduct(uPenParam*pow(hi,1)*jump(u), jump(curl(tau)))
    // //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad(curl(w))), jump(grad(v)))
    // //   // -innerProduct(uPenParam*pow(hi,3)*jump(grad(u)), jump(grad(curl(tau))))
    // //   // -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
    // //   // +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
    // //   // -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
    // //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))), jump(grad(q)))

    // //   , Khi
    // //   , macro
    // // );

    stokes3D.solve();

    // EXTRACT SOLUTION
    int nb_vort_dof = Uh.get_nb_dof();
    int nb_flux_dof = Vh.get_nb_dof();
    Rn_ data_wh     = stokes3D.rhs_(SubArray(nb_vort_dof, 0));
    Rn_ data_uh     = stokes3D.rhs_(SubArray(
        nb_flux_dof, nb_vort_dof)); // Rn_ data_uh = stokes.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));
    Rn_ data_ph     = stokes3D.rhs_(SubArray(
        Wh.get_nb_dof(),
        nb_vort_dof +
            nb_flux_dof)); // Rn_ data_ph = stokes.rhs_(SubArray(stokes_.get_nb_dof(),nb_vort_dof+nb_flux_dof));
    Fun_h wh(Uh, data_wh);
    Fun_h uh(Vh, data_uh);
    Fun_h ph(Wh, data_ph);

    //   // [Post process pressure]
    //   R meanP = integral(Khi,exactp,0);
    //   ExpressionFunFEM<Mesh> fem_p(ph,0,op_id);
    //   R meanPfem = integral(Khi,fem_p,0);
    //   // std::cout << meanP << std::endl;
    //   CutFEM<Mesh2> post(Wh);
    //   post.addLinear(
    //     innerProduct(1,q)
    //     , Khi
    //   );
    //   R area = post.rhs_.sum();
    //   ph.v -= meanPfem/area;
    //   ph.v += meanP/area;

    ExpressionFunFEM<Mesh> dx_uh0(uh, 0, op_dx);
    ExpressionFunFEM<Mesh> dy_uh1(uh, 1, op_dy);

    // [Paraview]
    {
        // Fun_h solw(Uh, fun_exact_w);
        // Fun_h solu(Vh, fun_exact_u); Fun_h soluErr(Vh, fun_exact_u);
        // Fun_h solp(Wh, fun_exact_p);
        // soluErr.v -= uh.v;
        // soluErr.v.map(fabs);
        // Fun_h divSolh(Wh, fun_div);
        // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

        // Paraview<Mesh> writer(Khi, "stokes_"+to_string(i)+".vtk");
        Paraview<Mesh> writer(Khi, "stokes3D_.vtk");
        writer.add(wh, "vorticity", 0, 3);
        writer.add(uh, "velocity", 0, 3);
        writer.add(ph, "pressure", 0, 1);
        writer.add(dx_uh0 + dy_uh1, "divergence");
        // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
        // writer.add(solp, "pressureExact" , 0, 1);
        // writer.add(solu, "velocityExact" , 0, 2);
        // writer.add(soluErr, "velocityError" , 0, 2);
        // writer.add(solh, "velocityError" , 0, 2);

        // writer.add(fabs(femDiv, "divergenceError");
    }
}
#endif