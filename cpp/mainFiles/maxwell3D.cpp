#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
// #include "../util/cputime.h"
#ifdef USE_MPI
#include "cfmpi.hpp"
#endif

#include "finiteElement.hpp"
#include "baseProblem.hpp"
#include "paraview.hpp"
#include "../num/matlab.hpp"

// #include "../num/gnuplot.hpp"

#define PROBLEM_UNFITTED_maxwell3D

#ifdef PROBLEM_UNFITTED_maxwell3D

using namespace globalVariable;

namespace Erik_Data_UNFITTED_maxwell3D {
R k = 1.;

R eps_r = 1.;

R3 shift(0.5, 0.5, 0.5);

R fun_levelSet(double *P, int i) {

    // return sqrt((P[0] - shift.x) * (P[0] - shift.x) + (P[1] - shift.y) * (P[1] - shift.y) +
    //             (P[2] - shift.z) * (P[2] - shift.z)) -
    //        0.35 + Epsilon;
    return (P[0] - shift.x) * (P[0] - shift.x) + (P[1] - shift.y) * (P[1] - shift.y) +
                (P[2] - shift.z) * (P[2] - shift.z) -
           0.35*0.35 + Epsilon;
}

R fun_rhs(double *P, int i) {
    R x = P[0], y = P[1], z = P[2];
    if (i == 0)
        return 2 * pi * pi * sin(y * pi) * sin(z * pi) - eps_r * (2 * x - 1) -
               eps_r * sin(y * pi) * sin(z * pi) * k * k;
    else if (i == 1)
        return 2 * pi * pi * sin(x * pi) * sin(z * pi) - eps_r * (2 * y - 1) -
               eps_r * sin(x * pi) * sin(z * pi) * k * k;
    else
        return 2 * pi * pi * sin(x * pi) * sin(y * pi) - eps_r * (2 * z - 1) -
               eps_r * sin(x * pi) * sin(y * pi) * k * k;
}

// R fun_boundary(double *P, int i) {
//     if (i == 0)
//         return 0.;
//     else if ()
// }

R fun_exact_u(double *P, int i, int dom) {
    R x = P[0], y = P[1], z = P[2];
    if (i == 0)
        return sin(pi * y) * sin(pi * z);
    else if (i == 1)
        return sin(pi * x) * sin(pi * z);
    else
        return sin(pi * x) * sin(pi * y);
}

R fun_exact_p(double *P, int i, int dom) {

    return (P[0] - shift.x) * (P[0] - shift.x) + (P[1] - shift.y) * (P[1] - shift.y) +
           (P[2] - shift.z) * (P[2] - shift.z) - 0.35 * 0.35;
}

R fun_kkk(double *P, int i) { return 0.5 * P[2]; }

} // namespace Erik_Data_UNFITTED_maxwell3D

using namespace Erik_Data_UNFITTED_maxwell3D;

int main(int argc, char **argv) {

    typedef TestFunction<Mesh3> FunTest;

    typedef FunFEM<Mesh3> Fun_h;

    typedef Mesh3 Mesh;

    typedef ActiveMeshT3 CutMesh;

    typedef FESpace3 Space;

    typedef CutFESpaceT3 CutSpace;

    const double cpubegin = CPUtime();

    MPIcf cfMPI(argc, argv);

    const int d = 3;

    int nx = 5;
    int ny = 5;
    int nz = 5;

    std::vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 4;
    for (int i = 0; i < iters; ++i) {

        Mesh3 Kh(nx, ny, nz, 0., 0., 0., 1., 1., 1.);

        const R hi = 1. / (nx - 1); // 1./(nx-1)

        //const R penaltyParam = 4e3; // 4e3, 8e2

        Space Uh_(Kh,
                  DataFE<Mesh>::Ned0); // Nedelec order 0 type 1

        Space Vh_(Kh, DataFE<Mesh>::RT0);
        Space Wh_(Kh, DataFE<Mesh>::P0);

        Lagrange3 VelocitySpace(2);
        Space Vel_h(Kh, VelocitySpace);

        // std::cout << fh0.v << std::endl;

        // Paraview<Mesh> writer(Kh, "maxwell3D_" + std::to_string(0) + ".vtk");

        // writer.add(fh0, "kkk", 0, 3);

        // FEM<Mesh3> maxwell3D_({&Uh_, &Vh_, &Wh_}); std::getchar();

        Space Lh(Kh, DataFE<Mesh>::P1);

        Fun_h levelSet(Lh, fun_levelSet);

        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        // [Remove exterior]

        ActiveMesh<Mesh> Khi(Kh);

        Khi.truncate(interface, 1);

        CutSpace Uh(Khi, Uh_);

        CutSpace Vh(Khi, Vh_);

        CutSpace Velh(Khi, Vel_h);

        CutSpace Wh(Khi, Wh_);

        // Interpolate data

        Fun_h fh(Velh, fun_rhs);

        Fun_h u0(Velh, fun_exact_u);

        // Init system matrix & assembly

        CutFEM<Mesh> maxwell3D(Uh);
        maxwell3D.add(Vh);
        maxwell3D.add(Wh);

        // CutFEM<Mesh> maxwell3D(Vh); maxwell3D.add(Wh);

        Normal n;

        /* Syntax:

        FunTest (fem space, #components, place in space)

        */

        FunTest w(Uh, 3, 0), tau(Uh, 3, 0);

        FunTest u(Vh, 3, 0), v(Vh, 3, 0), p(Wh, 1, 0), q(Wh, 1, 0);

        R mu = 1.;

        // maxwell3D.addBilinear( // w = curl u

        // innerProduct(curl(u), v)

        // , Khi

        // );

        // // [Bulk]

        // Eq 1
        maxwell3D.addBilinear( // w = curl u
            +innerProduct(1. / mu * w, tau) - innerProduct(u, curl(tau)), Khi);
        maxwell3D.addLinear(+innerProduct(cross(u0, n), tau), interface);

        // Eq 2
        maxwell3D.addBilinear( // mu Delta u + grad p
            +innerProduct(curl(w), v) - innerProduct(k * k * eps_r * u, v) + innerProduct(p, div(v)), Khi);
        maxwell3D.addLinear(+innerProduct(fh.exprList(), v), Khi);

        // Eq 3
        maxwell3D.addBilinear(-innerProduct(div(u), q), Khi);

        // // [Dirichlet Velocity BC]

        // const MeshParameter &itf_h(Parameter::measureIntegral);

        // maxwell3D.addBilinear( // int_Omg grad(p)*v = int_itf p v*t - int_Omg p div(v)

        // + innerProduct(p, v*n)

        // + innerProduct(1./hi*penaltyParam*u*n, v*n)

        // // - innerProduct(u*t, tau)

        // // + innerProduct(w, v*t)

        // // + innerProduct(1./hi*penaltyParam*u, v)

        // , interface

        // );

        // // [Sets uniqueness of the pressure]

        // // R meanP = integral(Khi,exactp,0);

        // // maxwell3D.addLagrangeMultiplier(

        // // innerProduct(1, p), 0

        // // , Khi

        // // );

        // // [Sets uniqueness of the pressure in another way such that divu = 0]

        // CutFEM<Mesh> lagr(Uh); lagr.add(Vh); lagr.add(Wh);

        // Rn zero_vec = lagr.rhs_;

        // lagr.addLinear(

        // innerProduct(1, p)

        // , Khi

        // );

        // Rn lag_row(lagr.rhs_); lagr.rhs_ = zero_vec;

        // lagr.addLinear(

        // innerProduct(1, v*n)

        // , interface

        // );

        // maxwell3D.addLagrangeVecToRowAndCol(lag_row,lagr.rhs_,0);

        // [Stabilization]

        double wPenParam = 1e0; // 1e1
        double uPenParam = 1e-2; // 1e-1 ~ 1/penParam (2e0 for (0,lamm,0))
        double pPenParam = 1e0;  // 1e0 (2e0 for (0,lamm,0))

        FunTest grad2un = grad(grad(u) * n) * n;

        FunTest grad2wn = grad(grad(w) * n) * n;

        // // FunTest grad2pn = grad(grad(p)*n)*n;

        // // FunTest grad2divun = grad(grad(div(u))*n)*n;

        maxwell3D.addFaceStabilization(

            /* "Primal" stab: (lw,0,la) */

            // innerProduct(uPenParam*pow(hi,1)*jump(w), jump(tau)) // [w in P1, continuous]
            // +innerProduct(wPenParam*pow(hi,3)*jump(grad(w)*n), jump(grad(tau)*n))
            // +innerProduct(uPenParam*pow(hi,5)*jump(grad2wn), jump(grad2wn))
            // +innerProduct(uPenParam*pow(hi,1)*jump(u), jump(v))
            // +innerProduct(uPenParam*pow(hi,3)*jump(grad(u)*n), jump(grad(v)*n))
            // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
            // -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
            // +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
            // -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
            // +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(u))) , jump(grad(q)))

            /* Mixed stab: (0,lm,0) */
            // innerProduct(uPenParam*pow(hi,1)*jump(w), jump(tau)) // [w in P1, continuous]
            // +innerProduct(wPenParam*pow(hi,3)*jump(grad(w)*n), jump(grad(tau)*n))
            // +innerProduct(wPenParam*pow(hi,5)*jump(grad2wn), jump(grad2wn))

                +innerProduct(uPenParam * pow(hi, 1) * jump(curl(w)), jump(v)) -
                innerProduct(uPenParam * pow(hi, 1) * jump(u), jump(curl(tau))) +
                innerProduct(uPenParam * pow(hi, 3) * jump(grad(curl(w))), jump(grad(v))) -
                innerProduct(uPenParam * pow(hi, 3) * jump(grad(u)), jump(grad(curl(tau)))) +
                innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) -
                innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) +
                innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) -
                innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(u))), jump(grad(q))),
            Khi
            //, macro

        );

        matlab::Export(maxwell3D.mat_[0], "mat" + std::to_string(i) + ".dat");

        maxwell3D.solve("mumps");

        // EXTRACT SOLUTION

        int nb_vort_dof = Uh.get_nb_dof();

        int nb_flux_dof = Vh.get_nb_dof();

        Rn_ data_wh = maxwell3D.rhs_(SubArray(nb_vort_dof, 0));

        Rn_ data_uh = maxwell3D.rhs_(SubArray(
            nb_flux_dof, nb_vort_dof)); // Rn_ data_uh = stokes.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));

        Rn_ data_ph = maxwell3D.rhs_(SubArray(Wh.get_nb_dof(), nb_vort_dof + nb_flux_dof)); //

        Fun_h wh(Uh, data_wh);

        Fun_h uh(Vh, data_uh);

        Fun_h ph(Wh, data_ph);

        // // [Post process pressure]

        // R meanP = integral(Khi,exactp,0);

        // ExpressionFunFEM<Mesh> fem_p(ph,0,op_id);

        // R meanPfem = integral(Khi,fem_p,0);

        // // std::cout << meanP << std::endl;

        // CutFEM<Mesh2> post(Wh);

        // post.addLinear(

        // innerProduct(1,q)

        // , Khi

        // );

        // R area = post.rhs_.sum();

        // ph.v -= meanPfem/area;

        // ph.v += meanP/area;

        ExpressionFunFEM<Mesh> dx_uh0(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> dy_uh1(uh, 1, op_dy);
        ExpressionFunFEM<Mesh> dz_uh2(uh, 2, op_dz);

        auto uh_0dx = dx(uh.expr(0));
        auto uh_1dy = dy(uh.expr(1));
        auto uh_2dz = dz(uh.expr(2));

        // [Paraview]

        {

            // Fun_h solw(Uh, fun_exact_w);

            Fun_h solu(Velh, fun_exact_u);
            Fun_h soluErr(Vh, fun_exact_u);

            Fun_h solp(Wh, fun_exact_p);

            soluErr.v -= uh.v;

            soluErr.v.map(fabs);

            // Fun_h divSolh(Wh, fun_div);

            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Khi, "maxwell_" + std::to_string(i) + ".vtk");

            writer.add(wh, "vorticity", 0, 1);

            writer.add(uh, "velocity", 0, 2);

            writer.add(ph, "pressure", 0, 1);

            // writer.add(dx_uh0+dy_uh1, "divergence");

            // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");

            writer.add(solp, "pressureExact", 0, 1);

            writer.add(solu, "velocityExact", 0, 2);

            writer.add(soluErr, "velocityError", 0, 2);

            R errU           = L2normCut(uh, fun_exact_u, 0, 3);
            R errP           = L2normCut(ph, fun_exact_p, 0, 1);
            double errDiv    = L2normCut(uh_0dx + uh_1dy + uh_2dz, Khi);
            double maxErrDiv = maxNormCut(uh_0dx + uh_1dy + uh_2dz, Khi);

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

            if (i == 4) {
                assert(0);
            } else {
                nx = 2 * nx - 1;
                ny = 2 * ny - 1;
                nz = 2 * nz - 1;
            }
        }
        std::cout << "\n"
                  << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ')
                  << "err p" << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ')
                  << "err u" << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
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
}

#endif