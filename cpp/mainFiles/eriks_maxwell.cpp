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

// f = 1/eps curl j => div f = 0 !
#define FITTED_WAVE_FORMULATION
// #define UNFITTED_WAVE_FORMULATION
// #define FITTED_KIKUCHI_FORMULATION
// #define UNFITTED_KIKUCHI_FORMULATION
// #define UNFITTED_3FIELD



#ifdef FITTED_WAVE_FORMULATION

    using namespace globalVariable;
    namespace Erik_Data_FITTED_maxwell3D { // f = 1/eps curl j => div f = 0 !
        R k = 1.;
        R eps = 1.;
        R mu = 1.;

        // Monks example
        // R fun_rhs(double *P, int i) {
        //     R x = P[0], y = P[1], z = P[2];
        //     if (i == 0)
        //         return pi*pi*sin(pi*z) - sin(pi*z);
        //     else if (i == 1)
        //         return 0;
        //     else
        //         return 0;
        // }
        // R fun_exact_u(double *P, int i) {
        //     R x = P[0], y = P[1], z = P[2];
        //     R cA = 1./sqrt(7./2);
        //     R ck = 1./(sqrt(13)/2);
        //     if (i == 0)
        //         return cA*cos(ck*(-3./2*x+y));
        //     else if (i == 1)
        //         return cA*3./2*cos(ck*(-3./2*x+y));
        //     else
        //         return cA*1./2*cos(ck*(-3./2*x+y));
        // }

        // Eriks example
        R fun_rhs(double *P, int i) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return pi*pi*sin(pi*z) - sin(pi*z);
            else if (i == 1)
                return 0;
            else
                return 0;
        }
        R fun_exact_u(double *P, int i) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return sin(pi*z);
            else if (i == 1)
                return 0;
            else
                return 0;
        }
        R fun_exact_curlu(double *P, int i) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return 0;
            else if (i == 1)
                return pi*cos(pi*z);
            else
                return 0;
        }

        R fun_0(double *P, int i) {
            return 0;
        }

    } 
    using namespace Erik_Data_FITTED_maxwell3D;
    int main(int argc, char **argv) { // (1.10) : curl(1/eps 1/mu curl(u)) - k^2 u = f, f = 1/eps curl j
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

        int iters = 1;
        for (int i = 0; i < iters; ++i) {

            // Mesh3 Kh(nx, ny, nz, 0., 0., 0., 1., 1., 1.);
            // Mesh3 Kh("../cpp/mainFiles/meshes/square_hole_"+std::to_string((i+1)*20)+".msh");
            // Mesh3 Kh("../cpp/mainFiles/meshes/square_hole_20.msh");
            //Mesh3 Kh("../cpp/mainFiles/meshes/mesh40.mesh");
            Mesh3 Kh("../cpp/mainFiles/meshes/transfinite");
            Kh.info();
            const R hi = 1. / (nx - 1);

            Space Uh(Kh, DataFE<Mesh>::Ned0); // Nedelec order 0 type 1
            // Space Wh_(Kh, DataFE<Mesh>::P1);

            Lagrange3 VelocitySpace(2);
            Space Velh(Kh, VelocitySpace);

            // Interpolate data
            Fun_h fh(Velh, fun_rhs);
            Fun_h u0(Velh, fun_exact_u);
            Fun_h curlu0(Velh, fun_exact_curlu);

            // Init system matrix & assembly
            FEM<Mesh> maxwell3D(Uh);
            // maxwell3D.add(Wh);

            /* Syntax:
            FunTest (fem space, #components, place in space)
            */
            FunTest u(Uh, 3, 0), v(Uh, 3, 0);
            // FunTest p(Wh, 1, 0), q(Wh, 1, 0);
            Normal n;

            // (1.10) : curl(1/eps 1/mu curl(u)) - k^2 u = f, f = 1/eps curl j
            R mui = 1./mu;
            R epsi = 1./eps;
            maxwell3D.addBilinear( 
                +innerProduct(epsi * mui * curl(u), curl(v))
                -innerProduct(k * k * u, v)
            , Kh);
            maxwell3D.addLinear(
                +innerProduct(fh.exprList(), v)
            , Kh);
            // BC
            // Essential
            // R pp = 1e2;
            // maxwell3D.addBilinear( 
            //     -innerProduct(1./mu * curl(u), cross(n, v))
            //     -innerProduct(1./mu * cross(n, u), curl(v))
            //     +innerProduct(cross(n, u), 1./hi*pp*cross(n, v))
            //     // +innerProduct(u, 1./hi*pp*v)
            // , Kh, INTEGRAL_BOUNDARY);
            // maxwell3D.addLinear(
            //     -innerProduct(cross(n, u0), 1./mu * curl(v))
            //     +innerProduct(cross(n, u0), 1./hi*pp*cross(n, v))
            //     // +innerProduct(u0.exprList(), 1./hi*pp*v)
            // , Kh, INTEGRAL_BOUNDARY);
            // Natural
            maxwell3D.addLinear(
                +innerProduct(cross(curlu0, n), 1./mu * v)
            , Kh, INTEGRAL_BOUNDARY);

            matlab::Export(maxwell3D.mat_[0], "mat" + std::to_string(i) + "Cut.dat");
            maxwell3D.solve("mumps");

            // EXTRACT SOLUTION

            int nb_electric_dof = Uh.get_nb_dof();
            // int nb_lagrange_dof = Wh.get_nb_dof();

            Rn_ data_uh = maxwell3D.rhs_(SubArray(nb_electric_dof, 0));
            // Rn_ data_ph = maxwell3D.rhs_(SubArray(nb_lagrange_dof, nb_electric_dof)); // Rn_ data_uh = stokes.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));

            Fun_h uh(Uh, data_uh);
            // Fun_h ph(Wh, data_ph);

            auto uh_0dx = dx(uh.expr(0));
            auto uh_1dy = dy(uh.expr(1));
            auto uh_2dz = dz(uh.expr(2));

            // [Paraview]
            {
                Fun_h solu(Velh, fun_exact_u);

                Fun_h soluErr(Uh, fun_exact_u);
                soluErr.v -= uh.v;
                soluErr.v.map(fabs);

                // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

                Paraview<Mesh> writer(Kh, "maxwell_" + std::to_string(i) + ".vtk");
                writer.add(uh, "electric_field", 0, 3);
                writer.add(uh_0dx + uh_1dy + uh_2dz, "divergence");
                // writer.add(uh_0dx, "divergence");
                writer.add(solu, "electric_field_exact", 0, 3);
                writer.add(soluErr, "electric_field_error", 0, 3);
            }

            R errU      = L2norm(uh, fun_exact_u, 0, 3);
            R errDiv    = L2norm(uh_0dx + uh_1dy + uh_2dz, fun_0, Kh);
            R maxErrDiv = maxNorm(uh_0dx + uh_1dy + uh_2dz, Kh);
            // R maxErrDiv = maxNorm(uh_0dx, Kh);

            h.push_back(hi);
            ul2.push_back(errU);
            divl2.push_back(errDiv);
            divmax.push_back(maxErrDiv);
            if (i == 0) {
                convu.push_back(0);
            } else {
                convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            }
            nx = 2 * nx - 1;
            ny = 2 * ny - 1;
            nz = 2 * nz - 1;
        }
        std::cout << "\n"
        << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ')
        << "err u" << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
        << "err divu"
        // << std::setw(15) << std::setfill(' ') << "conv divu"
        << std::setw(15) << std::setfill(' ')
        << "err maxdivu"
        // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
        << "\n"
        << std::endl;
        for (int i = 0; i < h.size(); ++i) {
            std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
            << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
            << divl2[i]
            // << std::setw(15) << std::setfill(' ') << convdivPr[i]
            << std::setw(15) << std::setfill(' ')
            << divmax[i]
            // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
            << std::endl;
        }
    }
#endif

#ifdef UNFITTED_WAVE_FORMULATION

    using namespace globalVariable;
    namespace Erik_Data_UNFITTED_maxwell3D {
        R k = 1.;
        R eps_r = 1.;
        R mu = 1.;

        R3 shift(0.5, 0.5, 0.5);

        R fun_levelSet(double *P, int i) {
            return (P[0] - shift.x) * (P[0] - shift.x) + (P[1] - shift.y) * (P[1] - shift.y) +
                (P[2] - shift.z) * (P[2] - shift.z) - 0.35 * 0.35 + Epsilon;
        }

        R fun_rhs(double *P, int i) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return 2 * pi * pi * sin(y * pi) * sin(z * pi) -
                    eps_r * sin(y * pi) * sin(z * pi) * k * k; //+ 2 * (P[0] - shift.x) - eps_r * (2 * x - 1);
            else if (i == 1)
                return 2 * pi * pi * sin(x * pi) * sin(z * pi) -
                    eps_r * sin(x * pi) * sin(z * pi) * k * k; //+ 2 * (P[1] - shift.y) - eps_r * (2 * y - 1);
            else
                return 2 * pi * pi * sin(x * pi) * sin(y * pi) -
                    eps_r * sin(x * pi) * sin(y * pi) * k * k; //+ 2 * (P[2] - shift.z) - eps_r * (2 * z - 1);
        }
        R fun_exact_u(double *P, int i, int dom) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return sin(pi * y) * sin(pi * z);
            else if (i == 1)
                return sin(pi * x) * sin(pi * z);
            else
                return sin(pi * x) * sin(pi * y);
        }
        R fun_exact_curlu(double *P, int i, int dom) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return pi*(cos(pi*y) - cos(pi*z))*sin(pi*x);
            else if (i == 1)
                return pi*(-cos(pi*x) + cos(pi*z))*sin(pi*y);
            else
                return pi*(cos(pi*x) - cos(pi*y))*sin(pi*z);
        }
        R fun_exact_p(double *P, int i, int dom) {
            return (P[0] - shift.x) * (P[0] - shift.x) + (P[1] - shift.y) * (P[1] - shift.y) +
                (P[2] - shift.z) * (P[2] - shift.z) - 0.35 * 0.35;
        }

        // Eriks example
        // R fun_rhs(double *P, int i, int dom) {
        //     R x = P[0], y = P[1], z = P[2];
        //     if (i == 0)
        //         return pi*pi*sin(pi*z) - sin(pi*z);
        //     else if (i == 1)
        //         return 0;
        //     else
        //         return 0;
        // }
        // R fun_exact_u(double *P, int i, int dom) {
        //     R x = P[0], y = P[1], z = P[2];
        //     if (i == 0)
        //         return sin(pi*z);
        //     else if (i == 1)
        //         return 0;
        //     else
        //         return 0;
        // }
        // R fun_exact_curlu(double *P, int i, int dom) {
        //     R x = P[0], y = P[1], z = P[2];
        //     if (i == 0)
        //         return 0;
        //     else if (i == 1)
        //         return pi*cos(pi*z);
        //     else
        //         return 0;
        // }

    } 
    using namespace Erik_Data_UNFITTED_maxwell3D;

    int main(int argc, char **argv) {
        typedef TestFunction<Mesh3> FunTest;
        typedef FunFEM<Mesh3> Fun_h;
        typedef Mesh3 Mesh;
        typedef ActiveMeshT3 CutMesh;
        typedef FESpace3 Space;
        typedef CutFESpaceT3 CutSpace;
        const double cpubegin = CPUtime();
        //MPIcf cfMPI(argc, argv);

        const int d = 3;

        int nx = 5;
        int ny = 5;
        int nz = 5;

        std::vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

        int iters = 3;
        for (int i = 0; i < iters; ++i) {

            Mesh3 Kh(nx, ny, nz, 0., 0., 0., 1., 1., 1.);
            const R hi = 1. / (nx - 1);

            Space Uh_(Kh, DataFE<Mesh>::Ned0); // Nedelec order 0 type 1

            Lagrange3 VelocitySpace(2);
            Space Vel_h(Kh, VelocitySpace);

            Space Lh(Kh, DataFE<Mesh>::P1);
            Fun_h levelSet(Lh, fun_levelSet);
            InterfaceLevelSet<Mesh> interface(Kh, levelSet);
            Normal n;

            // [Remove exterior]
            ActiveMesh<Mesh> Khi(Kh);
            Khi.truncate(interface, 1);
            Khi.info();

            CutSpace Velh(Khi, Vel_h);
            CutSpace Uh(Khi, Uh_);
            // CutSpace Wh(Khi, Wh_);

            // Interpolate data
            Fun_h fh(Velh, fun_rhs);
            Fun_h u0(Velh, fun_exact_u);
            Fun_h curlu0(Velh, fun_exact_u);

            // Init system matrix & assembly
            CutFEM<Mesh> maxwell3D(Uh);

            /* Syntax:
            FunTest (fem space, #components, place in space)
            */
            FunTest u(Uh, 3, 0), v(Uh, 3, 0);

            maxwell3D.addBilinear( 
                +innerProduct(1./mu * curl(u), curl(v))
                -innerProduct(k * k * eps_r * u, v)
            , Khi);
            maxwell3D.addLinear(
                +innerProduct(fh.exprList(), v)
            , Khi);
            R pp = 1e2;
            maxwell3D.addBilinear( 
                -innerProduct(1./mu * curl(u), cross(n, v))
                +innerProduct(1./mu * cross(n, u), curl(v))
                +innerProduct(cross(n, u), 1./hi*pp*cross(n, v))
            , interface);
            maxwell3D.addLinear(
                +innerProduct(cross(n, u0), 1./mu * curl(v))
                +innerProduct(cross(n, u0), 1./hi*pp*cross(n, v))
            , interface);
            // maxwell3D.addLinear(
            //     -innerProduct(cross(n, curlu0), 1./mu * v)
            // , interface);
            // // maxwell3D.addLinear(
            // //     +innerProduct(curlu0.exprList(), 1./mu * cross(n, v))
            // // , interface);

            // [Stabilization]
            double tau_a = 1e0;
            maxwell3D.addPatchStabilization(
            // W block
            - innerProduct(tau_a * pow(hi, 0) * jump(u), jump(v)) 
            , Khi);
            // maxwell3D.addFaceStabilization(
            //     + innerProduct(tau_a * pow(hi, 1) * jump(u), jump(v)) 
            //     + innerProduct(tau_a * pow(hi, 3) * jump(grad(u)*n), jump(grad(v)*n)) 
            // , Khi);

            matlab::Export(maxwell3D.mat_[0], "mat" + std::to_string(i) + "Cut.dat");
            maxwell3D.solve("umfpack");

            // EXTRACT SOLUTION

            int nb_electric_dof = Uh.get_nb_dof();
            Rn_ data_uh = maxwell3D.rhs_(SubArray(nb_electric_dof, 0));
            Fun_h uh(Uh, data_uh);

            auto uh_0dx = dx(uh.expr(0));
            auto uh_1dy = dy(uh.expr(1));
            auto uh_2dz = dz(uh.expr(2));

            // [Paraview]
            {
                Fun_h solu(Velh, fun_exact_u);
                Fun_h soluErr(Uh, fun_exact_u);

                soluErr.v -= uh.v;
                soluErr.v.map(fabs);

                // Fun_h divSolh(Wh, fun_div);
                // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

                Paraview<Mesh> writer(Khi, "maxwell_" + std::to_string(i) + ".vtk");
                writer.add(uh, "electric_field", 0, 3);
                // writer.add(ph, "pressure", 0, 1);
                // writer.add(dx_uh0+dy_uh1, "divergence");
                // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
                // writer.add(solp, "pressureExact", 0, 1);
                writer.add(solu, "electric_field_exact", 0, 3);
                writer.add(soluErr, "electric_field_error", 0, 3);
            }
            R errU           = L2normCut(uh, fun_exact_u, 0, 3);
            double errDiv    = L2normCut(uh_0dx + uh_1dy + uh_2dz, Khi);
            double maxErrDiv = maxNormCut(uh_0dx + uh_1dy + uh_2dz, Khi);

            ul2.push_back(errU);
            divl2.push_back(errDiv);
            divmax.push_back(maxErrDiv);
            h.push_back(hi);
            if (i == 0) {
                convu.push_back(0);
            } else {
                convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            }
            nx = 2 * nx - 1;
            ny = 2 * ny - 1;
            nz = 2 * nz - 1;
        }
        std::cout << "\n"
        << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ')
        << "err u" << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
        << "err divu"
        // << std::setw(15) << std::setfill(' ') << "conv divu"
        << std::setw(15) << std::setfill(' ')
        << "err maxdivu"
        // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
        << "\n"
        << std::endl;
        for (int i = 0; i < h.size(); ++i) {
            std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
            << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
            << divl2[i]
            // << std::setw(15) << std::setfill(' ') << convdivPr[i]
            << std::setw(15) << std::setfill(' ')
            << divmax[i]
            // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
            << std::endl;
        }
    }

#endif

#ifdef FITTED_KIKUCHI_FORMULATION

    using namespace globalVariable;
    namespace Erik_Data_FITTED_maxwell3D {
        R k = 1.;
        R eps_r = 1.;
        R mu = 1.;

        // Monks example
        // R fun_rhs(double *P, int i) {
        //     R x = P[0], y = P[1], z = P[2];
        //     if (i == 0)
        //         return pi*pi*sin(pi*z) - sin(pi*z);
        //     else if (i == 1)
        //         return 0;
        //     else
        //         return 0;
        // }
        // R fun_exact_u(double *P, int i) {
        //     R x = P[0], y = P[1], z = P[2];
        //     R cA = 1./sqrt(7./2);
        //     R ck = 1./(sqrt(13)/2);
        //     if (i == 0)
        //         return cA*cos(ck*(-3./2*x+y));
        //     else if (i == 1)
        //         return cA*3./2*cos(ck*(-3./2*x+y));
        //     else
        //         return cA*1./2*cos(ck*(-3./2*x+y));
        // }

        // Eriks example
        R fun_rhs(double *P, int i) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return pi*pi*sin(pi*z) - sin(pi*z) - (2*x-1);
            else if (i == 1)
                return 2*y-1;
            else
                return 0;
        }
        R fun_exact_u(double *P, int i) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return sin(pi*z);
            else if (i == 1)
                return 0;
            else
                return 0;
        }
        R fun_exact_curlu(double *P, int i) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return 0;
            else if (i == 1)
                return pi*cos(pi*z);
            else
                return 0;
        }
        R fun_exact_p(double *P, int i) {
            R x = P[0], y = P[1], z = P[2];
            return x*(x-1)-y*(y-1);
        }

        R fun_0(double *P, int i) {
            return 0;
        }

    } 
    using namespace Erik_Data_FITTED_maxwell3D;
    int main(int argc, char **argv) {
        typedef TestFunction<Mesh3> FunTest;
        typedef FunFEM<Mesh3> Fun_h;
        typedef Mesh3 Mesh;
        typedef ActiveMeshT3 CutMesh;
        typedef FESpace3 Space;
        typedef CutFESpaceT3 CutSpace;
        const double cpubegin = CPUtime();
        //MPIcf cfMPI(argc, argv);

        const int d = 3;

        int nx = 5;
        int ny = 5;
        int nz = 5;

        std::vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

        int iters = 3;
        for (int i = 0; i < iters; ++i) {

            Mesh3 Kh(nx, ny, nz, 0., 0., 0., 1., 1., 1.);
            Kh.info();
            const R hi = 1. / (nx - 1);

            Space Uh(Kh, DataFE<Mesh>::Ned0); // Nedelec order 0 type 1
            Space Wh(Kh, DataFE<Mesh>::P1);

            Lagrange3 VelocitySpace(2);
            Space Velh(Kh, VelocitySpace);

            // Interpolate data
            Fun_h fh(Velh, fun_rhs);
            Fun_h u0(Velh, fun_exact_u);
            Fun_h curlu0(Velh, fun_exact_curlu);

            // Init system matrix & assembly
            FEM<Mesh> maxwell3D(Uh);
            maxwell3D.add(Wh);

            /* Syntax:
            FunTest (fem space, #components, place in space)
            */
            FunTest u(Uh, 3, 0), v(Uh, 3, 0);
            FunTest p(Wh, 1, 0), q(Wh, 1, 0);
            Normal n;

            R pp = 1e2;
            maxwell3D.addBilinear( 
                +innerProduct(1./mu * curl(u), curl(v))
                -innerProduct(k * k * eps_r * u, v)

                -innerProduct(grad(p), v)
                -innerProduct(u, grad(q))
            , Kh);
            maxwell3D.addLinear(
                +innerProduct(fh.exprList(), v)
            , Kh);
            // BC
            // Essential
            maxwell3D.addBilinear( 
                -innerProduct(1./mu * curl(u), cross(n, v))
                -innerProduct(1./mu * cross(n, u), curl(v))
                +innerProduct(cross(n, u), 1./hi*pp*cross(n, v))
                // +innerProduct(u, 1./hi*pp*v)
            , Kh, INTEGRAL_BOUNDARY);
            maxwell3D.addLinear(
                -innerProduct(cross(n, u0), 1./mu * curl(v))
                +innerProduct(cross(n, u0), 1./hi*pp*cross(n, v))
                // +innerProduct(u0.exprList(), 1./hi*pp*v)
            , Kh, INTEGRAL_BOUNDARY);
            // Natural
            // maxwell3D.addLinear(
            //     -innerProduct(cross(n, curlu0), 1./mu * v)
            // , Kh, INTEGRAL_BOUNDARY);

            matlab::Export(maxwell3D.mat_[0], "mat" + std::to_string(i) + "Cut.dat");
            maxwell3D.solve("umfpack");

            // EXTRACT SOLUTION

            int nb_electric_dof = Uh.get_nb_dof();
            int nb_lagrange_dof = Wh.get_nb_dof();

            Rn_ data_uh = maxwell3D.rhs_(SubArray(nb_electric_dof, 0));
            Rn_ data_ph = maxwell3D.rhs_(SubArray(nb_lagrange_dof, nb_electric_dof)); // Rn_ data_uh = stokes.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));

            Fun_h uh(Uh, data_uh);
            Fun_h ph(Wh, data_ph);

            auto uh_0dx = dx(uh.expr(0));
            auto uh_1dy = dy(uh.expr(1));
            auto uh_2dz = dz(uh.expr(2));

            // [Paraview]
            {
                Fun_h solu(Velh, fun_exact_u);

                Fun_h soluErr(Uh, fun_exact_u);
                soluErr.v -= uh.v;
                soluErr.v.map(fabs);

                // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

                Paraview<Mesh> writer(Kh, "maxwell_" + std::to_string(i) + ".vtk");
                writer.add(uh, "electric_field", 0, 3);
                writer.add(uh_0dx + uh_1dy + uh_2dz, "divergence");
                // writer.add(uh_0dx, "divergence");
                writer.add(solu, "electric_field_exact", 0, 3);
                writer.add(soluErr, "electric_field_error", 0, 3);
            }

            R errP      = L2norm(ph, fun_exact_p, 0, 1);
            R errU      = L2norm(uh, fun_exact_u, 0, 3);
            R errDiv    = L2norm(uh_0dx + uh_1dy + uh_2dz, fun_0, Kh);
            R maxErrDiv = maxNorm(uh_0dx + uh_1dy + uh_2dz, Kh);
            // R maxErrDiv = maxNorm(uh_0dx, Kh);

            h.push_back(hi);
            pl2.push_back(errP);
            ul2.push_back(errU);
            divl2.push_back(errDiv);
            divmax.push_back(maxErrDiv);
            if (i == 0) {
                convp.push_back(0);
                convu.push_back(0);
            } else {
                convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
                convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            }
            nx = 2 * nx - 1;
            ny = 2 * ny - 1;
            nz = 2 * nz - 1;
        }
        std::cout << "\n"
        << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ')
        << "err p" << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ')
        << "err u" << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
        << "err divu"
        // << std::setw(15) << std::setfill(' ') << "conv divu"
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
            << std::setw(15) << std::setfill(' ')
            << divmax[i]
            // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
            << std::endl;
        }
    }
#endif

#ifdef UNFITTED_KIKUCHI_FORMULATION

    using namespace globalVariable;
    namespace Erik_Data_UNFITTED_maxwell3D {
        R k = 1.;
        R eps_r = 1.;
        R mu = 1.;

        R3 shift(0.5, 0.5, 0.5);

        R fun_levelSet(double *P, int i) {
            return (P[0] - shift.x) * (P[0] - shift.x) + (P[1] - shift.y) * (P[1] - shift.y) +
                (P[2] - shift.z) * (P[2] - shift.z) - 0.35 * 0.35 + Epsilon;
        }

        // R fun_rhs(double *P, int i) {
        //     R x = P[0], y = P[1], z = P[2];
        //     if (i == 0)
        //         return 2 * pi * pi * sin(y * pi) * sin(z * pi) -
        //             eps_r * sin(y * pi) * sin(z * pi) * k * k; //+ 2 * (P[0] - shift.x) - eps_r * (2 * x - 1);
        //     else if (i == 1)
        //         return 2 * pi * pi * sin(x * pi) * sin(z * pi) -
        //             eps_r * sin(x * pi) * sin(z * pi) * k * k; //+ 2 * (P[1] - shift.y) - eps_r * (2 * y - 1);
        //     else
        //         return 2 * pi * pi * sin(x * pi) * sin(y * pi) -
        //             eps_r * sin(x * pi) * sin(y * pi) * k * k; //+ 2 * (P[2] - shift.z) - eps_r * (2 * z - 1);
        // }
        // R fun_exact_u(double *P, int i, int dom) {
        //     R x = P[0], y = P[1], z = P[2];
        //     if (i == 0)
        //         return sin(pi * y) * sin(pi * z);
        //     else if (i == 1)
        //         return sin(pi * x) * sin(pi * z);
        //     else
        //         return sin(pi * x) * sin(pi * y);
        // }
        // R fun_exact_curlu(double *P, int i, int dom) {
        //     R x = P[0], y = P[1], z = P[2];
        //     if (i == 0)
        //         return pi*(cos(pi*y) - cos(pi*z))*sin(pi*x);
        //     else if (i == 1)
        //         return pi*(-cos(pi*x) + cos(pi*z))*sin(pi*y);
        //     else
        //         return pi*(cos(pi*x) - cos(pi*y))*sin(pi*z);
        // }
        // R fun_exact_p(double *P, int i, int dom) {
        //     return (P[0] - shift.x) * (P[0] - shift.x) + (P[1] - shift.y) * (P[1] - shift.y) +
        //         (P[2] - shift.z) * (P[2] - shift.z) - 0.35 * 0.35;
        // }

        // Eriks example
        R fun_rhs(double *P, int i, int dom) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return pi*pi*sin(pi*z) - sin(pi*z);
            else if (i == 1)
                return 0;
            else
                return 0;
        }
        R fun_exact_u(double *P, int i, int dom) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return sin(pi*z);
            else if (i == 1)
                return 0;
            else
                return 0;
        }
        R fun_exact_curlu(double *P, int i, int dom) {
            R x = P[0], y = P[1], z = P[2];
            if (i == 0)
                return 0;
            else if (i == 1)
                return pi*cos(pi*z);
            else
                return 0;
        }

    } 
    using namespace Erik_Data_UNFITTED_maxwell3D;

    int main(int argc, char **argv) {
        typedef TestFunction<Mesh3> FunTest;
        typedef FunFEM<Mesh3> Fun_h;
        typedef Mesh3 Mesh;
        typedef ActiveMeshT3 CutMesh;
        typedef FESpace3 Space;
        typedef CutFESpaceT3 CutSpace;
        const double cpubegin = CPUtime();
        //MPIcf cfMPI(argc, argv);

        const int d = 3;

        int nx = 5;
        int ny = 5;
        int nz = 5;

        std::vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

        int iters = 3;
        for (int i = 0; i < iters; ++i) {

            Mesh3 Kh(nx, ny, nz, 0., 0., 0., 1., 1., 1.);
            const R hi = 1. / (nx - 1);

            Space Uh_(Kh, DataFE<Mesh>::Ned0); // Nedelec order 0 type 1
            // Space Wh_(Kh, DataFE<Mesh>::P1);

            Lagrange3 VelocitySpace(2);
            Space Vel_h(Kh, VelocitySpace);

            Space Lh(Kh, DataFE<Mesh>::P1);
            Fun_h levelSet(Lh, fun_levelSet);
            InterfaceLevelSet<Mesh> interface(Kh, levelSet);
            Normal n;

            // [Remove exterior]
            ActiveMesh<Mesh> Khi(Kh);
            Khi.truncate(interface, 1);
            Khi.info();

            CutSpace Velh(Khi, Vel_h);
            CutSpace Uh(Khi, Uh_);
            // CutSpace Wh(Khi, Wh_);

            // Interpolate data
            Fun_h fh(Velh, fun_rhs);
            Fun_h u0(Velh, fun_exact_u);
            Fun_h curlu0(Velh, fun_exact_u);

            // Init system matrix & assembly
            CutFEM<Mesh> maxwell3D(Uh);
            // maxwell3D.add(Wh);

            /* Syntax:
            FunTest (fem space, #components, place in space)
            */
            FunTest u(Uh, 3, 0), v(Uh, 3, 0);
            // FunTest p(Wh, 1, 0), q(Wh, 1, 0);

            // Eq 1
            maxwell3D.addBilinear( 
                +innerProduct(1./mu * curl(u), curl(v))
                -innerProduct(k * k * eps_r * u, v)
            // +innerProduct(p, div(v))
            // -innerProduct(grad(p), v)
            , Khi);
            maxwell3D.addLinear(
                +innerProduct(fh.exprList(), v)
            , Khi);
            R pp = 1e2;
            // maxwell3D.addBilinear( 
            //     -innerProduct(1./mu * curl(u), cross(n, v))
            //     -innerProduct(1./mu * cross(n, u), curl(v))
            //     +innerProduct(cross(n, u), 1./hi*pp*cross(n, v))
            // , interface);
            // maxwell3D.addLinear(
            //     -innerProduct(cross(n, u0), 1./mu * curl(v))
            //     +innerProduct(cross(n, u0), 1./hi*pp*cross(n, v))
            // , interface);
            // maxwell3D.addLinear(
            //     -innerProduct(cross(n, curlu0), 1./mu * v)
            // , interface);
            maxwell3D.addLinear(
                +innerProduct(curlu0.exprList(), 1./mu * cross(n, v))
            , interface);

            // Eq 2
            // maxwell3D.addBilinear(
            //     // -innerProduct(div(u), q)
            //     +innerProduct(u, grad(q))
            //     , Khi);

            // [Stabilization]
            double tau_a = 1e0;
            double tau_b = 1e0;
            maxwell3D.addPatchStabilization(
            // W block
            + innerProduct(tau_a * pow(hi, 0) * jump(u), jump(v)) 
            // B blocks
            // + innerProduct(tau_b * pow(hi, 0) * jump(p), jump(div(v)))                     // -B^T block
            // - innerProduct(tau_b * pow(hi, 0) * jump(div(u)), jump(q))                     // B_0 block
            // - innerProduct(tau_b * pow(hi, 0) * jump(grad(p)), jump(v))                     // -B^T block
            // + innerProduct(tau_b * pow(hi, 0) * jump(u), jump(grad(q)))                     // B_0 block
            , Khi);
            // maxwell3D.addFaceStabilization(
            //     + innerProduct(tau_a * pow(hi, 1) * jump(u), jump(v)) 
            //     + innerProduct(tau_a * pow(hi, 3) * jump(grad(u)*n), jump(grad(v)*n)) 
            // , Khi);

            matlab::Export(maxwell3D.mat_[0], "mat" + std::to_string(i) + "Cut.dat");
            maxwell3D.solve("umfpack");

            // EXTRACT SOLUTION

            int nb_electric_dof = Uh.get_nb_dof();
            // int nb_lagrange_dof = Wh.get_nb_dof();

            Rn_ data_uh = maxwell3D.rhs_(SubArray(nb_electric_dof, 0));
            // Rn_ data_ph = maxwell3D.rhs_(SubArray(nb_lagrange_dof, nb_electric_dof)); // Rn_ data_uh = stokes.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));

            Fun_h uh(Uh, data_uh);
            // Fun_h ph(Wh, data_ph);

            auto uh_0dx = dx(uh.expr(0));
            auto uh_1dy = dy(uh.expr(1));
            auto uh_2dz = dz(uh.expr(2));

            // [Paraview]
            {

                // Fun_h solw(Uh, fun_exact_w);

                Fun_h solu(Velh, fun_exact_u);
                Fun_h soluErr(Uh, fun_exact_u);
                // Fun_h solp(Wh, fun_exact_p);

                soluErr.v -= uh.v;
                soluErr.v.map(fabs);

                // Fun_h divSolh(Wh, fun_div);
                // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

                Paraview<Mesh> writer(Khi, "maxwell_" + std::to_string(i) + ".vtk");

                writer.add(uh, "electric_field", 0, 3);
                // writer.add(ph, "pressure", 0, 1);
                // writer.add(dx_uh0+dy_uh1, "divergence");
                // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
                // writer.add(solp, "pressureExact", 0, 1);
                writer.add(solu, "electric_field_exact", 0, 3);
                writer.add(soluErr, "electric_field_error", 0, 3);

                R errU           = L2normCut(uh, fun_exact_u, 0, 3);
                R errP           = 0;
                // R errP           = L2normCut(ph, fun_exact_p, 0, 1);
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
                // nx = (int)(std::sqrt(2) * nx - 1);
                // ny = (int)(std::sqrt(2) * ny - 1);
                // nz = (int)(std::sqrt(2) * nz - 1);
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
            << std::setw(15) << std::setfill(' ')
            << divmax[i]
            // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
            << std::endl;
        }
    }

#endif

#ifdef UNFITTED_3FIELD

    using namespace globalVariable;

    namespace Erik_Data_UNFITTED_maxwell3D {
    R k = 1.;

    R eps_r = 1.;

    R3 shift(0.5, 0.5, 0.5);

    R fun_levelSet(double *P, int i) {
        return (P[0] - shift.x) * (P[0] - shift.x) + (P[1] - shift.y) * (P[1] - shift.y) +
            (P[2] - shift.z) * (P[2] - shift.z) - 0.35 * 0.35 + Epsilon;
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

        //MPIcf cfMPI(argc, argv);

        const int d = 3;

        int nx = 5;
        int ny = 5;
        int nz = 5;

        std::vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

        int iters = 4;
        for (int i = 0; i < iters; ++i) {
            Mesh3 Kh(nx, ny, nz, 0., 0., 0., 1., 1., 1.);
            const R hi = 1. / (nx - 1); // 1./(nx-1)

            Space Uh_(Kh,
                    DataFE<Mesh>::Ned0); // Nedelec order 0 type 1
            Space Vh_(Kh, DataFE<Mesh>::RT0);
            Space Wh_(Kh, DataFE<Mesh>::P0);

            Lagrange3 VelocitySpace(2);
            Space Vel_h(Kh, VelocitySpace);

            Space Lh(Kh, DataFE<Mesh>::P1);
            Fun_h levelSet(Lh, fun_levelSet);
            InterfaceLevelSet<Mesh> interface(Kh, levelSet);
            Normal n;

            // [Remove exterior]
            ActiveMesh<Mesh> Khi(Kh);
            Khi.truncate(interface, 1);

            CutSpace Velh(Khi, Vel_h);
            CutSpace Uh(Khi, Uh_);
            CutSpace Vh(Khi, Vh_);
            CutSpace Wh(Khi, Wh_);

            // Interpolate data
            Fun_h fh(Velh, fun_rhs);
            Fun_h u0(Velh, fun_exact_u);

            // Init system matrix & assembly
            CutFEM<Mesh> maxwell3D(Uh);
            maxwell3D.add(Vh);
            maxwell3D.add(Wh);

            /* Syntax:
            FunTest (fem space, #components, place in space)
            */

            FunTest w(Uh, 3, 0), tau(Uh, 3, 0);
            FunTest u(Vh, 3, 0), v(Vh, 3, 0), p(Wh, 1, 0), q(Wh, 1, 0);
            R mu = 1.;

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

            // [Stabilization]
            // order 1 with mumps: 1e-1, 2e0, 1e-1, 15
            double tau_w = 1e0;       // smaller tau_w seems to give larger condition number
            double tau_m = 1e0;
            double tau_a = 1e0;
            double tau_b = 1e0;

            // maxwell3D.addFaceStabilization(
            //     // WA

            //     // // W block
            //     // + innerProduct(tau_w * pow(hi, 1) * jump(w), jump(tau))
            //     // + innerProduct(tau_w * pow(hi, 3) * jump(grad(w)*n), jump(grad(tau)*n))
            //     // //+ innerProduct(tau_w * pow(hi, 5) * jump(grad(grad(w) * n)), jump(grad(grad(tau) * n)))

            //     // // A block
            //     // + innerProduct(tau_a * pow(hi, 1) * jump(u), jump(v))
            //     // + innerProduct(tau_a * pow(hi, 3) * jump(grad(u)*n), jump(grad(v)*n))
            //     // //+ innerProduct(tau_a * pow(hi, 5) * jump(grad(grad(u) * n)), jump(grad(grad(v) * n)))

            //     // // B blocks
            //     // + innerProduct(tau_b * pow(hi, 1) * jump(p), jump(div(v)))              // -B^T block
            //     // + innerProduct(tau_b * pow(hi, 3) * jump(grad(p)*n), jump(grad(div(v))*n))  // -B^T block
            //     // - innerProduct(tau_b * pow(hi, 1) * jump(div(u)), jump(q))              // B_0 block
            //     // - innerProduct(tau_b * pow(hi, 3) * jump(grad(div(u))*n), jump(grad(q)*n)) // B_0 block

            //     // WM

            //     // W block
            //     + innerProduct(tau_w * pow(hi, 1) * jump(w), jump(tau)) 
            //     + innerProduct(tau_w * pow(hi, 3) * jump(grad(w) * n), jump(grad(tau) * n))
            //     // M blocks
            //     + innerProduct(tau_m * pow(hi, 1) * jump(curl(w)), jump(v))                       // M block
            //     + innerProduct(tau_m * pow(hi, 3) * jump(grad(curl(w)) * n), jump(grad(v) * n))   // M block
            //     - innerProduct(tau_m * pow(hi, 1) * jump(u), jump(curl(tau)))                     // -M^T block
            //     - innerProduct(tau_m * pow(hi, 3) * jump(grad(u) * n), jump(grad(curl(tau)) * n)) // -M^T block
            //     // B blocks
            //     + innerProduct(tau_b * pow(hi, 1) * jump(p), jump(div(v)))                     // -B^T block
            //     + innerProduct(tau_b * pow(hi, 3) * jump(grad(p) * n), jump(grad(div(v)) * n)) // -B^T block
            //     - innerProduct(tau_b * pow(hi, 1) * jump(div(u)), jump(q))                     // B_0 block
            //     - innerProduct(tau_b * pow(hi, 3) * jump(grad(div(u)) * n), jump(grad(q) * n)) // B_0 block
            //     , Khi);

            maxwell3D.addPatchStabilization(
                // W block
                + innerProduct(tau_w  * jump(w), jump(tau)) 
                // M blocks
                + innerProduct(tau_m * jump(curl(w)), jump(v))                       // M block
                - innerProduct(tau_m * jump(u), jump(curl(tau)))                     // -M^T block
                // B blocks
                + innerProduct(tau_b * jump(p), jump(div(v)))                     // -B^T block
                - innerProduct(tau_b * jump(div(u)), jump(q))                     // B_0 block
                , Khi);//                                  , macro);

            matlab::Export(maxwell3D.mat_[0], "mat" + std::to_string(i) + ".dat");

            maxwell3D.solve("umfpack");

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

                writer.add(wh, "vorticity", 0, 3);

                writer.add(uh, "velocity", 0, 3);

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

                if (i == 10) {
                    assert(0);
                } else {
                    // nx = (int)(std::sqrt(2) * nx - 1);
                    // ny = (int)(std::sqrt(2) * ny - 1);
                    // nz = (int)(std::sqrt(2) * nz - 1);
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
