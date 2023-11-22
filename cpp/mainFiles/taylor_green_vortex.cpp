/**
 * @file navier_stokes_raising_drop_2D.cpp
 * @author Thomas Frachon
 * @brief Raising bubble using CutFEM. For details, one can see Example 2 in
 * "A cut finite element method for incompressible two-phase Navier-Stokes flows" by  Thomas Frachon and Sara Zahedi
 * @version 0.1
 * @date 2023-03-09
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "../tool.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;

const double nu  = 0.1;
const double rho = 1.;

double fun_rhs(R2 P, int i, const double t) { return 0.; }

double fun_u(R2 P, int i, const double t) {
    // double x = P[0], y = P[1];
    double x = P.x, y = P.y;
    if (i == 0)
        return std::sin(x) * std::cos(y) * std::exp(-2 * nu * t);
    else
        return (-1.) * std::cos(x) * std::sin(y) * std::exp(-2 * nu * t);
}

double fun_u_initial(R2 P, int i) {
    double x = P.x, y = P.y;
    if (i == 0)
        return std::sin(x) * std::cos(y);
    else
        return (-1.) * std::cos(x) * std::sin(y);
}

double fun_p(R2 P, int i, const double t) { return rho / 4 * (std::cos(2 * P.x) + std::cos(2 * P.y)); }

int main(int argc, char **argv) {

    // MPIcf cfMPI(argc, argv);
    const int thread_count = 1;
    Logger::initialize("log_navier_Stokes.txt");

    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/navier_stokes/taylor_green/data/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/navier_stokes/taylor_green/paraview/";

    double h = 0.2; // starting mesh size

    const size_t iterations = 1;
    for (int j = 0; j < iterations; ++j) {

        // Mesh
        const double lx = M_PI, ly = M_PI;
        const double x0 = 0., y0 = 0.;
        int nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
        Mesh2 Th(nx, ny, x0, y0, lx, ly);

        std::list<int> dirichlet{1, 2, 3, 4};

        // Th.info();

        // Time stepping
        int division_mesh_size     = 2;
        double final_time          = 0.5;
        double dT                  = h / division_mesh_size;
        int total_number_iteration = final_time / dT;
        dT                         = final_time / total_number_iteration;
        double time_step           = dT;
        double t0                  = 0.;

        Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);

        // Define the time quadrature formula
        const QuadratureFormular1d &qTime(*Lobatto(3));
        const Uint nb_quad_time   = qTime.n;
        const Uint ndf_time_slab  = Ih[0].NbDoF();
        const Uint last_quad_time = nb_quad_time - 1;

        ProblemOption optionProblem;
        optionProblem.solver_name_  = "umfpack";
        optionProblem.clear_matrix_ = true;
        std::vector<std::map<std::pair<int, int>, double>> mat_NL(thread_count);

        CutFEM<mesh_t> navier_stokes(qTime, thread_count, optionProblem);

        const double lambda_boundary = 10.;
        const double lambda_interior = 10.;

        // Finite element spaces
        Lagrange2 FEu(2); // for interpolating the exact solution
        space_t Vh_interpolation(Th, FEu);

        space_t Vh(Th, DataFE<mesh_t>::BDM1);
        space_t Ph(Th, DataFE<mesh_t>::P0);

        // std::vector<fct_t> vel(nb_quad_time);
        // for (int i = 0; i < nb_quad_time; ++i)
        //     vel[i].init(Vh, fun_boundary);

        if (iterations > 1) {
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "Iteration " << j + 1 << "/" << iterations << '\n';
        }

        std::cout << "h  = " << h << '\n';
        std::cout << "nx = " << nx << '\n';
        std::cout << "ny = " << ny << '\n';
        std::cout << "dT = " << dT << '\n';

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << '\n';

        int iter = 0;
        while (iter < total_number_iteration) {
            int current_iteration = iter;
            const TimeSlab &In(Ih[iter]);

            std::cout << " -------------------------------------------------------\n";
            std::cout << " -------------------------------------------------------\n";
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << '\n';
            std::cout << " Time      \t : \t" << current_iteration * time_step << '\n';
            std::cout << "dT = " << dT << '\n';

            // DEFINE TEST FUNCTIONS
            // ----------------------------------------------
            Normal n;
            Tangent t;
            // funtest_t u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);

            // uh^{k+1} = uh^k + du
            // ph^{k+1} = ph^k + dp
            funtest_t du(Vh, 2), dp(Ph, 1), v(Vh, 2), q(Ph, 1), dp1(Vh, 1, 0, 0);

            // DEFINE THE PROBLEM ON THE CORRECT SPACES
            // ----------------------------------------------
            navier_stokes.initSpace(Vh, In);
            navier_stokes.add(Ph, In);

            fct_t u_exact(Vh, In, fun_u);
            fct_t p_exact(Ph, In, fun_p);
            fct_t fh(Vh, In, fun_rhs); // rhs force
            fct_t gh(Vh, In, fun_u);   // Dirichlet boundary condition

            // std::cout << " Problem's DOF : \t" << navier_stokes.get_nb_dof() << std::endl;

            Rn data_init(navier_stokes.get_nb_dof(), 0.);
            // std::vector<double> data_init(navier_stokes.get_nb_dof());
            KN_<double> data_init_span(data_init);
            // std::span<double> data_init_span(data_init);
            navier_stokes.initialSolution(data_init_span);
            Rn data_all(data_init); //! Check this one!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // std::vector<double> data_all(data_init);

            int idxp0 = Vh.NbDoF() * In.NbDoF(); // index for when p starts in the array

            KN_<double> data_uh0(data_init(SubArray(Vh.NbDoF(), 0)));
            // std::span<double> data_uh0 = std::span<double>(data_init.data(), Vh.NbDoF()); // velocity in first
            // time DOF
            KN_<double> data_uh(data_init(SubArray(Vh.NbDoF() * In.NbDoF(), 0)));
            // std::span<double> data_uh =
            //     std::span<double>(data_all.data(), Vh.NbDoF() * In.NbDoF()); // velocity for all time DOFs

            KN_<double> data_ph(data_init(SubArray(Ph.NbDoF() * In.NbDoF(), idxp0)));
            // std::span<double> data_ph = std::span<double>(data_all.data() + idxp0, Ph.NbDoF() * In.NbDoF());

            fct_t u0(Vh, data_uh0);
            fct_t uh(Vh, In, data_uh);
            fct_t ph(Ph, In, data_ph);

            if (iter == 0) {
                interpolate(Vh, data_uh0, fun_u_initial);
            }

            // NEWTON ITERATION
            // ----------------------------------------------
            int newton_iterations = 0;
            // while (1) {
            //     double tt0 = MPIcf::Wtime();

            //     // Add terms that are in the residual
            //     if (newton_iterations == 0) {

            //         // time terms
            //         navier_stokes.addBilinear(innerProduct(dt(du), rho * v), Th, In);
            //         navier_stokes.addBilinear(innerProduct(rho * du, v), Th, 0, In);

            //         // a(t, du, vh) in In x Omega(t)
            //         navier_stokes.addBilinear(contractProduct(nu * grad(du), grad(v)), Th, In);
            //         navier_stokes.addBilinear(-innerProduct(nu * grad(du) * n, v) - innerProduct(du, nu * grad(v) *
            //         n) +
            //                                       innerProduct(lambda_boundary * du, v),
            //                                   Th, INTEGRAL_BOUNDARY, In);
            //         navier_stokes.addBilinear(-innerProduct(average(grad(du * t) * n, 0.5, 0.5), jump(v * t)) +
            //                                       innerProduct(jump(du * t), average(grad(v * t) * n, 0.5, 0.5)) +
            //                                       innerProduct(lambda_interior / h * (jump(du * t)), jump(v * t)),
            //                                   Th, INTEGRAL_INNER_EDGE_2D, In);

            //         // -b + b0
            //         navier_stokes.addBilinear(-innerProduct(dp, div(v)) + innerProduct(div(du), q), Th, In);
            //         navier_stokes.addBilinear(innerProduct(dp, v * n), Th, INTEGRAL_BOUNDARY, In);
            //     }

            //     //         for (int i = 0; i < nb_quad_time; ++i) { // computation of the curvature

            //     //             cutmesh_t cutTh(Kh);
            //     //             cutTh.createSurfaceMesh(*interface[i]);
            //     //             cutspace_t cutVh(cutTh, Vh1);
            //     //             Curvature<mesh_t> curvature(cutVh, interface[i]);

            //     //             auto tq    = qTime(i);
            //     //             double tid = (double)In.map(tq);

            //     //             auto data_H = curvature.solve();
            //     //             fct_t H(cutVh, data_H);
            //     //             stokes.addLinear(-sigma * innerProduct(H.exprList(), average(v, kappa2)), interface,
            //     In,
            //     //             i);

            //     //             // if (i == 0) {
            //     //             //     Paraview<mesh_t> writerS(cutTh, "raingDropExampleCurvature_" +
            //     //             std::to_string(ifig) + ".vtk");
            //     //             //     writerS.add(ls[i], "levelSet", 0, 1);
            //     //             //     writerS.add(H, "meanCurvature", 0, 2);
            //     //             // }
            //     //         }

            //     //         // if (iter == 0 && iterNewton == 0)
            //     //         //     globalVariable::verbose = 1;

            //     //         // if (iter == 0 && iterNewton == 0) {
            //     //         LOG_INFO << " Time assembly matrix A : \t" << MPIcf::Wtime() - tt0 << logger::endl;
            //     //         // }

            //     //
            //     // Add negative linear form
            //     navier_stokes.addLinear(-innerProduct(fh.exprList(), rho * v), Th, In);
            //     navier_stokes.addLinear(-innerProduct(u0.exprList(), rho * v), Th, 0, In);
            //     navier_stokes.addLinear(innerProduct(gh, nu * grad(v) * n) -
            //                                 innerProduct(gh.exprList(), lambda_boundary * v),
            //                             Th, INTEGRAL_BOUNDARY, In);

            //     navier_stokes.gather_map();
            //     navier_stokes.addMatMul(data_all); // multiply matrix with data (uh^k) and add result to rhs to
            //     assemble residual

            //     // Construct remaining terms of the Jacobian

            //     // Assemble Jacobian in the matrix mat_NL
            //     navier_stokes.set_map(mat_NL);
            //     mat_NL[0] = navier_stokes.mat_[0];

            //     funtest_t du1(Vh, 1, 0), du2(Vh, 1, 1), v1(Vh, 1, 0), v2(Vh, 1, 1);
            //     auto ux = uh.expr(0);
            //     auto uy = uh.expr(1);

            //     stokes.addBilinear(innerProduct(du1 * dx(ux) + du2 * dy(ux), rho * v1) +
            //                            innerProduct(du1 * dx(uy) + du2 * dy(uy), rho * v2) +
            //                            innerProduct(ux * dx(du1) + uy * dy(du1), rho * v1) +
            //                            innerProduct(ux * dx(du2) + uy * dy(du2), rho * v2),
            //                        Th, In);
            //     stokes.addLinear(innerProduct(ux * dx(ux) + uy * dy(ux), rho * v1) +
            //                          innerProduct(ux * dx(uy) + uy * dy(uy), rho * v2),
            //                      Th, In);

            //     stokes.addLagrangeMultiplier(innerProduct(1., dp1), 0., Th, In);

            //     // gather(mat_NL);
            //     // matlab::Export(mat_NL[0], "mat_1.dat");
            //     // return 0;
            //     stokes.solve(mat_NL[0], stokes.rhs_);

            //     //         LOG_INFO << " Time solver \t" << MPIcf::Wtime() - tt0 << logger::endl;

            //     //         std::span<double> dwu = std::span<double>(stokes.rhs_.data(), Wh.NbDoF() * In.NbDoF());
            //     //         const auto result =
            //     //             std::max_element(dwu.begin(), dwu.end(), [](int a, int b) { return std::abs(a) <
            //     //             std::abs(b); });
            //     //         double dist = *result;
            //     //         LOG_INFO << " Error || du ||_infty = " << dist << logger::endl;

            //     //         std::span<double> dw = std::span<double>(stokes.rhs_.data(), stokes.get_nb_dof());
            //     //         std::transform(data_all.begin(), data_all.end(), dw.begin(), data_all.begin(),
            //     //                        [](double a, double b) { return a - b; });

            //     //         stokes.rhs_.resize(stokes.get_nb_dof());
            //     //         stokes.rhs_ = 0.0;

            //     //         iterNewton += 1;

            //     //         if (iterNewton == 2 || dist < 1e-10) {
            //     //             std::span<double> data_all_span(data_all);
            //     //             stokes.saveSolution(data_all_span);
            //     //             stokes.cleanBuildInMatrix();
            //     //             stokes.set_map();

            //     //             if (iter == 0) {
            //     //                 LOG_INFO << " TIME NEWTON ITERATION : \t" << MPIcf::Wtime() - t0_newton <<
            //     //                 logger::endl;
            //     //             }

            //     //             break;
            //     //         }
            //     //     }

            //     //     // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
            //     //     // {
            //     //     //   Fun_h funX(Wh, fun_x);
            //     //     //   Fun_h fun1(Wh, fun_1);
            //     //     //   R tt0 = CPUtime();
            //     //     //   double areaBubble   = integral(fun1, 0, 1) ;
            //     //     //   double centerOfMass = integral(funX, 1, 1) / areaBubble ;
            //     //     //
            //     //     //   double Pb = integralSurf(fun1, 1);
            //     //     //   double ra = sqrt(areaBubble / M_PI);
            //     //     //   double Pa = 2*M_PI*ra;
            //     //     //   double circularity = Pa / Pb;
            //     //     //   double riseVelocity = integral(uh, 1, 1) / areaBubble;
            //     //     //
            //     //     //   double q    = integralSurf(us,  0, 0 , In.map(qTime[0]));
            //     //     //   double qend = integralSurf(us,  0, lastQuadTime,
            //     //     //   In.map(qTime[lastQuadTime])); double q0   = integralSurf(u0s, 0);
            //     //     //
            //     //     //
            //     //     //   if(iter == 0) initialConcentrationSurfactant = qend;
            //     //     //   LOG_INFO << "\n Features of the drop " << logger::endl;
            //     //     //   LOG_INFO << " Time                   ->    " <<
            //     //     //   GTime::current_time()+dT  << logger::endl; LOG_INFO << " Center Of Mass
            //     //     //   ->    " << centerOfMass << logger::endl; LOG_INFO << " Circularity ->
            //     //     //   " << circularity << logger::endl; LOG_INFO << " Rise velocity ->    "
            //     //     //   << riseVelocity << logger::endl; LOG_INFO << " Surfactant quantity
            //     //     //   init->    " << q0 << logger::endl; LOG_INFO << " Surfactant quantity ->
            //     //     //   " << q << logger::endl; LOG_INFO << " Surfactant quantity end ->   " <<
            //     //     //   qend << logger::endl; LOG_INFO << " |q_0 - q_end|           ->   " <<
            //     //     //   fabs(q - qend) << logger::endl; LOG_INFO << " Surfactant conservation
            //     //     //   ->   " << fabs(qend - initialConcentrationSurfactant) << logger::endl;
            //     //     //
            //     //     //
            //     //     //   outputData << GTime::current_time()+dT << "\t"
            //     //     //   << centerOfMass << "\t"
            //     //     //   << circularity << "\t"
            //     //     //   << riseVelocity << "\t"
            //     //     //   << areaBubble <<  "\t"
            //     //     //   << qend << "\t"
            //     //     //   << fabs(q - qend) << "\t"
            //     //     //   << fabs(qend - initialConcentrationSurfactant) << logger::endl;
            //     //     //
            //     //     // }
            //     //     //

            //     //     LOG_INFO << " Set Velocity " << logger::endl;
            //     //     for (int i = 0; i < nb_quad_time; ++i) {
            //     //         fct_t sol(Wh, In, data_uh);
            //     //         set_velocity(sol, vel[i], ls[i], In.map(qTime[i]));
            //     //     }

            //     //     //  -----------------------------------------------------
            //     //     //                     PLOTTING
            //     //     //  -----------------------------------------------------

            //     if (iterations == 1) {

            //         {
            //             Paraview<mesh_t> writer(Th, path_output_figures + "Th" + std::to_string(iter + 1) + ".vtk");
            //             // Paraview<mesh_t> writer(Thi, path_output_figures + "bulk_" + std::to_string(iter + 1) +
            //             // ".vtk");                writer.add(ls[0], "levelSet", 0, 1);

            //             writer.add(uh, "velocity", 0, 2);
            //             writer.add(ph, "pressure", 0, 1);

            //             writer.add(u_exact, "velocity_exact", 0, 2);
            //             writer.add(p_exact, "pressure_exact", 0, 1);
            //         }
            //     }

            //     //     bar++;
            //     iter += 1;
            //     //     globalVariable::verbose = 1;
            //     // }
            //     // bar.end();
            // }
            h *= 0.5;
        }
    }
}