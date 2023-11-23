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

const double mu  = 10.;
const double rho = .01;

double fun_levelset(R2 P) { return -1.; }

double fun_rhs(R2 P, int i, const double t) { return 0.; }

double fun_u(R2 P, int i, const double t) {
    // double x = P[0], y = P[1];
    double x = P.x, y = P.y;
    if (i == 0)
        return std::sin(x) * std::cos(y) * std::exp(-2 * mu * t);
    else
        return (-1.) * std::cos(x) * std::sin(y) * std::exp(-2 * mu * t);
}

double fun_u_d(R2 P, int i, int dom, const double t) {
    // double x = P[0], y = P[1];
    double x = P.x, y = P.y;
    if (i == 0)
        return std::sin(x) * std::cos(y) * std::exp(-2 * mu * t);
    else
        return (-1.) * std::cos(x) * std::sin(y) * std::exp(-2 * mu * t);
}

double fun_u_initial(R2 P, int i) {
    double x = P.x, y = P.y;
    if (i == 0)
        return std::sin(x) * std::cos(y);
    else
        return (-1.) * std::cos(x) * std::sin(y);
}

double fun_p(R2 P, int i, const double t) {
    return rho / 4 * (std::cos(2 * P.x) + std::cos(2 * P.y)) * std::exp(-4 * mu * t);
}

double fun_p_d(R2 P, int i, int dom, const double t) {
    return rho / 4 * (std::cos(2 * P.x) + std::cos(2 * P.y)) * std::exp(-4 * mu * t);
}

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
        const double lambda_interior = 100.;

        // Finite element spaces
        Lagrange2 FEu(2); // for interpolating the exact solution
        space_t Vh_interpolation(Th, FEu);

        // Taylor-Hood
        // space_t Vh(Th, FEu);
        // space_t Ph(Th, DataFE<mesh_t>::P1);

        // P0 x BDM1
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

            // uh^{k+1} = uh^k - du
            // ph^{k+1} = ph^k + dp
            funtest_t du(Vh, 2), dp(Ph, 1), v(Vh, 2), q(Ph, 1), dp1(Vh, 1, 0, 0);

            // DEFINE THE PROBLEM ON THE CORRECT SPACES
            // ----------------------------------------------
            // navier_stokes.initSpace(Vh, In);
            // navier_stokes.add(Ph, In);
            navier_stokes.initSpace(Vhn, In);
            navier_stokes.add(Phn, In);

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

                    // -b(v, p) + b(u, q) (or -b(v, p) + b0(u, q))
                    navier_stokes.addBilinear(-innerProduct(dp, div(v)) + innerProduct(div(du), q), Thi, In);
                    navier_stokes.addBilinear(innerProduct(dp, v * n) - innerProduct(du * n, q), Thi, INTEGRAL_BOUNDARY,
                                              In);
                    // navier_stokes.addBilinear(innerProduct(dp, v * n), Thi, INTEGRAL_BOUNDARY, In);
                }

                // Add -Lh(vh)
                navier_stokes.addLinear(-innerProduct(fh.exprList(), rho * v), Thi, In);
                navier_stokes.addLinear(-innerProduct(u0.exprList(), rho * v), Thi);
                navier_stokes.addLinear(innerProduct(gh.exprList(), mu * grad(v) * n) -
                                            innerProduct(gh.exprList(), lambda_boundary * v) +
                                            innerProduct(gh.exprList(), q * n),
                                        Thi, INTEGRAL_BOUNDARY, In);

                // navier_stokes.addLinear(innerProduct(gh.exprList(), mu * grad(v) * n) -
                //                             innerProduct(gh.exprList(), lambda_boundary * v),
                //                         Thi, INTEGRAL_BOUNDARY, In);

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

                navier_stokes.solve(mat_NL[0], navier_stokes.rhs_);

                // Compute norm of the difference in the succesive approximations
                std::span<double> dwu = std::span<double>(navier_stokes.rhs_.data(), Vhn.NbDoF() * In.NbDoF());
                const auto result =
                    std::max_element(dwu.begin(), dwu.end(), [](int a, int b) { return std::abs(a) < std::abs(b); });
                double dist = *result;
                std::cout << " Residual error: " << dist << "\n" << "\n";

            //     //         std::span<double> dw = std::span<double>(stokes.rhs_.data(), stokes.get_nb_dof());
            //     //         std::transform(data_all.begin(), data_all.end(), dw.begin(), data_all.begin(),
            //     //                        [](double a, double b) { return a - b; });

            //     //         stokes.rhs_.resize(stokes.get_nb_dof());
            //     //         stokes.rhs_ = 0.0;

            //     //         iterNewton += 1;

                if (dist < 1e-10) {
                    std::span<double> data_all_span(data_all);
                    navier_stokes.saveSolution(data_all_span);
                    navier_stokes.cleanBuildInMatrix();
                    navier_stokes.set_map();

                    break;
                }

                if (newton_iterations >= 3) {
                    std::cout << "Newton's method didn't converge in 3 iterations, breaking loop. \n";
                    break;
                }
            }

            // Compute L2 errors
            std::vector<double> sol_uh(Vhn.get_nb_dof());
            std::vector<double> sol_ph(Phn.get_nb_dof());

            std::vector<double> zeros_u(Vhn.get_nb_dof());
            std::vector<double> zeros_p(Phn.get_nb_dof());

            for (int n = 0; n < ndf_time_slab; ++n) {
                // get the DOFs of u corresponding to DOF n in time and sum with the previous n
                std::vector<double> u_dof_n(data_uh.begin() + n * Vhn.get_nb_dof(),
                                            data_uh.begin() + (n + 1) * Vhn.get_nb_dof());
                std::transform(sol_uh.begin(), sol_uh.end(), u_dof_n.begin(), sol_uh.begin(), std::plus<double>());

                // get the DOFs of p corresponding to DOF n in time and sum with the previous n
                std::vector<double> p_dof_n(data_ph.begin() + n * Phn.get_nb_dof(),
                                            data_ph.begin() + (n + 1) * Phn.get_nb_dof());
                std::transform(sol_ph.begin(), sol_ph.end(), p_dof_n.begin(), sol_ph.begin(), std::plus<double>());
            }

            fct_t fun_uh(Vhn, sol_uh);
            fct_t fun_ph(Phn, sol_ph);
            fct_t fun_zeros_u(Vhn, zeros_u);
            fct_t fun_zeros_p(Phn, zeros_p);

            error_uh = L2normCut(fun_uh, fun_u_d, current_time + dT, 0, 2);
            error_ph = L2normCut(fun_ph, fun_p_d, current_time + dT, 0, 1);

            u_norm = L2normCut(fun_zeros_u, fun_u_d, current_time + dT, 0, 2);
            p_norm = L2normCut(fun_zeros_p, fun_p_d, current_time + dT, 0, 1);

            std::cout << " ||u(T)-uh(T)||_2 = " << error_uh/u_norm << '\n';
            std::cout << " ||p(T)-ph(T)||_2 = " << error_ph/p_norm << '\n';

            errors_uh[j] = error_uh;
            errors_ph[j] = error_ph;

            // Plotting
            if (iterations == 1) {

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