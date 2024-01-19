
/**
 * @file navier_stokes_raising_drop_2D.cpp
 * @author Thomas Frachon
 * @brief Raising bubble using CutFEM. For details, one can see Example 2 in
 * "A cut finite element method for incompressible two-phase Navier-Stokes
 * flows" by  Thomas Frachon and Sara Zahedi
 * @version 0.1
 * @date 2023-03-09
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "../num/matlab.hpp"
#include "../tool.hpp"
#include "../FESpace/funfem_util.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;

const double mu  = 1.;
const double rho = 1.;

double fun_levelset(R2 P) { return -1.; }

double fun_rhs(R2 P, int i, const double t) { return 0.; }

double fun_u(R2 P, int i, const double t) {
    // double x = P[0], y = P[1];
    double x = P.x, y = P.y;
    if (i == 0)
        return std::sin(x) * std::cos(y) * std::exp(-2 * mu / rho * t);
    else
        return (-1.) * std::cos(x) * std::sin(y) * std::exp(-2 * mu / rho * t);
}

double fun_u_d(R2 P, int i, int dom, const double t) {
    // double x = P[0], y = P[1];
    double x = P.x, y = P.y;
    if (i == 0)
        return std::sin(x) * std::cos(y) * std::exp(-2 * mu / rho * t);
    else
        return (-1.) * std::cos(x) * std::sin(y) * std::exp(-2 * mu / rho * t);
}

double fun_u_initial(R2 P, int i) {
    double x = P.x, y = P.y;
    if (i == 0)
        return std::sin(x) * std::cos(y);
    else
        return (-1.) * std::cos(x) * std::sin(y);
}

double fun_p(R2 P, int i, const double t) {
    return rho / 4 * (std::cos(2 * P.x) + std::cos(2 * P.y)) * std::exp(-4 * mu / rho * t);
}

double fun_p_d(R2 P, int i, int dom, const double t) {
    return rho / 4 * (std::cos(2 * P.x) + std::cos(2 * P.y)) * std::exp(-4 * mu / rho * t);
}

#define epsnot // run with Eps(u) instead of grad(u)

int main(int argc, char **argv) {

    double h = 0.5; // starting mesh size

    const size_t iterations = 4;

    //MPIcf cfMPI(argc, argv);

    const std::string path_output_data    = "/NOBACKUP/smyrback/output_files/navier_stokes/taylor_green/data/";
    const std::string path_output_figures = "/NOBACKUP/smyrback/output_files/navier_stokes/taylor_green/paraview/";

    // if (MPIcf::IamMaster()) {
    //     std::filesystem::create_directories(path_output_data);
    //     std::filesystem::create_directories(path_output_figures);
    // }

    std::array<double, iterations> errors_uh, errors_ph, errors_div_uh, errors_uh_T, errors_ph_T, hs;
    for (int j = 0; j < iterations; ++j) {

        // Mesh
        const double lx = M_PI, ly = M_PI;
        const double x0 = 0., y0 = 0.;
        int nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
        Mesh2 Th(nx, ny, x0, y0, lx, ly);

        hs[j] = h;

        // Th.info();
        // Th.info();

        // Time stepping
        int division_mesh_size     = 2;
        double final_time          = .25;
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
        std::vector<std::map<std::pair<int, int>, double>> mat_NL(1);

        CutFEM<mesh_t> navier_stokes(qTime, optionProblem);

        // These parameters work fine
        // const double lambda_boundary = 100.;
        // const double lambda_interior = 10.;

        //const double lambda_boundary = 10.;
        const double lambda_boundary_tangent = 10. * mu / h;
        const double lambda_interior = 1. * mu / h;

        // Finite element spaces
        Lagrange2 FEu(2); // for interpolating the exact solution
        space_t Vh_interpolation(Th, FEu);

        // Taylor-Hood
        // space_t Vh(Th, FEu);
        // space_t Ph(Th, DataFE<mesh_t>::P1);

        // P0 x BDM1
        space_t Vh(Th, DataFE<mesh_t>::BDM1);
        space_t Ph(Th, DataFE<mesh_t>::P0);

        space_t Lh(Th, DataFE<mesh_t>::P1);
        TimeInterface<mesh_t> interface(qTime);
        double dt_levelset = dT / (nb_quad_time - 1);
        std::vector<fct_t> ls(nb_quad_time);

        for (int i = 0; i < nb_quad_time; ++i) {
            ls[i].init(Lh, fun_levelset);
        }

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

        std::cout << "number of time slabs \t : \t " << total_number_iteration << '\n';

        double error_uh = 0, error_ph = 0, error_div_uh = 0, u_norm = 0, p_norm = 0, error_I_uh = 0, error_I_ph = 0;
        std::vector<double> divergence_errors_t;
        int iter = 0;
        while (iter < total_number_iteration) {
            int current_iteration = iter;
            double current_time   = iter * time_step;

            const TimeSlab &In(Ih[iter]);

            for (int i = 0; i < nb_quad_time; ++i) {
                interface.init(i, Th, ls[i]);
            }

            cutmesh_t Thi(Th);
            Thi.truncate(interface, 1);

            cutspace_t Vhn(Thi, Vh);
            cutspace_t Phn(Thi, Ph);

            std::cout << " -------------------------------------------------------\n";
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << '\n';
            std::cout << " Time      \t : \t" << current_iteration * time_step << '\n';

            // DEFINE TEST FUNCTIONS
            // ----------------------------------------------
            Normal n;
            Tangent t;

            // uh^{k+1} = uh^k - du
            // ph^{k+1} = ph^k - dp
            // We will solve for (du, dp)
            funtest_t du(Vhn, 2), dp(Phn, 1), v(Vhn, 2), q(Phn, 1), dp1(Vhn, 1, 0, 0);

            // Create tensor-product spaces
            navier_stokes.initSpace(Vhn, In);
            navier_stokes.add(Phn, In);

            // Interpolate exact functions
            fct_t u_exact(Vhn, In, fun_u);
            fct_t p_exact(Phn, In, fun_p);
            fct_t fh(Vhn, In, fun_rhs); // rhs force
            fct_t gh(Vhn, In, fun_u);   // Dirichlet boundary condition

            // Initialize DOFs and data
            std::vector<double> data_init(navier_stokes.get_nb_dof());
            std::span<double> data_init_span(data_init);
            std::span<double> data_uh0 = std::span<double>(data_init.data(), Vhn.NbDoF()); // velocity in first time DOF

            if (iter == 0) {
                interpolate(Vhn, data_uh0, fun_u_initial);
            } else {
                navier_stokes.initialSolution(data_init_span);
            }

            std::vector<double> data_all(data_init);
            int idxp0 = Vhn.NbDoF() * In.NbDoF(); // index for when p starts in the array
            std::span<double> data_uh  = std::span<double>(data_all.data(), Vhn.NbDoF() * In.NbDoF()); // velocity for all time DOFs
            std::span<double> data_ph  = std::span<double>(data_all.data() + idxp0, Phn.NbDoF() * In.NbDoF());

            // Create FEM functions corresponding to the numerical solutions
            fct_t u0(Vhn, data_uh0);
            fct_t uh(Vhn, In, data_uh);
            fct_t ph(Phn, In, data_ph);

            auto boundary_dof = getBoundaryDof(Vhn, In);    // get DOFs in Vhn on the boundary
            setBoundaryDof(boundary_dof, gh, data_uh);      // set BC strongly on initial solution

            // Newton's method
            int newton_iterations = 0;
            while (1) {

                // Add terms that are in the residual
                if (newton_iterations == 0) {

                    // Time terms
                    navier_stokes.addBilinear(
                        + innerProduct(dt(du), rho * v)
                        + contractProduct(mu * grad(du), grad(v))
                        - innerProduct(dp, div(v)) 
                        + innerProduct(div(du), q)
                        , Thi
                        , In);

                    navier_stokes.addBilinear(
                        + innerProduct(rho * du, v)
                        , Thi
                        , 0
                        , In);

                    // A(t, u, v)
                    navier_stokes.addBilinear(
                        - innerProduct(mu * grad(du) * n, v)
                        - innerProduct(du, mu * grad(v) * n) 
                        + innerProduct(lambda_boundary_tangent * du * t, v * t)
                        + innerProduct(dp, v * n)
                        , Thi
                        , INTEGRAL_BOUNDARY
                        , In);

                    navier_stokes.addBilinear(
                        - innerProduct(average(grad(du * t) * n, 0.5, 0.5), mu * jump(v * t)) 
                        + innerProduct(jump(du * t), mu * average(grad(v * t) * n, 0.5, 0.5)) 
                        + innerProduct(lambda_interior * (jump(du * t)), jump(v * t))
                        , Thi
                        , INTEGRAL_INNER_EDGE_2D
                        , In);

                }

                // Add -Lh(vh)
                navier_stokes.addLinear(
                    - innerProduct(fh.exprList(), rho * v)
                    , Thi
                    , In);

                navier_stokes.addLinear(
                    - innerProduct(u0.exprList(), rho * v)
                    , Thi);


                navier_stokes.addLinear(
                    + innerProduct(gh.exprList(), mu * grad(v) * n) 
                    - innerProduct(gh.exprList(), lambda_boundary_tangent * v * t * t)
                    , Thi
                    , INTEGRAL_BOUNDARY
                    , In); 

                // set zero BC strongly for delta u, otherwise BC will accumulate with Newton iterations
                setBoundaryDof(navier_stokes.rhs_.size(), boundary_dof, 1., navier_stokes.mat_[0]);

                navier_stokes.gather_map();
                navier_stokes.addMatMul(data_all); // add B(uh^k, vh) to rhs

                // Construct remaining terms of the Jacobian

                // Assemble Jacobian in the matrix mat_NL
                navier_stokes.set_map(mat_NL);
                mat_NL[0] = navier_stokes.mat_[0];

                funtest_t du1(Vhn, 1, 0), du2(Vhn, 1, 1), v1(Vhn, 1, 0), v2(Vhn, 1, 1);
                auto ux = uh.expr(0);
                auto uy = uh.expr(1);

                // Linearized advection term
                navier_stokes.addBilinear(
                    + innerProduct(du1 * dx(ux) + du2 * dy(ux), rho * v1) 
                    + innerProduct(du1 * dx(uy) + du2 * dy(uy), rho * v2) 
                    + innerProduct(ux * dx(du1) + uy * dy(du1), rho * v1) 
                    + innerProduct(ux * dx(du2) + uy * dy(du2), rho * v2)
                    , Thi
                    , In);

                navier_stokes.addLinear(
                    + innerProduct(ux * dx(ux) + uy * dy(ux), rho * v1) 
                    + innerProduct(ux * dx(uy) + uy * dy(uy), rho * v2)
                    , Thi
                    , In);

                // set zero BC for update
                setBoundaryDof(boundary_dof, 0. , navier_stokes.rhs_);  
                setBoundaryDof(navier_stokes.rhs_.size(), boundary_dof, 1. , mat_NL[0]);  

                // Add Lagrange multipliers
                CutFEM<Mesh2> lagrange(qTime, optionProblem);
                lagrange.initSpace(Vhn, In);
                lagrange.add(Phn, In);
      
                Rn lag_row(lagrange.rhs_.size(), 0.);

                // Add multipliers in first and last time quadrature point
                for (int itq = 0; itq < 2; ++itq) {
                    std::fill(lagrange.rhs_.begin(),lagrange.rhs_.end(),0.);   
                                     lagrange.addLinear(innerProduct(1., q), Thi, itq * last_quad_time, In);
                    lag_row       = lagrange.rhs_;
                    std::fill(lagrange.rhs_.begin(),lagrange.rhs_.end(),0.); 
                                       lagrange.addLinear(innerProduct(1., v * n), Thi, INTEGRAL_BOUNDARY, itq * last_quad_time, In);
                    // lagrange.addLinear(innerProduct(1., q), Thi, itq * last_quad_time, In);
                    navier_stokes.addLagrangeVecToRowAndCol(lag_row, lagrange.rhs_, 0);
                }

                // Export matrix
                if (iter == total_number_iteration - 1) {
                    matlab::Export(mat_NL[0], path_output_data + "mat_" + std::to_string(j + 1) + ".dat");
                }

                navier_stokes.solve(mat_NL[0], navier_stokes.rhs_);


                // Compute norm of the difference in the succesive approximations
                std::span<double> dwu = std::span<double>(navier_stokes.rhs_.data(), Vhn.NbDoF() * In.NbDoF());
                // std::cout << dwu << "\n";
                // getchar();
                const auto result     = std::max_element(dwu.begin(), dwu.end());
                double dist           = *result;
                std::cout << " Residual error: " << dist << "\n"
                          << "\n";

                std::span<double> dw = std::span<double>(navier_stokes.rhs_.data(), navier_stokes.get_nb_dof());
                std::transform(data_all.begin(), data_all.end(), dw.begin(), data_all.begin(),
                               [](double a, double b) { return a - b; });

                navier_stokes.rhs_.resize(navier_stokes.get_nb_dof());
                       std::fill(navier_stokes.rhs_.begin(),navier_stokes.rhs_.end(),0.); 
                


                newton_iterations += 1;

                if (std::abs(dist) < 1e-10) {
                    std::span<double> data_all_span(data_all);
                    navier_stokes.saveSolution(data_all_span);
                    navier_stokes.cleanBuildInMatrix();
                    navier_stokes.set_map();

                    break;
                }

                if (newton_iterations >= 5) {
                    std::cout << "Newton's method didn't converge in 5 iterations, "
                                 "breaking loop. \n";
                    break;
                }
            }

            // Compute L2 errors
            std::vector<double> sol_uh(Vhn.get_nb_dof());
            std::vector<double> sol_ph(Phn.get_nb_dof());

            std::vector<double> zeros_u(Vhn.get_nb_dof());
            std::vector<double> zeros_p(Phn.get_nb_dof());

            for (int n = 0; n < ndf_time_slab; ++n) {
                // get the DOFs of u corresponding to DOF n in time and sum with the
                // previous n
                std::vector<double> u_dof_n(data_uh.begin() + n * Vhn.get_nb_dof(),
                                            data_uh.begin() + (n + 1) * Vhn.get_nb_dof());
                std::transform(sol_uh.begin(), sol_uh.end(), u_dof_n.begin(), sol_uh.begin(), std::plus<double>());

                // get the DOFs of p corresponding to DOF n in time and sum with the
                // previous n
                std::vector<double> p_dof_n(data_ph.begin() + n * Phn.get_nb_dof(),
                                            data_ph.begin() + (n + 1) * Phn.get_nb_dof());
                std::transform(sol_ph.begin(), sol_ph.end(), p_dof_n.begin(), sol_ph.begin(), std::plus<double>());
            }

            fct_t fun_uh(Vhn, sol_uh);
            fct_t fun_ph(Phn, sol_ph);
            fct_t fun_zeros_u(Vhn, zeros_u);
            fct_t fun_zeros_p(Phn, zeros_p);
            auto uh_0dx = dx(fun_uh.expr(0));
            auto uh_1dy = dy(fun_uh.expr(1));

            error_uh = L2normCut(fun_uh, fun_u_d, current_time + dT, 0, 2);
            error_ph = L2normCut(fun_ph, fun_p_d, current_time + dT, 0, 1);

            u_norm = L2normCut(fun_zeros_u, fun_u_d, current_time + dT, 0, 2);
            p_norm = L2normCut(fun_zeros_p, fun_p_d, current_time + dT, 0, 1);

            error_div_uh = maxNormCut(uh_0dx + uh_1dy, Thi);
            divergence_errors_t.push_back(error_div_uh);

            fct_t fun_uh_t(Vhn, In, data_uh);
            fct_t fun_ph_t(Phn, In, data_ph);

            error_I_uh += L2normCut_T(fun_uh_t, fun_u_d, Thi, In, qTime, 2);
            error_I_ph += L2normCut_T(fun_ph_t, fun_p_d, Thi, In, qTime, 1);

            // std::cout << " ||u(T)-uh(T)||_2 / ||u(T)||_2 = " << error_uh / u_norm
            // << '\n'; std::cout << " ||p(T)-ph(T)||_2 / ||p(T)||_2 = " << error_ph /
            // p_norm << '\n';

            // errors_uh[j] = error_uh / u_norm;
            // errors_ph[j] = error_ph / p_norm;

            std::cout << " ||u(T)-uh(T)||_2 = " << error_uh << '\n';
            std::cout << " ||p(T)-ph(T)||_2 = " << error_ph << '\n';

            std::cout << " int_{In} ||u(t)-uh(t)||^2_2 dt = " << error_I_uh << '\n';
            std::cout << " int_{In} ||p(t)-ph(t)||^2_2 dt = " << error_I_ph << '\n';

            errors_uh[j] = error_uh;
            errors_ph[j] = error_ph;

            std::cout << " ||div(uh(T))||_infty = " << error_div_uh << '\n';
            errors_div_uh[j] = error_div_uh;

            // Plotting
            if (iterations == 1) {

                Paraview<mesh_t> writerTh(Th, path_output_figures + "Th.vtk");
                Paraview<mesh_t> writer(Thi,
                                        path_output_figures + "navier_stokes_" + std::to_string(iter + 1) + ".vtk");
                writer.add(ls[0], "levelSet", 0, 1);
                // writer.add(uh, "velocity", 0, 2);
                writer.add(fun_uh, "velocity", 0, 2);
                // writer.add(ph, "pressure", 0, 1);
                writer.add(fun_ph, "pressure", 0, 1);
                writer.add(u_exact, "velocity_exact", 0, 2);
                writer.add(p_exact, "pressure_exact", 0, 1);

                writer.add(uh_0dx + uh_1dy, "divergence");

                fct_t u_exact_T(Vhn, fun_u, current_time + dT);
                fct_t p_exact_T(Phn, fun_p, current_time + dT);

                writer.add(fabs(fun_uh.expr() - u_exact_T.expr()), "velocity_error");
                writer.add(fabs(fun_ph.expr() - p_exact_T.expr()), "pressure_error");

                writer.writeActiveMesh(Thi, path_output_figures + "ActiveMesh" + std::to_string(iter + 1) + ".vtk");
            } else {
                Paraview<mesh_t> writerTh(Th, path_output_figures + "Th.vtk");
                Paraview<mesh_t> writer(Thi,
                                        path_output_figures + "navier_stokes_" + std::to_string(j) + ".vtk");
                writer.add(ls[0], "levelSet", 0, 1);
                // writer.add(uh, "velocity", 0, 2);
                writer.add(fun_uh, "velocity", 0, 2);
                // writer.add(ph, "pressure", 0, 1);
                writer.add(fun_ph, "pressure", 0, 1);
                writer.add(u_exact, "velocity_exact", 0, 2);
                writer.add(p_exact, "pressure_exact", 0, 1);

                writer.add(uh_0dx + uh_1dy, "divergence");

                fct_t u_exact_T(Vhn, fun_u, current_time + dT);
                fct_t p_exact_T(Phn, fun_p, current_time + dT);

                writer.add(fabs(fun_uh.expr() - u_exact_T.expr()), "velocity_error");
                writer.add(fabs(fun_ph.expr() - p_exact_T.expr()), "pressure_error");

                writer.writeActiveMesh(Thi, path_output_figures + "ActiveMesh" + std::to_string(j) + ".vtk");
            }

            iter += 1;
        }

        errors_uh_T[j] = std::sqrt(error_I_uh);
        errors_ph_T[j] = std::sqrt(error_I_ph);

        h *= 0.5;

        if (iterations == 1) {
            std::cout << "Divergence Errors Time = [";
            for (auto &err : divergence_errors_t) {

                std::cout << err;

                std::cout << ", ";
            }
            std::cout << "]\n";
        }
    }

    std::cout << std::setprecision(16);
    std::cout << '\n';
    std::cout << "Errors Velocity = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_uh.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';
    std::cout << "Errors Pressure = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_ph.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << "Divergence Errors = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_div_uh.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Errors [0, T] Velocity = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_uh_T.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Errors [0, T] Pressure = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors_ph_T.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << "h = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << hs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';
}
