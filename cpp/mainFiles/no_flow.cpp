


/**
 * @file no_flow.cpp
 * @author Sebastian Myrbäck
 * @brief Reproducing Ex 1.1 (no-flow Stokes eq) in Neilan et. al.
 * @version 0.1
 * @date 2024-03-04
 *
 * @copyright Copyright (c) 2024
 *
 */

// Dependencies
#include "../num/matlab.hpp"
#include "../num/gnuplot.hpp"
#include "../tool.hpp"
#include "../FESpace/funfem_util.hpp"
#include <filesystem>

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;

// User-defined parameters
const double nu  = 1.;
const double Ra = 1.;
double boundary_penalty = 1e1;
double tangential_penalty = 1e1;

namespace no_flow {

    double fun_levelset(R2 P) { return -1.; }

    double fun_rhs(R2 P, int i) { 
        const double y = P.y;
        if (i == 0)
            return 0.; 
        else 
            return Ra*(1-y+3*y*y);
    }

    double fun_u(double* P, int i) {
        return 0.;
    }

    double fun_u_d(R2 P, int i, int dom) {
        return 0.;
    }

    double fun_p(double* P, int i) {
        const double y = P[1];
        return Ra*(y*y*y - y*y/2 + y - 7./12);
    }

    double fun_p_d(double* P, int i, int dom) {
        const double y = P[1];
        return Ra*(y*y*y - y*y/2 + y - 7./12);
    }
}

using namespace no_flow;

#define SCOTT_VOGELIUS

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);

    const size_t mesh_refinements = 4;
    double h = 0.1;     // not used if mesh is imported
    const size_t thread_count = 1;

    ProblemOption optionProblem;
    optionProblem.solver_name_  = "mumps";
    optionProblem.clear_matrix_ = true;

    const std::string path_output_data    = "../output_files/stokes/no_flow/data/";
    const std::string path_output_figures = "../output_files/stokes/no_flow/paraview/";

    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_output_figures);
    }

    std::array<double, mesh_refinements> errors_uh, errors_ph, errors_div_uh, errors_grad_uh, hs, numb_elems;

    for (int j = 0; j < mesh_refinements; ++j) {

        // Mesh
        // const double lx = 1., ly = 1.;
        // const double x0 = 0., y0 = 0.;
        // int nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
        // mesh_t Th(nx, ny, x0, y0, lx, ly);
    
        mesh_t Th(("../cpp/meshes/no_flow_ps_" + std::to_string(j) + ".msh").c_str());
        //mesh_t Th(("../cpp/meshes/no_flow_ct_" + std::to_string(j) + ".msh").c_str());

        //int nx = (int)std::sqrt(Th.nbElements()); 
        // h = std::sqrt(2./Th.nbElements());
        h = Th[0].lenEdge(0);

        boundary_penalty = boundary_penalty * nu / h;

        // Finite element spaces
        Lagrange2 FEu(2); // for interpolating the exact solution
        space_t Vh_interpolation(Th, FEu);
        
        // P1-interpolated level set function
        space_t Lh(Th, DataFE<mesh_t>::P1);
        fct_t gamma(Lh, fun_levelset);
        InterfaceLevelSet<mesh_t> interface(Th, gamma);
        ActiveMesh<mesh_t> active_mesh(Th, interface); 
        
        gnuplot::save(Th, "no_flow" + std::to_string(j) + ".dat");
        //gnuplot::save(interface, "gamma" + std::to_string(j) + ".dat");

        // Scott-Vogelius element pair
    #ifdef SCOTT_VOGELIUS
        space_t V(Th, FEu);
        space_t P(Th, DataFE<mesh_t>::P1dc);
    #elif defined(TAYLOR_HOOD)
        space_t V(Th, FEu);
        space_t P(Th, DataFE<mesh_t>::P1);
    #elif defined(BDM1_P0)
        space_t V(Th, DataFE<mesh_t>::BDM1);
        space_t P(Th, DataFE<mesh_t>::P0);
    #endif

        // Cut finite element spaces
        cutspace_t Vh(active_mesh, V);
        cutspace_t Ph(active_mesh, P);

        // Create Stokes problem object
        CutFEM<mesh_t> no_flow(Vh, thread_count, optionProblem);
        no_flow.add(Ph);

        std::cout << "--------------------------------------------" << '\n';
        std::cout << "Iteration " << j + 1 << "/" << mesh_refinements << '\n';
        std::cout << "\n h  = " << h << '\n';
        std::cout << "N = " << Th.nbElements() << '\n';

        Normal n;
        Tangent t;

        funtest_t u(Vh, 2), p(Ph, 1), v(Vh, 2), q(Ph, 1);

        // Interpolate exact functions
        fct_t u_exact(Vh, fun_u);
        fct_t p_exact(Ph, fun_p);
        fct_t fh(Vh, fun_rhs);     // rhs force
        fct_t gh(Vh, fun_u);       // Dirichlet boundary condition

        // Add macro element ghost penalty
        //MacroElement<mesh_t> macro(Th, delta);
        
        
        // Bilinear form
        no_flow.addBilinear(
            + contractProduct(nu * grad(u), grad(v))
            - innerProduct(p, div(v)) 
            + innerProduct(div(u), q)
            , active_mesh);

        no_flow.addBilinear(
            - innerProduct(nu*grad(u)*n, v)
            - innerProduct(u, nu*grad(v)*n) 
            + innerProduct(boundary_penalty*u, v)
            + innerProduct(p, v*n)
            //+ innerProduct(u, q*n)      // make problem block-symmetric
            , active_mesh
            , INTEGRAL_BOUNDARY);

        #ifdef BDM1_P0
        no_flow.addBilinear(
            - innerProduct(average(nu * grad(u) * t * n, 0.5, 0.5), jump(v * t))
            - innerProduct(jump(u * t), average(nu * grad(v) * t * n, 0.5, 0.5))
            + innerProduct(tangential_penalty * (jump(u * t)), jump(v * t))
            , active_mesh
            , INTEGRAL_INNER_EDGE_2D);
        #endif

        // Linear form
        no_flow.addLinear(
            + innerProduct(fh.exprList(), v)
            , active_mesh);

        no_flow.addLinear(
            - innerProduct(gh.exprList(), nu*grad(v)*n) 
            + innerProduct(gh.exprList(), boundary_penalty*v)
            //+ innerProduct(gh.exprList(), q*n)  // compensate for added term
            , active_mesh
            , INTEGRAL_BOUNDARY); 

        // Add Lagrange multipliers

        // no_flow.addLagrangeMultiplier(
        //     + innerProduct(1., p)
        //     , 0.
        //     , active_mesh
        // );

        CutFEM<mesh_t> lagr(Vh);
        lagr.add(Ph);
        lagr.addLinear(innerProduct(1., p), active_mesh);
        //lagr.addLinear(innerProduct(1., div(v)), active_mesh);
        std::vector<double> lag_row(lagr.rhs_.begin(), lagr.rhs_.end());
        lagr.rhs_ = 0.;
        lagr.addLinear(innerProduct(1., v * n), active_mesh, INTEGRAL_BOUNDARY);
        //lagr.addLinear(innerProduct(1., p), active_mesh);

        no_flow.addLagrangeVecToRowAndCol(lag_row, lagr.rhsSpan(), 0.);

        // Export matrix
        matlab::Export(no_flow.mat_[0], path_output_data + "mat_" + std::to_string(j + 1) + ".dat");
    
        no_flow.solve(no_flow.mat_[0], no_flow.rhs_);

        // Extract numerical solution data
        std::span<double> data_uh = std::span<double>(no_flow.rhs_.data(), Vh.get_nb_dof());
        std::span<double> data_ph = std::span<double>(no_flow.rhs_.data() + Vh.get_nb_dof(), Ph.get_nb_dof());
        
        // Interpolate numerical solution onto FE spaces
        fct_t uh(Vh, data_uh);
        fct_t ph(Ph, data_ph);

        // Post process pressure
        // double meanP    = integral(Th, p_exact.expr(), 0);
        // double meanPfem = integral(Th, ph.expr(), 0);
        // FEM<mesh_t> post(Ph);
        // post.addLinear(innerProduct(1, q), Th);
        // double area = post.rhs_.sum();
        // ph.v -= meanPfem / area;
        // ph.v += meanP / area;

        std::vector<double> zeros_u(Vh.get_nb_dof());
        std::vector<double> zeros_p(Ph.get_nb_dof());

        auto uh_0dx = dx(uh.expr(0));
        auto uh_1dy = dy(uh.expr(1));
        auto uh_0dy = dy(uh.expr(0));
        auto uh_1dx = dx(uh.expr(1));

        auto u_0dx = dx(u_exact.expr(0));
        auto u_1dy = dy(u_exact.expr(1));
        auto u_0dy = dy(u_exact.expr(0));
        auto u_1dx = dx(u_exact.expr(1));

        // fct_t uh(Vh, data_u);
        // fct_t ph(Ph, data_p);
        fct_t fun_zeros_u(Vh, zeros_u);
        fct_t fun_zeros_p(Ph, zeros_p);
        fct_t u_error(Vh, fun_u);
        fct_t p_error(Ph, fun_p);
        

        double error_uh = L2normCut(uh, fun_u_d, 0, 2);
        double error_ph = L2normCut(ph, fun_p_d, 0, 1);
        double error_div_uh = maxNormCut(uh_0dx + uh_1dy, active_mesh);
        //double error_div_uh = maxNorm(uh_0dx + uh_1dy, Th);
        double error_grad_uh = std::sqrt(integral(active_mesh, (uh_0dx - u_0dx)*(uh_0dx - u_0dx) + (uh_0dy-u_0dy)*(uh_0dy-u_0dy) + (uh_1dx-u_1dx)*(uh_1dx-u_1dx) + (uh_1dy-u_1dy)*(uh_1dy-u_1dy)));

        // double u_norm = L2norm(fun_zeros_u, fun_u_d, 0, 2);
        // double p_norm = L2norm(fun_zeros_p, fun_p_d, 0, 1);

        u_error.v -= uh.v;
        p_error.v -= ph.v;

        // errors_uh[j] = error_uh / u_norm;
        // errors_ph[j] = error_ph / p_norm;

        std::cout << " ||u-uh||_2 = " << error_uh << '\n';
        std::cout << " ||p-ph||_2 = " << error_ph << '\n';
        std::cout << " ||grad(u-uh)||_2 = " << error_grad_uh << '\n';
        std::cout << " ||div(uh)||_infty = " << error_div_uh << '\n';

        errors_uh[j] = error_uh;
        errors_ph[j] = error_ph;
        errors_div_uh[j] = error_div_uh;
        errors_grad_uh[j] = error_grad_uh;
        hs[j] = h;
        numb_elems[j] = Th.nbElements();

        // Write solutions to Paraview
        Paraview<mesh_t> paraview(active_mesh, path_output_figures + "no_flow_" + std::to_string(j) + ".vtk");
        //paraview.add(gamma, "gamma", 0, 1);
        paraview.add(uh, "velocity", 0, 2);
        paraview.add(ph, "pressure", 0, 1);
        paraview.add(u_exact, "velocity_exact", 0, 2);
        paraview.add(p_exact, "pressure_exact", 0, 1);
        paraview.add(u_error, "u_error", 0, 2);
        paraview.add(p_error, "p_error", 0, 1);
        paraview.add(fabs(uh_0dx + uh_1dy), "divergence");
        paraview.add(fh, "rhs", 0, 2);
        //paraview.writeActiveMesh(active_mesh, path_output_figures + "active_mesh_" + std::to_string(j) + ".vtk");
        
        h *= 0.5;

    }

    std::cout << std::setprecision(16);
    std::cout << '\n';
    std::cout << "Errors Velocity = [";
    for (int i = 0; i < mesh_refinements; i++) {

        std::cout << errors_uh.at(i);
        if (i < mesh_refinements - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';
    std::cout << "Errors Pressure = [";
    for (int i = 0; i < mesh_refinements; i++) {

        std::cout << errors_ph.at(i);
        if (i < mesh_refinements - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "Gradient error = [";
    for (int i = 0; i < mesh_refinements; i++) {

        std::cout << errors_grad_uh.at(i);
        if (i < mesh_refinements - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << "Divergence Errors = [";
    for (int i = 0; i < mesh_refinements; i++) {

        std::cout << errors_div_uh.at(i);
        if (i < mesh_refinements - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';


    std::cout << "h = [";
    for (int i = 0; i < mesh_refinements; i++) {

        std::cout << hs.at(i);
        if (i < mesh_refinements - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';
    std::cout << "Number of mesh elements = [";
    for (int i = 0; i < mesh_refinements; i++) {

        std::cout << numb_elems.at(i);
        if (i < mesh_refinements - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    return 0;
}

