/**
 * @file stokes.cpp
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

using namespace globalVariable;

// User-defined parameters
double boundary_penalty = 4e3;
double tangential_penalty = 1.;
const double mu = 1.;
double tau_u = 1.;
double tau_p = 1.;
const double delta = 0.8;

double shift        = 0.5;
double radius = std::sqrt(0.2)-Epsilon;


namespace fictitious {

    const double K = 1e0;

    double fun_levelset(R2 P) { return std::sqrt((P.x - shift) * (P.x - shift) + (P.y - shift) * (P.y - shift)) - radius; }

    double fun_rhs(R2 P, int i, int dom) {
        R x = P.x;
        R y = P.y;
        if (i == 0)
            //return 40 * x * x * x - 40 * x * y * y - 32 * y + 16;
            return 40*K*x*(x*x - y*y) - 32*y + 16;

        else
            //return -40 * x * x * y + 32 * x + 40 * y * y * y - 16;
            return 32*x - 40*K*y*(x*x - y*y) - 16;

    }

    double fun_u(R2 P, int i, int dom) {
        double x = P.x;
        double y = P.y;
        if (i == 0)
            return 2 * (x * x - x + 0.25 + y * y - y) * (2 * y - 1);
        else
            return -2 * (x * x - x + 0.25 + y * y - y) * (2 * x - 1);
    }

    double fun_p(R2 P, int i, int dom) {
        double x = P.x;
        double y = P.y;
        //return 10 * (x * x - y * y) * (x * x - y * y);
        //return C * (x * x - y * y) * (x * x - y * y);

        double c = - 0.6702064327658224/(0.2*M_PI);     // normalizing constant
        return K*(10 * (x * x - y * y) * (x * x - y * y) + c);
    }

    double fun_zero(R2 P, int i, int dom) { return 0.; }
}

using namespace fictitious;

#define SCOTT_VOGELIUS

int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);

    const size_t mesh_refinements = 5;
    double h = 0.1;     // not used if mesh is imported
    
    std::string home(std::getenv("HOME"));
    std::string path_output_data(home + "/output_files/stokes/fictitious/data/");
    std::string path_output_figures(home + "/output_files/stokes/fictitious/paraview/");

    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_output_figures);
    }

    std::array<double, mesh_refinements> L2_errors_u, L2_errors_p, pwise_div, H1_errors_u, hs, numb_elems;
    std::array<double, mesh_refinements> mean_ps;

    for (int j = 0; j < mesh_refinements; ++j) {

        // Mesh
        // const double lx = 1., ly = 1.;
        // const double x0 = 0., y0 = 0.;
        //int nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
        //mesh_t Th(nx, ny, x0, y0, lx, ly);

        //mesh_t Th(("../cpp/meshes/no_flow_ps_" + std::to_string(j) + ".msh").c_str());
        //mesh_t Th(("../cpp/meshes/unit_square_" + std::to_string(j) + ".msh").c_str());
        mesh_t Th(("../cpp/mainFiles/meshes/no_flow_ct_" + std::to_string(j) + ".msh").c_str());



        //int nx = (int)std::sqrt(Th.nbElements()); 
        //h = std::sqrt(2./Th.nbElements());
        h = Th[0].lenEdge(0);

        //boundary_penalty = boundary_penalty * mu / h;
        //tangential_penalty = tangential_penalty * mu / h;
        
        // P1-interpolated level set function
        space_t Lh(Th, DataFE<mesh_t>::P1);
        fct_t gamma(Lh, fun_levelset);
        InterfaceLevelSet<mesh_t> interface(Th, gamma);

        gnuplot::save(Th, "stokes" + std::to_string(j) + ".dat");
        gnuplot::save(interface, "gamma" + std::to_string(j) + ".dat");

        ActiveMesh<mesh_t> active_mesh(Th); 
        active_mesh.truncate(interface, 1);

        MacroElement<mesh_t> macro(active_mesh, delta);


        // Finite element spaces
        Lagrange2 FE_velocity(2); 
        Lagrange2 FE_velocity_int(4); 
        // space_t V_interpolation(Th, FE_velocity_int);
        // space_t P_interpolation(Th, DataFE<mesh_t>::P4);   
        
        // Scott-Vogelius element pair
    #ifdef SCOTT_VOGELIUS
        space_t V(Th, FE_velocity);
        space_t P(Th, DataFE<mesh_t>::P1dc);
    #elif defined(TAYLOR_HOOD)
        space_t V(Th, FE_velocity);
        space_t P(Th, DataFE<mesh_t>::P1);
    #elif defined(BDM1_P0)
        space_t V(Th, DataFE<mesh_t>::BDM1);
        space_t P(Th, DataFE<mesh_t>::P0);
    #endif

        space_t V_interpolation(Th, FE_velocity);
        space_t P_interpolation(Th, DataFE<mesh_t>::P1dc);   

        // Cut finite element spaces
        cutspace_t Vh(active_mesh, V);
        cutspace_t Ph(active_mesh, P);

        // Create Stokes problem object
        //CutFEM<mesh_t> stokes(Vh, thread_count, optionProblem);
        CutFEM<mesh_t> stokes(Vh);
        stokes.add(Ph);

        std::cout << "--------------------------------------------" << '\n';
        std::cout << "Iteration " << j + 1 << "/" << mesh_refinements << '\n';
        std::cout << "\n h  = " << h << '\n';
        std::cout << "N = " << Th.nbElements() << '\n';

        Normal n;
        Tangent t;

        funtest_t u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);

        // Interpolate exact functions
        fct_t u_exact(V_interpolation, fun_u);
        fct_t p_exact(P_interpolation, fun_p);
        fct_t fh(V_interpolation, fun_rhs);     // rhs force
        fct_t gh(V_interpolation, fun_u);       // Dirichlet boundary condition
        
        // Bilinear form
        stokes.addBilinear(
            //+ contractProduct(mu * grad(u), grad(v))
            + contractProduct(2*mu * Eps(u), Eps(v))
            - innerProduct(p, div(v)) 
            + innerProduct(div(u), q)
            , active_mesh);

        stokes.addBilinear(
            //- innerProduct(mu*grad(u)*n, v)
            - innerProduct(2*mu*Eps(u)*n, v)
            // - innerProduct(u, mu*grad(v)*n) 
            - innerProduct(u, 2*mu*Eps(v)*n) 
            + innerProduct(boundary_penalty*mu/h*u, v)
            + innerProduct(p, v*n)
            //+ innerProduct(u, q*n)      // make problem block-symmetric
            , interface);

        #ifdef BDM1_P0
        stokes.addBilinear(
            - innerProduct(average(mu * grad(u) * t * n, 0.5, 0.5), jump(v * t))
            - innerProduct(jump(u * t), average(mu * grad(v) * t * n, 0.5, 0.5))
            + innerProduct(tangential_penalty * (jump(u * t)), jump(v * t))
            , active_mesh
            , INTEGRAL_INNER_EDGE_2D);
        #endif

        // Linear form
        stokes.addLinear(
            + innerProduct(fh.exprList(), v)
            , active_mesh);

        stokes.addLinear(
            //- innerProduct(gh.exprList(), mu*grad(v)*n) 
            - innerProduct(gh.exprList(), 2*mu*Eps(v)*n) 
            + innerProduct(gh.exprList(), boundary_penalty*mu/h*v)
            //+ innerProduct(gh.exprList(), q*n)  // compensate for added term
            , interface); 

        stokes.addPatchStabilization(
            + innerProduct(tau_u * std::pow(h, -2) * jump(u), jump(v)) 
            - innerProduct(tau_p * std::pow(h, 0) * jump(p), jump(div(v))) 
            + innerProduct(tau_p * std::pow(h, 0) * jump(div(u)), jump(q))
            , active_mesh
            //);
            , macro);



        // funtest_t grad2un = grad(grad(u) * n) * n;
        // stokes.addFaceStabilization(                                  // [h^(2k+1) h^(2k+1)]
        //     + innerProduct(tau_u * std::pow(h, 1) * jump(grad(u) * n), jump(grad(v) * n))
        //     + innerProduct(tau_u * std::pow(h, 3) * jump(grad2un), jump(grad2un))
        //     - innerProduct(tau_p * std::pow(h, 1) * jump(p), jump(div(v)))
        //     + innerProduct(tau_p * std::pow(h, 1) * jump(div(u)), jump(q)) 
        //     - innerProduct(tau_p * std::pow(h, 3) * jump(grad(p)), jump(grad(div(v)))) 
        //     + innerProduct(tau_p * std::pow(h, 3) * jump(grad(div(v))), jump(grad(q)))
        //     , active_mesh
        //     );
        //     //, macro);

        // Add Lagrange multipliers

        // stokes.addLagrangeMultiplier(
        //     + innerProduct(1., p)
        //     , 0.
        //     , active_mesh
        // );

        //double mean_p    = integral(active_mesh, p_exact, 0);

        CutFEM<mesh_t> lagr(Vh);
        lagr.add(Ph);
        lagr.addLinear(innerProduct(1., p), active_mesh);
        std::vector<double> lag_row(lagr.rhs_.begin(), lagr.rhs_.end());
        //lagr.rhs_ = 0.;
        std::fill(lagr.rhs_.begin(), lagr.rhs_.end(), 0.);
        lagr.addLinear(innerProduct(1., v * n), interface);
        //lagr.addLinear(innerProduct(1., p), active_mesh);

        stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhsSpan(), 0.);

        // Export matrix
        matlab::Export(stokes.mat_[0], path_output_data + "mat_" + std::to_string(j + 1) + ".dat");
    
        //stokes.solve(stokes.mat_[0], stokes.rhs_);
        stokes.solve("mumps");

        // Extract numerical solution data
        int nb_flux_dof           = Vh.get_nb_dof();
        std::span<double> data_uh = std::span<double>(stokes.rhs_.data(), nb_flux_dof);
        std::span<double> data_ph = std::span<double>(stokes.rhs_.data() + nb_flux_dof, Ph.get_nb_dof());
        
        // Interpolate numerical solution onto FE spaces
        fct_t uh(Vh, data_uh);
        fct_t ph(Ph, data_ph);

        // Post process pressure
        double meanP    = integral(active_mesh, p_exact, 0);
        mean_ps[j] = meanP;
        std::cout << "Mean exact pressure = " << meanP << '\n';
        double meanPfem = integral(active_mesh, ph.expr(), 0);
        std::cout << "Mean numerical pressure = " << meanPfem << '\n';
        // CutFEM<mesh_t> post(Ph);
        // post.addLinear(innerProduct(1, q), active_mesh);
        // double area = post.rhs_.sum();
        // ph.v -= meanPfem / area;
        // ph.v += meanP / area;

        auto uh_0dx = dx(uh.expr(0));
        auto uh_1dy = dy(uh.expr(1));
        auto uh_0dy = dy(uh.expr(0));
        auto uh_1dx = dx(uh.expr(1));

        auto u_0dx = dx(u_exact.expr(0));
        auto u_1dy = dy(u_exact.expr(1));
        auto u_0dy = dy(u_exact.expr(0));
        auto u_1dx = dx(u_exact.expr(1));

        double error_uh = L2normCut(uh, fun_u, 0, 2);
        double error_ph = L2normCut(ph, fun_p, 0, 1);
        double error_div_uh = maxNormCut(uh_0dx + uh_1dy, active_mesh);
        double error_grad_uh = std::sqrt(integral(active_mesh, (uh_0dx - u_0dx)*(uh_0dx - u_0dx) + (uh_0dy-u_0dy)*(uh_0dy-u_0dy) + (uh_1dx-u_1dx)*(uh_1dx-u_1dx) + (uh_1dy-u_1dy)*(uh_1dy-u_1dy)));

        fct_t u_error(Vh, fun_u);
        fct_t p_error(Ph, fun_p);
        
        u_error.v -= uh.v;
        p_error.v -= ph.v;

        std::cout << " ||u-uh||_2 = " << error_uh << '\n';
        std::cout << " ||p-ph||_2 = " << error_ph << '\n';
        std::cout << " ||grad(u-uh)||_2 = " << error_grad_uh << '\n';
        std::cout << " ||div(uh)||_infty = " << error_div_uh << '\n';


        L2_errors_u[j] = error_uh;
        L2_errors_p[j] = error_ph;
        pwise_div[j] = error_div_uh;
        H1_errors_u[j] = std::sqrt(error_uh*error_uh + error_grad_uh*error_grad_uh);
        hs[j] = h;
        numb_elems[j] = Th.nbElements();

        // Write solutions to Paraview
        Paraview<mesh_t> bg_mesh(Th, path_output_figures + "background_mesh_" + std::to_string(j) + ".vtk");
        Paraview<mesh_t> paraview(active_mesh, path_output_figures + "stokes_" + std::to_string(j) + ".vtk");
        paraview.add(gamma, "gamma", 0, 1);
        paraview.add(uh, "velocity", 0, 2);
        paraview.add(ph, "pressure", 0, 1);
        paraview.add(u_exact, "velocity_exact", 0, 2);
        paraview.add(p_exact, "pressure_exact", 0, 1);
        paraview.add(u_error, "u_error", 0, 2);
        paraview.add(p_error, "p_error", 0, 1);
        paraview.add(fabs(uh_0dx + uh_1dy), "divergence");
        paraview.add(fh, "rhs", 0, 2);
        
        h *= 0.5;

    }

    std::cout << std::setprecision(16);
    std::cout << '\n';
    std::cout << "Errors Velocity = [";
    for (int i = 0; i < mesh_refinements; i++) {

        std::cout << L2_errors_u.at(i);
        if (i < mesh_refinements - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';
    std::cout << "Errors Pressure = [";
    for (int i = 0; i < mesh_refinements; i++) {

        std::cout << L2_errors_p.at(i);
        if (i < mesh_refinements - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "H1 errors = [";
    for (int i = 0; i < mesh_refinements; i++) {

        std::cout << H1_errors_u.at(i);
        if (i < mesh_refinements - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << "Pointwise divergence = [";
    for (int i = 0; i < mesh_refinements; i++) {

        std::cout << pwise_div.at(i);
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


