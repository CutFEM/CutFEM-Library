/*
Here, we're solving the Poisson equation
-Δu = f in Ω
u = g on ∂ΩΩ
using CutFEM with an iterative solver.
The domain Ω is defined by the levelset function
φ(x,y) = sqrt((x-0.5)² + (y-0.5)²) - 0.25
which represents a circle of radius 0.25 centered at (0.5, 0.5).
The exact solution is u(x,y) = 2*sin(2πx)*sin(4πy),
and the corresponding right-hand side f and boundary condition g are derived from this exact solution.


Author: Sebastian Myrbäck, 
*/


#include "../cutfem.hpp"
#include <filesystem>
#include <iostream>
#include <iomanip>   
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>


using namespace globalVariable;

double fun_boundary(double* P, int elementComp) {
    return 2*std::sin(2*M_PI*P[0])*std::sin(4*M_PI*P[1]);
}
double fun_rhs(double* P, int elementComp) {
    return 40*std::pow(M_PI,2)*std::sin(2*M_PI*P[0])*std::sin(4*M_PI*P[1]);
}
double fun_exact(double* P, int elementComp, int domain) {
  return 2*std::sin(2*M_PI*P[0])*std::sin(4*M_PI*P[1]);
}
double fun_levelset(double* P, const int i) {
//   return -(std::sqrt((P[0]-0.5)*(P[0]-0.5) + (P[1]-0.5)*(P[1]-0.5)) - 0.25 - Epsilon);
return -((P[0]-0.5)*(P[0]-0.5) + (P[1]-0.5)*(P[1]-0.5) - 0.25*0.25 - Epsilon);
}


using mesh_t        = Mesh2;
using fespace_t     = GFESpace<mesh_t>;
using cut_fespace_t = CutFESpace<mesh_t>;
using funtest_t     = TestFunction<mesh_t>;
using fct_t         = FunFEM<mesh_t>;
using activemesh_t  = ActiveMesh<mesh_t>;
using interface_t   = InterfaceLevelSet<mesh_t>;
using paraview_t    = Paraview<mesh_t>;
using cutfem_t      = CutFEM<mesh_t>;

int main(int argc, char** argv ) {

    // Initialize MPI
    MPIcf cfMPI(argc,argv);

    ProblemOption option;
    // option.solver_name_ = "mumps";
    // option.solver_name_ = "eigen_cg";
    option.solver_name_ = "cg";
    option.order_space_element_quadrature_ = 5;
    option.it_maxit_ = 2000;
    option.it_tol_ = 1e-5;
    option.it_use_ic_ = false;
    option.verbose_ = 1;

    
    // nx, ny vertices in x and y direction
    int nx = 21, ny = 21;
    
    // side lengths
    const double lx = 1., ly = 1.;
    
    size_t n_refinements = 1;

    std::vector<double> h_values(n_refinements, 0.), L2_errors(n_refinements, 0.);
    std::vector<std::vector<double>> cg_histories(n_refinements);
    for (int i = 0; i < n_refinements; ++i) {

        double h =  lx / (nx - 1);
        std::cout << "h = " << h << std::endl;

        const double mu = 1.;
        double lambda = 40.*mu/h; 
        double ghost_penalty = 0.0002;

    
        mesh_t Th(nx, ny, 0., 0., lx, ly);
        
        fespace_t Lh(Th, DataFE<mesh_t>::P1);   // FE space for the levelset function
        fct_t levelset(Lh, fun_levelset);
        interface_t surface(Th, levelset);

        fespace_t V(Th, DataFE<mesh_t>::P1);   // FE space for u on the background mesh
        
        activemesh_t active_mesh(Th);     // Active mesh for both subdomains (first Omega_0, the outer domain, then Omega_1, the inner domain) 
        active_mesh.truncate(surface, -1);

        paraview_t writer_bg(active_mesh, "paraview_active" + std::to_string(i) + ".vtk");
        writer_bg.writeActiveMesh(active_mesh, "active_mesh.vtk");

        
        
        cut_fespace_t Vh(active_mesh, V);      
        
        // gnuplot::save(Th, "Th.dat");
        // gnuplot::save(interface, "interface.dat");
        // gnuplot::save(Wh,0,"Vh1.dat");
        // gnuplot::save(Wh,1,"Vh2.dat");

        cutfem_t poisson(Vh, option);

        std::cout << "NDOFS = " << poisson.get_nb_dof() << ", # active elements = " << active_mesh.get_nb_element() << "\n";
        Normal n;

        // Interpolate rhs and u on the boundary onto the FE space
        fct_t fh(Vh, fun_rhs);        
        fct_t gh(Vh, fun_boundary);
        funtest_t u(Vh, 1), v(Vh, 1);     // u,v have 1 component

        // assembly
        poisson.addBilinear(
            + innerProduct(mu*grad(u), grad(v))
            , active_mesh
        );
        
        // interface terms
        poisson.addBilinear(
            - innerProduct(mu*grad(u)*n, v)
            - innerProduct(u, mu*grad(v)*n)
            + innerProduct(lambda*u, v)
            , surface
        );

        poisson.addLinear(
            + innerProduct(fh.expr(), v)
            , active_mesh
        );
        
        poisson.addLinear(
            - innerProduct(gh.expr(), mu*grad(v)*n)   // added for symmetry
            + innerProduct(gh.expr(), lambda*v)       // penalty term
            , surface
        );

        // ghost-penalty stabilization
        // poisson.addPatchStabilization( 
        //     + innerProduct(ghost_penalty * mu * std::pow(h, -2) * jump(u), jump(v)) 
        //     , active_mesh
        // );

        poisson.addFaceStabilization(
            + innerProduct(ghost_penalty * mu * std::pow(h, 1) * jump(grad(u)*n), jump(grad(v)*n))
            , active_mesh
            );

        // std::vector<int> labels = {1,2,3,4};
        // BoundaryDirichlet bc(Vh, labels);
        // bc.apply_inhomogeneous(poisson.mat_);
        // bc.apply(poisson.rhs_, gh);
        // bc.finalize_inhomogeneous(poisson.mat_, poisson.rhs_);

        // auto dof_index_data  = poisson.get_dof_index_data(Vh, active_mesh);
        // auto dof_vertex_data = poisson.get_dof_vertex_data(Vh, active_mesh);

        // export matrix BEFORE solve (afterwards it gets reset)
        matlab::Export(poisson.mat_, "mat_gp1.dat");
        // matlab::Export(dof_index_data, "dof_index_data.dat");
        // matlab::Export(dof_vertex_data, "dof_vertex_data.dat");
        // poisson.solve("mumps");
        poisson.solve();
        cg_histories[i] = poisson.residual_history_;

        fct_t uh(Vh,  poisson.rhs_);    // solution to the system goes into poisson.rhs_
        double L2_error = L2normCut(uh, fun_exact, 0, 1);
        fct_t u_error(Vh, fun_exact);
        std::transform(u_error.v.begin(), u_error.v.end(), uh.v.begin(), u_error.v.begin(), std::minus<double>()); 
        
        std::cout << "L2 error = " << L2_error << "\n";

        // PRINT THE SOLUTION TO PARAVIEW
        // =====================================================
        paraview_t writer(active_mesh, "poisson_fictitious" + std::to_string(i) + ".vtk");
        writer.add(uh, "uh", 0, 1);
        writer.add(gh, "u_exact", 0, 1);
        writer.add(u_error, "u_error", 0, 1);
        writer.add(levelset, "levelset", 0, 1);

        nx = 2*nx - 1;
        ny = 2*ny - 1;

        h_values[i]  = h;
        L2_errors[i] = L2_error;
    }


    std::cout << "\nConvergence results:\n";
    std::cout << std::left
            << std::setw(12) << "Refinement"
            << std::setw(14) << "h"
            << std::setw(14) << "L2 Error"
            << std::setw(10) << "Rate"
            << "\n";

    for (size_t i = 0; i < n_refinements; ++i) {
        double rate = 0.0;
        if (i > 0) {
            rate = std::log(L2_errors[i-1]/L2_errors[i]) /
                std::log(h_values[i-1]/h_values[i]);
        }

        std::cout << std::left
                << std::setw(12) << i
                << std::setw(14) << std::setprecision(8)  << std::fixed << h_values[i]
                << std::setw(14) << std::setprecision(8)  << std::fixed << L2_errors[i]
                << std::setw(10) << std::setprecision(5)  << std::fixed << rate
                << "\n";
    }


    std::cout << "\nCG residual history (relative residual per iteration):\n";
    for (size_t i = 0; i < n_refinements; ++i) {
        std::cout << "  Refinement " << i << " (h=" << std::setprecision(4) << h_values[i]
                  << ", " << cg_histories[i].size() << " iters):\n";
        for (size_t k = 0; k < cg_histories[i].size(); ++k)
            std::cout << "    " << std::setw(4) << k
                      << "  " << std::setprecision(6) << std::scientific << cg_histories[i][k] << "\n";
    }

    return 0;
}
