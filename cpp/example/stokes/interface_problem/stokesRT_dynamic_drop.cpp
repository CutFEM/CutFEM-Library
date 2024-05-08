/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/

#include "../cutfem.hpp"
#include "matplotlibcpp.h"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;
using gamma_t    = InterfaceLevelSet<mesh_t>;
using paraview_t = Paraview<mesh_t>;

namespace plt = matplotlibcpp;

/// @brief Example of Stokes interface problem with zero velocity
/// @details We just want to illustrate that one get exact pressure jump
int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);

    // Paramers in problem
    constexpr double mu1   = 100.;
    constexpr double mu2   = 1.;
    constexpr double rad   = 2. / 3; // Radius of the drop
    constexpr double sigma = 10.;
    CutFEMParameter mu(mu1, mu2);

    double delta = 0.25;

    auto rhs = [](R2 P, int i) {
        double r2 = Norme2_2(P);
        double s  = std::exp(-r2);
        R2 R(4 * s * (r2 - 2) * P[1] + 3 * P[0] * P[0], -4 * s * (r2 - 2) * P[0]);
        return R[i];
    };
    auto falpha1 = [rad, mu1, mu2](double r) { return 1. / mu1 + (1. / mu2 - 1. / mu1) * std::exp(r * r - rad * rad); };
    auto falpha2 = [mu2](double r) { return 1. / mu2; };
    auto falpha  = [rad, &falpha1, &falpha2](double r) { return (r < rad) ? falpha2(r) : falpha1(r); };

    auto fun_velocity1 = [&falpha1](R2 P) {
        double r = Norme2(P);
        return falpha1(r) * exp(-r * r) * R2(-P[1], P[0]);
    };
    auto fun_velocity2 = [&falpha2](R2 P) {
        double r = Norme2(P);
        R2 R(-P[1], P[0]);
        return falpha2(r) * exp(-r * r) * R2(-P[1], P[0]);
    };
    auto fun_pressure1 = [](R2 P) { return pow(P[0], 3); };
    auto fun_pressure2 = [sigma](R2 P) { return pow(P[0], 3) + sigma * 3. / 2.; };

    auto u_exact = [&fun_velocity1, &fun_velocity2](R2 P, int ci, int domain) {
        return (domain == 0) ? fun_velocity1(P)[ci] : fun_velocity2(P)[ci];
    };
    auto p_exact = [&fun_pressure1, &fun_pressure2](R2 P, int ci, int domain) {
        return (domain == 0) ? fun_pressure1(P) : fun_pressure2(P);
    };
    auto Gamma = [rad](R2 P, int i) { return sqrt(P[0] * P[0] + P[1] * P[1]) - rad; };

    int n_iter      = 1;
    int nx          = 20;
    int num_threads = 1;
    int save_vtk    = 0;
    if (MPIcf::IamMaster()) {
        std::cout << " How many iteration? " << std::endl;
        std::cin >> n_iter;
        std::cout << " Save solution in vtk file for paraview? [0/1]?" << std::endl;
        std::cin >> save_vtk;
        std::cout << " How many threads? " << std::endl;
        std::cin >> num_threads;
    }
    MPIcf::Bcast(n_iter, MPIcf::Master(), 1);
    MPIcf::Bcast(num_threads, MPIcf::Master(), 1);

    omp_set_num_threads(num_threads);
    std::cout << "Number of threads set: " << num_threads << std::endl;
    std::cout << "OMP Number of threads: " << omp_get_max_threads() << std::endl;
    assert(num_threads > 0 && num_threads == omp_get_max_threads());

    auto t_start            = std::chrono::high_resolution_clock::now();
    globalVariable::verbose = 2;

    std::vector<double> h;
    std::vector<double> error_u;
    std::vector<double> error_p;
    std::vector<double> error_div;
    std::vector<double> max_div;

    for (int i = 0; i < n_iter; ++i) {
        mesh_t Kh(nx, nx, -1., -1., 2., 2.);
        Kh.info();

        // Finite element space for the levelset function
        space_t Lh(Kh, DataFE<Mesh2>::P1);
        // Interpolation of the levelset function
        fct_t levelSet(Lh, Gamma);
        // Interface
        gamma_t interface(Kh, levelSet);
        // Active mesh
        cutmesh_t Khi(Kh, interface);

        // Build spaces and cut spaces
        space_t Wh(Kh, DataFE<Mesh2>::BDM1);
        space_t Qh(Kh, DataFE<Mesh2>::P0);
        cutspace_t Vh(Khi, Wh);
        cutspace_t Ph(Khi, Qh);

        // Zero right hand side and velocity on the boundary
        fct_t gh(Vh, u_exact), fh(Vh, rhs);

        // Solve Stokes interface problem
        auto data_stokes = solver::cutfem::stokes::solve(Vh, Ph, interface, gh, fh, mu, sigma / rad, delta);

        // Extract solution
        std::span<double> data_uh{std::span(data_stokes.data(), Vh.get_nb_dof())};
        std::span<double> data_ph{std::span(data_stokes.data() + Vh.get_nb_dof(), Ph.get_nb_dof())};
        fct_t uh(Vh, data_uh);
        fct_t ph(Ph, data_ph);

        // Plot the solution
        if (MPIcf::IamMaster() && save_vtk) {
            paraview_t writer(Khi, "stokes_interface_dynamic_drop_" + std::to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
        }

        auto dxu = dx(uh.expr(0));
        auto dyu = dy(uh.expr(1));

        double errorL2_u       = L2normCut(uh, u_exact, 0, 2);
        double errorL2_p       = L2normCut(ph, p_exact, 0, 1);
        double errorL2_div     = L2normCut(fabs(dxu + dyu), Khi);
        double errorLinfty_div = maxNormCut(dxu + dyu, Khi);

        h.push_back(Khi[0].hElement());
        error_u.push_back(errorL2_u);
        error_p.push_back(errorL2_p);
        error_div.push_back(errorL2_div);
        max_div.push_back(errorLinfty_div);

        nx = 2 * nx - 1;
    }

    LOG_INFO << "L2 error u : \n " << error_u << logger::endl;
    LOG_INFO << "L2 error p : \n " << error_p << logger::endl;
    LOG_INFO << "L2 error div : \n " << error_div << logger::endl;
    LOG_INFO << "Max error div : \n " << max_div << logger::endl;

    auto t_end  = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(t_end - t_start).count();
    LOG_INFO << "Time: " << time << " sec." << logger::endl;

    if (MPIcf::IamMaster()) {

        std::vector<double> h2(h);
        std::transform(h.begin(), h.end(), h2.begin(), [](double x) { return x * x; });

        plt::figure();
        plt::title("Stokes: Static drop - Pressure");
        plt::loglog(h, error_p, {{"label", "$|| p_h - p_{ex}||_{L^2}$"}, {"marker", "+"}});
        plt::loglog(h, h, "k--", {{"label", "$h$"}});
        plt::legend();
        plt::tight_layout();

        std::transform(h.begin(), h.end(), h2.begin(), [](double x) { return 1e-1 * x * x; });
        plt::figure();
        plt::title("Stokes: Static drop - Velocity");
        plt::loglog(h, error_u, {{"label", "$|| u_h - u_{ex}||_{L^2}$"}, {"marker", "*"}});
        plt::loglog(h, h2, "k--", {{"label", "$h^2$"}});
        plt::xlabel("mesh size h");
        plt::legend();
        plt::tight_layout();

        plt::figure();
        plt::title("Stokes: Static drop - Divergence");
        plt::loglog(h, error_div, {{"label", "$|| div(u_h)||_{L^2}$"}, {"marker", "o"}});
        plt::loglog(h, max_div, {{"label", "max($div(u_h)$)"}, {"marker", "d"}});
        plt::xlabel("mesh size h");
        plt::legend();
        plt::tight_layout();
        plt::show();
    }
}