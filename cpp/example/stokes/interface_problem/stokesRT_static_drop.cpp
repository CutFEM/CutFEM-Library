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

#include "../../../tool.hpp"
#include "../../../matplotlib-cpp/matplotlibcpp.h"

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
    double rad   = 0.6; // Radius of the drop
    double sigma = 5.;
    // Viscosity in both subdomains
    CutFEMParameter mu(1., 1.);

    double delta = 0.25;

    auto p_exact = [sigma, rad](R2 P, int component, int domain) { return (domain == 0) ? 0 : sigma / (rad); };
    auto Gamma   = [rad](R2 P, int i) { return sqrt(P[0] * P[0] + P[1] * P[1]) - rad; };

    int n_iter = 1;
    int nx     = 20;
    std::cout << " How many iteration? " << std::endl;
    std::cin >> n_iter;
    std::cout << " Save solution in vtk file for paraview? [0/1]?" << std::endl;
    int save_vtk = 0;
    std::cin >> save_vtk;

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
        fct_t gh(Vh), fh(Vh);

        // Solve Stokes interface problem
        auto data_stokes = stokesSolver(Vh, Ph, interface, gh, fh, mu, sigma / rad, delta);

        // Extract solution
        std::span<double> data_uh{std::span(data_stokes.data(), Vh.get_nb_dof())};
        std::span<double> data_ph{std::span(data_stokes.data() + Vh.get_nb_dof(), Ph.get_nb_dof())};
        fct_t uh(Vh, data_uh);
        fct_t ph(Ph, data_ph);

        // Plot the solution
        if (MPIcf::IamMaster() && save_vtk) {
            paraview_t writer(Khi, "stokes_interface_static_drop_" + std::to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
        }

        auto dxu = dx(uh.expr(0));
        auto dyu = dy(uh.expr(1));

        double errorL2_u       = L2normCut(uh, 0, 2);
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

    plt::figure();
    plt::title("Stokes: Static drop");
    plt::plot(h, error_p, {{"label", "$|| p_h - p_{ex}||_{L^2}$"}, {"marker", "+"}});
    plt::plot(h, error_u, {{"label", "$|| u_h - u_{ex}||_{L^2}$"}, {"marker", "*"}});
    plt::plot(h, error_div, {{"label", "$|| div(u_h)||_{L^2}$"}, {"marker", "o"}});
    plt::plot(h, max_div, {{"label", "max($div(u_h)$)"}, {"marker", "d"}});
    plt::xlabel("mesh size h");
    plt::legend();
    plt::tight_layout();
    plt::show();
}
