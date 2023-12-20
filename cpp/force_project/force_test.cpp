#include "../force_project/util.hpp"
#include "../force_project/csvfile.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;
using gamma_t    = InterfaceLevelSet<mesh_t>;
using paraview_t = Paraview<mesh_t>;

using namespace force_project;

void run_advection(const mesh_t &Kh, double r_0, double xc, double yc, int number_of_time_step, double dt) {
    space_t Lh(Kh, DataFE<Mesh2>::P1);
    Lagrange2 FeU(1);
    space_t Uh(Kh, FeU);

    auto ls = [xc, yc, r_0](R2 P) { return sqrt((P.x - xc) * (P.x - xc) + (P.y - yc) * (P.y - yc)) - r_0; };
    fct_t levelSet(Lh, ls);
    fct_t velocity_field_0(Uh);
    fct_t velocity_field_1(Uh);

    // Zero right hand side and velocity on the boundary
    auto f_bc              = [](R2 P, int i, int dom) -> double { return (i == 0) * P.x - (i == 1) * P.y; };
    constexpr double sigma = 100.;
    CutFEMParameter mu(1., 100.);
    double delta = 1.;

    int ifig = 0;
    for (int i = 0; i < number_of_time_step; ++i) {

        std::cout << "\r";
        std::cout << "run advection: " << i + 1 << " / " << number_of_time_step;
        std::cout.flush();

        // Interface and active mesh
        gamma_t gamma(Kh, levelSet);
        cutmesh_t Khi(Kh, gamma);

        // 2a) Build spaces and cut spaces
        space_t Wh(Kh, DataFE<Mesh2>::BDM1);
        space_t Qh(Kh, DataFE<Mesh2>::P0);
        cutspace_t Vh(Khi, Wh);
        cutspace_t Ph(Khi, Qh);

        fct_t gh(Vh, f_bc), fh(Vh);

        // 2a) Compute mean curvature
        cutmesh_t Sh(Kh);
        Sh.createSurfaceMesh(gamma);
        cutspace_t cutVh(Sh, Uh);
        std::vector<double> data_H = solver::cutfem::curvature::solve(cutVh, gamma);
        fct_t H(cutVh, data_H);

        //  2) solve Stokes
        auto data_stokes = solver::cutfem::stokes::solve(Vh, Ph, gamma, gh, fh, mu, sigma, H, delta);

        // Extract solution
        std::span<double> data_uh{std::span(data_stokes.data(), Vh.get_nb_dof())};
        std::span<double> data_ph{std::span(data_stokes.data() + Vh.get_nb_dof(), Ph.get_nb_dof())};
        fct_t uh(Vh, data_uh);
        fct_t ph(Ph, data_ph);

        // Get velocity on the background mesh
        interpolateOnBackGroundMesh(velocity_field_1, uh, levelSet);

        if (i == 0) {
            velocity_field_0.v = velocity_field_1.v;
        }

        std::vector<double> data_ls =
            solver::fem::advection::streamline_diffusion::solve(levelSet, velocity_field_0, velocity_field_1, dt);

        double circularity = Circularity(Khi, gamma);

        LOG_INFO << " Circularity : " << circularity << logger::endl;

        levelSet.v         = data_ls;
        velocity_field_0.v = velocity_field_1.v;

        if (i % 2 == 0) {
            paraview_t writer(Kh, "advection_geometry" + std::to_string(ifig++) + ".vtk");
            writer.add(levelSet, "levelSet1", 0, 1);
            writer.add(velocity_field_1, "uh", 0, 2);
        }
    }
}

void run_force(const mesh_t &Kh, double r_0, double xc, double yc, int number_of_time_step, double dt) {

    space_t Lh(Kh, DataFE<Mesh2>::P1);
    Lagrange2 FeU(1);
    space_t Uh(Kh, FeU);

    // Zero right hand side and velocity on the boundary
    auto f_bc              = [](R2 P, int i, int dom) -> double { return (i == 0) * P.x - (i == 1) * P.y; };
    constexpr double sigma = 100.;
    CutFEMParameter mu(1., 100.);
    double delta = 1.;

    int ifig = 0;
    for (int i = 0; i < number_of_time_step; ++i) {

        std::cout << "\r";
        std::cout << "run force: " << i + 1 << " / " << number_of_time_step;
        std::cout.flush();

        auto ls = [xc, yc, r_0](R2 P) { return sqrt((P.x - xc) * (P.x - xc) + (P.y - yc) * (P.y - yc)) - r_0; };
        fct_t levelSet(Lh, ls);

        // Interface and active mesh
        gamma_t gamma(Kh, levelSet);
        cutmesh_t Khi(Kh, gamma);

        // 2a) Build spaces and cut spaces
        space_t Wh(Kh, DataFE<Mesh2>::BDM1);
        space_t Qh(Kh, DataFE<Mesh2>::P0);
        cutspace_t Vh(Khi, Wh);
        cutspace_t Ph(Khi, Qh);

        fct_t gh(Vh, f_bc), fh(Vh);

        //  2) solve Stokes
        auto data_stokes = solver::cutfem::stokes::solve(Vh, Ph, gamma, gh, fh, mu, sigma / r_0, delta);

        // Extract solution
        std::span<double> data_uh{std::span(data_stokes.data(), Vh.get_nb_dof())};
        fct_t uh(Vh, data_uh);

        R2 Fc = computeForce(uh, Khi);
        LOG_INFO << "Average bubble: Fc = ( " << Fc.x << " , " << Fc.y << " )" << logger::endl;

        xc = xc + dt * Fc.x;
        yc = yc + dt * Fc.y;
        if (i % 2 == 0) {
            paraview_t writer(Kh, "force_geometry" + std::to_string(ifig++) + ".vtk");
            writer.add(levelSet, "levelSet0", 0, 1);
        }
    }
}

int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);
    CutFEMLogger::initialize("../cpp/force_project/data/force_test_log.txt");

    MeshDimensions box;
    int nx = 50, ny = 50;
    int number_of_time_step = 120;

    double r_0 = 0.25;
    double xc  = 0.1;
    double yc  = 0.5;
    double dt  = 0.01;

    // Define Mesh
    mesh_t Kh(nx, ny, box.x_min, box.y_min, box.lx, box.ly);

    run_advection(Kh, r_0, xc, yc, number_of_time_step, dt);

    // run_force(Kh, r_0, xc, yc, number_of_time_step, dt);
}