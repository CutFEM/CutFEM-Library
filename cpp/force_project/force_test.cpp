




#include "../force_project/util.hpp"
#include "../force_project/csvfile.hpp"
#include "../example/stokes/interface_problem/stokes_solver.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;
using gamma_t    = InterfaceLevelSet<mesh_t>;
using paraview_t = Paraview<mesh_t>;

using namespace force_project;

int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);

    MeshDimensions box;
    int nx = 50, ny = 50;
    double r_0      = 0.25;
    double xc = 0.5;
    double yc = 0.5;


    // Define Mesh and Level Set Space
    mesh_t Kh(nx, ny, box.x_min, box.y_min, box.lx, box.ly);
    space_t Lh(Kh, DataFE<Mesh2>::P1);
    Lagrange2 FeU(1);
    space_t Uh(Kh, FeU);

    auto ls = [xc, yc, r_0](R2 P) { return sqrt((P.x - xc) * (P.x - xc) + (P.y - yc) * (P.y - yc)) - r_0; };
    fct_t levelSet(Lh, ls);
    fct_t velocity_field(Uh);


    // Zero right hand side and velocity on the boundary
    auto f_bc = [](R2 P, int i, int dom) -> double { return (i == 0) * P.x - (i == 1) * P.y; };
    constexpr double sigma = 10.;
    CutFEMParameter mu(1., 10.);
    double delta = 1.;
    double dt = 0.1;


    for( int i = 0; i< 1;++i) {

    // Interface and active mesh
    gamma_t interface(Kh, levelSet);
    cutmesh_t Khi(Kh, interface);

    // 2a) Build spaces and cut spaces
    space_t Wh(Kh, DataFE<Mesh2>::BDM1);
    space_t Qh(Kh, DataFE<Mesh2>::P0);
    cutspace_t Vh(Khi, Wh);
    cutspace_t Ph(Khi, Qh);

    fct_t gh(Vh, f_bc), fh(Vh);

    //  2) solve Stokes
    auto data_stokes = stokesSolver(Vh, Ph, interface, gh, fh, mu, sigma / r_0, delta);

    // Extract solution
    std::span<double> data_uh{std::span(data_stokes.data(), Vh.get_nb_dof())};
    std::span<double> data_ph{std::span(data_stokes.data() + Vh.get_nb_dof(), Ph.get_nb_dof())};
    fct_t uh(Vh, data_uh);
    fct_t ph(Ph, data_ph);


    // Get velocity on the background mesh
    interpolateOnBackGroundMesh(velocity_field, uh, levelSet);

    int ifig = 0;     
    paraview_t writer(Kh, "force_geometry" + std::to_string(ifig++) + ".vtk");
    writer.add(levelSet, "levelSet", 0,1);
    writer.add(velocity_field, "uh", 0,2);

    // levelSet = LevelSet::move(levelSet, uh, dt);




    }









}