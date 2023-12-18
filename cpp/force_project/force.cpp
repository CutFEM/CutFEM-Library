


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

    double t0 = CPUtime();
    auto t0_real_time = std::chrono::high_resolution_clock::now(); 

    CutFEMLogger::initialize("../cpp/force_project/data/force_log.txt");

    // 0) Read data
    // if (argc < 2) {
    //     std::cout << "Need an input data file" << std::endl;
    //     exit(EXIT_FAILURE);
    // }

    // std::string data_path = "../../libcutfem/cpp/force_project/data/" + std::string(argv[1]) + ".yaml";

    // YamlReader yaml_reader(data_path);
    // GeometryData geometry_data(yaml_reader);

    // some parameters for the level set functions
    MeshDimensions box;
    int nx = 100, ny = 100;
    int n_ls        = 1;
    double r_0      = 0.25;
    double delta_bc = 4. * box.lx / (nx - 1);
    double delta_l  = r_0 + 5. * delta_bc;
    double xc_min   = box.x_min + r_0 + delta_bc;
    double xc_max   = box.x_min + box.lx - r_0 - delta_bc;
    int Nx          = 2;
    double Lx       = xc_max - xc_min;
    double hx       = Lx / (Nx - 1);

    LOG_INFO << "delta_bc = " << delta_bc << logger::endl;
    LOG_INFO << "delta_l = " << delta_l << logger::endl;
    LOG_INFO << "xc_min = " << xc_min << logger::endl;
    LOG_INFO << "xc_max = " << xc_max << logger::endl;

    // 0b) Generate the level sets functions
    // std::vector<std::shared_ptr<LevelSetCircle>> levelSets_v =
    //     generateRandomLevelSets(xc_min, xc_max, n_ls, delta_l, r_0);

    // //  std::vector<std::shared_ptr<BaseLevelSetStruct>> levelSets_v;
    // // levelSets_v.push_back(std::make_shared<LevelSetCircle>(0., 0., r_0));

    // MultiLevelSet levelSets(levelSets_v);

    // 1) Generate geometry
    mesh_t Kh(nx, ny, box.x_min, box.y_min, box.lx, box.ly);

    // LevelSet
    space_t Lh(Kh, DataFE<Mesh2>::P1);

    double xc;
    double yc;
    std::cout << " Creating data for force problem" << std::endl;
     std::cout << " ";
    progress bar (" CPU time producing data", Nx*Nx, 1);
    try {
        csvfile csv("../cpp/force_project/data/test.csv", ",", std::ofstream::app);

        int ifig = 0, iter = 0;
        for (int ix = 0; ix < Nx; ++ix) {
            xc = xc_min + ix * hx;
            for (int jx = 0; jx < Nx; ++jx) {
                yc = xc_min + jx * hx;
                bar++;

                std::cout << "\r";;
                std::cout << "run : " << ++iter << " / " << Nx*Nx;
                std::cout.flush();

                LOG_INFO << "xc = " << xc << logger::endl;
                LOG_INFO << "yc = " << yc << logger::endl;

                auto ls = [xc, yc, r_0](R2 P) { return sqrt((P.x - xc) * (P.x - xc) + (P.y - yc) * (P.y - yc)) - r_0; };
                // fct_t levelSet(Lh, levelSets);
                fct_t levelSet(Lh, ls);
                // Interface and active mesh
                gamma_t interface(Kh, levelSet);
                cutmesh_t Khi(Kh, interface);

                // cutmesh_t Khi(Kh);
                // Khi.truncate(interface, 1);
                // cutspace_t Lh_cut(Khi, Lh);
                // fct_t signLevelSet(Lh_cut, sign_ls);
                // paraview_t writer(Khi, "Erik_geometry.vtk");

                // return 0;

                // 2a) Build spaces and cut spaces
                space_t Wh(Kh, DataFE<Mesh2>::BDM1);
                space_t Qh(Kh, DataFE<Mesh2>::P0);
                cutspace_t Vh(Khi, Wh);
                cutspace_t Ph(Khi, Qh);

                // Zero right hand side and velocity on the boundary
                auto f_bc = [](R2 P, int i, int dom) -> double { return (i == 0) * P.x - (i == 1) * P.y; };
                fct_t gh(Vh, f_bc), fh(Vh);
                constexpr double sigma = 10.;
                CutFEMParameter mu(1., 10.);
                double delta = 1.;

                // // 2) solve Stokes
                auto data_stokes = stokesSolver(Vh, Ph, interface, gh, fh, mu, sigma / r_0, delta);

                // Extract solution
                std::span<double> data_uh{std::span(data_stokes.data(), Vh.get_nb_dof())};
                std::span<double> data_ph{std::span(data_stokes.data() + Vh.get_nb_dof(), Ph.get_nb_dof())};
                fct_t uh(Vh, data_uh);
                fct_t ph(Ph, data_ph);

                auto uh_0dx      = dx(uh.expr(0));
                auto uh_1dy      = dy(uh.expr(1));
                double errDiv    = L2normCut(uh_0dx + uh_1dy, Khi);
                double maxErrDiv = maxNormCut(uh_0dx + uh_1dy, Khi);

                LOG_INFO << "L2 error div = " << errDiv << logger::endl;
                LOG_INFO << "Max error div = " << maxErrDiv << logger::endl;

                // Plot the solution
                // if (MPIcf::IamMaster()) {
                //     cutspace_t Lh_cut(Khi, Lh);
                //     // auto sign_ls = [](R2 P, int i, int d) { return 1. - 2 * d; };
                //     // fct_t signLevelSet(Lh_cut, sign_ls);

                //     paraview_t writer(Khi, "force_geometry" + std::to_string(ifig++) + ".vtk");
                //     // writer.add(signLevelSet, "signLevelSet", 0, 1);
                //     writer.add(uh, "velocity", 0, 2);
                //     writer.add(ph, "pressure", 0, 1);
                // }

                // 3) compute force
                // compute the mean of u.n on gamma
                // R2 Fc = computeForce(uh, interface);
                // LOG_INFO << "Average interface: Fc = ( " << Fc.x << " , " << Fc.y << " )" << logger::endl;
                R2 Fc = computeForce(uh, Khi);
                LOG_INFO << "Average bubble: Fc = ( " << Fc.x << " , " << Fc.y << " )" << logger::endl;
                // 4) save force

                // csv << levelSets_v[0]->xc << levelSets_v[0]->yc << Fc.x << Fc.y << endrow;
                csv << xc << yc << Fc.x << Fc.y << endrow;
            }
        }
    } catch (const std::exception &ex) {
        std::cout << "Exception was thrown: " << ex.what() << std::endl;
    }
    bar.end();
    double t1 = CPUtime();
    auto t1_real_time = std::chrono::high_resolution_clock::now(); 
    LOG_INFO << "CPU time: " << t1 - t0 << logger::endl;
    LOG_INFO << "Wall time : " << std::chrono::duration<double>(t1_real_time - t0_real_time).count() << " ms" << logger::endl;
}