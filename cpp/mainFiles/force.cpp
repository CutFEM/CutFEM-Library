#include "../tool.hpp"
#include <random>
#include <chrono>

#include "../force_project/util.hpp"
#include "../force_project/input_handler.hpp"
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

namespace force_project {
struct MeshDimensions {
    double x_min{-1.}, lx{2.};
    double y_min{-1.}, ly{2.};
    double z_min{-1.}, lz{2.};
};

struct BaseLevelSetStruct {
    BaseLevelSetStruct() {}
    virtual double operator()(R2 P) const = 0;

    virtual ~BaseLevelSetStruct() {}
};

struct LevelSetCircle : public BaseLevelSetStruct {

    LevelSetCircle(double xc, double yc, double r) : xc(xc), yc(yc), r(r) {}

    LevelSetCircle(LevelSetCircle &&)            = default;
    LevelSetCircle &operator=(LevelSetCircle &&) = default;

    double operator()(R2 P) const override { return sqrt((P.x - xc) * (P.x - xc) + (P.y - yc) * (P.y - yc)) - r; }

    double xc, yc, r;
};

struct MultiLevelSet {

    MultiLevelSet(std::vector<std::shared_ptr<LevelSetCircle>> levelSets) : levelSets(levelSets) {}

    double operator()(R2 P) const {
        double min = levelSets[0]->operator()(P);
        for (int i = 1; i < levelSets.size(); i++) {
            double val = levelSets[i]->operator()(P);
            if (val < min)
                min = val;
        }
        return min;
    }

    std::vector<std::shared_ptr<LevelSetCircle>> levelSets;
};

struct GeometryData {

    GeometryData(const YamlReader &input_data) {}

    MeshDimensions box;

    int nx, ny;
};

std::vector<std::shared_ptr<LevelSetCircle>> generateRandomLevelSets(double xc_min, double xc_max, int n_ls,
                                                                     double delta_l, double r_0) {

    std::vector<std::shared_ptr<LevelSetCircle>> levelSets_v;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(xc_min, xc_max);

    std::vector<double> xc_v;
    std::vector<double> yc_v;

    if (MPIcf::IamMaster()) {
        for (int k = 0; k < n_ls; ++k) {
            int count = 0;
            while (true) {
                double xc         = distribution(generator);
                double yc         = distribution(generator);
                bool is_satifying = checkIfValid(xc, yc, delta_l, levelSets_v);

                if (is_satifying) {
                    xc_v.push_back(xc);
                    yc_v.push_back(yc);
                    levelSets_v.push_back(std::make_shared<LevelSetCircle>(xc, yc, r_0));
                    break;
                }
                count += 1;
                if (count > 100) {
                    LOG_INFO << "Could not find a valid point" << logger::endl;
                    break;
                }
            }

            // double xc = 0.5 * std::cos(2. * k * M_PI / n_ls);
            // double yc = 0.5 * std::sin(2. * k * M_PI / n_ls);
        }
    }
    // broadcast results
    int n = xc_v.size();
    MPIcf::Bcast(n, MPIcf::Master(), 1);
    if (!MPIcf::IamMaster()) {
        xc_v.resize(n);
        yc_v.resize(n);
    }
    MPIcf::Bcast<double>(xc_v, MPIcf::Master());
    MPIcf::Bcast<double>(yc_v, MPIcf::Master());

    if (!MPIcf::IamMaster()) {
        for (int i = 0; i < xc_v.size(); ++i) {
            LOG_INFO << " Particle " << i << " - Center (" << xc_v[i] << ", " << yc_v[i] << ") , radius = " << r_0
                     << logger::endl;

            levelSets_v.push_back(std::make_shared<LevelSetCircle>(xc_v[i], yc_v[i], r_0));
        }
    }
    LOG_INFO << "Number of particles built: " << levelSets_v.size() << logger::endl;
    return levelSets_v;
}

struct DependencyParameter {

    double ratio_mu{10.};
    double radius_particle{0.5};
};

R2 computeForce(const fct_t &uh, const gamma_t &gamma) {

    Normal n;
    auto u_n = (uh * n);

    double fx = integral(u_n * n.x, gamma) / gamma.measure();
    double fy = integral(u_n * n.y, gamma) / gamma.measure();

    // double fx = integral(uh.expr(0) * n.x, gamma) / gamma.measure();
    // double fy = integral(uh.expr(1) * n.y, gamma) / gamma.measure();

    return R2(fx, fy);
}
R2 computeForce(const fct_t &uh, const cutmesh_t &Khi) {

    double area_bubble = integral(Khi, 1, 1, 0);

    double fx = integral(Khi, uh.expr(0), 1, 0) / area_bubble;
    double fy = integral(Khi, uh.expr(1), 1, 0) / area_bubble;
    return R2(fx, fy);
}

} // namespace force_project

using namespace force_project;

int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);

    double t0 = MPIcf::Wtime();

    CutFEMLogger::initialize("../../libcutfem/cpp/force_project/data/force_log.txt");

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
    int nx = 50, ny = 50;
    int n_ls        = 1;
    double r_0      = 0.25;
    double delta_bc = 2. * box.lx / (nx - 1);
    double delta_l  = r_0 + 5. * delta_bc;
    double xc_min   = box.x_min + r_0 + delta_bc;
    double xc_max   = box.x_min + box.lx - r_0 - delta_bc;
    int Nx          = 3;
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
    try {
        csvfile csv("../../libcutfem/cpp/force_project/data/test.csv", ",", std::ofstream::app);

        int ifig = 0;
        for (int ix = 0; ix < Nx; ++ix) {
            xc = xc_min + ix * hx;
            for (int jx = 0; jx < Nx; ++jx) {
                yc = xc_min + jx * hx;

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
    double t1 = MPIcf::Wtime();
    LOG_INFO << "Total time: " << t1 - t0 << logger::endl;
}