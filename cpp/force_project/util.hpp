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

#ifndef FORCE_PROJECT_UTIL_HPP
#define FORCE_PROJECT_UTIL_HPP

#include "../tool.hpp"
#include <random>
#include <chrono>

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;
using gamma_t    = InterfaceLevelSet<mesh_t>;

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

// struct GeometryData {

//     GeometryData(const YamlReader &input_data) {}

//     MeshDimensions box;

//     int nx, ny;
// };

// std::vector<std::shared_ptr<LevelSetCircle>> generateRandomLevelSets(double xc_min, double xc_max, int n_ls,
//                                                                      double delta_l, double r_0) {

//     std::vector<std::shared_ptr<LevelSetCircle>> levelSets_v;

//     unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//     std::default_random_engine generator(seed);
//     std::uniform_real_distribution<double> distribution(xc_min, xc_max);

//     std::vector<double> xc_v;
//     std::vector<double> yc_v;

//     if (MPIcf::IamMaster()) {
//         for (int k = 0; k < n_ls; ++k) {
//             int count = 0;
//             while (true) {
//                 double xc         = distribution(generator);
//                 double yc         = distribution(generator);
//                 bool is_satifying = checkIfValid(xc, yc, delta_l, levelSets_v);

//                 if (is_satifying) {
//                     xc_v.push_back(xc);
//                     yc_v.push_back(yc);
//                     levelSets_v.push_back(std::make_shared<LevelSetCircle>(xc, yc, r_0));
//                     break;
//                 }
//                 count += 1;
//                 if (count > 100) {
//                     LOG_INFO << "Could not find a valid point" << logger::endl;
//                     break;
//                 }
//             }

//             // double xc = 0.5 * std::cos(2. * k * M_PI / n_ls);
//             // double yc = 0.5 * std::sin(2. * k * M_PI / n_ls);
//         }
//     }
//     // broadcast results
//     int n = xc_v.size();
//     MPIcf::Bcast(n, MPIcf::Master(), 1);
//     if (!MPIcf::IamMaster()) {
//         xc_v.resize(n);
//         yc_v.resize(n);
//     }
//     MPIcf::Bcast<double>(xc_v, MPIcf::Master());
//     MPIcf::Bcast<double>(yc_v, MPIcf::Master());

//     if (!MPIcf::IamMaster()) {
//         for (int i = 0; i < xc_v.size(); ++i) {
//             LOG_INFO << " Particle " << i << " - Center (" << xc_v[i] << ", " << yc_v[i] << ") , radius = " << r_0
//                      << logger::endl;

//             levelSets_v.push_back(std::make_shared<LevelSetCircle>(xc_v[i], yc_v[i], r_0));
//         }
//     }
//     LOG_INFO << "Number of particles built: " << levelSets_v.size() << logger::endl;
//     return levelSets_v;
// }

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


template <typename T>
bool checkIfValid(double xc, double yc, double delta_l, const std::vector<std::shared_ptr<T>> &levelSets_v) {

    for (int i = 0; i < levelSets_v.size(); i++) {
        double val = levelSets_v[i]->operator()(R2(xc, yc));
        if (val < delta_l)
            return false;
    }
    return true;
}

} // namespace force_project

#endif