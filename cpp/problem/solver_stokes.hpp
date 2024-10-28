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

#ifndef CUTFEM_CPP_PROBLEM_SOLVER_STOKES_HPP
#define CUTFEM_CPP_PROBLEM_SOLVER_STOKES_HPP

namespace solver {

namespace cutfem {

namespace stokes {

/// @brief Solve the Stokes interface problem
/// @param Vh Velocity space
/// @param Ph Pressure space
/// @param interface The interface
/// @param gh boundary function
/// @param fh force function
/// @param mu viscosity
/// @param surface_tension surface tension force
/// @param delta macro element delta parameter
/// @return vector containing data of the velocity and the pressure
std::vector<double> solve(CutFESpace<Mesh2> &Vh, CutFESpace<Mesh2> &Ph, InterfaceLevelSet<Mesh2> &interface,
                          FunFEM<Mesh2> &gh, FunFEM<Mesh2> &fh, CutFEMParameter mu, double surface_tension,
                          double delta, std::list<int> dirichlet_labels, std::list<int> neumann_labels) {

    const auto &Khi           = Vh.get_mesh();
    double hi                 = Khi[0].hElement();
    double jump_penalty       = 10 / hi;
    double boundary_penalty   = 1. / hi;
    double tangential_penalty = 1. / hi;
    double kappa1             = 0.5;
    double kappa2             = 0.5;
    double ghost_penalty      = 1.;

    MacroElement<Mesh2> macro(Khi, delta);

    CutFEM<Mesh2> stokes(Vh);
    stokes.add(Ph);

    if (globalVariable::verbose > 0) {
        LOG_INFO << " Stokes has " << stokes.get_nb_dof() << " dofs" << logger::endl;
    }
    auto t0 = std::chrono::high_resolution_clock::now();

    Normal n;
    Tangent t;
    TestFunction<Mesh2> u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0), p1(Ph, 1, 0, 0);
    TestFunction<Mesh2> grad2un = grad(grad(u) * n) * n;
    stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                       Khi);

    stokes.addBilinear(-innerProduct(2 * mu * average(Eps(u) * n, kappa1, kappa2), jump(v)) +
                           innerProduct(jump(u), 2 * mu * average(Eps(v) * n, kappa1, kappa2)) +
                           innerProduct(average(p, kappa1, kappa2), jump(v * n)) +
                           innerProduct(jump_penalty * jump(u), jump(v)),
                       interface);

    if (Vh.basisFctType == BasisFctType::BDM1 || Vh.basisFctType == BasisFctType::RT0 ||
        Vh.basisFctType == BasisFctType::RT1) {
        stokes.addBilinear(-innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2), jump(v * t)) +
                               innerProduct(jump(u * t), average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
                               innerProduct(tangential_penalty * (jump(u * t)), jump(v * t)),
                           Khi, INTEGRAL_INNER_EDGE_2D);
    }

    if (!neumann_labels.empty()) {
        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) + innerProduct(p, v * n), Khi, INTEGRAL_BOUNDARY,
                           neumann_labels);
    }

    if (!dirichlet_labels.empty()) {
        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) + innerProduct(p, v * n) +
                               innerProduct(boundary_penalty * u, v) + innerProduct(u, 2 * mu * Eps(v) * n),
                           Khi, INTEGRAL_BOUNDARY, dirichlet_labels);
        stokes.addLinear(innerProduct(gh.exprList(), boundary_penalty * v) +
                             innerProduct(gh.exprList(), 2 * mu * Eps(v) * n),
                         Khi, INTEGRAL_BOUNDARY, dirichlet_labels);
    }

    stokes.addFaceStabilization(                                      // [h^(2k+1) h^(2k+1)]
        +innerProduct(ghost_penalty * pow(hi, -1) * jump(u), jump(v)) // [Method 1: Remove jump in vel]
            + innerProduct(ghost_penalty * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad2un), jump(grad2un))

            - innerProduct(ghost_penalty * pow(hi, 1) * jump(p), jump(div(v))) +
            innerProduct(ghost_penalty * pow(hi, 1) * jump(div(u)), jump(q)) -
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad(div(v))), jump(grad(q))),
        Khi);

    stokes.addLinear(innerProduct(fh.exprList(), v), Khi);

    stokes.addLinear(innerProduct(surface_tension, average(v * n, kappa2, kappa1)), interface);

    CutFEM<Mesh2> lagr(Vh);
    lagr.add(Ph);
    lagr.addLinear(innerProduct(1., p1), Khi);
    std::vector<double> lag_row(lagr.rhs_);
    std::fill(lagr.rhs_.begin(), lagr.rhs_.end(), 0.);
    lagr.addLinear(innerProduct(1, v * n), Khi, INTEGRAL_BOUNDARY);
    stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);

    auto t1 = std::chrono::high_resolution_clock::now();

    stokes.solve();

    auto t2 = std::chrono::high_resolution_clock::now();

    if (globalVariable::verbose > 0) {
        LOG_INFO << "Time to assembly: " << std::chrono::duration<double>(t1 - t0).count() << " sec." << logger::endl;
        LOG_INFO << "Time to solve: " << std::chrono::duration<double>(t2 - t1).count() << " sec." << logger::endl;
    }
    return stokes.rhs_;
}

/// @brief Solve the Stokes interface problem
/// @param Vh Velocity space
/// @param Ph Pressure space
/// @param interface The interface
/// @param gh boundary function
/// @param fh force function
/// @param mu viscosity
/// @param surface_tension surface tension force
/// @param delta macro element delta parameter
/// @return vector containing data of the velocity and the pressure
std::vector<double> solve(CutFESpace<Mesh2> &Vh, CutFESpace<Mesh2> &Ph, InterfaceLevelSet<Mesh2> &interface,
                          FunFEM<Mesh2> &gh, FunFEM<Mesh2> &fh, CutFEMParameter mu, double sigma,
                          const FunFEM<Mesh2> &H, double delta, std::list<int> dirichlet_labels,
                          std::list<int> neumann_labels) {

    const auto &Khi           = Vh.get_mesh();
    double hi                 = Khi[0].hElement();
    double jump_penalty       = 10 / hi;
    double boundary_penalty   = 1. / hi;
    double tangential_penalty = 1. / hi;
    double kappa1             = 0.5;
    double kappa2             = 0.5;
    double ghost_penalty      = 1.;

    MacroElement<Mesh2> macro(Khi, delta);

    CutFEM<Mesh2> stokes(Vh);
    stokes.add(Ph);

    LOG_INFO << " Stokes has " << stokes.get_nb_dof() << " dofs" << logger::endl;
    auto t0 = std::chrono::high_resolution_clock::now();

    Normal n;
    Tangent t;
    TestFunction<Mesh2> u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0), p1(Ph, 1, 0, 0);
    TestFunction<Mesh2> grad2un = grad(grad(u) * n) * n;

    stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                       Khi);

    stokes.addBilinear(-innerProduct(2 * mu * average(Eps(u) * n, kappa1, kappa2), jump(v)) +
                           innerProduct(jump(u), 2 * mu * average(Eps(v) * n, kappa1, kappa2)) +
                           innerProduct(average(p, kappa1, kappa2), jump(v * n)) +
                           innerProduct(jump_penalty * jump(u), jump(v)),
                       interface);

    if (Vh.basisFctType == BasisFctType::BDM1 || Vh.basisFctType == BasisFctType::RT0 ||
        Vh.basisFctType == BasisFctType::RT1) {
        stokes.addBilinear(-innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2), jump(v * t)) +
                               innerProduct(jump(u * t), average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
                               innerProduct(tangential_penalty * (jump(u * t)), jump(v * t)),
                           Khi, INTEGRAL_INNER_EDGE_2D);
    }

    if (!neumann_labels.empty()) {
        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v)
                           // +innerProduct(p, v * n)
                           ,
                           Khi, INTEGRAL_BOUNDARY, neumann_labels);
    }

    if (!dirichlet_labels.empty()) {
        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) + innerProduct(p, v * n) +
                               innerProduct(boundary_penalty * u, v) + innerProduct(u, 2 * mu * Eps(v) * n),
                           Khi, INTEGRAL_BOUNDARY, dirichlet_labels);

        stokes.addLinear(innerProduct(gh.exprList(), boundary_penalty * v) +
                             innerProduct(gh.exprList(), 2 * mu * Eps(v) * n),
                         Khi, INTEGRAL_BOUNDARY, dirichlet_labels);
    }

    stokes.addFaceStabilization(                                      // [h^(2k+1) h^(2k+1)]
        +innerProduct(ghost_penalty * pow(hi, -1) * jump(u), jump(v)) // [Method 1: Remove jump in vel]
            + innerProduct(ghost_penalty * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad2un), jump(grad2un))

            - innerProduct(ghost_penalty * pow(hi, 1) * jump(p), jump(div(v))) +
            innerProduct(ghost_penalty * pow(hi, 1) * jump(div(u)), jump(q)) -
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad(div(v))), jump(grad(q))),
        Khi);
    stokes.addLinear(innerProduct(fh.exprList(), v), Khi);

    stokes.addLinear(innerProduct(H.exprList(2), sigma * average(v, kappa2, kappa1)), interface);

    CutFEM<Mesh2> lagr(Vh);
    lagr.add(Ph);
    lagr.addLinear(innerProduct(1., p1), Khi);
    std::vector<double> lag_row(lagr.rhs_);
    std::fill(lagr.rhs_.begin(), lagr.rhs_.end(), 0.);
    lagr.addLinear(innerProduct(1, v * n), Khi, INTEGRAL_BOUNDARY);
    stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);

    auto t1 = std::chrono::high_resolution_clock::now();

    stokes.solve();

    auto t2 = std::chrono::high_resolution_clock::now();

    LOG_INFO << "Time to assembly: " << std::chrono::duration<double>(t1 - t0).count() << " sec." << logger::endl;
    LOG_INFO << "Time to solve: " << std::chrono::duration<double>(t2 - t1).count() << " sec." << logger::endl;

    return stokes.rhs_;
}

/// @brief Solve the Stokes interface problem
/// @param Vh Velocity space
/// @param Ph Pressure space
/// @param interface The interface
/// @param gh boundary function
/// @param fh force function
/// @param mu viscosity
/// @param surface_tension surface tension force
/// @param delta macro element delta parameter
/// @return vector containing data of the velocity and the pressure
std::vector<double> solve(CutFESpace<Mesh2> &Vh, CutFESpace<Mesh2> &Ph, InterfaceLevelSet<Mesh2> &interface,
                          FunFEM<Mesh2> &fh, CutFEMParameter mu, double sigma, const FunFEM<Mesh2> &H, double delta,
                          std::vector<std::pair<FunFEM<Mesh2> *, std::list<int> *>> &bc,
                          std::list<int> neumann_labels) {

    const auto &Khi           = Vh.get_mesh();
    double hi                 = Khi[0].hElement();
    double jump_penalty       = 10 / hi;
    double boundary_penalty   = 1. / hi;
    double tangential_penalty = 1. / hi;
    double kappa1             = 0.5;
    double kappa2             = 0.5;
    double ghost_penalty      = 1.;

    MacroElement<Mesh2> macro(Khi, delta);

    CutFEM<Mesh2> stokes(Vh);
    stokes.add(Ph);

    LOG_INFO << " Stokes has " << stokes.get_nb_dof() << " dofs" << logger::endl;
    auto t0 = std::chrono::high_resolution_clock::now();

    Normal n;
    Tangent t;
    TestFunction<Mesh2> u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0), p1(Ph, 1, 0, 0);
    TestFunction<Mesh2> grad2un = grad(grad(u) * n) * n;

    stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                       Khi);

    stokes.addBilinear(-innerProduct(2 * mu * average(Eps(u) * n, kappa1, kappa2), jump(v)) +
                           innerProduct(jump(u), 2 * mu * average(Eps(v) * n, kappa1, kappa2)) +
                           innerProduct(average(p, kappa1, kappa2), jump(v * n)) +
                           innerProduct(jump_penalty * jump(u), jump(v)),
                       interface);

    if (Vh.basisFctType == BasisFctType::BDM1 || Vh.basisFctType == BasisFctType::RT0 ||
        Vh.basisFctType == BasisFctType::RT1) {
        stokes.addBilinear(-innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2), jump(v * t)) +
                               innerProduct(jump(u * t), average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
                               innerProduct(tangential_penalty * (jump(u * t)), jump(v * t)),
                           Khi, INTEGRAL_INNER_EDGE_2D);
    }
    stokes.addFaceStabilization(                                      // [h^(2k+1) h^(2k+1)]
        +innerProduct(ghost_penalty * pow(hi, -1) * jump(u), jump(v)) // [Method 1: Remove jump in vel]
            + innerProduct(ghost_penalty * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad2un), jump(grad2un))

            - innerProduct(ghost_penalty * pow(hi, 1) * jump(p), jump(div(v))) +
            innerProduct(ghost_penalty * pow(hi, 1) * jump(div(u)), jump(q)) -
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad(div(v))), jump(grad(q))),
        Khi);
    stokes.addLinear(innerProduct(fh.exprList(), v), Khi);

    stokes.addLinear(innerProduct(H.exprList(2), sigma * average(v, kappa2, kappa1)), interface);

    if (!neumann_labels.empty()) {
        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) + innerProduct(p, v * n), Khi, INTEGRAL_BOUNDARY,
                           neumann_labels);
    }

    for (auto [gh, dirichlet_labels] : bc) {
        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) + innerProduct(p, v * n) +
                               innerProduct(boundary_penalty * u, v) + innerProduct(u, 2 * mu * Eps(v) * n),
                           Khi, INTEGRAL_BOUNDARY, *dirichlet_labels);
        stokes.addLinear(innerProduct(gh->exprList(), boundary_penalty * v) +
                             innerProduct(gh->exprList(), 2 * mu * Eps(v) * n),
                         Khi, INTEGRAL_BOUNDARY, *dirichlet_labels);
    }

    CutFEM<Mesh2> lagr(Vh);
    lagr.add(Ph);
    lagr.addLinear(innerProduct(1., p1), Khi);
    std::vector<double> lag_row(lagr.rhs_);
    std::fill(lagr.rhs_.begin(), lagr.rhs_.end(), 0.);
    for (auto [gh, dirichlet_labels] : bc) {
        lagr.addLinear(innerProduct(1, v * n), Khi, INTEGRAL_BOUNDARY, *dirichlet_labels);
    }
    if (!neumann_labels.empty()) {
        lagr.addLinear(innerProduct(1, v * n), Khi, INTEGRAL_BOUNDARY, neumann_labels);
    }

    stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);

    auto t1 = std::chrono::high_resolution_clock::now();

    stokes.solve();

    auto t2 = std::chrono::high_resolution_clock::now();

    LOG_INFO << "Time to assembly: " << std::chrono::duration<double>(t1 - t0).count() << " sec." << logger::endl;
    LOG_INFO << "Time to solve: " << std::chrono::duration<double>(t2 - t1).count() << " sec." << logger::endl;

    return stokes.rhs_;
}

std::vector<double> solve(CutFESpace<Mesh2> &Vh, CutFESpace<Mesh2> &Ph, InterfaceLevelSet<Mesh2> &interface,
                          FunFEM<Mesh2> &fh, CutFEMParameter mu, double surface_tension, double delta,
                          std::vector<std::pair<FunFEM<Mesh2> *, std::list<int> *>> &bc,
                          std::list<int> neumann_labels) {

    const auto &Khi           = Vh.get_mesh();
    double hi                 = Khi[0].hElement();
    double jump_penalty       = 10 / hi;
    double boundary_penalty   = 1. / hi;
    double tangential_penalty = 1. / hi;
    double kappa1             = 0.5;
    double kappa2             = 0.5;
    double ghost_penalty      = 1.;

    MacroElement<Mesh2> macro(Khi, delta);

    CutFEM<Mesh2> stokes(Vh);
    stokes.add(Ph);

    LOG_INFO << " Stokes has " << stokes.get_nb_dof() << " dofs" << logger::endl;
    auto t0 = std::chrono::high_resolution_clock::now();

    Normal n;
    Tangent t;
    TestFunction<Mesh2> u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0), p1(Ph, 1, 0, 0);
    TestFunction<Mesh2> grad2un = grad(grad(u) * n) * n;

    stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                       Khi);

    stokes.addBilinear(-innerProduct(2 * mu * average(Eps(u) * n, kappa1, kappa2), jump(v)) +
                           innerProduct(jump(u), 2 * mu * average(Eps(v) * n, kappa1, kappa2)) +
                           innerProduct(average(p, kappa1, kappa2), jump(v * n)) +
                           innerProduct(jump_penalty * jump(u), jump(v)),
                       interface);

    if (Vh.basisFctType == BasisFctType::BDM1 || Vh.basisFctType == BasisFctType::RT0 ||
        Vh.basisFctType == BasisFctType::RT1) {
        stokes.addBilinear(-innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2), jump(v * t)) +
                               innerProduct(jump(u * t), average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
                               innerProduct(tangential_penalty * (jump(u * t)), jump(v * t)),
                           Khi, INTEGRAL_INNER_EDGE_2D);
    }
    stokes.addFaceStabilization(                                      // [h^(2k+1) h^(2k+1)]
        +innerProduct(ghost_penalty * pow(hi, -1) * jump(u), jump(v)) // [Method 1: Remove jump in vel]
            + innerProduct(ghost_penalty * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad2un), jump(grad2un))

            - innerProduct(ghost_penalty * pow(hi, 1) * jump(p), jump(div(v))) +
            innerProduct(ghost_penalty * pow(hi, 1) * jump(div(u)), jump(q)) -
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
            innerProduct(ghost_penalty * pow(hi, 3) * jump(grad(div(v))), jump(grad(q))),
        Khi);
    stokes.addLinear(innerProduct(fh.exprList(), v), Khi);

    stokes.addLinear(innerProduct(surface_tension, average(v * n, kappa2, kappa1)), interface);

    if (!neumann_labels.empty()) {
        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) + innerProduct(p, v * n), Khi, INTEGRAL_BOUNDARY,
                           neumann_labels);
    }

    for (auto [gh, dirichlet_labels] : bc) {
        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) + innerProduct(p, v * n) +
                               innerProduct(boundary_penalty * u, v) + innerProduct(u, 2 * mu * Eps(v) * n),
                           Khi, INTEGRAL_BOUNDARY, *dirichlet_labels);
        stokes.addLinear(innerProduct(gh->exprList(), boundary_penalty * v) +
                             innerProduct(gh->exprList(), 2 * mu * Eps(v) * n),
                         Khi, INTEGRAL_BOUNDARY, *dirichlet_labels);
    }

    CutFEM<Mesh2> lagr(Vh);
    lagr.add(Ph);
    lagr.addLinear(innerProduct(1., p1), Khi);
    std::vector<double> lag_row(lagr.rhs_);
    std::fill(lagr.rhs_.begin(), lagr.rhs_.end(), 0.);
    for (auto [gh, dirichlet_labels] : bc) {
        lagr.addLinear(innerProduct(1, v * n), Khi, INTEGRAL_BOUNDARY, *dirichlet_labels);
    }
    if (!neumann_labels.empty()) {
        lagr.addLinear(innerProduct(1, v * n), Khi, INTEGRAL_BOUNDARY, neumann_labels);
    }

    stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);

    auto t1 = std::chrono::high_resolution_clock::now();

    stokes.solve();

    auto t2 = std::chrono::high_resolution_clock::now();

    LOG_INFO << "Time to assembly: " << std::chrono::duration<double>(t1 - t0).count() << " sec." << logger::endl;
    LOG_INFO << "Time to solve: " << std::chrono::duration<double>(t2 - t1).count() << " sec." << logger::endl;

    return stokes.rhs_;
}

} // namespace stokes
} // namespace cutfem

namespace fem {

namespace stokes {

/// @brief Solve the Stokes interface problem
/// @param Vh Velocity space
/// @param Ph Pressure space
/// @param interface The interface
/// @param gh boundary function
/// @param fh force function
/// @param mu viscosity
/// @param surface_tension surface tension force
/// @param delta macro element delta parameter
/// @return vector containing data of the velocity and the pressure
std::vector<double> solve(GFESpace<Mesh2> &Vh, GFESpace<Mesh2> &Ph, FunFEM<Mesh2> &gh, FunFEM<Mesh2> &fh, double mu,
                          std::list<int> dirichlet_labels, std::list<int> neumann_labels

) {

    const auto &Kh            = Vh.Th;
    double hi                 = Kh[0].hElement();
    double boundary_penalty   = 100. / hi;
    double tangential_penalty = 1. / hi;
    double kappa1             = 0.5;
    double kappa2             = 0.5;

    FEM<Mesh2> stokes(Vh);
    stokes.add(Ph);

    if (globalVariable::verbose > 0) {
        LOG_INFO << " Stokes has " << stokes.get_nb_dof() << " dofs" << logger::endl;
    }
    auto t0 = std::chrono::high_resolution_clock::now();

    Normal n;
    Tangent t;
    TestFunction<Mesh2> u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);
    TestFunction<Mesh2> grad2un = grad(grad(u) * n) * n;

    // Matrix assembly
    stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                       Kh);

    if (Vh.basisFctType == BasisFctType::BDM1 || Vh.basisFctType == BasisFctType::RT0 ||
        Vh.basisFctType == BasisFctType::RT1) {
        stokes.addBilinear(-innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2), jump(v * t)) +
                               innerProduct(jump(u * t), average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
                               innerProduct(tangential_penalty * (jump(u * t)), jump(v * t)),
                           Kh, INTEGRAL_INNER_EDGE_2D);
    }

    stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) + innerProduct(p, v * n) +
                           innerProduct(boundary_penalty * u, v) + innerProduct(u, 2 * mu * Eps(v) * n),
                       Kh, INTEGRAL_BOUNDARY, dirichlet_labels);
    stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) + innerProduct(p, v * n), Kh, INTEGRAL_BOUNDARY,
                       neumann_labels);

    stokes.addLinear(innerProduct(fh.exprList(), v), Kh);

    stokes.addLinear(innerProduct(gh.exprList(), boundary_penalty * v) +
                         innerProduct(gh.exprList(), 2 * mu * Eps(v) * n),
                     Kh, INTEGRAL_BOUNDARY, dirichlet_labels);

    // Lagrange multiplier
    FEM<Mesh2> lagr(Vh);
    lagr.add(Ph);
    lagr.addLinear(innerProduct(1., p), Kh);
    // stokes.addLagrangeVecToRowAndCol(lagr.rhs_, lagr.rhs_, 0);
    std::vector<double> lag_row(lagr.rhs_);
    std::fill(lagr.rhs_.begin(), lagr.rhs_.end(), 0.);
    lagr.addLinear(innerProduct(1, v * n), Kh, INTEGRAL_BOUNDARY);
    stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);

    auto t1 = std::chrono::high_resolution_clock::now();

    stokes.solve();

    auto t2 = std::chrono::high_resolution_clock::now();

    if (globalVariable::verbose > 0) {
        LOG_INFO << "Time to assembly: " << std::chrono::duration<double>(t1 - t0).count() << " sec." << logger::endl;
        LOG_INFO << "Time to solve: " << std::chrono::duration<double>(t2 - t1).count() << " sec." << logger::endl;
    }
    return stokes.rhs_;
}

/// @brief Solve the Stokes interface problem
/// @param Vh Velocity space
/// @param Ph Pressure space
/// @param interface The interface
/// @param gh boundary function
/// @param fh force function
/// @param mu viscosity
/// @param surface_tension surface tension force
/// @param delta macro element delta parameter
/// @return vector containing data of the velocity and the pressure
std::vector<double> solve(GFESpace<Mesh2> &Vh, GFESpace<Mesh2> &Ph, FunFEM<Mesh2> &fh, double mu,
                          std::vector<std::pair<FunFEM<Mesh2> *, std::list<int> *>> &bc, std::list<int> neumann_labels

) {

    const auto &Kh            = Vh.Th;
    double hi                 = Kh[0].hElement();
    double boundary_penalty   = 1. / hi;
    double tangential_penalty = 1. / hi;
    double kappa1             = 0.5;
    double kappa2             = 0.5;

    FEM<Mesh2> stokes(Vh);
    stokes.add(Ph);

    if (globalVariable::verbose > 0) {
        LOG_INFO << " Stokes has " << stokes.get_nb_dof() << " dofs" << logger::endl;
    }
    auto t0 = std::chrono::high_resolution_clock::now();

    Normal n;
    Tangent t;
    TestFunction<Mesh2> u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);
    TestFunction<Mesh2> grad2un = grad(grad(u) * n) * n;

    // Matrix assembly
    stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                       Kh);

    if (Vh.basisFctType == BasisFctType::BDM1 || Vh.basisFctType == BasisFctType::RT0 ||
        Vh.basisFctType == BasisFctType::RT1) {
        stokes.addBilinear(-innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2), jump(v * t)) +
                               innerProduct(jump(u * t), average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
                               innerProduct(tangential_penalty * (jump(u * t)), jump(v * t)),
                           Kh, INTEGRAL_INNER_EDGE_2D);
    }

    stokes.addLinear(innerProduct(fh.exprList(), v), Kh);

    if (!neumann_labels.empty()) {
        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) + innerProduct(p, v * n), Kh, INTEGRAL_BOUNDARY,
                           neumann_labels);
    }
    for (auto [gh, dirichlet_labels] : bc) {

        stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) + innerProduct(p, v * n) +
                               innerProduct(boundary_penalty * u, v) + innerProduct(u, 2 * mu * Eps(v) * n),
                           Kh, INTEGRAL_BOUNDARY, *dirichlet_labels);
        stokes.addLinear(innerProduct(gh->exprList(), boundary_penalty * v) +
                             innerProduct(gh->exprList(), 2 * mu * Eps(v) * n),
                         Kh, INTEGRAL_BOUNDARY, *dirichlet_labels);
    }

    // Lagrange multiplier
    FEM<Mesh2> lagr(Vh);
    lagr.add(Ph);
    lagr.addLinear(innerProduct(1., p), Kh);
    // stokes.addLagrangeVecToRowAndCol(lagr.rhs_, lagr.rhs_, 0);
    std::vector<double> lag_row(lagr.rhs_);
    std::fill(lagr.rhs_.begin(), lagr.rhs_.end(), 0.);
    lagr.addLinear(innerProduct(1, v * n), Kh, INTEGRAL_BOUNDARY);
    stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);

    auto t1 = std::chrono::high_resolution_clock::now();

    stokes.solve();

    auto t2 = std::chrono::high_resolution_clock::now();

    if (globalVariable::verbose > 0) {
        LOG_INFO << "Time to assembly: " << std::chrono::duration<double>(t1 - t0).count() << " sec." << logger::endl;
        LOG_INFO << "Time to solve: " << std::chrono::duration<double>(t2 - t1).count() << " sec." << logger::endl;
    }
    return stokes.rhs_;
}

} // namespace stokes

} // namespace fem

} // namespace solver

#endif