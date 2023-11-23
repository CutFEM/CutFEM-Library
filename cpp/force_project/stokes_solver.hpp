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

#ifndef FORCE_PROJECT_STOKES_SOLVER_HPP
#define FORCE_PROJECT_STOKES_SOLVER_HPP

namespace force_project {

std::vector<double> solveStokes(CutFESpace<Mesh2> &Vh, CutFESpace<Mesh2> &Ph, InterfaceLevelSet<Mesh2> &interface,
                                double radius) {

    const auto &Khi     = Vh.get_mesh();
    double hi           = Khi.Th.get_mesh_size();
    double penaltyParam = 1e1 / hi;
    double mu1          = 1.;
    double mu2          = 10.;
    double kappa1       = 0.5;
    double kappa2       = 0.5;
    double sigma        = 1.;
    double uPenParam    = 1e0;
    double pPenParam    = 1e0;

    auto f_bc = [](R2 P, int i, int dom) -> double { return (i == 0) * P.x - (i == 1) * P.y; };
    FunFEM<Mesh2> gh(Vh, f_bc);

    MacroElement<Mesh2> macro(Khi, 0.25);

    CutFEM<Mesh2> stokes(Vh);
    stokes.add(Ph);

    LOG_INFO << " Stokes has " << stokes.get_nb_dof() << " dofs" << logger::endl;

    Normal n;
    Tangent t;
    CutFEMParameter mu(mu1, mu2);
    TestFunction<Mesh2> u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0), p1(Ph, 1, 0, 0);
    TestFunction<Mesh2> grad2un = grad(grad(u) * n) * n;

    stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                       Khi);
    LOG_INFO << " Bilineair form on K done " << logger::endl;

    stokes.addBilinear(-innerProduct(average(2 * mu * Eps(u) * t * n, kappa1, kappa2), jump(v * t)) +
                           innerProduct(jump(u * t), average(2 * mu * Eps(v) * t * n, kappa1, kappa2)) +
                           innerProduct(1. / hi * (jump(u * t)), jump(v * t)),
                       Khi, INTEGRAL_INNER_EDGE_2D);
    LOG_INFO << " Bilineair form on edges done " << logger::endl;

    stokes.addBilinear(innerProduct(jump(u), -2 * mu * average(Eps(v) * n, kappa1, kappa2)) +
                           innerProduct(-2 * mu * average(Eps(u) * n, kappa1, kappa2), jump(v)) +
                           innerProduct(1. / hi / hi * jump(u), jump(v)) +
                           innerProduct(average(p, kappa1, kappa2), jump(v * n))
                       // - innerProduct(jump(u*n), average(q,0.5,0.5))
                       ,
                       interface);
    LOG_INFO << " Bilineair form on interface done " << logger::endl;

    stokes.addLinear(innerProduct(1. / radius, average(v * n, kappa2, kappa1)) * sigma, interface);
    LOG_INFO << " Linear form on interface done " << logger::endl;
    // stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v)                               // natural
    //                        + innerProduct(u, 2 * mu * Eps(v) * n) + innerProduct(p, v * n) // natural
    //                        + innerProduct(penaltyParam * u, v)
    //                    // - innerProduct(u*n, q)  // essential
    //                    ,
    //                    Khi, INTEGRAL_BOUNDARY, {1, 3});
    // stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) // natural
    //                                                          // + innerProduct(u, 2*mu*Eps(v)*n)
    //                                                          // + innerProduct(p, v*n)  // natural
    //                                                          // + innerProduct(penaltyParam*u, v)
    //                                                          // - innerProduct(u*n, q)  // essential
    //                    ,
    //                    Khi, INTEGRAL_BOUNDARY, {2, 4});
    // stokes.addLinear(innerProduct(1, v * n) // natural
    //                  ,
    //                  Khi, INTEGRAL_BOUNDARY, {4});
    // stokes.addLinear(innerProduct(-1, v * n) // natural
    //                  ,
    //                  Khi, INTEGRAL_BOUNDARY, {2});

    stokes.addBilinear(-innerProduct(2 * mu * Eps(u) * n, v) // natural
                           - innerProduct(u, 2 * mu * Eps(v) * n)
                           // + innerProduct(p, v*n)  // natural
                           + innerProduct(penaltyParam * u * t, v * t),
                       Khi, INTEGRAL_BOUNDARY);
    LOG_INFO << " Bilineair form on boundary done " << logger::endl;
    stokes.addLinear(-innerProduct(gh.exprList(2), 2 * mu * Eps(v) * n) + innerProduct(gh * t, penaltyParam * v * t),
                     Khi, INTEGRAL_BOUNDARY);
    LOG_INFO << " Linear form on boundary done " << logger::endl;

    stokes.setDirichlet(gh, Khi);

    stokes.addFaceStabilization(                                  // [h^(2k+1) h^(2k+1)]
        +innerProduct(uPenParam * pow(hi, -1) * jump(u), jump(v)) // [Method 1: Remove jump in vel]
            + innerProduct(uPenParam * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
            innerProduct(uPenParam * pow(hi, 3) * jump(grad2un), jump(grad2un))

            - innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
            innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
            innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
            innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(v))), jump(grad(q))),
        Khi, macro);
    LOG_INFO << " Face stabilization done " << logger::endl;

    CutFEM<Mesh2> lagr(Vh);
    lagr.add(Ph);
    lagr.addLinear(innerProduct(1., p), Khi);
    std::vector<double> lag_row(lagr.rhs_);
    std::fill(lagr.rhs_.begin(), lagr.rhs_.end(), 0.);
    lagr.addLinear(innerProduct(1, v * n), Khi, INTEGRAL_BOUNDARY);

    stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);

    // stokes.addLagrangeMultiplier(innerProduct(1., p1), 0., Khi);
    LOG_INFO << " Lagrange multiplier done " << logger::endl;

    stokes.solve();

    return stokes.rhs_;
}

} // namespace force_project

#endif