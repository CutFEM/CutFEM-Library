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

#ifndef CUTFEM_CPP_PROBLEM_SOLVER_ADVECTION_HPP
#define CUTFEM_CPP_PROBLEM_SOLVER_ADVECTION_HPP

#include "baseProblem.hpp"

template <typename mesh_t> class StreamLineDiffusion : public DefaultExpression {

    using fct_t = FunFEM<mesh_t>;

  public:
    StreamLineDiffusion(const fct_t &u, double hh, double ddt) : U(u), h(hh), dt(ddt) {}

    double evalOnBackMesh(const int k, int dom, const R *x, const R t, const R *normal) const final {
        double ux = U.eval(k, x, 0, op_id);
        double uy = U.eval(k, x, 1, op_id);
        return 2. / std::sqrt(1. / (dt * dt) + sqrt(ux * ux + uy * uy) / (h * h));
    }

  private:
    const fct_t &U;
    double h;
    double dt;
};

namespace solver {

namespace fem {

namespace advection {

namespace streamline_diffusion {

/// @brief Solve PDE advection in R2 using Crang Nicolson scheme for the time derivative
/// @param u_init Solution at time n
/// @param betap velocity field at time n
/// @param beta velocity at time n+1
/// @param dt time step
/// @return vector containing the dof of the solution
std::vector<double> solve(const FunFEM<Mesh2> &u_init, const FunFEM<Mesh2> &betap, const FunFEM<Mesh2> &beta,
                          double dt) {

    const auto &Vh = u_init.getSpace();
    const auto &Th = Vh.Th;
    FEM<Mesh2> problem(Vh);
    TestFunction<Mesh2> u(Vh, 1), v(Vh, 1);

    double h = Th[0].hElement();
    auto U   = beta.exprList();
    auto Up  = betap.exprList();
    auto u0  = u_init.expr();

    auto tauSD = std::make_shared<StreamLineDiffusion<Mesh2>>(beta, h, dt);

    problem.addBilinear((u + (U * (0.5 * dt * grad(u))), v + U * (tauSD * grad(v))), Th);
    problem.addLinear(innerProduct(u0 - 0.5 * dt * (Up[0] * dx(u0) + Up[1] * dy(u0)), v + U * (tauSD * grad(v))), Th);

    problem.solve();

    return problem.rhs_;
}

} // namespace streamline_diffusion
} // namespace advection
} // namespace fem
} // namespace solver

#endif