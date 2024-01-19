

#include "../tool.hpp"
#include <random>
#include <chrono>

#include "../force_project/util.hpp"
// #include "../force_project/input_handler.hpp"
// #include "../force_project/stokes_solver.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using paraview_t = Paraview<mesh_t>;

int main(int argc, char **argv) {

    double t0 = MPIcf::Wtime();

    mesh_t Th("../../mesh/kiln_100.msh");

    paraview_t writer(Th, "kiln_bed.vtk");

    Lagrange2 FEu(2);
    space_t Vh(Th, FEu);
    space_t Ph(Th, DataFE<Mesh2>::P1);

    // auto fw = [](R2 P, int i) { return (i == 0) * P.y - (i == 1) * P.x; };
    // auto fs = [](R2 P, int i) { return (i == 0) * 2. / globalVariable::Pi + (i == 1) * 0.; };
    double speed = 0.3246312408709453;
    auto fw      = [speed](R2 P, int i) -> double { return speed / norm2(P) * ((i == 0) * P[1] - (i == 1) * P[0]); };
    auto fs      = [speed](R2 P, int i) -> double { return (i == 0) * speed; };

    fct_t f_w(Vh, fw);
    fct_t f_s(Vh, fs);

    FEM<mesh_t> stokes(Vh);
    stokes.add(Ph);

    TestFunction<Mesh2> u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);
    Normal n;
    double mu      = 1.;
    double lambdaB = 1e5;

    // stokes.addBilinear(contractProduct(2 * mu * Eps(u), Eps(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
    //                    Th);

    // stokes.addBilinear(innerProduct(lambdaB * u, v) + innerProduct(p, v * n) - innerProduct(u * n, q) -
    //                        innerProduct(2. * mu * Eps(u) * n, v) - innerProduct(u, 2. * mu * Eps(v) * n),
    //                    Th, INTEGRAL_BOUNDARY, {1, 2});

    // stokes.addLinear(innerProduct(f_w.exprList(2), lambdaB * v) - innerProduct(f_w.exprList(2), 2. * mu * Eps(v) * n)
    // -
    //                      innerProduct(f_w.exprList(2), p * n),
    //                  Th, INTEGRAL_BOUNDARY, {1});

    // stokes.addLinear(innerProduct(f_s.exprList(2), lambdaB * v) - innerProduct(f_s.exprList(2), 2. * mu * Eps(v) * n)
    // -
    //                      innerProduct(f_s.exprList(2), p * n),
    //                  Th, INTEGRAL_BOUNDARY, {2});

    stokes.addBilinear(contractProduct(2 * mu * grad(u), grad(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                       Th);

    stokes.addBilinear(innerProduct(lambdaB * u, v) + innerProduct(p, v * n) - innerProduct(u * n, q) -
                           innerProduct(2. * mu * grad(u) * n, v) - innerProduct(u, 2. * mu * grad(v) * n),
                       Th, INTEGRAL_BOUNDARY, {1, 2});

    stokes.addLinear(innerProduct(f_w.exprList(2), lambdaB * v) - innerProduct(f_w.exprList(2), 2. * mu * grad(v) * n) -
                         innerProduct(f_w.exprList(2), p * n),
                     Th, INTEGRAL_BOUNDARY, {1});

    stokes.addLinear(innerProduct(f_s.exprList(2), lambdaB * v) - innerProduct(f_s.exprList(2), 2. * mu * grad(v) * n) -
                         innerProduct(f_s.exprList(2), p * n),
                     Th, INTEGRAL_BOUNDARY, {2});

    stokes.addLagrangeMultiplier(innerProduct(1., q), 0., Th);

    double t1 = MPIcf::Wtime();
    std::cout << "Assembly time: " << t1 - t0 << std::endl;

    stokes.solve("umfpack");

    double t2 = MPIcf::Wtime();
    std::cout << "Solver time: " << t2 - t1 << std::endl;

    // Extract solution
    std::span<double> data_uh{std::span(stokes.rhs_.data(), Vh.get_nb_dof())};
    std::span<double> data_ph{std::span(stokes.rhs_.data() + Vh.get_nb_dof(), Ph.get_nb_dof())};
    fct_t uh(Vh, data_uh);
    fct_t ph(Ph, data_ph);

    writer.add(uh, "velocity", 0, 2);
    writer.add(ph, "pressure", 0, 1);

    double t3 = MPIcf::Wtime();
    std::cout << "Total time: " << t3 - t2 << std::endl;

    return 0;
}