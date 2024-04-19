#include "cpp/cutfem.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;
using MatMap     = std::map<std::pair<int, int>, double>;

double c0 = 0.5;
double fun_levelSet(const R2 P, const int i) { return -P.x - P.y + c0; }
double fun_boundary(const R2 P, int elementComp, double t) { return sin(std::numbers::pi * (P.x + P.y - 4 * t)); }
double fun_solution(const R2 P, int elementComp, int domain, double t) {
    if (domain == 0)
        return sin(std::numbers::pi * (P.x + P.y - 4 * t));
    else
        return 4. / 3 * sin(4. / 3 * std::numbers::pi * (P.x + P.y - 3 * t - c0 / 4));
}

void assembly(const space_t &Wh, const Interface<mesh_t> &interface, MatMap &Ah, MatMap &Mh) {

    auto t0 = std::chrono::high_resolution_clock::now();
    const cutmesh_t &Khi(Wh.get_mesh());
    CutFEM<Mesh2> problem(Wh);
    CutFEM_R2 beta({R2(3, 1), R2(2, 1)});
    CutFEMParameter lambda(0., 1.);

    const MeshParameter &h(Parameter::h);
    CutFEMParameter lambdaE(3, 2); // max B = sqrt(10)/sqrt(5)

    double lambdaB = 3;
    double Cstab   = 5e-1;
    double Cstabt  = 5e-1;

    Normal n;
    funtest_t u(Wh, 1), v(Wh, 1);

    funtest_t Hessu = grad(grad(u)), Hessv = grad(grad(v));
    // BUILDING A
    // =====================================================
    problem.set_map(Ah);
    problem.addBilinear(innerProduct(beta * u, grad(v)), Khi);

    problem.addBilinear(-jump(innerProduct(beta * u, v * n)) - innerProduct(jump(beta * u * n), jump(lambda * v)),
                        interface);
    //[ u \beta \cdot n]=3u_1-2u_2=0
    problem.addFaceStabilization(-innerProduct(jump(u), Cstab * jump(v)) -
                                     innerProduct((h ^ 2) * jump(grad(u)), Cstab * jump(grad(v)))
                                 // - contractProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
                                 ,
                                 Khi);

    // F(u)_e = {B.u} - lambda_e/2 [u]
    problem.addBilinear(-innerProduct(average(beta * u * n), jump(v)) - innerProduct(0.5 * lambdaE * jump(u), jump(v)),
                        Khi, INTEGRAL_INNER_EDGE_2D);

    problem.addBilinear(-innerProduct(beta * u * n, v), Khi, INTEGRAL_BOUNDARY, {2, 3} // label other boundary
    );
    problem.addBilinear(-innerProduct(u, beta * (0.5 * v) * n) - innerProduct(u, lambdaB * 0.5 * v), Khi,
                        INTEGRAL_BOUNDARY, {1, 4} // label left boundary
    );

    // BUILDING (M + S)
    // =====================================================
    problem.set_map(Mh);
    problem.addBilinear(innerProduct(u, v), Khi);
    problem.addFaceStabilization(innerProduct(h * jump(u), Cstab * jump(v)) +
                                     innerProduct((h ^ 3) * jump(grad(u)), Cstab * jump(grad(v))) +
                                     contractProduct((h ^ 5) * jump(Hessu), Cstab * jump(Hessv)),
                                 Khi);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << " Time assembly \t" << std::chrono::duration<double>(t1 - t0).count() << std::endl;
}

void solve_problem(const space_t &Wh, const Interface<mesh_t> &interface, std::span<double> u0, std::span<double> uh,
                   MatMap &Ah, MatMap &Mh, double tn) {

    double t0 = MPIcf::Wtime();
    const cutmesh_t &Khi(Wh.get_mesh());

    ProblemOption optionProblem;
    optionProblem.solver_name_  = "umfpack";
    optionProblem.clear_matrix_ = false;
    CutFEM<Mesh2> problem(Wh, optionProblem);

    CutFEM_R2 beta({R2(3, 1), R2(2, 1)});
    double lambdaB = 3;

    fct_t gh(Wh, fun_solution, tn);
    auto gx = gh.expr();

    Normal n;
    funtest_t u(Wh, 1), v(Wh, 1);

    // MULTIPLYING A * u_0  => rhs
    // =====================================================
    // problem.mat = Ah;
    // problem.addMatMul(u0);
    int N = problem.get_nb_dof();
    multiply(N, N, Ah, u0, problem.rhs_);

    problem.addLinear(-innerProduct(gx, beta * (0.5 * v) * n) + innerProduct(gx, lambdaB * 0.5 * v), Khi,
                      INTEGRAL_BOUNDARY, {1, 4} // label left boundary
    );
    // std::cout << " Time building problem \t" << MPIcf::Wtime() - t0 << std::endl;

    // SOLVING  (M+S)E'(t) = rhs
    // =====================================================
    problem.solve(Mh, problem.rhs_);

    std::copy(problem.rhs_.begin(), problem.rhs_.end(), uh.begin());
}

int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);
    auto t0 = std::chrono::high_resolution_clock::now();

    // OUTPUT FILE
    // =====================================================
    std::ofstream outputData("output_test.txt");

    // DEFINITION OF THE MESH and SPACE
    // =====================================================
    int nx = 50;
    int ny = 50;
    mesh_t Th(nx, ny, -1., -1., 2., 2.); // [-1,1]*[-1,1]
    space_t Vh(Th, DataFE<mesh_t>::P2dc);

    int order = Vh.polynomialOrder;

    // DEFINITION OF SPACE AND TIME PARAMETERS
    // =====================================================
    double tid      = 0;
    double meshSize = 2. / nx;
    double dt       = meshSize / 3 / sqrt(10) / std::pow(2., order);
    double tend     = 0.5;
    int niteration  = tend / dt;
    dt              = tend / niteration;
    double errSum   = 0;
    double qu0      = 0;
    std::cout << "Mesh size h = \t" << meshSize << std::endl;
    std::cout << "Time step dt = \t" << dt << std::endl;

    // DEFINITION OF THE LEVELSET
    // =====================================================
    space_t Lh(Th, DataFE<Mesh2>::P1);
    fct_t levelSet(Lh, fun_levelSet);

    // CONSTRUCTION INTERFACE AND CUTSPACE
    // =====================================================
    InterfaceLevelSet<mesh_t> interface(Th, levelSet);
    cutmesh_t Khi(Th, interface);
    cutspace_t Wh(Khi, Vh);

    // DECLARATION OF THE VECTOR CONTAINING THE solution
    // =====================================================
    std::vector<double> u0(Wh.NbDoF(), 0.);
    interpolate(Wh, u0, fun_solution, 0.);
    std::vector<double> uh(u0);

    // ASSEMBLY THE CONSTANT PART
    // ==================================================
    MatMap Ah, Mh;
    assembly(Wh, interface, Ah, Mh);

    // RESOLUTION OF THE PROBLEM_MIXED_DARCY
    // ==================================================
    std::vector<double> u1(Wh.NbDoF(), 0.);
    int ifig = 1;
    for (int i = 0; i < niteration; ++i) {

        std::fill(u1.begin(), u1.end(), 0.);
        // THIRD ORDER RK
        // =================================================
        solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
        Scheme::RK3::step1(u0.begin(), u0.end(), uh.begin(), dt, u1.begin());
        // u1 += dt * uh;
        // u2 += 3. / 4 * u0 + 1. / 4 * u1;
        solve_problem(Wh, interface, u1, uh, Ah, Mh, tid + dt);
        Scheme::RK3::step2(u0.begin(), u0.end(), u1.begin(), uh.begin(), dt, u1.begin());

        // u2 += 1. / 4 * dt * uh;
        solve_problem(Wh, interface, u1, uh, Ah, Mh, tid + 0.5 * dt);
        Scheme::RK3::step3(u0.begin(), u0.end(), u1.begin(), uh.begin(), dt, uh.begin());
        // uh *= 2. / 3 * dt;
        // uh += 1. / 3 * u0 + 2. / 3 * u2;

        u0 = uh;
        tid += dt;

        // COMPUTATION OF THE L2 ERROR
        // =================================================
        interpolate(Wh, u1, fun_solution, tid);

        // u1 -= uh;
        std::transform(u1.begin(), u1.end(), uh.begin(), u1.begin(), std::minus<double>());

        fct_t femErrh(Wh, u1);
        auto femErr = femErrh.expr();

        double errU = sqrt(integral(Khi, femErr * femErr));
        errSum += errU;

        // PLOT THE SOLUTION
        // ==================================================
        if (MPIcf::IamMaster() && i % 10 == 0 || i + 1 == niteration) {
            fct_t femSolh(Wh, uh);

            // fixe for paraview number format
            std::transform(uh.begin(), uh.end(), uh.begin(), [](double x) { return (std::abs(x) < 1e-16) ? 0. : x; });
            fct_t sol(Wh, uh);

            Paraview<mesh_t> writer(Khi, "hyperbolicLawExample1_" + std::to_string(ifig++) + ".vtk");
            writer.add(sol, "uh", 0, 1);

            fct_t solex(Wh, fun_solution, tid);
            writer.add(solex, "uex", 0, 1);
        }
        std::cout << "Iteration " << i << " / " << niteration << " \t time = " << tid << " , || u-uex ||_2 = " << errU
                  << std::endl;
    }

    std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << " Time computation \t" << std::chrono::duration<double>(t1 - t0).count() << std::endl;
    return 0;
}
