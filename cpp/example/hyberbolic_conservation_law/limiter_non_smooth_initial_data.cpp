#include "../tool.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;
using MatMap     = std::map<std::pair<int, int>, R>;

using namespace globalVariable;

// double fun_levelSet(R2 P) { return -(P[1] + 0.5 * P[0] - 0.73); }
double fun_levelSet(R2 P) { return -(P[1] + 0.5 * P[0] - 1.73); }

// double fun_initial(R2 P, int elementComp, int domain) { return sin(pi * (P[0] + P[1])); }
// double fun_boundary(R2 P, int elementComp, int domain, double t) { return sin(pi * (P[0] + P[1] - 4 * t)); }

double fun_initial(R2 P, int elementComp, int domain) {
    double xs = 0.5, ys = 0.5;
    double r = 0.3;
    double v = (P[0] - xs) * (P[0] - xs) + (P[1] - ys) * (P[1] - ys);
    if (v < r * r)
        return exp(-pow(P[0] - xs, 2) / (0.5 * r) - pow(P[1] - ys, 2) / (0.5 * r));
    else
        return 0.;
}
double fun_boundary(R2 P, int elementComp, int domain, double t) { return 0.; }

void assembly(const space_t &Wh, const Interface<mesh_t> &interface, MatMap &Ah, MatMap &Mh) {

    double t0 = getTime();
    const ActiveMesh<mesh_t> &Khi(Wh.get_mesh());
    CutFEM<mesh_t> problem(Wh);
    CutFEM_R2 beta({R2(3, 1), R2(3, 1)});
    CutFEMParameter lambda(0, 1.);
    CutFEMParameter lambdaE(sqrt(10), sqrt(10)); // max B = sqrt(10)/sqrt(5)
    const MeshParameter &h(Parameter::h);

    double lambdaB = sqrt(10);
    double Cstab   = 1e-2;
    double Cstabt  = 1e-2;

    Normal n;
    funtest_t u(Wh, 1), v(Wh, 1);

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
    problem.addFaceStabilization(
        innerProduct(h * jump(u), Cstab * jump(v)) + innerProduct((h ^ 3) * jump(grad(u)), Cstab * jump(grad(v))), Khi);

    std::cout << " Time assembly \t" << getTime() - t0 << std::endl;
}

void solve_problem(const space_t &Wh, const Interface<mesh_t> &interface, const std::span<double> u0,
                   std::span<double> uh, MatMap &Ah, MatMap &Mh, double tn) {

    double t0 = getTime();

    const ActiveMesh<mesh_t> &Khi(Wh.get_mesh());

    ProblemOption optionProblem;
    optionProblem.solver_name_  = "umfpack";
    optionProblem.clear_matrix_ = false;
    CutFEM<mesh_t> problem(Wh, optionProblem);

    CutFEM_R2 beta({R2(3, 1), R2(3, 1)});
    double lambdaB = sqrt(10);

    Normal n;
    funtest_t u(Wh, 1), v(Wh, 1);
    fct_t gh(Wh, fun_boundary, tn);

    // MULTIPLYING A * u_0  => rhs
    // =====================================================
    int N = problem.get_nb_dof();
    multiply(N, N, Ah, u0, problem.rhs_);

    problem.addLinear(-innerProduct(gh.expr(), beta * (0.5 * v) * n) + innerProduct(gh.expr(), lambdaB * 0.5 * v), Khi,
                      INTEGRAL_BOUNDARY, {1, 4} // label left boundary
    );

    // SOLVING  (M+S)E'(t) = rhs
    // =====================================================
    problem.solve(Mh, problem.rhs_);

    std::copy(problem.rhs_.begin(), problem.rhs_.end(), uh.begin());
}

int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);
    const double cpubegin = getTime();

    // OUTPUT FILE
    // =====================================================
    std::ofstream outputData("output_testP1.txt");

    // DEFINITION OF THE MESH and SPACE
    // =====================================================
    int nx = 20;
    int ny = 20;
    mesh_t Th(nx, ny, 0., 0., 2., 2.);

    space_t Vh(Th, DataFE<mesh_t>::P1dc);

    // DEFINITION OF SPACE AND TIME PARAMETERS
    // =====================================================
    double tid      = 0;
    double meshSize = 2. / (nx + 1);
    double dt       = meshSize / 3 / sqrt(10) * 0.5; // 0.3 * meshSize / 3 ;  // h /10
    double tend     = 0.3;
    int niteration  = tend / dt;
    dt              = tend / niteration;
    double errSum   = 0;

    std::cout << "mesh_t size h = \t" << meshSize << std::endl;
    std::cout << "Time step dt = \t" << dt << std::endl;

    // DEFINITION OF THE LEVELSET
    // =====================================================
    space_t Lh(Th, DataFE<mesh_t>::P1);
    fct_t levelSet(Lh, fun_levelSet);

    Lagrange2 FE_beta(1);
    space_t Uh(Th, FE_beta);

    // CONSTRUCTION INTERFACE AND CUTSPACE
    // =====================================================
    InterfaceLevelSet<mesh_t> interface(Th, levelSet);
    cutmesh_t Khi(Th, interface);
    cutspace_t Wh(Khi, Vh);
    cutspace_t Qh(Khi, Uh);

    MacroElement<mesh_t> macro(Khi, 0.5);

    // DECLARATION OF THE VECTOR CONTAINING THE solution
    // =====================================================
    std::vector<double> u0(Wh.NbDoF(), 0.);
    interpolate(Wh, u0, fun_initial);
    std::vector<double> uh(u0);
    std::vector<double> uh_tild(u0);
    fct_t fun_u0(Wh, u0);

    // Plot the macro elements
    // {
    //     Paraview<mesh_t> writer(Khi, "limiter_smooth_solution_0.vtk");
    //     writer.add(fun_u0, "uhLimiter", 0, 1);
    //     writer.add(fun_u0, "uhNoLimiter", 0, 1);
    //     writer.add(levelSet, "levelSet", 0, 1);
    //     writer.writeMacroInnerEdge(macro, 0, "limiter_smooth_solution_macro_inner_edge1.vtk");
    //     writer.writeMacroOutterEdge(macro, 0, "limiter_smooth_solution_macro_outter_edge1.vtk");
    //     writer.writeMacroInnerEdge(macro, 1, "limiter_smooth_solution_macro_inner_edge2.vtk");
    //     writer.writeMacroOutterEdge(macro, 1, "limiter_smooth_solution_macro_outter_edge2.vtk");
    // }

    double qu0           = integral(Khi, fun_u0, 0);
    const auto minmax_u0 = std::minmax_element(uh.begin(), uh.end());
    double min_u0        = *minmax_u0.first;
    double max_u0        = *minmax_u0.second;

    // ASSEMBLY THE CONSTANT PART
    // ==================================================
    MatMap Ah, Mh;
    assembly(Wh, interface, Ah, Mh);
    // RESOLUTION OF THE PROBLEM_MIXED_DARCY
    // ==================================================

    int ifig = 1;
    for (int i = 0; i < niteration; ++i) {
        std::cout << " ------------------------------ " << std::endl;
        std::cout << "Iteration " << i + 1 << " / " << niteration << " \t time = " << tid << std::endl;

        // THIRD ORDER RK
        // =================================================

        std::vector<double> u1(Wh.NbDoF(), 0.);
        fct_t fun_u1(Wh, u1);
        fct_t fun_uh(Wh, uh);
        std::map<int, double> u_mean;

        if (Wh.basisFctType == BasisFctType::P1dc) {
            solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
            Scheme::RK3::step1(u0.begin(), u0.end(), uh.begin(), dt, u1.begin());
            std::vector<double> u1_tild = limiter::CutFEM::applyBoundPreservingLimiter(fun_u1, min_u0, max_u0, macro);

            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + dt);
            Scheme::RK3::step2(u0.begin(), u0.end(), u1_tild.begin(), uh.begin(), dt, u1.begin());
            u1_tild = limiter::CutFEM::applyBoundPreservingLimiter(fun_u1, min_u0, max_u0, macro);

            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + 0.5 * dt);
            Scheme::RK3::step3(u0.begin(), u0.end(), u1_tild.begin(), uh.begin(), dt, uh.begin());
            uh_tild = limiter::CutFEM::applyBoundPreservingLimiter(fun_uh, min_u0, max_u0, macro);

        } else if (Wh.basisFctType == BasisFctType::P0) {

            solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
            Scheme::RK3::step1(u0.begin(), u0.end(), uh.begin(), dt, u1.begin());
            std::vector<double> u1_tild = limiter::CutFEM::extendToMacro(fun_u1, u_mean, macro);

            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + dt);
            Scheme::RK3::step2(u0.begin(), u0.end(), u1_tild.begin(), uh.begin(), dt, u1.begin());
            u1_tild = limiter::CutFEM::extendToMacro(fun_u1, u_mean, macro);

            solve_problem(Wh, interface, u1_tild, uh, Ah, Mh, tid + 0.5 * dt);
            Scheme::RK3::step3(u0.begin(), u0.end(), u1_tild.begin(), uh.begin(), dt, uh.begin());
            uh_tild = limiter::CutFEM::extendToMacro(fun_uh, u_mean, macro);
        } else {
            assert(0);
        }

        u0 = uh_tild;
        tid += dt;

        // COMPUTATION OF THE L2 ERROR
        // =================================================
        // fct_t fun_uh(Wh, uh);
        fct_t fun_uh_tild(Wh, uh_tild);
        double qu             = integral(Khi, fun_uh_tild, 0);
        auto [min_uh, max_uh] = limiter::CutFEM::findMinAndMaxValue(fun_uh);
        auto [min_u1, max_u1] = limiter::CutFEM::findMinAndMaxValue(fun_uh_tild);

        if (min_u0 <= min_u1 + Epsilon && max_u1 <= max_u0 + Epsilon) {

            std::cout << " Maximum principle satified! " << std::endl;
            std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
            std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
            std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;
        } else {
            std::cout << " Maximum principle not satified! " << std::endl;
            std::cout << "[m, M] = [ " << min_u0 << " , " << max_u0 << " ]" << std::endl;
            std::cout << min_uh << " < u_{h,M}  < " << max_uh << std::endl;
            std::cout << min_u1 << " < u1_{h,M} < " << max_u1 << std::endl;
        }

        // PLOT THE SOLUTION
        // ==================================================
        if (MPIcf::IamMaster() && i % 5 == 0 || i + 1 == niteration) {

            Paraview<mesh_t> writer(Khi, "limiter_non_smooth_solution100_" + std::to_string(ifig++) + ".vtk");
            writer.add(fun_uh, "uhNoLimiter", 0, 1);
            writer.add(fun_uh_tild, "uhLimiter", 0, 1);
        }

        std::cout << std::setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu - qu0) << std::setprecision(6)
                  << std::endl;
        outputData << i << "\t" << tid << "\t" << std::setprecision(16) << qu << "\t" << std::setprecision(16)
                   << fabs(qu - qu0) << "\t" << std::setprecision(5) << std::endl;
    }

    std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
    std::cout << " Time computation \t" << getTime() - cpubegin << std::endl;
    return 0;
}
