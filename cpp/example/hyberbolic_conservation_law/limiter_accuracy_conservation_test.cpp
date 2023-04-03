#include "../tool.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;
using MatMap     = std::map<std::pair<int, int>, double>;

using namespace globalVariable;

double fun_levelSet(R2 P) { return P[0] * P[0] + P[1] * P[1] - 1; }
double fun_initial(R2 P, int elementComp, int domain) {
    double r0 = 0.4;
    double a  = 0.3;
    return 0.5 * (1 - tanh(((P[0] - r0) * (P[0] - r0) + P[1] * P[1]) / (a * a) - 1));
}
double fun_theta(R2 P) {
    double r = Norme2(P);
    int s_y  = util::fsign(P[1]);
    if (P[0] < -Epsilon)
        return (s_y * pi + atan(P[1] / P[0]));
    else if (P[0] > Epsilon)
        return atan(P[1] / (P[0]));
    else
        // return s_y * pi / 2;
        return atan(P[1] / (P[0] + Epsilon));
}

double fun_solution(R2 P, int elementComp, int domain, double t) {
    double r     = Norme2(P);
    double theta = fun_theta(P);
    R2 Q(r * cos(theta - 2 * pi * t), r * sin(theta - 2 * pi * t));
    return fun_initial(Q, 0, 0);
}
double fun_boundary(R2 P, int elementComp, double t) { return 0.; }
double fun_velocity(R2 P, int elementComp, int domain) { return (elementComp == 0) ? -2 * pi * P[1] : 2 * pi * P[0]; }

void assembly(const space_t &Wh, const Interface<mesh_t> &interface, MatMap &Ah, MatMap &Mh, const fct_t &beta) {

    double t0 = getTime();
    const ActiveMesh<mesh_t> &Khi(Wh.get_mesh());
    CutFEM<mesh_t> problem(Wh);
    CutFEMParameter lambda(0, 1.);
    double lambdaE = sqrt(2) * pi;

    const MeshParameter &h(Parameter::h);

    double Cstab  = 1e-2;
    double Cstabt = 1e-2;

    Normal n;
    funtest_t u(Wh, 1), v(Wh, 1);

    // BUILDING A
    // =====================================================
    problem.set_map(Ah);
    double tt = getTime();

    problem.addBilinear(innerProduct(beta.exprList(2) * u, grad(v)), Khi);

    // problem.addBilinear(
    //   - jump(innerProduct(beta*u,v*n))
    //   - innerProduct(jump(beta*u*n), jump(lambda*v))
    //   , interface
    // );
    //[ u \beta \cdot n]=3u_1-2u_2=0
    problem.addFaceStabilization(-innerProduct(jump(u), Cstab * jump(v)) -
                                     innerProduct((h ^ 2) * jump(grad(u)), Cstab * jump(grad(v)))
                                 // - contractProduct((h^4)*jump(Hessu), Cstab*jump(Hessv))
                                 ,
                                 Khi);

    // F(u)_e = {B.u} - lambda_e/2 [u]
    problem.addBilinear(-innerProduct(average(beta.exprList(2) * u * n), jump(v)) -
                            innerProduct(0.5 * lambdaE * jump(u), jump(v)),
                        Khi, INTEGRAL_INNER_EDGE_2D);

    // problem.addBilinear(-innerProduct(beta.exprList(2) * u * n, v), interface
    //                     // , Khi
    //                     // , boundary
    //                     // , {2, 3}          // label other boundary
    // );
    // problem.addBilinear(-innerProduct(u, beta.exprList(2) * (0.5 * v) * n) -
    //                         innerProduct(u, lambdaB * (0.5 * v)),
    //                     interface
    //                     // , Khi
    //                     // , INTEGRAL_BOUNDARY
    //                     // , {1, 4}          // label left boundary
    // );
    // problem.addBilinear(-innerProduct(beta.exprList(2) * u, v * n) -
    //                         innerProduct(beta.exprList(2) * u * n, lambdaB *
    //                         v),
    //                     interface);

    // BUILDING (M + S)
    // =====================================================
    problem.set_map(Mh);
    problem.addBilinear(innerProduct(u, v), Khi);
    problem.addFaceStabilization(innerProduct(h * jump(u), Cstab * jump(v)) +
                                     innerProduct((h ^ 3) * jump(grad(u)), Cstab * jump(grad(v)))
                                 // + contractProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
                                 ,
                                 Khi);

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

    // CutFEM_R2 beta({R2(1,1), R2(1,1)});
    // double lambdaB = sqrt(10);

    Normal n;
    funtest_t u(Wh, 1), v(Wh, 1);
    // fct_t gh(Wh, fun_boundary, tn);
    // Expression2 gx(gh, 0, op_id);

    // MULTIPLYING A * u_0  => rhs
    // =====================================================
    // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
    // mAh.addMatMul(u0, problem.rhs_);
    int N = problem.get_nb_dof();
    multiply(N, N, Ah, u0, problem.rhs_);

    // problem.addLinear(
    //    - innerProduct(gx, beta*  (0.5*v)*n)
    //    + innerProduct(gx, lambdaB*0.5*v)
    //    , Khi
    //    , boundary
    //   , {1,4}          // label left boundary
    // );

    // SOLVING  (M+S)E'(t) = rhs
    // =====================================================
    problem.solve(Mh, problem.rhs_);

    std::copy(problem.rhs_.begin(), problem.rhs_.end(), uh.begin());
}

int main(int argc, char **argv) {

#ifdef USE_MPI
    MPIcf cfMPI(argc, argv);
#endif
    const double cpubegin = getTime();

    // OUTPUT FILE
    // =====================================================
    std::ofstream outputData("output_smooth2_nx20_P0.txt");

    // DEFINITION OF THE MESH and SPACE
    // ====================================================
    int nx = 40;                                       // 160;
    int ny = 40;                                       // 160;
    mesh_t Th(nx, ny, -1.0075, -1.0075, 2.015, 2.015); // [-1,1]*[-1,1]
    space_t Vh(Th, DataFE<mesh_t>::P1dc);

    // DEFINITION OF SPACE AND TIME PARAMETERS
    // =====================================================
    double tid      = 0;
    double meshSize = 2. / nx;
    double dt       = meshSize / 5 / sqrt(10) * 0.5; // 0.3 * meshSize / 3 ;  // h /10
    double tend     = 1.;
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
    ActiveMesh<mesh_t> Khi(Th);
    Khi.truncate(interface, 1);
    cutspace_t Wh(Khi, Vh);
    cutspace_t Qh(Khi, Uh);
    fct_t beta(Qh, fun_velocity);

    MacroElement<mesh_t> macro(Khi, 0.5);

    // DECLARATION OF THE VECTOR CONTAINING THE solution
    // =====================================================
    std::vector<double> u0(Wh.NbDoF(), 0.);
    interpolate(Wh, u0, fun_initial);
    std::vector<double> uh(u0);
    std::vector<double> uh_tild(u0);
    fct_t fun_u0(Wh, u0);

    // Plot the macro elements
    {
        Paraview<mesh_t> writer(Khi, "limiter_accuracy_test_0.vtk");
        writer.add(fun_u0, "uhLimiter", 0, 1);
        writer.add(levelSet, "levelSet", 0, 1);
        writer.writeMacroInnerEdge(macro, 0, "pei_macro_inner_edge1.vtk");
        writer.writeMacroOutterEdge(macro, 0, "pei_macro_outter_edge1.vtk");
    }

    double qu0           = integral(Khi, fun_u0, 0);
    const auto minmax_u0 = std::minmax_element(uh.begin(), uh.end());
    double min_u0        = *minmax_u0.first;
    double max_u0        = *minmax_u0.second;

    // ASSEMBLY THE CONSTANT PART
    // ==================================================
    MatMap Ah, Mh;
    assembly(Wh, interface, Ah, Mh, beta);

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

        // COMPUTATION OF THE L2 ERROR
        // =================================================
        std::vector<double> usol(Wh.NbDoF(), 0.);
        interpolate(Wh, usol, fun_solution, tid);
        std::vector<double> uerr(usol);
        std::transform(uh_tild.begin(), uh_tild.end(), uerr.begin(), uerr.begin(),
                       [](double a, double b) { return a - b; });

        fct_t femErrh(Wh, uerr);
        fct_t fun_ex(Wh, usol);

        double errU = sqrt(integral(Khi, (femErrh.expr() ^ 2)));
        errSum += errU;

        // PLOT THE SOLUTION
        // ==================================================
#ifdef USE_MPI
        if (MPIcf::IamMaster() && i % 10000 == 0 || i + 1 == niteration) {
#else
        if (i % 1 == 0 || i + 1 == niteration) {
#endif

            Paraview<mesh_t> writer(Khi, "test_accuracyP0_" + std::to_string(ifig++) + ".vtk");
            writer.add(fun_uh, "uhNoLimiter", 0, 1);
            writer.add(fun_uh_tild, "uhLimiter", 0, 1);
            writer.add(femErrh, "error", 0, 1);
            writer.add(fun_ex, "u_exact", 0, 1);
            // writer.add(fun_thet , "yx", 0, 1);
        }
        // if (i == 10)
        //    return 0;
        // std::cout << std::setprecision(16) << "q(u) = " << qu << std::endl;
        std::cout << " || u-uex ||_2 = " << errU << std::endl;
        std::cout << std::setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu - qu0) << std::setprecision(6)
                  << std::endl;

        if (i == 0) {
            outputData << "it \t tid \t errU \t q \t\t errQ \t\t minU \t\t maxU " << std::endl;
        }
        outputData << i << "\t" << tid << "\t" << std::setprecision(16) << errU << "\t" << std::setprecision(16) << qu
                   << "\t" << std::setprecision(16) << fabs(qu - qu0) << "\t" << std::setprecision(16) << min_u1 << "\t"
                   << std::setprecision(16) << max_u1 << "\t" << std::setprecision(5) << std::endl;
    }

    std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
    std::cout << " Time computation \t" << getTime() - cpubegin << std::endl;
    return 0;
}
