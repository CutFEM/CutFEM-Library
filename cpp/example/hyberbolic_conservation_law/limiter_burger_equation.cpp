#include "../cutfem.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;
using MatMap     = std::map<std::pair<int, int>, double>;

using namespace globalVariable;

double fun_levelSet(R2 P, ) { return -(P[1] - 0.5 * P[0] - 1. / 4); }
double fun_initial(R2 P, int elementComp, int domain) { return 1 + 0.5 * sin(pi * (P[0] + P[1])); }
double fun_scalar(double x) { return 1 + 0.5 * sin(pi * x); }
double fun_solution(R2 P, int elementComp, int domain, double t) {
    double x0 = P[0] + P[1];
    double xy = x0;
    for (int i = 0; i <= 50; ++i) {

        double x1 = xy - 2 * t * fun_scalar(x0);
        if (fabs(x1 - x0) < 1e-10)
            return fun_scalar(x1);
        x0 = x1;
        if (i == 100) {
            std::cout << " Burger Newton not converged \t" << fabs(xy - 2 * t * fun_scalar(x0) - x0) << std::endl;

            getchar();
        }
    }
}

void assembly(const fespace_t &Wh, const Interface<mesh_t> &interface, MatMap &Ah, MatMap &Mh) {

    double t0 = getTime();
    const ActiveMesh<mesh_t> &Khi(Wh.get_mesh());

    CutFEM<Mesh2> problem(Wh);
    CutFEMParameter lambda(0, 1.);
    CutFEMParameter lambdaE(sqrt(10), sqrt(10)); // max B = sqrt(10)/sqrt(5)
    // CutFEMParameter lambdaE(3., 3.);  // max B = sqrt(10)/sqrt(5)

    const MeshParameter &h(Parameter::h);
    double Cstab  = 1e-2;
    double Cstabt = 1e-2;
    Normal n;
    FunTest u(Wh, 1), v(Wh, 1);
    // double lambdaB = 2.;//sqrt(10);

    // BUILDING A
    // =====================================================
    problem.set_map(Ah);

    // F(u)_e = {B.u} - lambda_e/2 [u]
    problem.addBilinear(-innerProduct(0.5 * lambdaE * jump(u), jump(v))
                        // - innerProduct(average(beta*u*n), jump(v))
                        ,
                        Khi, INTEGRAL_INNER_EDGE_2D);

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

    // problem.addBilinear(
    //   - innerProduct(u, beta*(0.5*v)*n)
    //   - innerProduct(u, lambdaB*0.5*v)
    //   , Khi
    //   , INTEGRAL_BOUNDARY
    //   , {1, 4}          // label left boundary
    // );

    // BUILDING (M + S)
    // =====================================================
    problem.set_map(Mh);
    problem.addBilinear(innerProduct(u, v), Khi);
    problem.addFaceStabilization(
        innerProduct(h * jump(u), Cstab * jump(v)) + innerProduct((h ^ 3) * jump(grad(u)), Cstab * jump(grad(v))), Khi);

    std::cout << " Time assembly \t" << getTime() - t0 << std::endl;
}
void solve_problem(const fespace_t &Wh, const Interface<mesh_t> &interface, Rn &u0, Rn &uh, MatMap &Ah, MatMap &Mh,
                   double tn) {

    double t0 = getTime();

    const ActiveMesh<mesh_t> &Khi(Wh.get_mesh());
    ProblemOption optionProblem;
    optionProblem.solver_name_  = "umfpack";
    optionProblem.clear_matrix_ = false;
    CutFEM<Mesh2> problem(Wh, optionProblem);

    // CutFEM_R2 beta({R2(1,1), R2(1,1)});
    double lambdaB = sqrt(10);
    // CutFEMParameter lambdaE(sqrt(10), sqrt(10));  // max B = sqrt(10)/sqrt(5)
    // double lambdaB = 3.;//sqrt(10);
    CutFEMParameter lambda(0, 1.);

    Normal N;
    FunTest v(Wh, 1);
    Fun_h fun_Un(Wh, u0);
    Expression Un(fun_Un, 0, op_id);

    Fun_h gh(Wh, fun_solution, tn);
    Expression gx(gh, 0, op_id);

    CutFEM_R2 beta({R2(3, 1), R2(3, 1)});

    // MULTIPLYING A * u_0  => rhs
    // =====================================================
    MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
    mAh.addMatMul(u0, problem.rhs_);

    // matlab::Export(problem.rhs_, "rhs0.dat");
    // problem.rhs_ = 0.;
    // CONSTRUCT THE RHS
    double tt = getTime();

    problem.addLinear(innerProduct(fluxFctX(Un), dx(v)) + innerProduct(fluxFctY(Un), dy(v)), Khi);

    problem.addLinear(
        // - innerProduct(0.5*lambdaE*jump(Un)  , jump(v))
        -innerProduct(average(fluxFctN(Un)), jump(v)), Khi, INTEGRAL_INNER_EDGE_2D);

    problem.addLinear(-jump(innerProduct(fluxFctN(Un), v)) - innerProduct(jump(fluxFctN(Un)), jump(lambda * v)),
                      interface);
    // matlab::Export(problem.rhs_, "rhs1.dat");

    // BOUNDARY CONDITION IN RHS
    problem.addLinear(-innerProduct(fluxFctN(Un), 0.5 * v) - innerProduct(Un, lambdaB * 0.5 * v), Khi, INTEGRAL_BOUNDARY
                      // , {1,4}          // boundary in
    );
    problem.addLinear(-innerProduct(fluxFctN(gx), 0.5 * v) + innerProduct(gx, lambdaB * 0.5 * v), Khi, INTEGRAL_BOUNDARY
                      // , {1,4}          // label left boundary
    );

    // problem.addLinear(
    //    - innerProduct(fluxFctN(Un), v)
    //   , Khi
    //   , INTEGRAL_BOUNDARY
    //   , {2,3}        // boundary out
    // );

    // matlab::Export(problem.rhs_, "rhs1.dat");
    // std::cout << "hey " << std::endl;
    // getchar();

    // SOLVING  (M+S)E'(t) = rhs
    // =====================================================
    tt = getTime();
    problem.solve(Mh, problem.rhs_);
    // std::cout << "Time solver \t" << getTime()-tt << std::endl;

    uh           = problem.rhs_;
    problem.rhs_ = 0.;
}
