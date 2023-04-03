#ifdef EXAMPLE1
double c0 = 0.5;
R fun_levelSet(const R2 P, const int i) { return -P.x - P.y + c0; }
R fun_boundary(const R2 P, int elementComp, double t) { return sin(pi * (P.x + P.y - 4 * t)); }
R fun_solution(const R2 P, int elementComp, int domain, double t) {
    // return 2*sin(pi*(P.x+P.y-3*t));
    if (domain == 0)
        return sin(pi * (P.x + P.y - 4 * t));
    else
        return 4. / 3 * sin(4. / 3 * pi * (P.x + P.y - 3 * t - c0 / 4));
}

void assembly(const Space &Wh, const Interface<Mesh> &interface, MatMap &Ah, MatMap &Mh) {

    double t0 = MPIcf::Wtime();
    const ActiveMesh<Mesh> &Khi(Wh.get_mesh());
    CutFEM<Mesh2> problem(Wh);
    CutFEM_R2 beta({R2(3, 1), R2(2, 1)});
    CutFEMParameter lambda(0., 1.);

    const MeshParameter &h(Parameter::h);
    CutFEMParameter lambdaE(3, 2); // max B = sqrt(10)/sqrt(5)

    double lambdaB = 3;
    double Cstab   = 5e-1;
    double Cstabt  = 5e-1;

    Normal n;
    TestFunction2 u(Wh, 1), v(Wh, 1);

    TestFunction2 Hessu = grad(grad(u)), Hessv = grad(grad(v));
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
                        Khi, innerFacet);

    problem.addBilinear(-innerProduct(beta * u * n, v), Khi, boundary, {2, 3} // label other boundary
    );
    problem.addBilinear(-innerProduct(u, beta * (0.5 * v) * n) - innerProduct(u, lambdaB * 0.5 * v), Khi, boundary,
                        {1, 4} // label left boundary
    );

    // BUILDING (M + S)
    // =====================================================
    problem.set_map(Mh);
    problem.addBilinear(innerProduct(u, v), Khi);
    problem.addFaceStabilization(innerProduct(h * jump(u), Cstab * jump(v)) +
                                     innerProduct((h ^ 3) * jump(grad(u)), Cstab * jump(grad(v))) +
                                     contractProduct((h ^ 5) * jump(Hessu), Cstab * jump(Hessv)),
                                 Khi);

    std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}

void solve_problem(const Space &Wh, const Interface<Mesh> &interface, const Rn &u0, Rn &uh, MatMap &Ah, MatMap &Mh,
                   double tn) {

    double t0 = MPIcf::Wtime();
    const ActiveMesh<Mesh> &Khi(Wh.get_mesh());

    ProblemOption optionProblem;
    optionProblem.clear_matrix_ = false;
    CutFEM<Mesh2> problem(Wh, optionProblem);

    CutFEM_R2 beta({R2(3, 1), R2(2, 1)});
    double lambdaB = 3;

    Fun_h gh(Wh, fun_solution, tn);
    Expression2 gx(gh, 0, op_id);

    Normal n;
    FunTest u(Wh, 1), v(Wh, 1);

    // MULTIPLYING A * u_0  => rhs
    // =====================================================
    // problem.mat = Ah;
    // problem.addMatMul(u0);
    int N = problem.get_nb_dof();
    multiply(N, N, Ah, u0, problem.rhs_);

    problem.addLinear(-innerProduct(gx, beta * (0.5 * v) * n) + innerProduct(gx, lambdaB * 0.5 * v), Khi, boundary,
                      {1, 4} // label left boundary
    );
    // std::cout << " Time building problem \t" << MPIcf::Wtime() - t0 << std::endl;

    // SOLVING  (M+S)E'(t) = rhs
    // =====================================================
    problem.solve(Mh, problem.rhs_);

    uh = problem.rhs_;
}

int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);
    const double cpubegin = MPIcf::Wtime();

    // OUTPUT FILE
    // =====================================================
    std::ofstream outputData("output_test.txt");

    // DEFINITION OF THE MESH and SPACE
    // =====================================================
    int nx = 50;
    int ny = 50;
    Mesh Th(nx, ny, -1., -1., 2., 2.); // [-1,1]*[-1,1]
    Space Vh(Th, DataFE<Mesh2>::P1dc);

    // DEFINITION OF SPACE AND TIME PARAMETERS
    // =====================================================
    double tid      = 0;
    double meshSize = 2. / nx;
    double dt       = meshSize / 3 / sqrt(10) * 0.5; // 0.3 * meshSize / 3 ;  // h /10
    double tend     = 0.5;
    int niteration  = tend / dt;
    dt              = tend / niteration;
    double errSum   = 0;
    double qu0      = 0;
    std::cout << "Mesh size h = \t" << meshSize << std::endl;
    std::cout << "Time step dt = \t" << dt << std::endl;

    // DEFINITION OF THE LEVELSET
    // =====================================================
    Space Lh(Th, DataFE<Mesh2>::P1);
    Fun_h levelSet(Lh, fun_levelSet);

    // CONSTRUCTION INTERFACE AND CUTSPACE
    // =====================================================
    InterfaceLevelSet<Mesh> interface(Th, levelSet);
    ActiveMesh<Mesh> Khi(Th, interface);
    CutSpace Wh(Khi, Vh);

    // DECLARATION OF THE VECTOR CONTAINING THE solution
    // =====================================================
    Rn u0(Wh.NbDoF(), 0.);
    interpolate(Wh, u0, fun_solution, 0.);
    Fun_h Un(Wh, u0);
    Rn uh(u0);

    // Fun2_h femSolh(Wh, uh);

    // if(MPIcf::IamMaster()) {
    //   Fun2_h solex(Wh, fun_solution, tid);
    //   Fun2_h sol(Wh, uh);
    //   Paraview2 writer(Wh, levelSet, "testDetector.vtk");
    //   writer.add(sol  , "uh" , 0, 1);
    //   // writer.add(levelSet, "levelSet", 0, 1);
    // }

    // ASSEMBLY THE CONSTANT PART
    // ==================================================
    MatMap Ah, Mh;
    assembly(Wh, interface, Ah, Mh);

    // RESOLUTION OF THE PROBLEM_MIXED_DARCY
    // ==================================================
    int ifig = 1;
    for (int i = 0; i < niteration; ++i) {

        // EULER METHOD
        // =================================================
        // solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
        // std::cout << uh << std::endl;
        //
        // uh *= dt;
        // uh += u0;

        Rn u1(u0);
        Rn u2(Wh.NbDoF(), 0.);
        // THIRD ORDER RK
        // =================================================
        solve_problem(Wh, interface, u0, uh, Ah, Mh, tid);
        u1 += dt * uh;
        u2 += 3. / 4 * u0 + 1. / 4 * u1;
        solve_problem(Wh, interface, u1, uh, Ah, Mh, tid + dt);
        u2 += 1. / 4 * dt * uh;
        solve_problem(Wh, interface, u2, uh, Ah, Mh, tid + 0.5 * dt);
        uh *= 2. / 3 * dt;
        uh += 1. / 3 * u0 + 2. / 3 * u2;

        u0 = uh;
        tid += dt;

        // COMPUTATION OF THE L2 ERROR
        // =================================================
        interpolate(Wh, u1, fun_solution, tid);
        u1 -= uh;
        Fun_h femErrh(Wh, u1);
        Expression2 femErr(femErrh, 0, op_id);

        R errU = sqrt(integral(Khi, femErr * femErr));
        errSum += errU;

        // Expression2 femSol(femSolh, 0, op_id);
        // R qu = integral(Khi, femSol);

        // PLOT THE SOLUTION
        // ==================================================
        if (MPIcf::IamMaster() && i % 10 == 0 || i + 1 == niteration) {
            Fun2_h femSolh(Wh, uh);
            // Fun2_h solex(Wh, fun_solution, tid);
            for (int j = 0; j < uh.size(); ++j) {
                if (fabs(uh(j)) < 1e-16)
                    uh(j) = 0.;
            }
            Fun2_h sol(Wh, uh);
            Paraview<Mesh> writer(Khi, "hyperbolicLawExample1_" + to_string(ifig++) + ".vtk");
            writer.add(sol, "uh", 0, 1);
            // writer.add(solex, "uex", 0, 1);
        }
        std::cout << "Iteration " << i << " / " << niteration << " \t time = " << tid << " , || u-uex ||_2 = " << errU
                  << std::endl;
        // outputData << i << "\t"
        //            << tid << "\t"
        //            << setprecision(16) << qu << "\t"
        //            << setprecision(16) << fabs(qu-qu0) << "\t"
        //            << setprecision(5) << std::endl;

        // return 0;
    }

    std::cout << "Error  - sum || u ||_2 / N = " << errSum / niteration << std::endl;
    std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
    return 0;
}

#endif
