
double c0 = 0.25;
R fun_levelSet(const R2 P, const int i) { return -P.x - P.y + c0; }
R fun_initial(const R2 P, int elementComp, int domain) {

    if (domain == 1)
        return 0;
    double xs = -0.3, ys = -0.3;
    double r = 0.3;
    double v = (P.x - xs) * (P.x - xs) + (P.y - ys) * (P.y - ys);
    if (v < r * r)
        return 1; // exp(-pow(P.x - xs,2)/r - pow(P.y - ys,2)/r);
    else
        return 0;
}

void assembly(const Space &Wh, const Interface<Mesh> &interface, MatMap &Ah, MatMap &Mh) {

    double t0 = MPIcf::Wtime();
    const ActiveMesh<Mesh> &Khi(Wh.get_mesh());
    CutFEM<Mesh2> problem(Wh);
    CutFEM_R2 beta({R2(3, 1), R2(1, 2)});
    CutFEMParameter lambda(0, 1.);
    CutFEMParameter lambdaE(sqrt(10), sqrt(5)); // max B = sqrt(10)/sqrt(5)
    const MeshParameter &h(Parameter::h);
    // double h = 2./40;

    double lambdaB = sqrt(10);
    double Cstab   = 1e-2;
    double Cstabt  = 1e-2;

    Normal n;
    FunTest u(Wh, 1), v(Wh, 1);

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
                                     innerProduct((h ^ 3) * jump(grad(u)), Cstab * jump(grad(v)))
                                 // + contractProduct((h^5)*jump(Hessu), Cstab*jump(Hessv))
                                 ,
                                 Khi);

    std::cout << " Time assembly \t" << MPIcf::Wtime() - t0 << std::endl;
}
void solve_problem(const Space &Wh, const Interface<Mesh> &interface, const Rn &u0, Rn &uh, MatMap &Ah, MatMap &Mh,
                   double tn) {

    double t0 = MPIcf::Wtime();

    const ActiveMesh<Mesh> &Khi(Wh.get_mesh());

    ProblemOption optionProblem;
    optionProblem.solver_name_  = "umfpack";
    optionProblem.clear_matrix_ = false;
    CutFEM<Mesh2> problem(Wh, optionProblem);

    CutFEM_R2 beta({R2(3, 1), R2(1, 2)});
    double lambdaB = sqrt(10);

    Normal n;
    FunTest u(Wh, 1), v(Wh, 1);

    // MULTIPLYING A * u_0  => rhs
    // =====================================================
    // MatriceMap<double> mAh(problem.nb_dof_, problem.nb_dof_, Ah);
    // mAh.addMatMul(u0, problem.rhs_);
    int N = problem.get_nb_dof();
    multiply(N, N, Ah, u0, uh);

    // SOLVING  (M+S)E'(t) = rhs
    // =====================================================
    problem.solve(Mh, uh); // problem.rhs_);

    // uh = problem.rhs_;
}

int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);
    const double cpubegin = MPIcf::Wtime();

    // OUTPUT FILE
    // =====================================================
    std::ofstream outputData("output_test.txt");

    // DEFINITION OF THE MESH and SPACE
    // =====================================================
    int nx = 200;
    int ny = 200;
    Mesh Th(nx, ny, -1., -1., 2., 2.); // [-1,1]*[-1,1]
    Space Vh(Th, DataFE<Mesh2>::P1dc);

    // DEFINITION OF SPACE AND TIME PARAMETERS
    // =====================================================
    double tid      = 0;
    double meshSize = 2. / nx;
    double dt       = meshSize / 5 / sqrt(10) * 0.5; // 0.3 * meshSize / 3 ;  // h /10
    double tend     = 0.4;
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
    interpolate(Wh, u0, fun_initial);
    Fun_h Un(Wh, u0);
    Rn uh(u0);

    Fun_h fun_uh(Wh, uh);
    // Expression2 femSol(femSolh, 0, op_id);
    qu0 = integral(Khi, fun_uh, 0);

    // ASSEMBLY THE CONSTANT PART
    // ==================================================
    MatMap Ah, Mh;
    assembly(Wh, interface, Ah, Mh);

    // RESOLUTION OF THE PROBLEM_MIXED_DARCY
    // ==================================================
    int ifig = 1;
    for (int i = 0; i < niteration; ++i) {

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
        // Fun2_h femSolh(Wh, uh);
        // Expression2 femSol(femSolh, 0, op_id);
        R qu = integral(Khi, fun_uh, 0);

        // PLOT THE SOLUTION
        // ==================================================
        // if(MPIcf::IamMaster() && i%5 == 0 || i+1 == niteration) {
        //   // Fun2_h solex(Wh, fun_solution, tid);
        //   for(int j=0;j<uh.size();++j) {
        //     if(fabs(uh(j)) < 1e-16 ) uh(j) = 0.;
        //   }
        //   // Fun2_h sol(Wh, uh);
        //   Paraview<Mesh> writer(Khi, "conservation_"+to_string(ifig++)+".vtk");
        //   writer.add(fun_uh, "uh", 0, 1);
        //   // writer.add(solex, "uex", 0, 1);
        // }
        std::cout << "Iteration " << i + 1 << " / " << niteration << " \t time = " << tid << std::endl;

        std::cout << setprecision(16) << "q(u) = " << qu << std::endl;
        std::cout << setprecision(16) << "|q(u) - q(u0)| = " << fabs(qu - qu0) << setprecision(6) << std::endl;
        outputData << i << "\t" << tid << "\t" << setprecision(16) << qu << "\t" << setprecision(16) << fabs(qu - qu0)
                   << "\t" << setprecision(5) << std::endl;
    }

    std::cout << " Time computation \t" << MPIcf::Wtime() - cpubegin << std::endl;
    return 0;
}
