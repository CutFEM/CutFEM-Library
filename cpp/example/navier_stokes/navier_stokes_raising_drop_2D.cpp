/**
 * @file navier_stokes_raising_drop_2D.cpp
 * @author Thomas Frachon
 * @brief Raising bubble using CutFEM. For details, one can see Example 2 in
 * "A cut finite element method for incompressible two-phase Navier-Stokes flows" by  Thomas Frachon and Sara Zahedi
 * @version 0.1
 * @date 2023-03-09
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "../tool.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;

/// @brief Levelset function
/// @param P
/// @return
double fun_levelSet(R2 P) { return sqrt((P.x - 0.5) * (P.x - 0.5) + (P.y - 0.5) * (P.y - 0.5)) - 0.25; }

double fun_rhs(R2 P, int i) { return (i < 1) ? 0. : -0.98; }

double fun_boundary(R2 P, int i) { return 0.; }

double fun_x(R2 P, int i) { return P[i]; }

double fun_1(R2 P, int i) { return 1; }

/// @brief Define penalty parameter on the boundary
class LambdaBoundary : public VirtualParameter {
  public:
    const CutFEMParameter &mu_;
    double G_, H_;
    LambdaBoundary(const CutFEMParameter &mu, double G, double H) : mu_(mu), G_(G), H_(H) {}
    double evaluate(int domain, double h, double meas, double measK, double measCut) const {
        double gamma = meas / h;
        double alpha = measK / h / h;
        double val   = mu_.evaluate(domain, h, meas, measK, measCut);
        return val / h * (G_ + H_ * gamma / alpha);
    }
};

/// @brief Define penalty paramter on interface
class LambdaGamma : public VirtualParameter {
  public:
    const double mu1_;
    const double mu2_;
    LambdaGamma(const double m1, const double m2) : mu1_(m1), mu2_(m2) {}
    double evaluate(int domain, double h, double meas, double measK, double measCut) const {
        double gamma  = meas / h;
        double alphaK = measK / h / h;
        return (0.5 * mu1_ + 0.5 * mu2_) * (100 + 10 * gamma) / (alphaK);
    }
};

/// @brief Define weights for {}
class WeightKappa : public VirtualParameter {
  public:
    const double mu1_;
    const double mu2_;
    WeightKappa(const double m1, const double m2) : mu1_(m1), mu2_(m2) {}
    double mu(int i) const { return (i == 0) ? mu1_ : mu2_; }
    double evaluate(int dom0, double h, double meas, double measK, double measCut0) const {
        int dom1        = (dom0 == 0);
        double measCut1 = meas - measCut0;
        double alphaK0  = 0.5;
        // measCut0 / h / h;
        double alphaK1  = 0.5;
        // measCut1 / h / h;
        return (mu(dom1) * alphaK0) / (mu(dom0) * alphaK1 + mu(dom1) * alphaK0);
    }
};

/// @brief Define  weights for <> operator
class Weight2Kappa : public VirtualParameter {
  public:
    const double mu1_;
    const double mu2_;
    Weight2Kappa(const double m1, const double m2) : mu1_(m1), mu2_(m2) {}
    double mu(int i) const { return (i == 0) ? mu1_ : mu2_; }
    double evaluate(int dom0, double h, double meas, double measK, double measCut0) const {
        int dom1        = (dom0 == 0);
        double measCut1 = meas - measCut0;
        double alphaK0  = 0.5;
        // measCut0 / h / h;
        double alphaK1  = 0.5;
        // measCut1 / h / h;
        return (mu(dom0) * alphaK1) / (mu(dom0) * alphaK1 + mu(dom1) * alphaK0);
    }
};

int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);
    // int thread_count_tmp = 1;
    // cout << "Threads: ";
    // cin >> thread_count_tmp;
    // MPIcf::Bcast(thread_count_tmp, MPIcf::Master(), 1);
    // const int thread_count = thread_count_tmp;
    // omp_set_num_threads(thread_count);
    int thread_count = 1;
    double cpubegin  = MPIcf::Wtime();

    CutFEMLogger::initialize("log_navier_Stokes.txt");

    // MESH DEFINITION
    // ---------------------------------------------
    int nx = 40;
    int ny = 80;
    Mesh2 Kh(nx, ny, 0., 0., 1., 2.);
    double mesh_size = (1. / (nx - 1));

    std::list<int> dirichlet{1, 3};
    std::list<int> neumann{2, 4};

    LOG_INFO << " ------------------------------------" << logger::endl;
    LOG_INFO << " Background mesh " << logger::endl;
    Kh.info();

    // TIME DEFINITION
    // ---------------------------------------------
    int division_mesh_size     = 2;
    double final_time          = 1;
    double dT                  = mesh_size / division_mesh_size;
    int total_number_iteration = final_time / dT;
    dT                         = final_time / total_number_iteration;
    double time_step           = dT;
    double t0                  = 0.;

    Mesh1 Qh(total_number_iteration + 1, t0, final_time);
    FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);

    // Define the time quadrature formula
    const QuadratureFormular1d &qTime(*Lobatto(3));
    const Uint nb_quad_time   = qTime.n;
    const Uint ndf_time_slab  = Ih[0].NbDoF();
    const Uint last_quad_time = nb_quad_time - 1;

    // PROBLEM AND PARAMETER DEFINITION
    // ----------------------------------------------
    ProblemOption optionProblem;
    optionProblem.solver_name_  = "mumps";
    optionProblem.clear_matrix_ = true;
    std::vector<std::map<std::pair<int, int>, double>> mat_NL(thread_count);

    CutFEM<mesh_t> stokes(qTime, thread_count, optionProblem);

    CutFEMParameter mu(10., 1.);
    CutFEMParameter rho(1000., 100.);
    CutFEMParameter invmu(0.1, 1.);
    LambdaBoundary lambdaB(mu, 10, 100);
    LambdaGamma lambdaG(mu(0), mu(1));
    WeightKappa kappa(mu(0), mu(1));
    Weight2Kappa kappa2(mu(0), mu(1));
    const double sigma = 24.5;

    // WRITE THE PROBLEM PARAMETERS
    // ----------------------------------------------
    LOG_INFO << " ------------------------------------" << logger::endl;
    LOG_INFO << " Discretization of the problem : \n"
             << " h = " << mesh_size << "\n t_begin = " << t0 << " and t_final = " << final_time << "\n dt = " << dT
             << "\n number of iteration = " << total_number_iteration << logger::endl;
    LOG_INFO << " ------------------------------------" << logger::endl;
    LOG_INFO << " Parameters of the problem : \n"
             << " mu_1 = " << mu(0) << " and mu_2 = " << mu(1) << " \n"
             << " rho_1 = " << rho(0) << " and rho_2 = " << rho(1) << " \n"
             << " the surface tension sigma = " << sigma << logger::endl;
    LOG_INFO << " ------------------------------------ " << logger::endl;

    // SPACE DEFINITION
    // ---------------------------------------------
    Lagrange2 FEu(2);
    space_t Vh(Kh, FEu);
    space_t Lh(Kh, DataFE<mesh_t>::P1);
    Lagrange2 FEcurv(1);
    space_t Vh1(Kh, FEcurv);

    const int frequency_plotting = 1;

    // INTERPOLATION OF THE VELOCITY (BOUNDARY CONDITION)
    // ----------------------------------------------
    LOG_INFO << " ------------------------------------" << logger::endl;
    LOG_INFO << " Interpolate the velocity " << logger::endl;
    std::vector<fct_t> vel(nb_quad_time);
    for (int i = 0; i < nb_quad_time; ++i)
        vel[i].init(Vh, fun_boundary);

    // DECLARATION OF INTERFACE AND LEVELSET
    // ----------------------------------------------
    LOG_INFO << " ------------------------------------" << logger::endl;
    LOG_INFO << " Create interface and levelSet " << logger::endl;
    TimeInterface<mesh_t> interface(qTime);
    double dt_levelSet = dT / (nb_quad_time - 1);
    std::vector<fct_t> ls_k(nb_quad_time), ls(nb_quad_time);

    for (int i = 0; i < nb_quad_time; ++i)
        ls[i].init(Lh, fun_levelSet);

    //   projection(ls_k[0], ls[nbTime-1]);
    //   levelSet.setStrongBC({2,4});

    //   CReinitialization<Mesh> reinitialization;
    //   reinitialization.number_iteration = 4;
    //   reinitialization.epsilon_diffusion = 1e-3;
    //   reinitialization.dt = dT/8;
    //   reinitialization.ON_OFF = "ON";
    //   reinitialization.ODE_method = "Euler";
    //   reinitialization.mass_correction = "ON";
    //   reinitialization.precision_correction = 1e-6;
    //   reinitialization.max_iteration_correction = 10;
    //   reinitialization.info();
    //   const int frequencyReinitialization = 5;

    //   list<int> dirichlet = {1,2,3,4,20,40};

    // INTERPOLATE RHS AND BOUNDARY CONDITION
    // ----------------------------------------------
    LOG_INFO << " ------------------------------------" << logger::endl;
    LOG_INFO << " Interpolate boundary and rhs " << logger::endl;
    fct_t fh(Vh, fun_rhs);
    fct_t gh(Vh, fun_boundary);

    LOG_INFO << " ------------------------------------" << logger::endl;
    LOG_INFO << " ------------------------------------" << logger::endl;

    int iter = 0, ifig = 0;
    progress bar(" Navier-Stokes solver ", total_number_iteration, 2);
    while (iter < total_number_iteration) {
        int current_iteration = iter;
        const TimeSlab &In(Ih[iter]);

        // MOVE THE LEVELSET AND CREATE INTERFACES
        // ----------------------------------------------
        {
            LOG_INFO << "\n ------------------------------------" << logger::endl;
            LOG_INFO << " Move the levelSet and create interfaces " << logger::endl;

            double t0               = MPIcf::Wtime();
            globalVariable::verbose = 0;

            // ls.begin()->swap(ls[last_quad_time]);
            swap(ls[0], ls[last_quad_time]);

            for (int i = 0; i < nb_quad_time; ++i) {

                // here we want to reinitializa the levelSet.
                // !!! WE CANNOT DO ON THE FIRST LEVELSET BECAUSE IT HAS TO MATCH
                // THE LAST STEP if(iter%frequencyReinitialization == 0 && i == 1 &&
                // iter > 0) {
                //   LOG_INFO << " ------- reinitialization LevelSet -------" <<
                //   logger::endl; reinitialization.perform(ls_k[i], ls[i]);
                // }

                interface.init(i, Kh, ls[i]);
                if (i < nb_quad_time - 1) {
                    LevelSet::move(ls[i], vel[last_quad_time], vel[last_quad_time], dt_levelSet, ls[i + 1]);
                }
            }
            LOG_INFO << " Time moving interface : \t" << MPIcf::Wtime() - t0 << logger::endl;
        }

        // CREATE ACTIVE MESH AND CUTSPACES
        // ----------------------------------------------
        double t0_cmesh = MPIcf::Wtime();
        cutmesh_t Kh_i(Kh, interface);
        cutspace_t Wh(Kh_i, Vh);
        cutspace_t Ph(Kh_i, Lh);

        if (iter == 0) {
            LOG_INFO << " Cut Mesh " << logger::endl;
            Kh_i.info();
            LOG_INFO << " Cut Space for velocity " << logger::endl;
            Wh.info();
            LOG_INFO << " Cut Space for pressure " << logger::endl;
            Ph.info();
            LOG_INFO << " TIME BUILDING CUT MESH/SPACE : \t" << MPIcf::Wtime() - t0_cmesh << logger::endl;
        }

        // DEFINE TEST FUNCTIONS
        // ----------------------------------------------
        Normal n;
        funtest_t du(Wh, 2), dp(Ph, 1), v(Wh, 2), q(Ph, 1), dp1(Wh, 1, 0, 0);
        funtest_t Eun  = (Eps(du) * n);
        funtest_t D2nu = grad(grad(du) * n) * n, D2nv = grad(grad(v) * n) * n;

        // DEFINE THE PROBLEM ON THE CORRECT SPACES
        // ----------------------------------------------
        stokes.initSpace(Wh, In);
        stokes.add(Ph, In);

        LOG_INFO << " Problem's DOF : \t" << stokes.get_nb_dof() << logger::endl;

        // INITIALIZE VECTORS USED
        // -------------------------------------
        double t0_init_sol = MPIcf::Wtime();

        std::vector<double> data_init(stokes.get_nb_dof());
        std::span<double> data_init_span(data_init);
        stokes.initialSolution(data_init_span);
        std::vector<double> data_all(data_init);

        int idxp0 = Wh.NbDoF() * In.NbDoF();

        std::span<double> data_uh0 = std::span<double>(data_init.data(), Wh.NbDoF());
        std::span<double> data_uh  = std::span<double>(data_all.data(), Wh.NbDoF() * In.NbDoF());
        std::span<double> data_ph  = std::span<double>(data_all.data() + idxp0, Ph.NbDoF() * In.NbDoF());

        fct_t u0(Wh, data_uh0);
        fct_t uh(Wh, In, data_uh);
        fct_t ph(Ph, In, data_ph);

        if (iter == 0) {
            interpolate(Wh, data_uh0, fun_boundary);
        }
        LOG_INFO << " Time data initialization : \t" << MPIcf::Wtime() - t0_init_sol << logger::endl;

        // NEWTON ITERATION
        // ----------------------------------------------
        int iterNewton   = 0;
        double t0_newton = MPIcf::Wtime();
        LOG_INFO << " NEWTON ITERATION " << logger::endl;
        while (1) {
            double tt0 = MPIcf::Wtime();
            if (iterNewton == 0) {
                stokes.addBilinear(innerProduct(dt(du), rho * v) + contractProduct(2 * mu * Eps(du), Eps(v)) -
                                       innerProduct(dp, div(v)) + innerProduct(div(du), q),
                                   Kh_i, In);
                stokes.addBilinear(innerProduct(jump(du), -2 * mu * average(Eps(v) * n, kappa)) +
                                       innerProduct(-2 * mu * average(Eps(du) * n, kappa), jump(v)) +
                                       innerProduct(lambdaG * jump(du), jump(v)) +
                                       innerProduct(average(dp, kappa), jump(v * n)) -
                                       innerProduct(jump(du * n), average(q, kappa)),
                                   interface, In);

                stokes.addBilinear(innerProduct(lambdaB * du, v) + innerProduct(dp, v * n) - innerProduct(du * n, q) -
                                       innerProduct(2. * mu * Eps(du) * n, v) - innerProduct(du, 2. * mu * Eps(v) * n),
                                   Kh_i, INTEGRAL_BOUNDARY, In, dirichlet);

                stokes.addBilinear(innerProduct(lambdaB * du * n, v * n) + innerProduct(dp, v * n) -
                                       innerProduct(du * n, q) - innerProduct(2. * mu * Eps(du) * n * n, v * n) -
                                       innerProduct(du * n, 2. * mu * Eps(v) * n * n),
                                   Kh_i, INTEGRAL_BOUNDARY, In, neumann);

                double h1 = mesh_size;
                double h3 = pow(mesh_size, 3);
                stokes.addFaceStabilization(1e-1 * h1 * innerProduct(rho * jump(grad(du) * n), mu * jump(grad(v) * n)) +
                                                1e-1 * h3 * innerProduct(rho * jump(D2nu), mu * jump(D2nv)) +
                                                1e-1 * h3 * innerProduct(jump(grad(dp) * n), invmu * jump(grad(q) * n)),
                                            Kh_i, In);

                // impose initial condition
                stokes.addBilinear(innerProduct(rho * du, v), Kh_i);
            }

            for (int i = 0; i < nb_quad_time; ++i) { // computation of the curvature

                cutmesh_t cutTh(Kh);
                cutTh.createSurfaceMesh(*interface[i]);
                cutspace_t cutVh(cutTh, Vh1);
                Curvature<mesh_t> curvature(cutVh, interface[i]);

                auto tq    = qTime(i);
                double tid = (double)In.map(tq);

                auto data_H = curvature.solve();
                fct_t H(cutVh, data_H);
                stokes.addLinear(-sigma * innerProduct(H.exprList(), average(v, kappa2)), interface, In, i);

                // if (i == 0) {
                //     Paraview<mesh_t> writerS(cutTh, "raingDropExampleCurvature_" + std::to_string(ifig) + ".vtk");
                //     writerS.add(ls[i], "levelSet", 0, 1);
                //     writerS.add(H, "meanCurvature", 0, 2);
                // }
            }

            // if (iter == 0 && iterNewton == 0)
            //     globalVariable::verbose = 1;

            // if (iter == 0 && iterNewton == 0) {
            LOG_INFO << " Time assembly matrix A : \t" << MPIcf::Wtime() - tt0 << logger::endl;
            // }

            //
            // CONSTRUCTION LINEAR PART
            //
            tt0 = MPIcf::Wtime();
            stokes.addLinear(-innerProduct(fh.exprList(), rho * v), Kh_i, In);
            stokes.addLinear(-innerProduct(u0.exprList(), rho * v), Kh_i);
            stokes.addLinear(-innerProduct(gh.exprList(), lambdaB * v) +
                                 innerProduct(gh.exprList(), 2. * mu * Eps(v) * n) + innerProduct(gh.exprList(), q * n),
                             Kh_i, INTEGRAL_BOUNDARY, In);

            LOG_INFO << " Time assembly rhs L : \t" << MPIcf::Wtime() - tt0 << logger::endl;

            tt0 = MPIcf::Wtime();
            stokes.gather_map();
            stokes.addMatMul(data_all);
            LOG_INFO << " Time matrix multiplication : \t" << MPIcf::Wtime() - tt0 << logger::endl;

            //
            // CONSTRUCTION NON LINEAR PART
            //
            tt0 = MPIcf::Wtime();
            stokes.set_map(mat_NL);
            mat_NL[0] = stokes.mat_[0];

            funtest_t du1(Wh, 1, 0), du2(Wh, 1, 1), v1(Wh, 1, 0), v2(Wh, 1, 1);
            auto ux = uh.expr(0);
            auto uy = uh.expr(1);

            stokes.addBilinear(innerProduct(du1 * dx(ux) + du2 * dy(ux), rho * v1) +
                                   innerProduct(du1 * dx(uy) + du2 * dy(uy), rho * v2) +
                                   innerProduct(ux * dx(du1) + uy * dy(du1), rho * v1) +
                                   innerProduct(ux * dx(du2) + uy * dy(du2), rho * v2),
                               Kh_i, In);
            stokes.addLinear(innerProduct(ux * dx(ux) + uy * dy(ux), rho * v1) +
                                 innerProduct(ux * dx(uy) + uy * dy(uy), rho * v2),
                             Kh_i, In);

            stokes.addLagrangeMultiplier(innerProduct(1., dp1), 0., Kh_i, In);
            // stokes.addLagrangeMultiplier(
            //   innerProduct(1.,dp1), 0.
            //   , Kh_i
            //   // , 0
            //   // , In
            // );

            LOG_INFO << " Time assembly non linear matrix \t" << MPIcf::Wtime() - tt0 << logger::endl;

            // gather(mat_NL);
            // matlab::Export(mat_NL[0], "mat_1.dat");
            // return 0;
            stokes.solve(mat_NL[0], stokes.rhs_);

            LOG_INFO << " Time solver \t" << MPIcf::Wtime() - tt0 << logger::endl;

            std::span<double> dwu = std::span<double>(stokes.rhs_.data(), Wh.NbDoF() * In.NbDoF());
            const auto result =
                std::max_element(dwu.begin(), dwu.end(), [](int a, int b) { return std::abs(a) < std::abs(b); });
            double dist = *result;
            LOG_INFO << " Error || du ||_infty = " << dist << logger::endl;

            std::span<double> dw = std::span<double>(stokes.rhs_.data(), stokes.get_nb_dof());
            std::transform(data_all.begin(), data_all.end(), dw.begin(), data_all.begin(),
                           [](double a, double b) { return a - b; });

            stokes.rhs_.resize(stokes.get_nb_dof());
            // stokes.rhs_ = 0.0;
            std::fill(stokes.rhs_.begin(), stokes.rhs_.end(), 0.0);

            iterNewton += 1;

            if (iterNewton == 2 || dist < 1e-10) {
                std::span<double> data_all_span(data_all);
                stokes.saveSolution(data_all_span);
                stokes.cleanBuildInMatrix();
                stokes.set_map();

                if (iter == 0) {
                    LOG_INFO << " TIME NEWTON ITERATION : \t" << MPIcf::Wtime() - t0_newton << logger::endl;
                }

                break;
            }
        }

        // COMPUTATION OF THE CARACTERISTICS OF THE DROPS
        // {
        //   Fun_h funX(Wh, fun_x);
        //   Fun_h fun1(Wh, fun_1);
        //   R tt0 = CPUtime();
        //   double areaBubble   = integral(fun1, 0, 1) ;
        //   double centerOfMass = integral(funX, 1, 1) / areaBubble ;
        //
        //   double Pb = integralSurf(fun1, 1);
        //   double ra = sqrt(areaBubble / M_PI);
        //   double Pa = 2*M_PI*ra;
        //   double circularity = Pa / Pb;
        //   double riseVelocity = integral(uh, 1, 1) / areaBubble;
        //
        //   double q    = integralSurf(us,  0, 0 , In.map(qTime[0]));
        //   double qend = integralSurf(us,  0, lastQuadTime,
        //   In.map(qTime[lastQuadTime])); double q0   = integralSurf(u0s, 0);
        //
        //
        //   if(iter == 0) initialConcentrationSurfactant = qend;
        //   LOG_INFO << "\n Features of the drop " << logger::endl;
        //   LOG_INFO << " Time                   ->    " <<
        //   GTime::current_time()+dT  << logger::endl; LOG_INFO << " Center Of Mass
        //   ->    " << centerOfMass << logger::endl; LOG_INFO << " Circularity ->
        //   " << circularity << logger::endl; LOG_INFO << " Rise velocity ->    "
        //   << riseVelocity << logger::endl; LOG_INFO << " Surfactant quantity
        //   init->    " << q0 << logger::endl; LOG_INFO << " Surfactant quantity ->
        //   " << q << logger::endl; LOG_INFO << " Surfactant quantity end ->   " <<
        //   qend << logger::endl; LOG_INFO << " |q_0 - q_end|           ->   " <<
        //   fabs(q - qend) << logger::endl; LOG_INFO << " Surfactant conservation
        //   ->   " << fabs(qend - initialConcentrationSurfactant) << logger::endl;
        //
        //
        //   outputData << GTime::current_time()+dT << "\t"
        //   << centerOfMass << "\t"
        //   << circularity << "\t"
        //   << riseVelocity << "\t"
        //   << areaBubble <<  "\t"
        //   << qend << "\t"
        //   << fabs(q - qend) << "\t"
        //   << fabs(qend - initialConcentrationSurfactant) << logger::endl;
        //
        // }
        //

        LOG_INFO << " Set Velocity " << logger::endl;
        for (int i = 0; i < nb_quad_time; ++i) {
            fct_t sol(Wh, In, data_uh);
            set_velocity(sol, vel[i], ls[i], In.map(qTime[i]));
        }

        //  -----------------------------------------------------
        //                     PLOTTING
        //  -----------------------------------------------------

        if (MPIcf::IamMaster() && (iter % frequency_plotting == 0 || iter + 1 == total_number_iteration)) {

            {
                std::string filename = "raisingDropExample_" + std::to_string(ifig) + ".vtk";
                LOG_INFO << " Plotting -> " << filename << logger::endl;
                Paraview<mesh_t> writer(Kh_i, filename);
                writer.add(ls[0], "levelSet", 0, 1);
                writer.add(uh, "velocity", 0, 2);
                writer.add(ph, "pressure", 0, 1);
            }
            ifig++;
        }

        bar++;
        iter += 1;
        globalVariable::verbose = 1;
    }
    bar.end();

    LOG_INFO << "\n\n ------------------------------------" << logger::endl;
    LOG_INFO << " -----------------------------------  " << logger::endl;
    LOG_INFO << " Computation time : \t " << MPIcf::Wtime() - cpubegin << logger::endl;
}
