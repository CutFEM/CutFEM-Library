/*
   We consider a stationary surface problem.
*/

// Dependencies
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <array>
#include <iostream>
#include <experimental/filesystem>
#ifdef USE_MPI
#include "cfmpi.hpp"
#endif
#include "finiteElement.hpp"

// #include "GenericMapping.hpp"
#include "levelSet.hpp"
#include "baseProblem.hpp"
#include "../time_stuff.hpp"
// #include "projection.hpp"
// #include "../util/redirectOutput.hpp"
#include "../num/matlab.hpp"
#include "../num/gnuplot.hpp"
#include "paraview.hpp"

// #include "../algoim/quadrature_general.hpp"
#include "../problem/AlgoimIntegration.hpp"

using namespace globalVariable;

namespace LaplaceBeltrami {

// RHS for surface variable
double fun_rhs0(const R2 P, int elementComp) {
    R x = P.x, y = P.y;
    return (9 * y * (3 * x * x - y * y)) / (x * x + y * y);
}

double fun_testbasisfunction(const R2 P, int elementComp) {
    R x = P.x, y = P.y;
    return 0.5 * x * y;
}

// Exact solution on surface
double fun_uSurface(const R2 P, int elementComp) { return 3.0 * P.x * P.x * P.y - pow(P.y, 3); }

// Exact solution on surface for time-step t
double fun_uSurfaceT(const R2 P, int elementComp, double t) { return 3.0 * P.x * P.x * P.y - pow(P.y, 3); }

// Level-set function
double fun_levelSet(const R2 P, const int i) { return (P.x) * (P.x) + (P.y) * (P.y) - 1.0 - Epsilon; }

double fun_one(const R2 P, const int i) { return 1.; }

// Velocity Field
R fun_velocity(const R2 P, const int i) { return 0.; }
template <int N> struct Levelset {

    double t;

    // level set function
    template <typename T> T operator()(const algoim::uvector<T, N> &x) const { return x(0) * x(0) + x(1) * x(1) - 1; }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2.0 * x(0), 2.0 * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(const R2 P) {
        R normalize = 1. / sqrt(4. * P.x * P.x + 4. * P.y * P.y);
        return R2(2.0 * P.x * normalize, 2.0 * P.y * normalize);
    }
};

} // namespace LaplaceBeltrami

namespace ConvectionDiffusionConstVel {

// RHS for surface variable
double fun_rhs0(const R2 P, int elementComp) {
    R x = P.x, y = P.y;
    return -(3. * (-3 * pow(x, 4) + 2. * pow(x, 3) * y + 9. * x * x * y * y - 45. * x * x * y - 6. * x * pow(y, 3) +
                   15. * pow(y, 3))) /
           (5 * (x * x + y * y));
}

// Exact solution on surface
double fun_uSurface(const R2 P, int elementComp) { return 3.0 * P.x * P.x * P.y - pow(P.y, 3); }

// Exact solution on surface for time-step t
double fun_uSurfaceT(const R2 P, int elementComp, double t) { return 3.0 * P.x * P.x * P.y - pow(P.y, 3); }

// Level-set function
double fun_levelSet(const R2 P, const int i) { return ((P.x) * (P.x) + (P.y) * (P.y) - 1.0); }

// Velocity Field
R fun_velocity(const R2 P, const int i) {
    if (i == 0)
        return 0.4;
    else
        return 0.6;
}

double fun_one(const R2 P, const int i) { return 1.; }

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename T> T operator()(const algoim::uvector<T, N> &x) const { return x(0) * x(0) + x(1) * x(1) - 1; }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2.0 * x(0), 2.0 * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(const R2 P) {
        R normalize = 1. / sqrt(4. * P.x * P.x + 4. * P.y * P.y);
        return R2(2.0 * P.x * normalize, 2.0 * P.y * normalize);
    }
};

} // namespace ConvectionDiffusionConstVel

namespace ConvectionDiffusion {

// RHS for surface variable
double fun_rhs0(const R2 P, int elementComp) {
    R x = P.x, y = P.y;
    return 3 * (-pow(x, 5) + 2 * pow(x, 3) * y * y + 9 * x * x * y + 3 * x * pow(y, 4) - 3 * pow(y, 3)) /
           (x * x + y * y);
}

// Exact solution on surface
double fun_uSurface(const R2 P, int elementComp) { return 3.0 * P.x * P.x * P.y - pow(P.y, 3); }

// Exact solution on surface for time-step t
double fun_uSurfaceT(const R2 P, int elementComp, double t) { return 3.0 * P.x * P.x * P.y - pow(P.y, 3); }

// Level-set function
double fun_levelSet(const R2 P, const int i) { return ((P.x) * (P.x) + (P.y) * (P.y) - 1.0); }

// Velocity Field
R fun_velocity(const R2 P, const int i) {
    if (i == 0)
        return P.y;
    else
        return -P.x;
}

double fun_one(const R2 P, const int i) { return 1.; }

template <int N> struct Levelset {

    double t;

    // level set function
    template <typename T> T operator()(const algoim::uvector<T, N> &x) const { return x(0) * x(0) + x(1) * x(1) - 1; }

    // gradient of level set function
    template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

        return algoim::uvector<T, N>(2.0 * x(0), 2.0 * x(1));
    }

    // normal = grad(phi)/norm(grad(phi))
    R2 normal(const R2 P) {
        R normalize = 1. / sqrt(4. * P.x * P.x + 4. * P.y * P.y);
        return R2(2.0 * P.x * normalize, 2.0 * P.y * normalize);
    }
};

} // namespace ConvectionDiffusion

using namespace ConvectionDiffusion;

#define algoim     // option: "algoim", "quad", "triangle"
#define fem
#define use_h

// Setup two-dimensional class types
const int d = 2;
#if defined(triangle)
typedef Mesh2 Mesh;
#else
typedef MeshQuad2 Mesh;
#endif
typedef GFESpace<Mesh> FESpace;
typedef CutFESpace<Mesh> CutSpace;
typedef TestFunction<Mesh> FunTest;
typedef FunFEM<Mesh> Fun_h;

#if defined(fem)
int main(int argc, char **argv) {

    verbose = 0; // 2 for more info

    MPIcf cfMPI(argc, argv);

    std::cout << std::setprecision(16);

    // Mesh settings and data objects
    const size_t iterations = 5; // number of mesh refinements   (set to 1 to
                                 // run only once and plot to paraview)
    int nx = 20, ny = 15;        // starting mesh size (only apply if use_n is defined)
    double h = 0.1;              // starting mesh size

    std::array<double, iterations> errors; // array to hold L2 errors vs h
    std::array<double, iterations> gamma_length_h;
    std::array<double, iterations> nxs; // array to hold mesh sizes
    std::array<double, iterations> nys; // array to hold mesh sizes
    std::array<double, iterations> hs;  // array to hold mesh sizes

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

        // Mesh

        double x0 = -1.5 + 0.0003, y0 = -1.5 + 0.0003, lx = 3., ly = 3.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, x0, y0, lx, ly);

        std::cout << "Number of elements in background mesh: " << Th.nbElements() << "\n";

        // Paths to store data
        const std::string path_output_data = "../output_files/surface/algoim/laplace_beltrami/data/";
        const std::string path_figures     = "../output_files/surface/algoim/laplace_beltrami/paraview/";

        // Create directory if not already existent
        if (MPIcf::IamMaster()) {
            std::filesystem::create_directories(path_output_data);
            std::filesystem::create_directories(path_figures);
        }

        // Data file to hold problem data
        std::ofstream output_data(path_output_data + "data.dat", std::ofstream::out);

        hs.at(j)  = h;
        nxs.at(j) = nx;
        nys.at(j) = ny;

        std::cout << "h = " << h << "\n";
        std::cout << "nx = " << nx << "\n";
        std::cout << "ny = " << ny << "\n";

        // CG stabilization parameters
        double tau1 = 0.01;

        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space
        GFESpace<Mesh> Vh(Th, DataFE<Mesh>::P1); // continuous basis functions

        // Velocity field
#if defined(algoim) || defined(quad)
        LagrangeQuad2 FEvelocity(1);
#else
        Lagrange2 FEvelocity(1);
#endif
        FESpace VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // std::cout << vel.v << "\n";
        // return 0;

        // Level-set function
        GFESpace<Mesh> Lh(Th, DataFE<Mesh>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Th, levelSet);
        Levelset<2> phi;

        // Create active meshes
        ActiveMesh<Mesh> ThGamma(Th);
        ThGamma.createSurfaceMesh(interface);

        std::cout << "Number of elements in active mesh: " << ThGamma.get_nb_element() << "\n";
        CutSpace Wh(ThGamma, Vh);

#if defined(algoim)
        AlgoimCutFEM<Mesh, Levelset<2>> surfactant(Wh, phi);
#else
        CutFEM<Mesh> surfactant(Wh);
#endif

        double errL2 = 0., intGamma = 0., intGammaSaye = 0.;

        std::vector<double> gamma_length, gamma_length_saye;

        Fun_h funrhs(Wh, fun_rhs0);
        Fun_h funone(Wh, fun_one);
        // Fun_h funtestbasisfunction(Wh, fun_testbasisfunction);

        Normal n;

        // gnuplot::save(Th);
        // gnuplot::save<Mesh, Levelset<2>>(interface, phi);
        // gnuplot::save<Mesh>(interface);

        // Test and Trial functions
        FunTest u(Wh, 1), v(Wh, 1);

        // Diffusion
        surfactant.addBilinear(+innerProduct(gradS(u), gradS(v)), interface);

        // Convection
        surfactant.addBilinear(+innerProduct((vel.exprList() * gradS(u)), v), interface);

        // Stabilization
        double stab_surf_face = tau1;
        surfactant.addFaceStabilization(+innerProduct(stab_surf_face * jump(grad(u) * n), jump(grad(v) * n)), ThGamma);

        // Add RHS on surface
        surfactant.addLinear(+innerProduct(funrhs.expr(), v), interface);

        // Add Lagrange multiplier
        surfactant.addLagrangeMultiplier(innerProduct(1., v), 0., interface);

        matlab::Export(surfactant.mat_[0], path_output_data + "mat_" + std::to_string(j) + ".dat");

        // Solve linear system
        surfactant.solve("mumps");

        KN_<double> dw(surfactant.rhs_(SubArray(surfactant.get_nb_dof(), 0)));
        Fun_h us(Wh, dw);

#if defined(algoim)
        errL2    = L2_norm_surface(us, fun_uSurface, interface, phi, 0, 1);
        intGamma = integral_algoim<Levelset<2>, Fun_h>(funone, interface, 0, phi);
#else
        errL2    = L2normSurf(us, fun_uSurface, interface, 0, 1);
        intGamma = integral(funone, interface, 0);
#endif
        errors.at(j)         = errL2;
        gamma_length_h.at(j) = fabs(intGamma - 2 * pi);

        if (iterations == 1) {
            Paraview<Mesh> ThB(Th, path_figures + "Th.vtk");
            Paraview<Mesh> writer(ThGamma, path_figures + "surfactant.vtk");

            Fun_h uS_ex(Wh, fun_uSurface);
            ExpressionFunFEM<Mesh> vx(vel, 0, op_id);
            ExpressionFunFEM<Mesh> vy(vel, 1, op_id);
            writer.add(us, "surfactant", 0, 1);
            writer.add(uS_ex, "surfactant_exact", 0, 1);
            writer.add(fabs(us.expr() - uS_ex.expr()), "surfactant_error");
            ThB.add(vx, "velx");
            ThB.add(vy, "vely");
            // writer.add(funtestbasisfunction, "test_basis_function", 0, 1);
            writer.add(levelSet, "levelSet", 0, 1);
            writer.writeActiveMesh(ThGamma, path_figures + "ActiveMesh.vtk");
            writer.writeFaceStab(ThGamma, 0, path_figures + "Edges.vtk");
#if defined(algoim)
            writer.writeAlgoimQuadrature(ThGamma, phi, 2, path_figures + "AlgoimQuadrature.vtk");
#endif
        }

        // Refine mesh
#ifdef use_n
        nx *= 2;
        ny *= 2;
#elif defined(use_h)
        h *= 0.5;
#endif
    }

    std::cout << "\n";
    std::cout << "Errors = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << errors.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]"
              << "\n";

    std::cout << "\n";

    std::cout << "\n";
    std::cout << "Length Gamma = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << gamma_length_h.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]"
              << "\n";

    std::cout << "\n";

    std::cout << "h = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << hs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]"
              << "\n";
    std::cout << "nx = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << nxs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]"
              << "\n";

    return 0;
}

#elif defined(integration)
int main(int argc, char **argv) {

    verbose = 0; // 2 for more info

    MPIcf cfMPI(argc, argv);

    std::cout << std::setprecision(16);

    // Mesh settings and data objects
    const size_t iterations = 5; // number of mesh refinements   (set to 1 to
                                 // run only once and plot to paraview)
    int nx = 20, ny = 15;        // starting mesh size (only apply if use_n is defined)
    double h  = 0.1;             // starting mesh size
    // double h  = 0.00625/sqrt(2);             // starting mesh size
    double dT = 0.0625;

    int total_number_iteration;
    double time_step;
    double t0 = 0.;

    std::array<double, iterations> errors; // array to hold L2 errors vs h
    std::array<double, iterations> gamma_length_h;
    std::array<double, iterations> gamma_length_h_saye;
    std::array<double, iterations> nxs; // array to hold mesh sizes
    std::array<double, iterations> nys; // array to hold mesh sizes
    std::array<double, iterations> hs;  // array to hold mesh sizes
    std::array<double, iterations> dts;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

#if defined(frachon1) or defined(zahedi1)
        const double lx = 4., ly = 4.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -2, -2, lx, ly);
#elif defined(zahedi2) || defined(example1)
        const double lx = 1., ly = 1.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, 0., 0., lx, ly);
#elif defined(shi1) || defined(shi2) || defined(shi3) || defined(deckelnick2toshi1)
        const double lx = 8., ly = 6.;

#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif
        Mesh Th(nx, ny, -3, -3, lx, ly);

        // std::string f = "../mesh/square_seb_"+std::to_string(j+1)+".msh";
        // Mesh Th(f.c_str());
#elif defined(deckelnick) || defined(deckelnick2)
        const double lx = 4.8, ly = 4.8;
        // const double lx = 8., ly = 6.;
#ifdef use_h
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
#elif defined(use_n)
        h = lx / (nx - 1);
#endif

        Mesh Th(nx, ny, -2.4, -2.4, lx, ly);
        // Mesh Th(nx, ny, -3, -3, lx, ly);
#endif

        // Parameters
        double tfinal          = .5; // Final time

#ifdef use_t
        total_number_iteration = int(tfinal / dT);
#else
        int divisionMeshSize = 2;

        // int divisionMeshSize = 2*3*pi;
        // int divisionMeshSize = 18;

        double dT = h / divisionMeshSize;

        total_number_iteration = int(tfinal / dT);
#endif
        dT                     = tfinal / total_number_iteration;
        time_step              = dT;

        hs.at(j)  = h;
        nxs.at(j) = nx;
        nys.at(j) = ny;
        dts.at(j) = dT;

        if (iterations > 1) {
            std::cout << "--------------------------------------------"
                      << "\n";
            std::cout << "--------------------------------------------"
                      << "\n";
            std::cout << "--------------------------------------------"
                      << "\n";
            std::cout << "--------------------------------------------"
                      << "\n";
            std::cout << "Iteration " << j + 1 << "/" << iterations << "\n";
        }
        std::cout << "h = " << h << "\n";
        std::cout << "nx = " << nx << "\n";
        std::cout << "ny = " << ny << "\n";
        std::cout << "dT = " << dT << "\n";

        double D = 1.;

        // CG stabilization parameters
        double tau1 = 1., tau2 = .1;

        // Background FE Space, Time FE Space & Space-Time Space
        // 2D Domain space
        GFESpace<Mesh> Vh(Th, DataFE<Mesh>::P1); // continuous basis functions

        // 1D Time mesh
        double final_time = total_number_iteration * time_step;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);
        // 1D Time space
        FESpace1 Ih(Qh, DataFE<Mesh1>::P1Poly);
        // Quadrature data
        const QuadratureFormular1d &qTime(*Lobatto(3));

        const Uint nbTime       = qTime.n;
        const Uint ndfTime      = Ih[0].NbDoF();
        const Uint lastQuadTime = nbTime - 1;

        // Velocity field
#if defined(shi1) || defined(shi2) || defined(shi3) || defined(deckelnick2)
        Lagrange2 FEvelocity(0);
#elif defined(example1) || defined(deckelnick)
        Lagrange2 FEvelocity(1);
#else
        Lagrange2 FEvelocity(2);
#endif
        GFESpace<Mesh> VelVh(Th, FEvelocity);
        Fun_h vel(VelVh, fun_velocity);

        // Set up level-set function
        GFESpace<Mesh> Lh(Th, DataFE<Mesh>::P1);
        double dt_levelSet = dT / (nbTime - 1);
        std::vector<Fun_h> ls(nbTime);
#if defined(levelsetexact)
        for (int i = 0; i < nbTime; i++)
            ls[i].init(Lh, fun_levelSet, 0.);
#elif defined(levelsetsolve)
        for (int i = 0; i < nbTime; i++)
            ls[i].init(Lh, fun_levelSet);
#endif

        // Declare time dependent interface
        TimeInterface<Mesh> interface(qTime);

        Levelset<2> phi;

        AlgoimCutFEM<Mesh, Levelset<2>> surfactant(qTime, phi);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << "\n";

        int iter = 0;

        double intGamma = 0., intGammaSaye = 0.; // hold integrals of rhs and Neumann bcs

        std::vector<double> gamma_length, gamma_length_saye;

        // Iterate over time-slabs
        while (iter < total_number_iteration) {

            int current_iteration = iter;
            double tid            = iter * time_step;

            const TimeSlab &In(Ih[iter]);

            std::cout << " ----------------------------------------------------"
                         "--------- "
                      << "\n";
            std::cout << " ----------------------------------------------------"
                         "--------- "
                      << "\n";
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << "\n";
            std::cout << " Time      \t : \t" << current_iteration * time_step << "\n";

            swap(ls[0], ls[lastQuadTime]);
            // computation of the interface in each of the three quadrature
            // points
            for (int i = 0; i < nbTime; ++i) {

#if defined(levelsetexact) && not defined(example2)
                R tt = In.Pt(R1(qTime(i).x));
                ls[i].init(Lh, fun_levelSet, tt);
#endif
                interface.init(i, Th, ls[i]);

#if defined(levelsetsolve) || defined(example2)
                // We solve for the level-set using Crank-Nicholson in time
                if (i < lastQuadTime) {
                    LevelSet::move_2D(ls[i], vel, vel, dt_levelSet, ls[i + 1]);
                }
#endif
            }

            // Create active meshes
            ActiveMesh<Mesh> ThGamma(Th);
            ThGamma.createSurfaceMesh(interface);

            Fun_h funone(Vh, In, fun_one);
            Fun_h funx(Vh, In, fun_x);
            auto fundx   = dx(funx.expr());
            // intGamma = integral(funone, interface(1), 0);
            intGamma     = integral(funone, In, interface, 0);
            // Levelset<2> phi;
            intGammaSaye = integral_saye(funone, In, interface, 0, phi);

            gamma_length.push_back(intGamma / dT);

            gamma_length_saye.push_back(intGammaSaye / dT);

            // 6.27844864122449

            // Levelset<2> phi2;
            // phi2.t = tid;

            // algoim::QuadratureRule<2> q2 = algoim::quadGen<2>(phi2, algoim::HyperRectangle<double, 2>(-2.4, 2.4), 2,
            // -1, 4);

            // std::cout << "[";
            // for (auto &node : q2.nodes) {
            //     std::cout << "[" << node.x(0) << ", " << node.x(1) << "]," << "\n";
            // }
            // std::cout << "]" << "\n";

            // getchar();

            iter++;
        }

        gamma_length_h.at(j)      = intGamma;
        gamma_length_h_saye.at(j) = intGammaSaye / dT;

        std::cout << "\n";
        std::cout << "Length of Gamma(t) vs t = [";
        for (auto &length : gamma_length) {

            std::cout << length;

            std::cout << ", ";
        }
        std::cout << "]"
                  << "\n";

        std::cout << "\n";
        std::cout << "Length of Gamma(t) Saye vs t = [";
        for (auto &length : gamma_length_saye) {

            std::cout << length;

            std::cout << ", ";
        }
        std::cout << "]"
                  << "\n";

        // Refine mesh

#ifdef use_n
        nx *= 2;
        ny *= 2;
#elif defined(use_t)
        dT *= 0.5;
#elif defined(use_h)
        // h *= sqrt(0.5);     //! CHANGE BACK
        h *= 0.5;
#endif
    }

    std::cout << "\n";
    std::cout << "Length of Gamma(t) = [";
    for (auto &length : gamma_length_h) {

        std::cout << length;

        std::cout << ", ";
    }
    std::cout << "]"
              << "\n";

    std::cout << "\n";
    std::cout << "Length of Gamma(t) Saye = [";
    for (auto &length : gamma_length_h_saye) {

        std::cout << length;

        std::cout << ", ";
    }
    std::cout << "]"
              << "\n";

    return 0;
}

#endif