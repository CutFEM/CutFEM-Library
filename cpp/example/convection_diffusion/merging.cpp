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

/**
 * @brief Space-time CutFEM for simulating the convection-diffusion equation on colliding circles/spheres.
 * @note 

 *  Problem:
    Find u in Omega(t) such that

    dt(u) + div(beta*u - D*grad(u)) = f,    in Omega(t).
                        n*D*grad(u) = 0,    on Gamma(t).

 *  Two schemes are implemented:
 *  Non-conservative scheme,
 *  Conservative scheme.

 *  Two types of stabilization is possible: full stabilization or macroelement stabilization.
 *  Read more about the methods in our paper: https://www.sciencedirect.com/science/article/pii/S0045782524005012
*/

// Dependencies
#include "../cutfem.hpp"
#include "../problem/AlgoimIntegration.hpp"
#include "../num/matlab.hpp"
#include <string>


namespace Merging2D {

const double R0       = .5;
const double D        = .1;
const double T        = 1.5;

double fun_one(double *P, const int i) { return 1.; }


double fun_levelSet(double *P, const int i, const double t) {

    const double dist1 = P[0]*P[0] + (P[1] - (t-0.75))*(P[1] - (t-0.75));
    const double dist2 = P[0]*P[0] + (P[1] - (0.75-t))*(P[1] - (0.75-t));

    return std::min(dist1, dist2) - R0 * R0;
}

R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;

}

// Velocity field
R fun_velocity(double *P, const int i, const double t) {
    
    if (((P[1] > 0) && (t <= T/2)) || ((P[1] <= 0) && (t > T/2))) 
        return (i==0) ? 0 : -1;
    
    else if (((P[1] <= 0) && (t <= T/2)) || ((P[1] > 0) && (t > T/2))) 
        return (i==0) ? 0 : 1;
    else {
        std::cerr << "ERROR\n";
        return 0;
    }
        
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i, const int domain) { 
    if (P[1] > 0)
        return 1;
    else
        return -1;
}

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    return 0.;
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    return 0.;
}
} // namespace Merging2D


namespace Merging3D {

const double R0       = .5;
const double D        = .1;
const double T        = 1.5;

double fun_one(double *P, const int i) { return 1.; }


double fun_levelSet(double *P, const int i, const double t) {

    const double dist1 = P[0]*P[0] + P[1]*P[1] + (P[2] - (t-0.75))*(P[2] - (t-0.75));
    const double dist2 = P[0]*P[0] + P[1]*P[1] + (P[2] - (0.75-t))*(P[2] - (0.75-t));

    return std::min(dist1, dist2) - R0 * R0;
}

R fun_rhsBulk(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;
}

// The rhs Neumann boundary condition
R fun_neumann_Gamma(double *P, const int i, const R t) {
    R x = P[0], y = P[1];

    return 0.;

}

// Velocity field
R fun_velocity(double *P, const int i, const double t) {
    
    if (((P[2] > 0) && (t <= T/2)) || ((P[2] <= 0) && (t > T/2))) 
        return (i==2) ? -1 : 0;
    
    else if (((P[2] <= 0) && (t <= T/2)) || ((P[2] > 0) && (t > T/2))) 
        return (i==2) ? 1 : 0;
    else {
        std::cerr << "ERROR\n";
        return 0;
    }
        
}

// Initial solution bulk
R fun_uBulkInit(double *P, const int i, const int domain) { 
    if (P[2] > 0)
        return 1;
    else
        return -1;
}

// Exact solution bulk
R fun_uBulk(double *P, const int i, const R t) {
    return 0.;
}

R fun_uBulkD(double *P, const int i, const int d, const R t) {
    return 0.;
}
} // namespace Merging3D


// Define example, method and stabilization
#define merging2D         // example (merging2D/merging3D)
#define conservative  // method (conservative/non_conservative)
#define fullstab          // stabilization (fullstab/macro)
#define K 1               // polynomial order in time (1)
#define M 1               // polynomial order in space (1)

#if defined(merging2D)
using mesh_t        = Mesh2;
#elif defined(merging3D)
using mesh_t        = Mesh3;
#endif
using fun_test_t    = TestFunction<mesh_t>;
using fct_t         = FunFEM<mesh_t>;
using activemesh_t  = ActiveMesh<mesh_t>;
using fespace_t     = GFESpace<mesh_t>;
using cut_fespace_t = CutFESpace<mesh_t>;

std::vector<const GTypeOfFE<mesh_t> *> FE_space = {&DataFE<mesh_t>::P1, &DataFE<mesh_t>::P2};
std::vector<const GTypeOfFE<Mesh1> *> FE_time   = {&DataFE<Mesh1>::P0Poly, &DataFE<Mesh1>::P1Poly, &DataFE<Mesh1>::P2Poly,
                                                   &DataFE<Mesh1>::P3Poly};



int main(int argc, char **argv) {

    std::string example, method, stabilization;

// Do not touch
#if defined(merging2D)
    using namespace Merging2D;
#elif defined(merging3D)
    using namespace Merging3D;
#else
#error "No example defined"
#endif

#if defined(non_conservative)
    method = "non_conservative";
#elif defined(conservative)
    method              = "conservative";
#else
#error "No method defined"
#endif

#if defined(fullstab)
    stabilization = "fullstab";
#elif defined(macro)
    #if defined(merging3D)
    std::cout << "Macro stabilization not implemented for 3D" << '\n';
    exit(1);
    #endif
    stabilization       = "macro";
#else
#error "No stabilization defined"
#endif

    MPIcf cfMPI(argc, argv);

    const int k = K;
    const int m = M;

    assert(k >= 0 && k <= 3);
    assert(m >= 1 && m <= 3);

    const size_t iterations = 1;   // number of mesh refinements
    double h                = 0.05625; // starting mesh size
    int nx, ny, nz;

    const double cfl_number = 1./3;
    int total_number_iteration;
    double dT       = 0.1;
    
    const double t0 = 0., tfinal = T;

    // Time integration quadrature
    const size_t quadrature_order_time = 3;
    const QuadratureFormular1d &qTime(*Lobatto(quadrature_order_time)); // specify order of quadrature in time
    const Uint nbTime       = qTime.n;
    const Uint lastQuadTime = nbTime - 1;

    // Space integration quadrature
    const std::string solver("umfpack");

    CutFEM<mesh_t> convdiff(qTime);

    // Method parameters
    const double tau   = 0.1;  // stabilization constant
    const double delta = 0.3; // macro parameter

    assert(tau > 0.);
    assert(delta > 0. && delta <= 1.);
    
    // Arrays to hold data
    std::array<double, iterations> hs, dts, global_conservation_errors;

    // Iterate over mesh sizes
    for (int j = 0; j < iterations; ++j) {

        #if defined(merging2D)
        const double x0 = -0.6 - globalVariable::Epsilon, y0 = -1.35 - globalVariable::Epsilon;
        const double lx = 1.2, ly = 2.7;
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1;
        mesh_t Th(nx, ny, x0, y0, lx, ly);
        #elif defined(merging3D)
        const double x0 = -0.6 - globalVariable::Epsilon;
        const double y0 = -0.6 - globalVariable::Epsilon;
        const double z0 = -1.35 - globalVariable::Epsilon;
        const double lx = 1.2, ly = 1.2, lz = 2.7;
        nx = (int)(lx / h) + 1, ny = (int)(ly / h) + 1, nz = (int)(lz / h) + 1;
        mesh_t Th(nx, ny, nz, x0, y0, z0, lx, ly, lz);
        #endif

        dT                     = cfl_number * h;
        total_number_iteration = int(tfinal / dT);
        dT                     = tfinal / total_number_iteration;

        hs.at(j)  = h;
        dts.at(j) = dT;

        if (iterations > 1) {
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "--------------------------------------------" << '\n';
            std::cout << "Iteration " << j + 1 << "/" << iterations << '\n';
        }

        std::cout << "h  = " << h << '\n';
        std::cout << "dT = " << dT << '\n';

        fespace_t Vh(Th, *FE_space[m-1]);
        // 1D Time mesh
        double final_time = total_number_iteration * dT;
        Mesh1 Qh(total_number_iteration + 1, t0, final_time);

        // 1D Time space
        FESpace1 Ih(Qh, *FE_time[k]);
        const Uint ndof_time_slab = Ih[0].NbDoF();

        // Velocity field
        #if defined(merging2D)
        Lagrange2 FEvelocity(1);
        #elif defined(merging3D)
        Lagrange3 FEvelocity(1);
        #endif
        fespace_t VelVh(Th, FEvelocity);
        std::vector<fct_t> vel(nbTime);

        // Declare time dependent interface
        TimeInterface<mesh_t> interface(qTime);

        // Visualization of the level set
        fespace_t Lh(Th, DataFE<mesh_t>::P1);
        double dt_levelSet = dT / (nbTime - 1);
        std::vector<fct_t> ls(nbTime);
        for (int i = 0; i < nbTime; i++)
            ls[i].init(Lh, fun_levelSet, 0.);

        std::cout << "Number of time slabs \t : \t " << total_number_iteration << '\n';

        int iter                  = 0;
        double mass_last_previous = 0., mass_initial = 0., mass_last = 0.;
        double intF = 0, intF_total = 0;
        double global_conservation_error = 0;
        std::vector<double> global_conservation_errors_t;

        // Iterate over time-slabs
        while (iter < total_number_iteration) {

            int current_iteration = iter;
            double current_time   = iter * dT;
            const TimeSlab &In(Ih[iter]);

            std::cout << " -------------------------------------------------------\n";
            std::cout << " -------------------------------------------------------\n";
            std::cout << " Iteration \t : \t" << iter + 1 << "/" << total_number_iteration << '\n';
            std::cout << " Time      \t : \t" << current_iteration * dT << '\n';
            std::cout << "dT = " << dT << '\n';

            
            // Initialization of the interface in each quadrature point
            for (int i = 0; i < nbTime; ++i) {

                R tt  = In.Pt(R1(qTime(i).x));

                vel[i].init(VelVh, fun_velocity, tt);

                ls[i].init(Lh, fun_levelSet, tt);

                interface.init(i, Th, ls[i]);

            }

            
            // Create active meshes
            activemesh_t Thi(Th);

            Thi.truncate(interface, 1);

            //  Cut FE space
            cut_fespace_t Wh(Thi, Vh);

            // Initialize the convection-diffusion problem
            convdiff.initSpace(Wh, In);

            // Objects needed for the weak form
            Normal n;

            // Test and Trial functions
            fun_test_t u(Wh, 1), v(Wh, 1);

            std::vector<double> data_init(convdiff.get_nb_dof(), 0.); // initial data total
            std::span<double> data_uh0 = std::span<double>(data_init.data(), Wh.NbDoF());

            if (iter == 0) {
                interpolate(Wh, data_uh0, fun_uBulkInit);
            } else {
                convdiff.initialSolution(data_init);
            }

            std::vector<double> data_all(data_init);
            std::span<double> data_uh = std::span<double>(data_all.data(), Wh.NbDoF() * In.NbDoF());

            fct_t fh(Wh, In, fun_rhsBulk);
            fct_t uh0(Wh, data_uh0);
            fct_t uh(Wh, In, data_uh);

// Variational formulation
#if defined(non_conservative)
            convdiff.addBilinear(+innerProduct(u, v), Thi, 0, In);
            // Impose initial condition
            convdiff.addLinear(+innerProduct(uh0.expr(), v), Thi, 0, In);
            convdiff.addBilinear(+innerProduct(dt(u), v) + innerProduct(D * grad(u), grad(v)), Thi, In);

            for (int i = 0; i < nbTime; ++i) {
                convdiff.addBilinear(+innerProduct((vel[i].exprList() * grad(u)), v), Thi, In, i);
            }

#elif defined(conservative)
            convdiff.addBilinear(+innerProduct(u, v), Thi, (int)lastQuadTime, In);

            // Impose initial condition
            if (iter == 0) {
                // convdiff.addLinearExact(fun_uBulk, +innerProduct(1, v), Thi, 0, In);
                convdiff.addLinear(fun_uBulkInit, +innerProduct(1, v), Thi, 0, In);
            } else {
                convdiff.addLinear(+innerProduct(uh0.expr(), v), Thi, 0, In);
            }

            // convdiff.addLinear(+innerProduct(uh0.expr(), v), Thi, 0, In);

            convdiff.addBilinear(-innerProduct(u, dt(v)) + innerProduct(D * grad(u), grad(v)), Thi, In);

            for (int i = 0; i < nbTime; ++i) {
                convdiff.addBilinear(-innerProduct(u, (vel[i].exprList() * grad(v))), Thi, In, i);
            }
#endif
            // Source function
            // convdiff.addLinearExact(fun_rhsBulk, +innerProduct(1, v), Thi, In);
            // convdiff.addLinear(+innerProduct(fh.expr(), v), Thi, In);

            // Stabilization

#if defined(fullstab)

            // convdiff.addFaceStabilization(
            //     + innerProduct(h * tau * jump(grad(u) * n), jump(grad(v) * n)) 
            //     // + innerProduct(h * h * h * tau * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n))
            //     , Thi, In);

            convdiff.addPatchStabilization(+innerProduct(tau / h / h * jump(u), jump(v)), Thi, In);

#elif defined(macro)

            MacroElementPartition<mesh_t> TimeMacro(Thi, delta);

            // convdiff.addFaceStabilization(
            //     + innerProduct(h * tau * jump(grad(u) * n), jump(grad(v) * n)) 
            //     // + innerProduct(h * h * h * tau2 * jump(grad(grad(u) * n) * n), jump(grad(grad(v) * n) * n)),
            //     , Thi, In, TimeMacro);

            convdiff.addPatchStabilization(+innerProduct(tau / h / h * jump(u), jump(v)), Thi, In, TimeMacro);

#endif

            // if (iter == total_number_iteration - 1) {
            //     matlab::Export(convdiff.mat_[0], "mat_" + std::to_string(j + 1) + ".dat");
            // }

            // Solve linear system
            convdiff.solve(solver);

            std::span<double> rhs = std::span<double>(convdiff.rhs_.data(), convdiff.get_nb_dof());
            data_all.assign(rhs.begin(), rhs.end());
            convdiff.saveSolution(data_all);

            // Compute error of numerical solution
            std::vector<double> sol(Wh.get_nb_dof());

            for (int n = 0; n < ndof_time_slab; n++) {
                std::vector<double> u_dof_n(data_uh.begin() + n * Wh.get_nb_dof(),
                                            data_uh.begin() + (n + 1) * Wh.get_nb_dof());
                std::transform(sol.begin(), sol.end(), u_dof_n.begin(), sol.begin(),
                               std::plus<double>()); // sum up all dofs
            }

            fct_t funuh(Wh, sol);

            // Compute conservation error
            intF = integral(Thi, In, fh, 0, qTime);
            intF_total += intF;

            // mass_last = integral_algoim(funuh, Thi, phi, In, qTime, lastQuadTime, quadrature_order_space);
            mass_last          = integral(Thi, funuh, 0, lastQuadTime);

            if (iter == 0) {
                mass_initial       = integral(Thi, fun_uBulkInit, 0, 0);
                mass_last_previous = mass_initial;
            }

            global_conservation_error = (mass_last - mass_initial);

            std::cout << "global_conservation_error: " << global_conservation_error << "\n";

            mass_last_previous = mass_last; // set current last to previous last for next time slab

            global_conservation_errors[j] = global_conservation_error;
            global_conservation_errors_t.push_back(global_conservation_error);


            if ((iterations == 1) && (h > 0.001)) {
                //Fun_h sol_h(Wh, sol);
                Paraview<mesh_t> writerTh(Th, "Th.vtk");
                #if defined(merging2D)
                Paraview<mesh_t> writer(Thi, "merging_circles_" + std::to_string(iter + 1) + ".vtk");
                #elif defined(merging3D)
                Paraview<mesh_t> writer(Thi, "merging_spheres_" + std::to_string(iter + 1) + ".vtk");
                #endif
                writer.add(uh0, "bulk_0", 0, 1);
                writer.add(funuh, "bulk_N", 0, 1);

                // fct_t uBex(Wh, fun_uBulk, current_time);
                // fct_t uBex_N(Wh, fun_uBulk, current_time+dT);
                
                // writer.add(ls[0], "levelSet0", 0, 1);
                // writer.add(ls[1], "levelSet1", 0, 1);
                // writer.add(ls[2], "levelSet2", 0, 1);
                writer.writeActiveMesh(Thi, "ActiveMesh" + std::to_string(iter + 1) + ".vtk");
            }

            iter++;
        }


        if (iterations == 1) {
            std::cout << "\n";
            std::cout << "Global conservation errors = [";
            for (auto &err : global_conservation_errors_t) {

                std::cout << err;

                std::cout << ", ";
            }
            std::cout << "]\n";
        }

        h *= 0.5;
        // dT *= 0.5;
    }

    std::cout << '\n';
    std::cout << "Global Conservation Errors = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << global_conservation_errors.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "h = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << hs.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';
    std::cout << '\n';

    std::cout << "dT = [";
    for (int i = 0; i < iterations; i++) {

        std::cout << dts.at(i);
        if (i < iterations - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << '\n';

    return 0;
}
