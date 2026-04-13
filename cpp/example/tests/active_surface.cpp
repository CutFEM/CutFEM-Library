/*
active_mesh_surfaces file is part of CutFEM-Library.

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

// Dependencies
#include "../cutfem.hpp"

#include <cmath>
#include <filesystem>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <iomanip>
// #include <gmsh.h> // Dependency for gmsh function

// Type Aliases for Readability
using mesh_t           = Mesh3;
using fun_test_t       = TestFunction<mesh_t>;
using fct_t            = FunFEM<mesh_t>;
using activemesh_t     = ActiveMesh<mesh_t>;
using interface_t      = Interface<mesh_t>;
using time_interface_t = TimeInterface<mesh_t>;
using space_t          = GFESpace<mesh_t>;
using cutspace_t       = CutFESpace<mesh_t>;
using fem_t            = FEM<mesh_t>;
using cutfem_t         = CutFEM<mesh_t>;
using lagrange_t       = Lagrange3;
using paraview_t       = Paraview<mesh_t>;
using Rd_t             = mesh_t::Rd;


struct ReInit {
    // Reinitialization parameters
    int number_iteration             = 4;
    double epsilon_diffusion         = 1e-3;
    double dt                        = 1e-3;

    // Mass correction parameters
    bool mass_correction             = true;
    double precision_correction      = 1e-6;
    int max_iteration_correction     = 10;

    void info() const {
        std::cout << "\n--- Reinitialization Settings ---\n";
        std::cout << " number_iteration         -> " << number_iteration << "\n"
                  << " epsilon_diffusion        -> " << epsilon_diffusion << "\n"
                  << " dt                       -> " << dt << "\n"
                  << " mass_correction          -> " << (mass_correction ? "ON" : "OFF") << "\n"
                  << " precision_correction     -> " << precision_correction << "\n"
                  << " max_iteration_correction -> " << max_iteration_correction << "\n"
                  << "-------------------------------\n" << std::endl;
   }
};


using namespace globalVariable; // to access some globally defined constants


// Numerical examples namespace
namespace Test1 {

    const std::string example = "8april_reinit_interface_refinement";

    // --- Physical and Numerical Constants ---
    const double D  = 1.0; // Diffusion coefficient
    const double xi = 5.;
    const double eta0 = 0.01;
    const double eta1 = 1.0;
    const double u_max_init = 0.;

    // Penalty Parameters
    const double boundary_penalty_stokes    = 100.0;
    const double interface_penalty_stokes   = 100.0;
    const double ghost_penalty_stokes       = 1e-2;
    const double ghost_penalty_conc         = 1e-1;
    const double ghost_penalty_surf_tension = 1e-1;
    const double ghost_penalty_curv         = 1.0;

    // Time and Geometry Constants
    const double T      = 4.0;
    // const double T      = 0.5;
    const double radius = 1.0;

    // Background Domain Definition
    const double X0 = -1.5;
    const double Y0 = -1.5;
    const double Z0 = -1.5;
    const double DX = 3.0;
    const double DY = 3.0;
    const double DZ = 3.0;

    // NEW STATE: Elongated in the x-direction to allow the polarized volume to "swim"
    // const double X0 = -1.5;
    // const double Y0 = -3.;
    // const double Z0 = -3.;
    // const double DX = 6.0;  // Increased to give 4.5 units of space in the +x direction
    // const double DY = 6.0;
    // const double DZ = 6.0;

    // Random Perturbation Generator (Static, initialized once)
    static std::mt19937_64 rng(std::random_device{}());
    static std::uniform_real_distribution<double> dist(-1e-2, 1e-2);


    // --- Function Definitions ---

    double fun_levelset(double* P, int i) {
        const double x = P[0], y = P[1], z = P[2];
        return std::sqrt(x*x + y*y + z*z) - radius + 1e-14;
    }
    double rhs(double *P, int component, double t) {
        return 0.0;
    }

    double cinitial(double* P, int component) {
        // OLD STATE: Uncorrelated grid-scale random noise
        // double perturbation = dist(rng);
        
        // NEW STATE: Smooth spatial gradient along the x-axis
        // double perturbation = 1e-4 * P[0];
        double perturbation = 0.2 * P[0];
        return 1.0 + perturbation;
    }

    double c0initial(double* P, int component) {
        return 1.0;
    }

    double uinitial(double *P, int component) {
        return 0.0;
    }

    double uboundary(double *P, int component) {
        return 0.0;
    }

} // namespace Test1

// Define the test case for the preprocessor check
#define TEST1
#define USE_CURVATURE

// Select the test case and apply 'using namespace' globally
#if defined(TEST1)
    using namespace Test1;
#else
#error "No example defined"
#endif

// Class to manage the implicit geometry and level set function
class ImplicitGeometry {

public:
    // Constructors
    ImplicitGeometry(const mesh_t &Th, const QuadratureFormular1d &q_time)
        : Lh_(Th, DataFE<mesh_t>::P1), levelsets_(q_time.n), interface_(q_time), n_(q_time.n) {}

    ImplicitGeometry(const mesh_t &Th, const QuadratureFormular1d &q_time, ReInit &reinit)
        : Lh_(Th, DataFE<mesh_t>::P1), levelsets_(q_time.n), interface_(q_time), n_(q_time.n), reinit_settings_(reinit) {}

    // Methods
    void initialize(std::function<double(double *, int)> phi) {
        for (auto& ls : levelsets_) {
            ls.init(Lh_, phi);
        }
    }

    fct_t update(fct_t &vel, double dt, bool reinit = false, bool shift_barycenter = false) {
        const int domain_inner = 1;
        constexpr int D_mesh = mesh_t::D; // Use D_mesh to avoid confusion with D from Test1 namespace
        fct_t applied_shift = vel;
        std::fill(applied_shift.v.begin(), applied_shift.v.end(), 0.0);

        // Shift last level set to the first position for the new time slab
        std::copy(levelsets_[n_-1].array().begin(), levelsets_[n_-1].array().end(), levelsets_[0].array().begin());

        // Loop over the quadrature points in [t_{n-1}, t_n]
        for (int i = 0; i < n_; ++i) {

            // Perform reinitialization at the start of the time slab (i=1 corresponds to the first intermediate step)
            if (reinit && (i == 1)) {
                InterfaceLevelSet<mesh_t> interface0(Lh_.Th, levelsets_[i]);
                const double vol0 = compute_volume(interface0, domain_inner, i);
                std::cout << "Volume before reinitialization: " << vol0 << "\n";

                // Perform reinitialization (changes levelsets_[i])
                reinitialize(interface0, levelsets_[i]);

                InterfaceLevelSet<mesh_t> interface1(Lh_.Th, levelsets_[i]);
                const double vol1 = compute_volume(interface1, domain_inner, i);
                std::cout << "Volume after reinitialization: " << vol1 << "\n";

                if (reinit_settings_.mass_correction) {
                    // Correct mass, returns the shift amount, delta
                    correct(levelsets_[i], vol0, i);

                    InterfaceLevelSet<mesh_t> interface2(Lh_.Th, levelsets_[i]);
                    double vol3 = compute_volume(interface2, domain_inner, i);
                    std::cout << "Volume after mass correction: " << vol3 << "\n";
                }
            }

            interface_.init(i, Lh_.Th, levelsets_[i]);

            if (shift_barycenter) {
                const double vol = compute_volume(*interface_[i], domain_inner, i);
                std::array<double, 3> vel_ref{0.0, 0.0, 0.0};
                activemesh_t Thi(Lh_.Th);
                const int dom_truncate = 1; // truncates the part where phi > 0
                Thi.truncate(*interface_[i], dom_truncate);

                // Compute volume-averaged velocity
                for (int d = 0; d < D_mesh; ++d) {
                    vel_ref[d] = integral(Thi, vel, d) / vol;
                }

                if (i < n_ - 1) {
                    // Apply velocity shift for advection step
                    fct_t vel_shifted = vel;
                    std::size_t nb_dofs_per_comp = vel_shifted.size() / D_mesh;
                    std::size_t idx = 0;
                    for (std::size_t j = 0; j < nb_dofs_per_comp; ++j) {
                        for (int d = 0; d < D_mesh; ++d) {
                            vel_shifted.array()[idx++] -= vel_ref[d];
                        }
                    }

                    // Advect level set
                    LevelSet::move(levelsets_[i], vel_shifted, vel_shifted, dt, levelsets_[i + 1]);
                } else { // last time quadrature point
                    // Keep the frame velocity as a full vector field so shifting can be done outside.
                    std::size_t idx = 0;
                    const std::size_t nb_dofs_per_comp = applied_shift.size() / D_mesh;
                    for (std::size_t j = 0; j < nb_dofs_per_comp; ++j) {
                        for (int d = 0; d < D_mesh; ++d) {
                            applied_shift.array()[idx++] = vel_ref[d];
                        }
                    }
                }
            } else {
                // Advect level set without velocity shift
                if (i < n_ - 1) {
                    LevelSet::move(levelsets_[i], vel, vel, dt, levelsets_[i + 1]);
                }
            }
        }

        return applied_shift;
    }

    double compute_volume(const interface_t &gamma, int domain, int itq) const {
        activemesh_t Thi(Lh_.Th);
        // domain 1 = inner (positive level set sign), so truncate with sign 1
        const int dom_truncate = (domain == 1) ? 1 : -1;
        Thi.truncate(gamma, dom_truncate);
        return integral(Thi, 1.0);
    }


    void reinitialize(const interface_t &gamma, fct_t &phi) {
        const double dt = reinit_settings_.dt;
        const double epsilon_diffusion = reinit_settings_.epsilon_diffusion;

        using FElement_t = space_t::FElement;
        using QF_t       = FElement_t::QF;

        // Using order 5 quadrature for reinitialization
        const QF_t &qf(*QF_Simplex<Rd_t>(5));

        fem_t reinitialization(Lh_);
        fct_t levelset_init(Lh_);   // "old" levelset
        levelset_init.init(phi);

        // KNM for basis function and derivatives
        KNMK<double> fu(Lh_[0].NbDoF(), 1, op_Dall);

        int iter = 0;
        while (iter < reinit_settings_.number_iteration) {
            std::fill(reinitialization.rhs_.begin(), reinitialization.rhs_.end(), 0.0);

            // Assemble
            for (int k = Lh_.first_element(); k < Lh_.last_element(); k += Lh_.next_element()) {
                const FElement_t &FK(Lh_[k]);
                const int kb = Lh_.Th(FK.T);
                const Partition cutK = gamma.get_partition(kb);
                const double h = FK.T.hElement();

                // Loop over cut element partition parts
                int loop_over_both_sides = 0;
                for (auto it = cutK.element_begin(loop_over_both_sides); it != cutK.element_end(loop_over_both_sides); ++it) {
                    const double meas   = cutK.measure(it);
                    const double signU0 = cutK.whatSign(it);
                    // Streamline diffusion parameter Tsd
                    // const double Tsd    = 1.0 / (std::sqrt(1.0 / dt / dt + 1.0 / h / h));
                    const double Tsd    = 0.5 / (std::sqrt(1.0 / dt / dt + 1.0 / h / h));

                    for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                        GQuadraturePoint<Rd_t> ip(qf[ipq]);
                        Rd_t mip = cutK.mapToPhysicalElement(it, ip);
                        const double Cint = meas * ip.getWeight();
                        FK.BF(FK.T.toKref(mip), fu);

                        double phix = phi.eval(k, mip, 0, op_id);
                        Rd_t dphi;
                        dphi.x = phi.eval(k, mip, 0, op_dx);
                        dphi.y = phi.eval(k, mip, 0, op_dy);
                        dphi.z = phi.eval(k, mip, 0, op_dz);

                        // Normalize gradient, handle near-zero norm
                        double norm_dphi  = dphi.norm();
                        if (std::fabs(norm_dphi) < 1e-14) {
                            norm_dphi = 1e-9;
                        }

                        const Rd_t w = dphi / norm_dphi * signU0;
                        const double ff = signU0;

                        // Explicit Euler (LHS and RHS assembly)
                        for (int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
                            // Loop for the bilinear form (LHS)
                            for (int j = FK.dfcbegin(0); j < FK.dfcend(0); ++j) {
                                reinitialization(FK.loc2glb(i), FK.loc2glb(j)) += Cint * (
                                    (1.0 / dt) * fu(j, 0, op_id) * (fu(i, 0, op_id) + Tsd * (w.x * fu(i, 0, op_dx) + w.y * fu(i, 0, op_dy) + w.z * fu(i, 0, op_dz))) +
                                    epsilon_diffusion * (fu(i, 0, op_dx) * fu(j, 0, op_dx) + fu(i, 0, op_dy) * fu(j, 0, op_dy) + fu(i, 0, op_dz) * fu(j, 0, op_dz))
                                );
                            }
                            // Linear form (RHS)
                            reinitialization(FK.loc2glb(i)) += Cint * (
                                ((1.0 / dt) * phix - (w, dphi) + ff) * (fu(i, 0, op_id) + Tsd * (w.x * fu(i, 0, op_dx) + w.y * fu(i, 0, op_dy) + w.z * fu(i, 0, op_dz)))
                            );
                        }
                    }
                }
            }

            // Solve
            reinitialization.solve("mumps");
            phi.init(reinitialization.rhs_);
            iter++;
        }
    }


    double correct(fct_t &phi, double initial_volume, int itq) {
        int domain_inner = 1;
        const double dV = reinit_settings_.precision_correction;

        // Compute volume before correction (vol0) and initial volume difference (diff0)
        InterfaceLevelSet<mesh_t> newInterface0(Lh_.Th, phi);
        double vol0 = compute_volume(newInterface0, domain_inner, itq);
        const double diffInitial = vol0 - initial_volume;
        double diff0 = diffInitial;

        if (std::fabs(diff0) < dV) {
            return 0.0;
        }

        double dz0   = 0.0;
        double dz1   = diff0; // Approximation for first step
        double diff1 = diff0;
        int nn       = 0;

        // Bisection search to bracket the zero
        while (diff1 * diff0 > 0) {
            std::transform(phi.v.begin(), phi.v.end(), phi.v.begin(), [dz1](double value) { return value + dz1; });
            InterfaceLevelSet<mesh_t> newInterface1(Lh_.Th, phi);
            double vol1 = compute_volume(newInterface1, domain_inner, itq);
            diff1       = vol1 - initial_volume;
            nn += 1;
            if (std::fabs(diff1) < dV) {
                // FIX: If we hit the tolerance during the bracketing phase, 
                // we leave phi shifted and return.
                return dz1; 
            }
        }
        
        // Undo the transformation to reset phi for the precise root finding
        dz1 = nn * dz1;
        std::transform(phi.v.begin(), phi.v.end(), phi.v.begin(), [dz1](double value) { return value - dz1; });

        // Root finding using bisection method
        for (int i = 0; i < reinit_settings_.max_iteration_correction; ++i) {
            double dz2 = 0.5 * (dz0 + dz1);
            
            // Apply shift
            std::transform(phi.v.begin(), phi.v.end(), phi.v.begin(), [dz2](double value) { return value + dz2; });

            InterfaceLevelSet<mesh_t> newInterface2(Lh_.Th, phi);
            double vol2  = compute_volume(newInterface2, domain_inner, itq);
            double diff2 = vol2 - initial_volume;

            if (std::fabs(diff2) < dV) {
                // FIX: We reached the required precision! 
                // DO NOT undo the shift. Leave phi modified and return.
                std::cout << " Mass correction error \t" << diff2 << std::endl;
                return dz2;
            }

            if (i == reinit_settings_.max_iteration_correction - 1 && std::fabs(diffInitial) > std::fabs(diff2)) {
                // FIX: Max iterations reached, but this shift is better than what we started with.
                // DO NOT undo the shift. Keep the best state.
                std::cout
                    << " Correction of the levelset function did not reach the required precision ("
                    << std::setprecision(10) << dV << ")" << std::endl;
                std::cout << " Mass correction error \t" << diff2 << std::endl;
                return dz2;
            }
            
            // This specific shift wasn't accurate enough and we haven't hit the max iterations.
            // Now we DO undo the shift so we can try the next bisection interval.
            std::transform(phi.v.begin(), phi.v.end(), phi.v.begin(), [dz2](double value) { return value - dz2; });

            // Update bisection interval
            if (diff0 * diff2 < 0) {
                dz1   = dz2;
                diff1 = diff2;
            } else {
                dz0   = dz2;
                diff0 = diff2;
            }
        }
        
        std::cout << " Correction of levelSet failed!! Returning -1.\n";
        return -1.0;
    }

    // Accessors for level set data
    fct_t &levelset(size_t i) { return levelsets_.at(i); }
    std::vector<fct_t> &levelsets() { return levelsets_; }
    const fct_t &levelset(size_t i) const { return levelsets_.at(i); }
    const std::vector<fct_t> &levelsets() const { return levelsets_; }

    // Accessors for interface data
    time_interface_t &time_interface() { return interface_; }
    InterfaceLevelSet<mesh_t> &interface(size_t i) { return *static_cast<InterfaceLevelSet<mesh_t>*>(interface_[i]); }
    const time_interface_t &time_interface() const { return interface_; }
    const InterfaceLevelSet<mesh_t> &interface(size_t i) const { return *static_cast<InterfaceLevelSet<mesh_t>*>(interface_[i]); }

    fct_t normal_phi(const space_t &V_normal, const int itq) {
        // L2 projection of grad(phi)/norm(grad(phi)) to get the normal on the background mesh
        fem_t l2_projection(V_normal);
        fun_test_t n_proj(V_normal, 3), v_proj(V_normal, 3);
        fun_test_t v_projx(V_normal, 1, 0), v_projy(V_normal, 1, 1), v_projz(V_normal, 1, 2);

        // Bilinear form: (n_proj, v_proj)_L2
        l2_projection.addBilinear(
            + innerProduct(n_proj, v_proj),
            Lh_.Th
        );

        auto phi = levelset(itq);
        // Correct way to write inv_norm_grad_phi:
        auto inv_norm_grad_phi = pow(sqrt(dx(phi.expr())*dx(phi.expr()) + dy(phi.expr())*dy(phi.expr()) + dz(phi.expr())*dz(phi.expr())), -1);

        // Linear form: - (grad(phi)/|grad(phi)|, v_proj)_L2
        l2_projection.addLinear(
            - innerProduct(dx(phi.expr()), inv_norm_grad_phi * v_projx)
            - innerProduct(dy(phi.expr()), inv_norm_grad_phi * v_projy)
            - innerProduct(dz(phi.expr()), inv_norm_grad_phi * v_projz),
            Lh_.Th
        );

        l2_projection.solve("mumps");
        data_normal = l2_projection.rhs_;
        return fct_t(V_normal, data_normal);
    }


    fct_t tangential_projection(fct_t& vel, fct_t &normal, size_t itq) {
        // Create surface mesh and corresponding cut space
        surface_mesh_ = std::make_unique<activemesh_t>(Lh_.Th);
        surface_mesh_->createSurfaceMesh(*interface_[itq]);
        cutVelVh_ = std::make_unique<cutspace_t>(*surface_mesh_, *vel.Vh);

        ProblemOption option;
        option.order_space_element_quadrature_ = 8;
        cutfem_t projection(*cutVelVh_, option);

        fun_test_t u(*cutVelVh_, 3), v(*cutVelVh_, 3);
        fun_test_t vx(*cutVelVh_, 1, 0), vy(*cutVelVh_, 1, 1), vz(*cutVelVh_, 1, 2);

        const double h = surface_mesh_->Th.get_mesh_size();

        // Bilinear form: L2 projection
        projection.addBilinear(
            + innerProduct(u, v),
            *interface_[itq]
        );

        // Patch Stabilization
        projection.addPatchStabilization(
            + innerProduct(0.01 * std::pow(h, 0) * jump(u), jump(v)),
            *surface_mesh_
        );

        // Element/Volume stabilization (using BaseFEM::addBilinear)
        projection.BaseFEM::addBilinear(
            + innerProduct(0.01 * std::pow(h, 2) * (normal.exprList() * grad(u)), (normal.exprList() * grad(v))),
            *surface_mesh_
        );

        // Linear form: Tangential part of vel on the interface
        // This calculates ( (I - n*n^T) * vel, v )_{\Gamma}
        projection.addLinear(
            // x-component: (vel_x - nx*(nx*vel_x + ny*vel_y + nz*vel_z)) * vx
            + innerProduct(vel.expr(0), vx)
            - innerProduct(normal.expr(0) * normal.expr(0) * vel.expr(0), vx)
            - innerProduct(normal.expr(0) * normal.expr(1) * vel.expr(1), vx)
            - innerProduct(normal.expr(0) * normal.expr(2) * vel.expr(2), vx)

            // y-component: (vel_y - ny*(nx*vel_x + ny*vel_y + nz*vel_z)) * vy
            - innerProduct(normal.expr(1) * normal.expr(0) * vel.expr(0), vy)
            + innerProduct(vel.expr(1), vy)
            - innerProduct(normal.expr(1) * normal.expr(1) * vel.expr(1), vy)
            - innerProduct(normal.expr(1) * normal.expr(2) * vel.expr(2), vy)

            // z-component: (vel_z - nz*(nx*vel_x + ny*vel_y + nz*vel_z)) * vz
            - innerProduct(normal.expr(2) * normal.expr(0) * vel.expr(0), vz)
            - innerProduct(normal.expr(2) * normal.expr(1) * vel.expr(1), vz)
            + innerProduct(vel.expr(2), vz) // Corrected from vel.expr(1) in original
            - innerProduct(normal.expr(2) * normal.expr(2) * vel.expr(2), vz),
            *interface_[itq]
        );

        projection.solve("mumps");
        data_tangential_projection = projection.rhs_;
        return fct_t(*cutVelVh_, data_tangential_projection);
    }

    size_t n() const { return n_; }
    const activemesh_t &get_surface_mesh() const { return *surface_mesh_; }

private:
    space_t Lh_; // space for level set function
    std::vector<fct_t> levelsets_;
    time_interface_t interface_;
    std::vector<double> data_normal;
    std::vector<double> data_tangential_projection;
    std::unique_ptr<activemesh_t> surface_mesh_;
    std::unique_ptr<cutspace_t> cutVelVh_;

    ReInit reinit_settings_;
    // std::vector<double> data; // Not used in this scope

    const size_t n_;
};


// Class for solving the concentration advection-diffusion problem
class SolveConcentration {

public:
    explicit SolveConcentration(const QuadratureFormular1d &q_time) {
        problem_object_ = std::make_unique<cutfem_t>(q_time);
    }

    void initialize(const space_t &V, const TimeSlab &In, const time_interface_t &interface) {
        active_mesh_ = std::make_unique<activemesh_t>(V.Th);
        active_mesh_->createSurfaceMesh(interface);

        Vh_ = std::make_unique<cutspace_t>(*active_mesh_, V);
        problem_object_->initSpace(*Vh_, In);
    }

    // Overload 1: Initialize data by calling initialSolution (likely uses zeros)
    void initialize_data() {
        initial_data_.resize(problem_object_->get_nb_dof(), 0.0);
        problem_object_->initialSolution(initial_data_);
    }

    // Overload 2: Initialize data by interpolation with a function
    void initialize_data(std::function<double(double*, int)> c0) {
        initial_data_.resize(problem_object_->get_nb_dof(), 0.0);
        std::span<double> span_data(initial_data_.data(), Vh_->NbDoF());
        interpolate(*Vh_, span_data, c0);
    }

    void solve(const ImplicitGeometry &geometry, const TimeSlab &In, fct_t &vel, fct_t &sol) {
        fct_t ch0(*Vh_, initial_data_);

        fun_test_t c(*Vh_, 1), r(*Vh_, 1);

        const double h = active_mesh_->Th.get_mesh_size();

        Normal n;

        problem_object_->addBilinear(+innerProduct(c, r), geometry.interface(0), In, 0);
        problem_object_->addLinear(+innerProduct(ch0.expr(), r), geometry.interface(0), In, 0);

        problem_object_->addBilinear(
            + innerProduct(dt(c), r)
            + innerProduct((vel.exprList() * grad(c)), r)
            + innerProduct(divS(vel) * c, r) 
            + innerProduct(D * gradS(c), gradS(r)), 
            geometry.time_interface(),
            In
        );

        // Stabilization
        // problem_object_->addPatchStabilization(
        //     + innerProduct(ghost_penalty_conc * std::pow(h, -3) * jump(c), jump(r)),
        //     *active_mesh_, In
        // );

        problem_object_->addFaceStabilization(
            + innerProduct(ghost_penalty_conc * jump(grad(c) * n), jump(grad(r) * n)),
            *active_mesh_, In);

        // Solve linear system
        problem_object_->solve(solver_name);

        // Map solution back to the background mesh
        problem_object_->saveSolutionBackMesh(problem_object_->rhs_, sol);
    }

private:
    std::unique_ptr<cutfem_t> problem_object_;
    std::unique_ptr<activemesh_t> active_mesh_;
    std::unique_ptr<cutspace_t> Vh_;

    std::vector<double> initial_data_;

    const std::string solver_name = "mumps";
};


// Class for solving the Stokes interface problem
class SolveStokes {
public:
    SolveStokes(const space_t &V, const space_t &Q, const interface_t &interface)
        : interface_(interface)
    {
        active_mesh_  = std::make_unique<activemesh_t>(V.Th, interface_);
        surface_mesh_ = std::make_unique<activemesh_t>(V.Th);
        surface_mesh_->createSurfaceMesh(interface_);

        Vh_ = std::make_unique<cutspace_t>(*active_mesh_, V);
        Qh_ = std::make_unique<cutspace_t>(*active_mesh_, Q);

        problem_object_ = std::make_unique<cutfem_t>(*Vh_);
        problem_object_->add(*Qh_); // Add pressure space
    }

    void solve(const fct_t &Hn, const fct_t &ch, const fct_t &c0h, const fct_t &u_boundary) {
        CutFEMParameter eta(2*eta0, 2*eta1);
    
        Normal n;
        // const double h = active_mesh_->Th.get_largest_mesh_size();
        const double h = active_mesh_->Th.get_mesh_size();

        // Test and Trial functions
        fun_test_t u(*Vh_, 3), v(*Vh_, 3);
        fun_test_t p(*Qh_, 1), q(*Qh_, 1);

        fun_test_t v0(*Vh_, 1, 0), v1(*Vh_, 1, 1), v2(*Vh_, 1, 2);

        auto t_solve_start = now();

        problem_object_->addBilinear(
            + contractProduct(Eps(u), eta * Eps(v))
            - innerProduct(p, div(v))
            + innerProduct(div(u), q)
            , *active_mesh_);


        // Outer boundary terms
        problem_object_->addBilinear(
            - innerProduct(eta * Eps(u) * n, v)
            + innerProduct(u, eta * Eps(v) * n)
            + innerProduct(boundary_penalty_stokes / h * eta * u, v)
            + innerProduct(p, v*n)
            - innerProduct(u*n, q)
            , active_mesh_->Th
            , INTEGRAL_BOUNDARY
        );

        problem_object_->addLinear(
            + innerProduct(u_boundary.exprList(), eta*Eps(v)*n)
            + innerProduct(u_boundary.exprList(), boundary_penalty_stokes / h * eta * v)
            - innerProduct(u_boundary*n, q)
            , active_mesh_->Th
            , INTEGRAL_BOUNDARY);

        // Interface terms
        problem_object_->addBilinear(
          - innerProduct(average(eta * Eps(u) * n), jump(v))
          + innerProduct(jump(u), average(eta * Eps(v) * n)) // added for anti-symmetry
          + innerProduct(average(p), jump(v * n))
          - innerProduct(jump(u * n), average(q))       // added for symmetry
          + innerProduct(interface_penalty_stokes / h * jump(u), jump(v))   // interface penalty for stability
          , interface_);

        // Concentration-dependent force (Marangoni/Surface tension force on RHS)
        // f(c) = 2*c^2/(c_0^2 + c^2)
        auto fh = 2.0 * ch.expr() * ch.expr() / (ch.expr() * ch.expr() + c0h.expr() * c0h.expr());
        // df/dc = 4*c*c_0^2 / (c_0^2 + c^2)^2
        auto dfhdc = 4.0 * ch.expr() * c0h.expr() * c0h.expr() / (
            (ch.expr() * ch.expr() + c0h.expr() * c0h.expr()) * (ch.expr() * ch.expr() + c0h.expr() * c0h.expr())
        );

        problem_object_->addLinear(
            // Surface tension term: xi * ( Hn*f(c) * n + grad_Gamma(f(c)) )
            // First term: xi * ( Hn * f(c), v)_Gamma
            // + xi*innerProduct(2./Test1::radius, fh*average(v*n))
            + xi * innerProduct(Hn.exprList(), fh * average(v))
            // Second term: xi * ( grad_Gamma(f(c)), v)_Gamma
            + xi * innerProduct(dxS(ch, 0), dfhdc * average(v0))
            + xi * innerProduct(dyS(ch, 0), dfhdc * average(v1))
            + xi * innerProduct(dzS(ch, 0), dfhdc * average(v2)),
            interface_
        );

        // Ghost Penalty Stabilization (Volume)
        // problem_object_->addPatchStabilization(
        //     + innerProduct(ghost_penalty_stokes * eta * std::pow(h, -2) * jump(u), jump(v))
        //     + innerProduct(ghost_penalty_stokes * std::pow(h, 0) * jump(p), jump(q)),
        //     *active_mesh_
        // );
        problem_object_->addFaceStabilization(
            + innerProduct(ghost_penalty_stokes * eta * std::pow(h, 1) * jump(grad(u)*n), jump(grad(v)*n)) 
            + innerProduct(ghost_penalty_stokes * eta * std::pow(h, 3) * jump(grad(grad(u)*n)*n), jump(grad(grad(v)*n)*n)) 
            + innerProduct(ghost_penalty_stokes * std::pow(h, 1) * jump(p), jump(q)) 
            + innerProduct(ghost_penalty_stokes * std::pow(h, 3) * jump(grad(p)*n), jump(grad(q)*n)) 
            , *active_mesh_
            );

        // Enforce mean zero pressure condition (Lagrange Multiplier)
        problem_object_->addLagrangeMultiplier(
            + innerProduct(1.0, q),
            0,
            *active_mesh_
        );

        auto t_solve_end = now();
        std::cout << "  * Time for assembling: " << seconds(t_solve_start, t_solve_end) << "\n";

        t_solve_start = now();
        problem_object_->solve("mumps");
        t_solve_end = now();

        std::cout << "  * Time for solving: " << seconds(t_solve_start, t_solve_end) << "\n";
    }

        void solve(const fct_t &Sh, const fct_t &u_boundary) {
        CutFEMParameter eta(2*eta0, 2*eta1);
    
        Normal n;
        // const double h = active_mesh_->Th.get_largest_mesh_size();
        const double h = active_mesh_->Th.get_mesh_size();

        // Test and Trial functions
        fun_test_t u(*Vh_, 3), v(*Vh_, 3);
        fun_test_t p(*Qh_, 1), q(*Qh_, 1);

        auto t_solve_start = now();

        problem_object_->addBilinear(
            + contractProduct(Eps(u), eta * Eps(v))
            - innerProduct(p, div(v))
            + innerProduct(div(u), q)
            , *active_mesh_);


        // Outer boundary terms
        problem_object_->addBilinear(
            - innerProduct(eta * Eps(u) * n, v)
            + innerProduct(u, eta * Eps(v) * n)
            + innerProduct(boundary_penalty_stokes / h * eta * u, v)
            + innerProduct(p, v*n)
            - innerProduct(u*n, q)
            , active_mesh_->Th
            , INTEGRAL_BOUNDARY
        );

        problem_object_->addLinear(
            + innerProduct(u_boundary.exprList(), eta*Eps(v)*n)
            + innerProduct(u_boundary.exprList(), boundary_penalty_stokes / h * eta * v)
            - innerProduct(u_boundary*n, q)
            , active_mesh_->Th
            , INTEGRAL_BOUNDARY);

        // Interface terms
        problem_object_->addBilinear(
          - innerProduct(average(eta * Eps(u) * n), jump(v))
          + innerProduct(jump(u), average(eta * Eps(v) * n)) // added for anti-symmetry
          + innerProduct(average(p), jump(v * n))
          - innerProduct(jump(u * n), average(q))       // added for symmetry
          + innerProduct(interface_penalty_stokes / h * jump(u), jump(v))   // interface penalty for stability
          , interface_);

        problem_object_->addLinear(
            + innerProduct(Sh.exprList(), average(v))
            , interface_
        );

        // Ghost Penalty Stabilization (Volume)
        // problem_object_->addPatchStabilization(
        //     + innerProduct(ghost_penalty_stokes * eta * std::pow(h, -2) * jump(u), jump(v))
        //     + innerProduct(ghost_penalty_stokes * std::pow(h, 0) * jump(p), jump(q)),
        //     *active_mesh_
        // );
        problem_object_->addFaceStabilization(
            + innerProduct(ghost_penalty_stokes * eta * std::pow(h, 1) * jump(grad(u)*n), jump(grad(v)*n)) 
            + innerProduct(ghost_penalty_stokes * eta * std::pow(h, 3) * jump(grad(grad(u)*n)*n), jump(grad(grad(v)*n)*n)) 
            + innerProduct(ghost_penalty_stokes * std::pow(h, 1) * jump(p), jump(q)) 
            + innerProduct(ghost_penalty_stokes * std::pow(h, 3) * jump(grad(p)*n), jump(grad(q)*n)) 
            , *active_mesh_
            );

        // Enforce mean zero pressure condition (Lagrange Multiplier)
        problem_object_->addLagrangeMultiplier(
            + innerProduct(1.0, q),
            0,
            *active_mesh_
        );

        auto t_solve_end = now();
        std::cout << "  * Time for assembling: " << seconds(t_solve_start, t_solve_end) << "\n";

        t_solve_start = now();
        problem_object_->solve("mumps");
        t_solve_end = now();

        std::cout << "  * Time for solving: " << seconds(t_solve_start, t_solve_end) << "\n";
    }


    void info() const {
        std::cout << "\n--- Stokes Problem Info ---\n";
        std::cout << "Number of DOFs: " << problem_object_->get_nb_dof() << "\n";
        std::cout << "xi = " << xi << ", eta0 = " << eta0 << ", eta1 = " << eta1 << "\n \n";
    }

    void getSolutionBackMesh(fct_t &vel, fct_t &levelset) {
        fct_t sol(*Vh_, problem_object_->rhs_);
        // Interpolate solution (which is defined on the cut elements) back to the background mesh
        interpolateOnBackGroundMesh(vel, sol, levelset);
    }

private:
    const interface_t &interface_;
    std::unique_ptr<activemesh_t> active_mesh_;
    std::unique_ptr<activemesh_t> surface_mesh_;
    std::unique_ptr<cutspace_t> Vh_;
    std::unique_ptr<cutspace_t> Qh_;
    std::unique_ptr<cutfem_t> problem_object_;
};


// Class for computing the mean curvature vector H*n
class Curvature {
public:
    Curvature(const space_t &V, const interface_t &interface)
        : interface_(interface)
    {
        surface_mesh_ = std::make_unique<activemesh_t>(V.Th);
        surface_mesh_->createSurfaceMesh(interface_);
        Vh_ = std::make_unique<cutspace_t>(*surface_mesh_, V);
    }

    fct_t solve() {
        constexpr int D_mesh = mesh_t::Rd::d; // Use D_mesh to avoid conflict
        cutfem_t problem(*Vh_);
        const double h = surface_mesh_->Th.get_mesh_size();

        TestFunction<mesh_t> H(*Vh_, D_mesh), v(*Vh_, D_mesh);
        Normal n;
        Rnm Id(D_mesh, D_mesh); // Identity matrix
        Id = 0.0;
        for (int i = 0; i < D_mesh; ++i)
            Id(i, i) = 1.0;

        problem.addBilinear(
            + innerProduct(H, v)
            // + 1e-1*std::pow(h, 0)*innerProduct(grad(H) * n, grad(v) * n)   // surface stabilization, try h^2
            , interface_
        );

        problem.addLinear(
            - contractProduct(Id, gradS(v)), // gradS is space-only gradient
            interface_
        );

        // Face Stabilization (Discontinuous Galerkin type)
        problem.addFaceStabilization(
            + innerProduct(ghost_penalty_curv * std::pow(h, 0) * jump(grad(H) * n), jump(grad(v) * n)),
            *surface_mesh_
        );
        problem.solve("mumps"); // Assuming "mumps" is the default solver
        data = problem.rhs_;
        return fct_t(*Vh_, data);
    }

    const activemesh_t &mesh() const { return *surface_mesh_; }

private:
    const interface_t &interface_;
    std::unique_ptr<activemesh_t> surface_mesh_;
    std::unique_ptr<cutspace_t> Vh_;
    std::vector<double> data;
};


class SurfaceTension {
public:
    SurfaceTension(const space_t &V, const interface_t &interface)
        : interface_(interface)
    {
        surface_mesh_ = std::make_unique<activemesh_t>(V.Th);
        surface_mesh_->createSurfaceMesh(interface_);
        Vh_ = std::make_unique<cutspace_t>(*surface_mesh_, V);
    }

    template <typename Fct>
    fct_t solve(const Fct &fc) {
        constexpr int D_mesh = mesh_t::Rd::d; // Use D_mesh to avoid conflict
        cutfem_t problem(*Vh_);
        const double h = surface_mesh_->Th.get_mesh_size();

        TestFunction<mesh_t> S(*Vh_, D_mesh), w(*Vh_, D_mesh);
        Normal n;
        Rnm Id(D_mesh, D_mesh); // Identity matrix
        Id = 0.0;
        for (int i = 0; i < D_mesh; ++i)
            Id(i, i) = 1.0;

        problem.addBilinear(
            + innerProduct(S, w)
            // + 1e-1*std::pow(h, 0)*innerProduct(grad(S) * n, grad(w) * n)   // surface stabilization
            , interface_
        );

        problem.addLinear(
            - contractProduct(Id, fc*gradS(w)), // gradS is space-only gradient
            interface_
        );

        // Face Stabilization (Discontinuous Galerkin type)
        problem.addFaceStabilization(
            + innerProduct(ghost_penalty_surf_tension * std::pow(h, 0) * jump(grad(S) * n), jump(grad(w) * n)),
            *surface_mesh_
        );
        problem.solve("mumps"); // Assuming "mumps" is the default solver
        data = problem.rhs_;
        return fct_t(*Vh_, data);
    }

    const activemesh_t &mesh() const { return *surface_mesh_; }

private:
    const interface_t &interface_;
    std::unique_ptr<activemesh_t> surface_mesh_;
    std::unique_ptr<cutspace_t> Vh_;
    std::vector<double> data;
};

/*
void generate_gmsh(int argc, char** argv, const double h_min, const double h_max) {
    gmsh::initialize(argc, argv);

    gmsh::model::add("cube");
    try {
        gmsh::model::occ::addBox(-3., -3., -3., 6., 6., 6., 1);
    } catch (...) {
        gmsh::logger::write("Could not create OpenCASCADE shapes: bye!");
        return;
    }

    gmsh::model::occ::synchronize();


    // THOMAS CODE THAT DOESN'T WORK WITH MUMPS

    // gmsh::model::geo::addPoint(0., 0., 0., 0.5, 1001);

    // gmsh::model::geo::synchronize();

    // gmsh::model::mesh::field::add("Distance", 1);
    // gmsh::model::mesh::field::setNumbers(1, "PointsList", {static_cast<double>(1001)});
    // gmsh::model::mesh::field::add("Threshold", 2);
    // gmsh::model::mesh::field::setNumber(2, "InField", 1);
    // // gmsh::model::mesh::field::setNumber(2, "SizeMin", 0.1);
    // // gmsh::model::mesh::field::setNumber(2, "SizeMax", 1.);
    // gmsh::model::mesh::field::setNumber(2, "SizeMin", h_min);
    // gmsh::model::mesh::field::setNumber(2, "SizeMax", h_max);
    // gmsh::model::mesh::field::setNumber(2, "DistMin", 2);
    // gmsh::model::mesh::field::setNumber(2, "DistMax", 3.5);

    // gmsh::model::mesh::field::setAsBackgroundMesh(2);
    // gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    // gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    // gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);


    // CHAT GPT'S ALTERNATIVE THAT WORKS
    // --- Sphere-driven sizing (Ball field) ---
    const int f = gmsh::model::mesh::field::add("Ball");
    gmsh::model::mesh::field::setNumber(f, "XCenter", 0.0);
    gmsh::model::mesh::field::setNumber(f, "YCenter", 0.0);
    gmsh::model::mesh::field::setNumber(f, "ZCenter", 0.0);
    gmsh::model::mesh::field::setNumber(f, "Radius", 1.5);  // sphere radius
    gmsh::model::mesh::field::setNumber(f, "VIn",   h_min);  // size inside radius
    gmsh::model::mesh::field::setNumber(f, "VOut",  h_max);  // size outside radius


    // Another box refined region
    // const int f = gmsh::model::mesh::field::add("Box");
    // gmsh::model::mesh::field::setNumber(f, "VIn",  h_min); // fine size inside
    // gmsh::model::mesh::field::setNumber(f, "VOut", h_max); // coarse size outside

    // // Refined box bounds (spans all x, half-size in y,z)
    // gmsh::model::mesh::field::setNumber(f, "XMin", -11.0);
    // gmsh::model::mesh::field::setNumber(f, "XMax",  3.0);
    // gmsh::model::mesh::field::setNumber(f, "YMin", -3.5); // half of 7
    // gmsh::model::mesh::field::setNumber(f, "YMax",  3.5);
    // gmsh::model::mesh::field::setNumber(f, "ZMin", -3.5);
    // gmsh::model::mesh::field::setNumber(f, "ZMax",  3.5);

    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    // Let the field dominate
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

    gmsh::model::mesh::generate(3);

    gmsh::write("cube.mesh");

    gmsh::finalize();
}


void refine_interface(int argc, char **argv, const double h_fine, const double h_coarse, const double band_width = 0.3) {
    gmsh::initialize(argc, argv);

    gmsh::model::add("cube_interface");
    try {
        gmsh::model::occ::addBox(X0, Y0, Z0, DX, DY, DZ, 1);
    } catch (...) {
        gmsh::logger::write("Could not create OpenCASCADE shapes: bye!");
        return;
    }

    gmsh::model::occ::synchronize();

    // Compute the absolute distance from the unit sphere surface: |sqrt(x²+y²+z²) - radius|
    // This equals zero on the interface and increases away from it in both directions.
    const int f_dist = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f_dist, "F",
        "Abs(Sqrt(x*x + y*y + z*z) - " + std::to_string(radius) + ")");

    // Threshold field: h_fine at the interface (distance=0), h_coarse beyond band_width
    const int f_thresh = gmsh::model::mesh::field::add("Threshold");
    gmsh::model::mesh::field::setNumber(f_thresh, "InField",  f_dist);
    gmsh::model::mesh::field::setNumber(f_thresh, "SizeMin",  h_fine);
    gmsh::model::mesh::field::setNumber(f_thresh, "SizeMax",  h_coarse);
    gmsh::model::mesh::field::setNumber(f_thresh, "DistMin",  0.0);
    gmsh::model::mesh::field::setNumber(f_thresh, "DistMax",  band_width);

    gmsh::model::mesh::field::setAsBackgroundMesh(f_thresh);

    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints",         0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature",      0);

    gmsh::model::mesh::generate(3);

    gmsh::write("cube_interface.mesh");

    gmsh::finalize();
}
*/

// --- Main Program ---
int main(int argc, char **argv) {
    MPIcf cfMPI(argc, argv);

    // omp_set_num_threads(4);

    // Setup output paths
    const std::string base_output_path = "../output_files/active_surfaces/" + example + "/";
    // const std::string base_output_path = "/cfs/klemming/scratch/s/smyrback/output_files/" + example + "/";
    std::string path_output_data(base_output_path + "data/");
    std::string path_output_figures(base_output_path + "paraview/");

    if (MPIcf::IamMaster()) {
        std::filesystem::create_directories(path_output_data);
        std::filesystem::create_directories(path_output_figures);
    }

    std::ofstream output_data(path_output_data + "output.dat", std::ofstream::out);


    bool reshift_center = true;
    bool reinit_on = true;

    // Create background mesh
    double h = 0.15;
    const int nx = static_cast<int>(DX / h) + 1;
    const int ny = static_cast<int>(DY / h) + 1;
    const int nz = static_cast<int>(DZ / h) + 1;
    mesh_t Th(nx, ny, nz, X0, Y0, Z0, DX, DY, DZ);
    // generate_gmsh(argc, argv, h, 2*h);
    // mesh_t Th("cube.mesh", MeshFormat::mesh_gmsh);
    // refine_interface(argc, argv, 0.08, 0.2);
    // mesh_t Th("cube_interface.mesh", MeshFormat::mesh_gmsh);
    // double h = Th.get_mesh_size();

    // h = Th.get_mesh_size();

    // Discretize in time (Calculate stable time step dT)
    // dT = 0.95*sqrt((eta0+eta1)*h^3 / (8*pi*xi))
    // double dT = 0.95 * std::sqrt((eta0 + eta1) * h * h * h / (8.0 * M_PI * xi));
    double dT = h/10;

    // Time integration quadrature
    constexpr size_t quadrature_order_time = 5;
    const QuadratureFormular1d &q_time(*Lobatto(quadrature_order_time));
    const int n_quad_pts_time    = q_time.n;
    const int last_quad_pt_time  = n_quad_pts_time - 1;

    int num_time_steps = static_cast<int>(T / dT);
    dT = T / num_time_steps;

    Mesh1 TimeDiscretization(num_time_steps + 1, 0.0, T);
    FESpace1 Ih(TimeDiscretization, DataFE<Mesh1>::P1Poly);
    // const int ndf_time_slab = Ih[0].NbDoF(); // Not used
    double dt_levelset = dT / (n_quad_pts_time - 1);

    // Set reinitialization parameters
    ReInit reinitialization;
    reinitialization.number_iteration = 5;
    reinitialization.epsilon_diffusion = 1e-2;
    reinitialization.dt = dT / 8.0;
    reinitialization.mass_correction = true;
    reinitialization.precision_correction = 1e-3;
    reinitialization.max_iteration_correction = 10;
    reinitialization.info();

    const int frequency_reinit = 3; 
    const int frequency_paraview = 1; 

    // --- Finite Element Spaces ---
    lagrange_t FEvelocity(2);
    lagrange_t FEcurv(1);
    lagrange_t FEnormal(1);

    space_t VelV(Th, FEvelocity);       // Velocity space (P2, vector)
    space_t CurvV(Th, FEcurv);          // Curvature space (P1, vector/scalar)
    space_t NormalV(Th, FEnormal);      // Normal vector space (P1, vector)
    space_t V(Th, DataFE<mesh_t>::P1);  // Concentration space (P1, scalar)
    // space_t V_int(Th, DataFE<mesh_t>::P2); // Not used
    space_t Q(Th, DataFE<mesh_t>::P1);  // Pressure space (P1, scalar)

    // --- Function Objects (initialization) ---
    fct_t vel(VelV, uinitial);
    fct_t u_boundary(VelV, uboundary);
    fct_t ch(V, cinitial);      // Concentration
    fct_t c0h(V, c0initial);    // Reference concentration

    // Output initialization
    if (MPIcf::IamMaster()) {
        std::cout << "--------------------------------------------" << '\n';
        std::cout << "--- Simulation Setup ---" << '\n';
        std::cout << "h: " << h << ", dT: " << dT << ", n_elements: " << Th.nbElements() << ", n_times_steps: " << num_time_steps << ", T = " << T <<  "\n";
        std::cout << "\n";

        output_data << "h: " << h << ", DX = " << DX << ", DY = " << DY << ", DZ = " << DZ << ", dT: " << dT << ", n_elements: " << Th.nbElements() << ", n_times_steps: " << num_time_steps << ", T = " << T << ", xi = " << xi << ", eta0 = " << eta0 << ", eta1 = " << eta1 << "\n";

        #ifdef USE_CURVATURE
        output_data << "ghost_penalty_curv = " << ghost_penalty_curv << "\n";
        #else 
        output_data << "ghost_penalty_surf_tension = " << ghost_penalty_surf_tension << "\n";
        #endif

        if (reinit_on) {
            output_data << "REINIT: number_iteration: " << reinitialization.number_iteration << ", epsilon_diffusion: " << reinitialization.epsilon_diffusion << ", frequency_reinit: " << frequency_reinit << "\n \n";
        }

        
        output_data << std::right
                    << std::setw(25) << "time step"
                    << std::setw(25) << "time"
                    << std::setw(25) << "mass error"
                    << std::setw(25) << "diff conc" << "\n";

        output_data.flush();
    }

    // --- Geometry and Solver Setup ---
    ImplicitGeometry geometry(Th, q_time, reinitialization);
    geometry.initialize(fun_levelset);

    // Compute initial mass
    InterfaceLevelSet<mesh_t> interface_initial(Th, geometry.levelset(0));
    double mass_0 = integral(ch, interface_initial, 0);

    SolveConcentration concentration(q_time);

    // --- Time Loop ---
    std::vector<double> diff_c(num_time_steps, 0.0);
    std::vector<double> mass_error(num_time_steps, 0.0);
    double max_vel = u_max_init;
    size_t iter = 0, iter_paraview = 0;
    while (iter < num_time_steps) {
        auto t0 = now();

        double current_time   = iter * dT;

        if (MPIcf::IamMaster()) {
            std::cout << "--- Time Step " << iter + 1 << "/" << num_time_steps << ", t = "
                      << std::setw(4) << std::fixed << std::setprecision(3) << current_time + dT << " ---\n";
            std::cout << "---------------------------\n";
        }

        const TimeSlab &In(Ih[iter]);

        // 1. Update geometry (Level Set Advection)
        auto t_start = now();
        

        fct_t vel_shifted(VelV, 0.0);

        if (reinit_on) {
            bool do_reinit = ((iter % frequency_reinit == 0) && (iter > 0));
            if (std::fabs(max_vel) > Epsilon)
                dt_levelset = std::min(h/max_vel, dT) / (n_quad_pts_time - 1);    // use smallest value of advective CFL condition or surface tension CFL condition
            vel_shifted = geometry.update(vel, dt_levelset, do_reinit, reshift_center);
        } else {
            bool do_reinit      = false;
            vel_shifted = geometry.update(vel, dt_levelset, do_reinit, reshift_center);
        }

        if (reshift_center) {
            // 1. Accumulate the total frame velocity into the boundary condition
            // Mathematical equivalent: u_boundary_{new} = u_boundary_{old} - v_avg
            std::transform(
                u_boundary.v.begin(),
                u_boundary.v.end(),
                vel_shifted.v.begin(),
                u_boundary.v.begin(),
                [](double ub, double vs) { return ub - vs; }
            );

            // 2. Shift the internal fluid velocity to match the new frame
            // Mathematical equivalent: v_{new} = v_{old} - v_avg
            std::transform(
                vel.v.begin(),
                vel.v.end(),
                vel_shifted.v.begin(),
                vel.v.begin(),
                [](double u, double vs) { return u - vs; }
            );
        }

        auto t_end = now();

        if (MPIcf::IamMaster()) {
            std::cout << "Time to update geometry: " << seconds(t_start, t_end) << "s with dt_levelset = " << dt_levelset << "\n";
            std::cout << "    advective CFL dt = " << h/max_vel << ", surface tension CFL dt = " << dT << "\n" ;
        }
            

        // 2. Solve time-dependent problem for concentration
        t_start = now();
        concentration.initialize(V, In, geometry.time_interface());

        if (iter == 0)
            concentration.initialize_data(cinitial);
        else
            concentration.initialize_data(); // Initialize with data from the previous step

        concentration.solve(geometry, In, vel, ch);

        // Compute mass and concentration difference on the surface
        const double mass_n = integral(ch, geometry.interface(last_quad_pt_time), 0);
        const double max_c = max_val_surface(ch.expr(), geometry.interface(last_quad_pt_time));
        const double min_c = min_val_surface(ch.expr(), geometry.interface(last_quad_pt_time));

        diff_c[iter] = max_c - min_c;
        mass_error[iter] = std::fabs(mass_n - mass_0);
        t_end = now();

        if (MPIcf::IamMaster())
            std::cout << "Time to solve for concentration: " << seconds(t_start, t_end) << "s\n";


        #ifdef USE_CURVATURE
        // 3. Compute Curvature
        t_start = now();
        Curvature curvature(CurvV, geometry.interface(last_quad_pt_time));
        fct_t Hn = curvature.solve(); // Hn is the mean curvature vector (H*n)
        t_end = now();

        if (MPIcf::IamMaster())
            std::cout << "Time to compute curvature: " << seconds(t_start, t_end) << "\n";

        // 4. Solve Stokes interface problem
        if (MPIcf::IamMaster())
            std::cout << "Solving Stokes...\n";
        t_start = now();

        SolveStokes stokes(VelV, Q, geometry.interface(last_quad_pt_time));
        stokes.info();
        stokes.solve(Hn, ch, c0h, u_boundary);
        #else
        
        // 3. Compute Surface tension
        t_start = now();
        SurfaceTension surface_tension(CurvV, geometry.interface(last_quad_pt_time));

        auto fh = 2.0 * xi * ch.expr() * ch.expr() / (ch.expr() * ch.expr() + c0h.expr() * c0h.expr());
        fct_t Sh = surface_tension.solve(fh); 
        t_end = now();

        if (MPIcf::IamMaster())
            std::cout << "Time to compute surface tension: " << seconds(t_start, t_end) << "\n";

        // 4. Solve Stokes interface problem
        if (MPIcf::IamMaster())
            std::cout << "Solving Stokes...\n";
        t_start = now();

        SolveStokes stokes(VelV, Q, geometry.interface(last_quad_pt_time));
        stokes.info();
        stokes.solve(Sh, u_boundary);
        #endif
        

        stokes.getSolutionBackMesh(vel, geometry.levelset(last_quad_pt_time));

        // update maximal velocity
        const double max_ux = max_val(vel.expr(0), Th);
        const double max_uy = max_val(vel.expr(1), Th);
        const double max_uz = max_val(vel.expr(2), Th);
        max_vel = std::sqrt(max_ux*max_ux + max_uy*max_uy + max_uz*max_uz);

        t_end = now();

        if (MPIcf::IamMaster())
            std::cout << "Total time Stokes solve: " << seconds(t_start, t_end) << "\n";


        // 5. Compute Normal and Tangential Projection
        auto normal_phi = geometry.normal_phi(NormalV, last_quad_pt_time);
        auto vel_gamma = geometry.tangential_projection(vel, normal_phi, last_quad_pt_time);

        if (MPIcf::IamMaster()) {
            std::cout << "Mass conservation error = " << std::setw(4) << std::scientific << std::setprecision(2) << mass_error[iter] << "\n";
            std::cout << "Max - min c = " << diff_c[iter] << "\n";

            output_data << std::fixed << std::setprecision(16) << std::right
                        << std::setw(25) << iter + 1
                        << std::setw(25) << current_time + dT
                        << std::setw(25) << mass_error[iter]
                        << std::setw(25) << diff_c[iter] << '\n';
            output_data.flush();
        }

        auto tn = now();

        if (MPIcf::IamMaster())
            std::cout << "TOTAL TIME FOR FULL TIME STEP: " << std::fixed << std::setprecision(4) << seconds(t0, tn) << "s\n\n";

        // 6. Paraview Output
        if (MPIcf::IamMaster() && (iter % frequency_paraview == 0 || iter == num_time_steps - 1)) {

            // Output solution fields on the background mesh
            paraview_t writer(Th, path_output_figures + "concentration" + std::to_string(iter_paraview + 1) + ".vtk");
            writer.add(ch, "ch", 0, 1);

            if (reshift_center) {
                // --- Calculate Laboratory Frame Velocity ---
                // v_lab = v_comoving - u_boundary = v_comoving + v_avg (since u_boundary was shifted by -v_avg)
                fct_t vel_lab = vel; // Copy the FE space and structure
                std::transform(
                    vel.v.begin(), 
                    vel.v.end(), 
                    u_boundary.v.begin(), 
                    vel_lab.v.begin(), 
                    [](double v_c, double u_b) { return v_c - u_b; }
                );
                writer.add(vel, "vel_comoving", 0, 3); // Wind-tunnel frame
                writer.add(vel_lab, "vel_lab", 0, 3);  // Physical laboratory frame
            } else {
                writer.add(vel, "vel", 0, 3);
            }
            
            // writer.add(2.0 * ch.expr() * ch.expr() / (ch.expr() * ch.expr() + c0h.expr() * c0h.expr()), "fh");
            // writer.add(geometry.levelset(0), "phi_0", 0, 1);
            writer.add(geometry.levelset(last_quad_pt_time), "phi_n", 0, 1);

            #ifdef USE_CURVATURE
            // Output curvature on the surface mesh
            paraview_t para_curv(curvature.mesh(), path_output_figures + "curvature_" + std::to_string(iter_paraview + 1) + ".vtk");
            para_curv.add(Hn, "Hn", 0, 3);
            para_curv.add(geometry.levelset(last_quad_pt_time), "levelset_n", 0, 1);
            #else
            paraview_t para_surf(surface_tension.mesh(), path_output_figures + "surface_tension_" + std::to_string(iter + 1) + ".vtk");
            para_surf.add(Sh, "Sh", 0, 3);
            para_surf.add(geometry.levelset(last_quad_pt_time), "levelset_n", 0, 1);
            #endif
            // // Output velocity on the surface mesh
            // paraview_t para_surface(geometry.get_surface_mesh(), path_output_figures + "surface_" + std::to_string(iter + 1) + ".vtk");
            // para_surface.add(geometry.levelset(last_quad_pt_time), "levelset_n", 0, 1);
            // para_surface.add(vel_gamma, "vel_gamma", 0, 3);

            // // Output normal on the background mesh
            // paraview_t para_bg(Th, path_output_figures + "background_" + std::to_string(iter + 1) + ".vtk");
            // para_bg.add(geometry.levelset(last_quad_pt_time), "phi_n", 0, 1);
            // para_bg.add(normal_phi, "normal_phi_n", 0, 3);

            iter_paraview++;
        }

        iter++;
        std::cout << "---------------------------\n";
    }

    // --- Final Output Summary ---
    if (MPIcf::IamMaster()) {
        std::cout << std::setprecision(16);
        std::cout << "\n--- Final Results ---\n";

        std::cout << "diff_c = [";
        for (size_t i = 0; i < num_time_steps; i++) {
            std::cout << diff_c.at(i) << (i < num_time_steps - 1 ? ", " : "");
        }
        std::cout << "]\n \n";

        std::cout << "mass_error = [";
        for (size_t i = 0; i < num_time_steps; i++) {
            std::cout << mass_error.at(i) << (i < num_time_steps - 1 ? ", " : "");
        }
        std::cout << "]\n";
        std::cout << "---------------------\n";
    }

    return 0;
}
