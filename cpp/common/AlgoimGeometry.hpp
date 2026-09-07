#ifndef ALGOIM_GEOMETRY_HPP
#define ALGOIM_GEOMETRY_HPP

// ---------------------------------------------------------------------------
//  Algoim cut-geometry toolkit shared by the active-surface drivers.
//
//  Every algoim driver (axisymmetric / planar 2D quad meshes and the 3D hexa
//  mesh) needs the same handful of geometry utilities built on the algoim
//  high-order cut quadrature:
//
//    * SurfaceSample / collect_surface_samples / closest_point_on_interface
//        - a brute-force closest-point fallback from surface quadrature points.
//    * HocpClosestPointMap
//        - the algoim high-order closest-point (HOCP) map on a structured
//          Cartesian background grid (algoim/hocp.hpp), used for the
//          signed-distance-preserving closest-point extension velocity and for
//          reinitialization.
//    * build_extension_velocity
//        - beta(x) = vel(cp(x)); advecting phi with beta keeps |grad phi| = 1
//          because beta is constant along the interface normals (Neiva 2026).
//    * volume_integral_domain_algoim / integral_stored_surface /
//      stored_surface_measure
//        - diagnostic integrals that reuse the interface's STORED per-node cut
//          rule (setUseStoredInterfaceRule), so measured mass/area/volume match
//          the conserved quantities the solver assembles with.
//
//  These used to be copy-pasted (and drift out of sync) at the top of every
//  driver.  They are collected here dimension-generically: everything is
//  templated on the mesh type `M` and level-set/function type `L`, with the
//  space dimension deduced from `M::Rd::d`.  The 2D and 3D behaviour is
//  bit-for-bit the same as the previous hand-written per-driver copies.
//
//  Usage in a driver (after its own `using mesh_t = MeshQuad2;` etc. aliases):
//
//      #include "common/AlgoimGeometry.hpp"
//      using namespace algoim_geometry;   // brings the helpers + config into scope
//
//  Assumes `cutfem.hpp` has already been included by the translation unit.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "AlgoimInterface.hpp"
#include "../algoim/hocp.hpp"

namespace algoim_geometry {

// ---------------------------------------------------------------------------
//  Quadrature configuration.  Shared defaults for the algoim cut rules; drivers
//  referenced these as file-scope constants, so they keep the same names here.
// ---------------------------------------------------------------------------
inline constexpr int algoim_surface_quad_deg    = 5;
inline constexpr int algoim_volume_quad_deg     = 5;
inline constexpr int hocp_closest_point_degree  = 3;

// ProblemOption preconfigured with the algoim cut-rule degrees.  Pass a
// positive `order_space_element_quadrature` to also raise the uncut element
// quadrature order (used by the MMS drivers).
inline ProblemOption algoim_options(const int order_space_element_quadrature = -1) {
    ProblemOption option = defaultProblemOption;
    if (order_space_element_quadrature > 0)
        option.order_space_element_quadrature_ = order_space_element_quadrature;
    option.algoim_surface_quad_deg_ = algoim_surface_quad_deg;
    option.algoim_vol_quad_deg_     = algoim_volume_quad_deg;
    option.algoim_bernstein_deg_    = std::max(algoim_surface_quad_deg, algoim_volume_quad_deg);
    return option;
}

// ---------------------------------------------------------------------------
//  Small dimension-generic conversions between the library point type `Rd`
//  (R2 / R3) and an algoim uvector<double, N>.
// ---------------------------------------------------------------------------
template <typename Rd>
inline algoim::uvector<double, Rd::d> to_uvector(const Rd &P) {
    algoim::uvector<double, Rd::d> v;
    for (int a = 0; a < Rd::d; ++a)
        v(a) = P[a];
    return v;
}

template <typename Rd, typename UV>
inline Rd to_point(const UV &v) {
    if constexpr (Rd::d == 2)
        return Rd(v(0), v(1));
    else
        return Rd(v(0), v(1), v(2));
}

// ---------------------------------------------------------------------------
//  Brute-force closest-point fallback: sample the surface quadrature points of
//  the interface and, for a query point, return the nearest sample (and the
//  background element it lives on).  Every MPI rank needs the complete sample
//  cloud because extension interpolation is performed on the complete,
//  replicated background space.  Restricting this loop with first_element()/
//  next_element() would give each rank a different fallback geometry.
// ---------------------------------------------------------------------------
template <typename Rd>
struct SurfaceSample {
    Rd point;
    int kb = -1;
};

template <typename M, typename L>
std::vector<SurfaceSample<typename M::Rd>>
collect_surface_samples(const AlgoimInterface<M, L> &gamma, L &phi) {
    using Rd = typename M::Rd;
    std::vector<SurfaceSample<Rd>> samples;
    (void)phi;

    for (int iface = 0; iface < static_cast<int>(gamma.size()); ++iface) {
        const int kb = static_cast<int>(gamma.idxElementOfFace(iface));
        const AlgoimQuadratureRule<M> *quad_rule = gamma.get_cut_quadrature(kb);
        if (quad_rule == nullptr)
            continue;

        samples.reserve(samples.size() + quad_rule->points.size());
        for (size_t ipq = 0; ipq < quad_rule->points.size(); ++ipq)
            samples.push_back({quad_rule->points[ipq], kb});
    }

    return samples;
}

template <typename Rd>
Rd closest_point_on_interface(const std::vector<SurfaceSample<Rd>> &samples,
                              const Rd &P,
                              int &kback) {
    double best = std::numeric_limits<double>::max();
    Rd cp = P;
    kback = -1;

    for (const SurfaceSample<Rd> &sample : samples) {
        const Rd d = P - sample.point;
        const double dist2 = (d, d);
        if (dist2 < best) {
            best  = dist2;
            cp    = sample.point;
            kback = sample.kb;
        }
    }
    return cp;
}

// ---------------------------------------------------------------------------
//  Algoim high-order closest-point (HOCP) map on the structured Cartesian
//  background grid.  Reconstructs the grid geometry from the mesh vertices,
//  builds cell polynomials of the level set, and answers closest-point queries
//  via algoim::ComputeHighOrderCP.  Dimension-generic (N = M::Rd::d).
// ---------------------------------------------------------------------------
template <typename M, typename L>
class HocpClosestPointMap {
    using mesh_t = M;
    using Rd     = typename M::Rd;
    static constexpr int N = Rd::d;

    using Vec      = algoim::uvector<double, N>;
    using IVec     = algoim::uvector<int, N>;
    using Poly     = typename algoim::StencilPoly<N, hocp_closest_point_degree>::T_Poly;
    using CellPoly = algoim::detail::CellPoly<N, Poly>;

    struct GridPhi {
        const HocpClosestPointMap *map = nullptr;
        L *phi = nullptr;

        double operator()(IVec idx) const {
            for (int a = 0; a < N; ++a)
                idx(a) = std::clamp(idx(a), 0, map->extent_(a) - 1);

            const Vec x = map->xmin_ + idx * map->dx_;
            const Rd point = to_point<Rd>(x);
            const int kb = map->locate_element(point);
            const int k = phi->idxElementFromBackMesh(kb, 0);
            return phi->eval(k, point, 0, op_id);
        }
    };

  public:
    HocpClosestPointMap(const mesh_t &Th, L &phi) : Th_(Th) {
        if (!infer_structured_grid())
            return;

        GridPhi grid_phi{this, &phi};
        algoim::detail::createCellPolynomials<N, GridPhi, Poly>(extent_, grid_phi, dx_, false, cells_);
        if (cells_.empty())
            return;

        algoim::detail::samplePolynomials<N, Poly>(cells_, 2, dx_, xmin_, points_, pointcells_);
        if (points_.empty())
            return;

        kdtree_ = std::make_unique<algoim::KDTree<double, N>>(points_);
        double diagonal_sqr = 0.0;
        for (int a = 0; a < N; ++a)
            diagonal_sqr += algoim::util::sqr(xmax_(a) - xmin_(a));
        const double band_radius_sqr = diagonal_sqr;
        const double overlap_radius  = 0.5 * algoim::max(dx_);
        const double cp_tol_sqr =
            algoim::util::sqr(std::max(1.0e-14, std::pow(algoim::max(dx_), Poly::order)));

        hocp_ = std::make_unique<algoim::ComputeHighOrderCP<N, Poly>>(
            band_radius_sqr, overlap_radius, cp_tol_sqr, cells_, *kdtree_, points_, pointcells_, dx_, xmin_);
    }

    bool ready() const { return static_cast<bool>(hocp_); }

    // `status`, when given, receives an algoim::ClosestPointStatus value naming
    // why a query was rejected.  Diagnostic only: the numerical path and the
    // return value are identical whether or not it is requested.
    bool closest_point(const Rd &query, Rd &cp, int &kb, int *status = nullptr) const {
        if (status)
            *status = algoim::cp_ok;
        if (!hocp_) {
            if (status)
                *status = algoim::cp_no_seed_in_band;
            return false;
        }

        Vec x = to_uvector(query);
        Vec cp_vec;
        if (!hocp_->compute(x, cp_vec, nullptr, status))
            return false;

        for (int a = 0; a < N; ++a)
            cp_vec(a) = std::clamp(cp_vec(a), xmin_(a), xmax_(a));

        cp = to_point<Rd>(cp_vec);
        kb = locate_element(cp);
        return true;
    }

    int seed_count() const { return static_cast<int>(points_.size()); }
    int cell_polynomial_count() const { return static_cast<int>(cells_.size()); }

    // ---- Diagnostics for the sampling path the cell polynomials are built from ----
    // createCellPolynomials sees the level set ONLY through GridPhi, i.e. through an
    // element lookup followed by a FunFEM evaluation.  If that evaluation is wrong, every
    // polynomial, seed and closest point is built on a corrupted field while the level-set
    // coefficients themselves look perfect.  These accessors let a driver replay exactly
    // that sampling and compare it against a reference, without duplicating the geometry.
    int grid_extent(int a) const { return extent_(a); }

    Rd grid_point(const std::array<int, N> &idx) const {
        Vec x = xmin_;
        for (int a = 0; a < N; ++a)
            x(a) += std::clamp(idx[a], 0, extent_(a) - 1) * dx_(a);
        return to_point<Rd>(x);
    }

    // The structured grid inferred from the mesh vertices.  Everything the map does is
    // built on these six numbers; if they are wrong, every sampled value, cell polynomial
    // and seed is wrong even though the level-set coefficients are perfect.
    Rd grid_origin() const { return to_point<Rd>(xmin_); }
    Rd grid_upper() const { return to_point<Rd>(xmax_); }
    double grid_spacing(int a) const { return h_(a); }
    double grid_sample_spacing(int a) const { return dx_(a); }
    int grid_cells(int a) const { return ncells_(a); }

    // The element GridPhi pairs with a grid point.  Exposed so a driver can check the
    // structured-numbering assumption (index = ix + nx*iy + nx*ny*iz) against the mesh.
    int grid_element(const std::array<int, N> &idx) const {
        return locate_element(grid_point(idx));
    }

    double sample_grid_levelset(const std::array<int, N> &idx, L &phi) const {
        IVec i;
        for (int a = 0; a < N; ++a)
            i(a) = idx[a];
        GridPhi grid_phi{this, &phi};
        return grid_phi(i);
    }

    // Names for the diagnostic status codes, for driver-side reporting.
    static const char *status_name(int status) {
        switch (status) {
        case algoim::cp_ok: return "ok";
        case algoim::cp_no_seed_in_band: return "no_seed_in_band";
        case algoim::cp_newton_left_ball: return "newton_left_overlap_ball";
        case algoim::cp_newton_max_steps: return "newton_max_steps";
        case algoim::cp_nonfinite: return "nonfinite_closest_point";
        case algoim::cp_residual_above_tolerance: return "residual_above_tolerance";
        default: return "unknown";
        }
    }

  private:
    static std::vector<double> unique_coordinates(std::vector<double> values) {
        std::sort(values.begin(), values.end());
        const double span = values.empty() ? 1.0 : std::max(1.0, values.back() - values.front());
        const double tol = 1.0e-12 * span;

        std::vector<double> unique;
        unique.reserve(values.size());
        for (const double value : values) {
            if (unique.empty() || std::abs(value - unique.back()) > tol)
                unique.push_back(value);
        }
        return unique;
    }

    bool infer_structured_grid() {
        if (Th_.nbVertices() == 0 || Th_.nbElements() == 0)
            return false;

        std::array<std::vector<double>, N> coords;
        for (int a = 0; a < N; ++a)
            coords[a].reserve(static_cast<size_t>(Th_.nbVertices()));

        for (int iv = 0; iv < Th_.nbVertices(); ++iv)
            for (int a = 0; a < N; ++a)
                coords[a].push_back(Th_(iv)[a]);

        int total_cells = 1;
        for (int a = 0; a < N; ++a) {
            coords[a] = unique_coordinates(std::move(coords[a]));
            if (coords[a].size() < 2)
                return false;
            ncells_(a) = static_cast<int>(coords[a].size()) - 1;
            total_cells *= ncells_(a);
        }
        if (total_cells != Th_.nbElements())
            return false;

        for (int a = 0; a < N; ++a) {
            xmin_(a) = coords[a].front();
            xmax_(a) = coords[a].back();
            h_(a)    = (xmax_(a) - xmin_(a)) / static_cast<double>(ncells_(a));
            if (h_(a) <= 0.0)
                return false;
            dx_(a)     = 0.5 * h_(a);
            extent_(a) = 2 * ncells_(a) + 1;
        }
        return true;
    }

    int locate_element(const Rd &point) const {
        int index = 0;
        int stride = 1;
        for (int a = 0; a < N; ++a) {
            int i = static_cast<int>(std::floor((point[a] - xmin_(a)) / h_(a)));
            i = std::clamp(i, 0, ncells_(a) - 1);
            index += i * stride;
            stride *= ncells_(a);
        }
        return index;
    }

    const mesh_t &Th_;
    Vec xmin_;
    Vec xmax_;
    Vec dx_;
    Vec h_;
    IVec extent_;
    IVec ncells_;
    std::vector<CellPoly> cells_;
    std::vector<Vec> points_;
    std::vector<int> pointcells_;
    std::unique_ptr<algoim::KDTree<double, N>> kdtree_;
    std::unique_ptr<algoim::ComputeHighOrderCP<N, Poly>> hocp_;
};

// ---------------------------------------------------------------------------
//  Closest-point extension velocity beta(x) = vel(cp(x)) on the space Vh.
//
//  A closest-point map is only needed in a tube around Gamma.  Asking for it
//  throughout a large box introduces avoidable queries at the medial axis,
//  where a closest point is non-unique.  A finite band therefore uses
//
//      beta = u(cp(x))                                      |phi| <= core*h,
//      beta = w u(cp(x)) + (1-w) u(x)              core*h < |phi| < outer*h,
//      beta = u(x)                                         |phi| >= outer*h,
//
//  with the C2 quintic smootherstep w.  The sampled field is interpolated into
//  the continuous finite-element space Vh, so the final beta_h remains C0
//  across element boundaries.  Failed HOCP queries can either reject the field
//  or use the nearest point in the complete stored Algoim surface-rule cloud;
//  the resilient policy finally uses u(x) if that cloud is unavailable.
// ---------------------------------------------------------------------------
enum class ExtensionFailurePolicy {
    Reject,
    SurfaceSamplesThenRawVelocity
};

struct ExtensionVelocityOptions {
    // Widths are multiples of gamma.get_mesh().get_mesh_size().  Infinity gives
    // the historical whole-domain closest-point extension.
    double pure_extension_band_cells = std::numeric_limits<double>::infinity();
    double blend_outer_band_cells = std::numeric_limits<double>::infinity();
    ExtensionFailurePolicy failure_policy =
        ExtensionFailurePolicy::SurfaceSamplesThenRawVelocity;
};

template <typename Rd>
struct ExtensionVelocityDiagnostics {
    static constexpr int nb_status = algoim::cp_residual_above_tolerance + 1;
    static constexpr int nb_bands = 10;

    bool hocp_map_ready = false;
    int hocp_seed_count = 0;
    int hocp_cell_polynomial_count = 0;
    unsigned long surface_sample_count = 0;
    double mesh_size = 0.0;
    double pure_extension_band_width = std::numeric_limits<double>::infinity();
    double blend_outer_band_width = std::numeric_limits<double>::infinity();

    // Counts are scalar interpolation-point visits (component 0), not unique
    // global degrees of freedom; shared element nodes are intentionally visited
    // more than once by interpolate().
    unsigned long interpolation_point_visits = 0;
    unsigned long pure_extension_point_visits = 0;
    unsigned long blend_point_visits = 0;
    unsigned long raw_velocity_point_visits = 0;
    unsigned long cut_cell_point_visits = 0;
    unsigned long cut_cell_nonpure_point_visits = 0;
    unsigned long hocp_queries = 0;
    unsigned long hocp_successes = 0;
    unsigned long hocp_failures = 0;
    unsigned long surface_sample_fallbacks = 0;
    unsigned long raw_velocity_recoveries = 0;
    unsigned long nonfinite_extension_component_recoveries = 0;
    unsigned long nonfinite_raw_components = 0;
    unsigned long mpi_count_mismatch_fields = 0;

    std::array<unsigned long, nb_status> failures_by_status{};
    std::array<unsigned long, nb_bands> queries_by_band{};
    std::array<unsigned long, nb_bands> failures_by_band{};

    bool have_first_failure = false;
    int first_failure_rank = -1;
    Rd first_failure_point;
    double first_failure_phi = std::numeric_limits<double>::quiet_NaN();
    int first_failure_status = algoim::cp_ok;
    double nearest_failure_abs_phi = std::numeric_limits<double>::infinity();
    double max_surface_fallback_distance = 0.0;
    double max_surface_fallback_sdf_mismatch = 0.0;
    double max_extension_raw_component_delta = 0.0;
    double max_fallback_raw_component_delta = 0.0;
};

inline double extension_blend_weight(const double abs_phi,
                                     const double core_width,
                                     const double outer_width) {
    if (!std::isfinite(outer_width) || abs_phi <= core_width)
        return 1.0;
    if (abs_phi >= outer_width)
        return 0.0;
    const double s = (outer_width - abs_phi) / (outer_width - core_width);
    return s * s * s * (s * (6.0 * s - 15.0) + 10.0);
}

inline int extension_distance_band(const double abs_phi, const double mesh_size,
                                   const int nb_bands) {
    if (!(mesh_size > 0.0) || !std::isfinite(abs_phi))
        return nb_bands - 1;
    return std::clamp(static_cast<int>(std::floor(abs_phi / mesh_size)),
                      0, nb_bands - 1);
}

template <typename M, typename L>
ExtensionVelocityDiagnostics<typename M::Rd>
build_extension_velocity(const GFESpace<M> &Vh,
                         const AlgoimInterface<M, L> &gamma,
                         L &phi,
                         const L &vel,
                         L &beta,
                         const ExtensionVelocityOptions &options) {
    using Rd = typename M::Rd;
    using Diagnostics = ExtensionVelocityDiagnostics<Rd>;

    const bool whole_domain =
        !std::isfinite(options.pure_extension_band_cells)
        && !std::isfinite(options.blend_outer_band_cells);
    if (!whole_domain
        && (!(options.pure_extension_band_cells >= 0.0)
            || !(options.blend_outer_band_cells
                 > options.pure_extension_band_cells)
            || !std::isfinite(options.blend_outer_band_cells))) {
        throw std::invalid_argument(
            "closest-point extension bands require 0 <= core < outer < infinity");
    }

    const HocpClosestPointMap<M, L> hocp_cp(gamma.get_mesh(), phi);
    if (!hocp_cp.ready() && options.failure_policy == ExtensionFailurePolicy::Reject)
        throw std::runtime_error(
            "HOCP closest-point extension map is unavailable; no fallback allowed");

    const bool resilient = options.failure_policy
        == ExtensionFailurePolicy::SurfaceSamplesThenRawVelocity;
    const std::vector<SurfaceSample<Rd>> samples = resilient
        ? collect_surface_samples(gamma, phi)
        : std::vector<SurfaceSample<Rd>>{};

    Diagnostics local;
    local.hocp_map_ready = hocp_cp.ready();
    local.hocp_seed_count = hocp_cp.seed_count();
    local.hocp_cell_polynomial_count = hocp_cp.cell_polynomial_count();
    local.surface_sample_count = static_cast<unsigned long>(samples.size());
    local.mesh_size = gamma.get_mesh().get_mesh_size();
    local.pure_extension_band_width =
        options.pure_extension_band_cells * local.mesh_size;
    local.blend_outer_band_width =
        options.blend_outer_band_cells * local.mesh_size;

    static bool reported_closest_point_backend = false;
    if (!reported_closest_point_backend && MPIcf::IamMaster()) {
        std::cout << "closest-point extension = ";
        if (whole_domain) {
            std::cout << "whole-domain Algoim HOCP";
        } else {
            std::cout << "banded Algoim HOCP (pure through "
                      << options.pure_extension_band_cells
                      << " h, C2 blend to bulk velocity at "
                      << options.blend_outer_band_cells << " h)";
        }
        std::cout << "; recovery = "
                  << (resilient
                      ? "complete stored-surface sample cloud, then raw bulk velocity"
                      : "reject field")
                  << "; map = " << hocp_cp.seed_count() << " seeds; samples = "
                  << samples.size() << '\n';
        reported_closest_point_backend = true;
    }

    enum class ExtensionSource { RawVelocity, Hocp, SurfaceSample };
    Rd cached_query;
    Rd cached_closest_point;
    int cached_closest_element = -1;
    double cached_weight = 0.0;
    ExtensionSource cached_source = ExtensionSource::RawVelocity;

    auto fun_ext = [&](int elem, std::span<double> P, int comp) -> double {
        if (comp == 0) {
            for (int a = 0; a < Rd::d; ++a)
                cached_query[a] = P[a];
            cached_closest_element = -1;
            cached_source = ExtensionSource::RawVelocity;

            const double phi_query = phi.eval(elem, cached_query, 0, op_id);
            const double abs_phi = std::fabs(phi_query);
            cached_weight = extension_blend_weight(
                abs_phi, local.pure_extension_band_width,
                local.blend_outer_band_width);

            ++local.interpolation_point_visits;
            if (cached_weight == 1.0)
                ++local.pure_extension_point_visits;
            else if (cached_weight == 0.0)
                ++local.raw_velocity_point_visits;
            else
                ++local.blend_point_visits;
            if (gamma.isCut(elem)) {
                ++local.cut_cell_point_visits;
                if (cached_weight != 1.0)
                    ++local.cut_cell_nonpure_point_visits;
            }

            if (cached_weight > 0.0) {
                ++local.hocp_queries;
                const int band = extension_distance_band(
                    abs_phi, local.mesh_size, Diagnostics::nb_bands);
                ++local.queries_by_band[band];

                int status = algoim::cp_ok;
                if (hocp_cp.closest_point(cached_query, cached_closest_point,
                                          cached_closest_element, &status)) {
                    cached_source = ExtensionSource::Hocp;
                    ++local.hocp_successes;
                } else {
                    ++local.hocp_failures;
                    ++local.failures_by_status[std::clamp(
                        status, 0, Diagnostics::nb_status - 1)];
                    ++local.failures_by_band[band];
                    local.nearest_failure_abs_phi = std::min(
                        local.nearest_failure_abs_phi, abs_phi);
                    if (!local.have_first_failure) {
                        local.have_first_failure = true;
                        local.first_failure_point = cached_query;
                        local.first_failure_phi = phi_query;
                        local.first_failure_status = status;
                    }

                    if (resilient && !samples.empty()) {
                        cached_closest_point = closest_point_on_interface(
                            samples, cached_query, cached_closest_element);
                        if (cached_closest_element >= 0) {
                            cached_source = ExtensionSource::SurfaceSample;
                            ++local.surface_sample_fallbacks;
                            const Rd displacement =
                                cached_query - cached_closest_point;
                            const double distance = std::sqrt(
                                (displacement, displacement));
                            local.max_surface_fallback_distance = std::max(
                                local.max_surface_fallback_distance, distance);
                            local.max_surface_fallback_sdf_mismatch = std::max(
                                local.max_surface_fallback_sdf_mismatch,
                                std::fabs(distance - abs_phi));
                        }
                    }
                    if (cached_source == ExtensionSource::RawVelocity)
                        ++local.raw_velocity_recoveries;
                }
            }
        }

        const double raw_value = vel.eval(elem, cached_query, comp, op_id);
        if (!std::isfinite(raw_value))
            ++local.nonfinite_raw_components;
        if (cached_source == ExtensionSource::RawVelocity)
            return raw_value;

        const double extension_value = vel.eval(
            cached_closest_element, cached_closest_point, comp, op_id);
        if (!std::isfinite(extension_value)) {
            ++local.nonfinite_extension_component_recoveries;
            return raw_value;
        }

        const double delta = std::fabs(extension_value - raw_value);
        local.max_extension_raw_component_delta = std::max(
            local.max_extension_raw_component_delta, delta);
        if (cached_source == ExtensionSource::SurfaceSample)
            local.max_fallback_raw_component_delta = std::max(
                local.max_fallback_raw_component_delta, delta);
        return cached_weight * extension_value
             + (1.0 - cached_weight) * raw_value;
    };
    interpolate(Vh, beta.array(), fun_ext);

    Diagnostics global = local;
    const auto reduce_max = [](const unsigned long local_value) {
        unsigned long global_value = 0;
        MPIcf::AllReduce(local_value, global_value, MPI_MAX);
        return global_value;
    };
    const auto reduce_max_double = [](const double local_value) {
        double global_value = 0.0;
        MPIcf::AllReduce(local_value, global_value, MPI_MAX);
        return global_value;
    };

    int map_ready = local.hocp_map_ready ? 1 : 0;
    int global_map_ready = 0;
    MPIcf::AllReduce(map_ready, global_map_ready, MPI_MIN);
    global.hocp_map_ready = global_map_ready != 0;
    MPIcf::AllReduce(local.hocp_seed_count, global.hocp_seed_count, MPI_MAX);
    MPIcf::AllReduce(local.hocp_cell_polynomial_count,
                     global.hocp_cell_polynomial_count, MPI_MAX);
    global.surface_sample_count = reduce_max(local.surface_sample_count);
    global.interpolation_point_visits = reduce_max(local.interpolation_point_visits);
    global.pure_extension_point_visits = reduce_max(local.pure_extension_point_visits);
    global.blend_point_visits = reduce_max(local.blend_point_visits);
    global.raw_velocity_point_visits = reduce_max(local.raw_velocity_point_visits);
    global.cut_cell_point_visits = reduce_max(local.cut_cell_point_visits);
    global.cut_cell_nonpure_point_visits = reduce_max(
        local.cut_cell_nonpure_point_visits);
    global.hocp_queries = reduce_max(local.hocp_queries);
    global.hocp_successes = reduce_max(local.hocp_successes);
    global.hocp_failures = reduce_max(local.hocp_failures);
    global.surface_sample_fallbacks = reduce_max(local.surface_sample_fallbacks);
    global.raw_velocity_recoveries = reduce_max(local.raw_velocity_recoveries);
    global.nonfinite_extension_component_recoveries = reduce_max(
        local.nonfinite_extension_component_recoveries);
    global.nonfinite_raw_components = reduce_max(local.nonfinite_raw_components);
    for (int i = 0; i < Diagnostics::nb_status; ++i)
        global.failures_by_status[i] = reduce_max(local.failures_by_status[i]);
    for (int i = 0; i < Diagnostics::nb_bands; ++i) {
        global.queries_by_band[i] = reduce_max(local.queries_by_band[i]);
        global.failures_by_band[i] = reduce_max(local.failures_by_band[i]);
    }
    double nearest_failure = local.nearest_failure_abs_phi;
    MPIcf::AllReduce(nearest_failure, global.nearest_failure_abs_phi, MPI_MIN);
    global.max_surface_fallback_distance = reduce_max_double(
        local.max_surface_fallback_distance);
    global.max_surface_fallback_sdf_mismatch = reduce_max_double(
        local.max_surface_fallback_sdf_mismatch);
    global.max_extension_raw_component_delta = reduce_max_double(
        local.max_extension_raw_component_delta);
    global.max_fallback_raw_component_delta = reduce_max_double(
        local.max_fallback_raw_component_delta);

    // Geometry and interpolation are replicated, so these counters must agree
    // across ranks.  Keep the maxima above for the forensic values and record
    // how many key fields disagree instead of silently assuming replication.
    const auto differs_across_ranks = [](const unsigned long local_value,
                                         const unsigned long global_maximum) {
        unsigned long global_minimum = 0;
        MPIcf::AllReduce(local_value, global_minimum, MPI_MIN);
        return global_minimum != global_maximum;
    };
    global.mpi_count_mismatch_fields += differs_across_ranks(
        static_cast<unsigned long>(local.hocp_seed_count),
        static_cast<unsigned long>(global.hocp_seed_count));
    global.mpi_count_mismatch_fields += differs_across_ranks(
        static_cast<unsigned long>(local.hocp_cell_polynomial_count),
        static_cast<unsigned long>(global.hocp_cell_polynomial_count));
    global.mpi_count_mismatch_fields += differs_across_ranks(
        local.surface_sample_count, global.surface_sample_count);
    global.mpi_count_mismatch_fields += differs_across_ranks(
        local.interpolation_point_visits, global.interpolation_point_visits);
    global.mpi_count_mismatch_fields += differs_across_ranks(
        local.pure_extension_point_visits, global.pure_extension_point_visits);
    global.mpi_count_mismatch_fields += differs_across_ranks(
        local.blend_point_visits, global.blend_point_visits);
    global.mpi_count_mismatch_fields += differs_across_ranks(
        local.raw_velocity_point_visits, global.raw_velocity_point_visits);
    global.mpi_count_mismatch_fields += differs_across_ranks(
        local.cut_cell_nonpure_point_visits,
        global.cut_cell_nonpure_point_visits);
    global.mpi_count_mismatch_fields += differs_across_ranks(
        local.hocp_queries, global.hocp_queries);
    global.mpi_count_mismatch_fields += differs_across_ranks(
        local.hocp_failures, global.hocp_failures);
    global.mpi_count_mismatch_fields += differs_across_ranks(
        local.surface_sample_fallbacks, global.surface_sample_fallbacks);
    global.mpi_count_mismatch_fields += differs_across_ranks(
        local.raw_velocity_recoveries, global.raw_velocity_recoveries);

    int first_failure_rank = local.have_first_failure
        ? MPIcf::my_rank() : MPIcf::size();
    MPIcf::AllReduce(first_failure_rank, global.first_failure_rank, MPI_MIN);
    global.have_first_failure = global.first_failure_rank < MPIcf::size();
    if (global.have_first_failure) {
        if (MPIcf::my_rank() != global.first_failure_rank) {
            global.first_failure_status = algoim::cp_ok;
            global.first_failure_phi = std::numeric_limits<double>::quiet_NaN();
        }
        MPIcf::Bcast(global.first_failure_status, global.first_failure_rank, 1);
        MPIcf::Bcast(global.first_failure_phi, global.first_failure_rank, 1);
        for (int a = 0; a < Rd::d; ++a)
            MPIcf::Bcast(global.first_failure_point[a],
                         global.first_failure_rank, 1);
    }

    if (options.failure_policy == ExtensionFailurePolicy::Reject
        && global.hocp_failures != 0) {
        throw std::runtime_error(
            "HOCP closest-point extension failed at "
            + std::to_string(global.hocp_failures)
            + " scalar interpolation-point visits on at least one MPI rank; "
              "no fallback applied");
    }
    return global;
}

// Backward-compatible whole-domain entry point used by the axisymmetric
// driver.  `false` preserves the strict behavior of the former boolean API.
template <typename M, typename L>
ExtensionVelocityDiagnostics<typename M::Rd>
build_extension_velocity(const GFESpace<M> &Vh,
                         const AlgoimInterface<M, L> &gamma,
                         L &phi,
                         const L &vel,
                         L &beta,
                         const bool allow_surface_sample_fallback = true) {
    ExtensionVelocityOptions options;
    options.failure_policy = allow_surface_sample_fallback
        ? ExtensionFailurePolicy::SurfaceSamplesThenRawVelocity
        : ExtensionFailurePolicy::Reject;
    return build_extension_velocity(Vh, gamma, phi, vel, beta, options);
}

// ---------------------------------------------------------------------------
//  Volume integral over one phase using freshly generated algoim volume rules.
//  Two overloads: an expression-weighted integral and a pure measure.
// ---------------------------------------------------------------------------
template <typename M, typename L>
double volume_integral_domain_algoim(const FunFEM<M> &integrand,
                                     const M &Th,
                                     const AlgoimInterface<M, L> &interface,
                                     L &phi,
                                     const int domain,
                                     ProblemOption option = algoim_options()) {
    ActiveMesh<M> active_mesh(Th, interface);

    double local_value = 0.0;
    for (int k = active_mesh.first_element(); k < active_mesh.last_element(); k += active_mesh.next_element()) {
        if (active_mesh.get_domain_element(k) != domain || active_mesh.isInactive(k, 0))
            continue;

        const int kb = active_mesh.idxElementInBackMesh(k);
        phi.setElementFromBackMesh(kb, 0);
        const auto quad_rule = quadGenVol(active_mesh[k], phi, option, domain);

        for (size_t ipq = 0; ipq < quad_rule.points.size(); ++ipq)
            local_value += quad_rule.weights[ipq] * integrand.eval(kb, quad_rule.points[ipq], 0, op_id);
    }

    double global_value = 0.0;
    MPIcf::AllReduce(local_value, global_value, MPI_SUM);
    return global_value;
}

template <typename M, typename L>
double volume_integral_domain_algoim(const M &Th,
                                     const AlgoimInterface<M, L> &interface,
                                     L &phi,
                                     const int domain,
                                     ProblemOption option = algoim_options()) {
    ActiveMesh<M> active_mesh(Th, interface);

    double local_value = 0.0;
    for (int k = active_mesh.first_element(); k < active_mesh.last_element(); k += active_mesh.next_element()) {
        if (active_mesh.get_domain_element(k) != domain || active_mesh.isInactive(k, 0))
            continue;

        const int kb = active_mesh.idxElementInBackMesh(k);
        phi.setElementFromBackMesh(kb, 0);
        const auto quad_rule = quadGenVol(active_mesh[k], phi, option, domain);

        for (const double weight : quad_rule.weights)
            local_value += weight;
    }

    double global_value = 0.0;
    MPIcf::AllReduce(local_value, global_value, MPI_SUM);
    return global_value;
}

// ---------------------------------------------------------------------------
//  Surface integrals over an interface using its STORED per-node algoim cut
//  rule (setUseStoredInterfaceRule).  `integral_stored_surface` weights an
//  expression (t forwarded for space-time fields); `stored_surface_measure`
//  is the plain area/length.
// ---------------------------------------------------------------------------
template <typename M, typename L>
double integral_stored_surface(const std::shared_ptr<const ExpressionVirtual> &fh,
                               const AlgoimInterface<M, L> &interface,
                               const double t = 0.0) {
    double local_value = 0.0;
    for (int iface = interface.first_element(); iface < interface.last_element();
         iface += interface.next_element()) {
        const int kb = static_cast<int>(interface.idxElementOfFace(iface));
        const AlgoimQuadratureRule<M> *rule = interface.get_cut_quadrature(kb);
        if (rule == nullptr)
            continue;
        for (size_t ipq = 0; ipq < rule->points.size(); ++ipq)
            local_value += rule->weights[ipq]
                         * fh->evalOnBackMesh(kb, 0, rule->points[ipq], t, rule->normals[ipq]);
    }
    double global_value = 0.0;
    MPIcf::AllReduce(local_value, global_value, MPI_SUM);
    return global_value;
}

template <typename M, typename L>
double stored_surface_measure(const AlgoimInterface<M, L> &interface) {
    double local_value = 0.0;
    for (int iface = interface.first_element(); iface < interface.last_element();
         iface += interface.next_element()) {
        const int kb = static_cast<int>(interface.idxElementOfFace(iface));
        const AlgoimQuadratureRule<M> *rule = interface.get_cut_quadrature(kb);
        if (rule == nullptr)
            continue;
        for (const double weight : rule->weights)
            local_value += weight;
    }
    double global_value = 0.0;
    MPIcf::AllReduce(local_value, global_value, MPI_SUM);
    return global_value;
}

} // namespace algoim_geometry

#endif // ALGOIM_GEOMETRY_HPP
