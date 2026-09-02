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
//  background element it lives on).  Used only when the HOCP map cannot be
//  built (e.g. an unstructured background mesh).
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

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
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
            const int kb = map->locate_element(to_point<Rd>(x));
            phi->setElementFromBackMesh(kb, 0);
            return (*phi)(x);
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

    bool closest_point(const Rd &query, Rd &cp, int &kb) const {
        if (!hocp_)
            return false;

        Vec x = to_uvector(query);
        Vec cp_vec;
        if (!hocp_->compute(x, cp_vec))
            return false;

        for (int a = 0; a < N; ++a)
            cp_vec(a) = std::clamp(cp_vec(a), xmin_(a), xmax_(a));

        cp = to_point<Rd>(cp_vec);
        kb = locate_element(cp);
        return true;
    }

    int seed_count() const { return static_cast<int>(points_.size()); }

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
//  `vel` is the (single-valued) background-mesh velocity; reading it at the
//  closest point cp(x) on Gamma yields the surface velocity that the interface
//  terms are built on.  By default it prefers the HOCP map and falls back to
//  brute-force surface samples.  Pass allow_surface_sample_fallback=false for
//  a strict HOCP-only path that rejects the whole beta field on any failed
//  query.  The backend is reported once (master rank).
// ---------------------------------------------------------------------------
template <typename M, typename L>
void build_extension_velocity(const GFESpace<M> &Vh,
                              const AlgoimInterface<M, L> &gamma,
                              L &phi,
                              const L &vel,
                              L &beta,
                              const bool allow_surface_sample_fallback = true) {
    using Rd = typename M::Rd;

    const HocpClosestPointMap<M, L> hocp_cp(gamma.get_mesh(), phi);
    if (!hocp_cp.ready() && !allow_surface_sample_fallback)
        throw std::runtime_error(
            "HOCP closest-point extension map is unavailable; no fallback allowed");

    // Existing drivers retain the geometric fallback by default. Production
    // drivers that require one unambiguous geometry path can disable it; their
    // temporary beta field is then discarded if any HOCP query fails.
    const std::vector<SurfaceSample<Rd>> samples = allow_surface_sample_fallback
        ? collect_surface_samples(gamma, phi)
        : std::vector<SurfaceSample<Rd>>{};

    static bool reported_closest_point_backend = false;
    if (!reported_closest_point_backend && MPIcf::IamMaster()) {
        std::cout << "closest-point extension backend = "
                  << (hocp_cp.ready() ? "Algoim HOC closest point" : "surface quadrature sample fallback");
        if (hocp_cp.ready())
            std::cout << " (" << hocp_cp.seed_count() << " seeds)";
        if (!allow_surface_sample_fallback)
            std::cout << " [strict: no fallback]";
        std::cout << "\n";
        reported_closest_point_backend = true;
    }

    std::size_t failed_queries = 0;
    auto fun_ext = [&](int /*t*/, std::span<double> P, int comp) -> double {
        Rd Pq;
        for (int a = 0; a < Rd::d; ++a)
            Pq[a] = P[a];
        int kb = -1;
        Rd cp;
        if (!hocp_cp.closest_point(Pq, cp, kb)) {
            ++failed_queries;
            if (!allow_surface_sample_fallback)
                return 0.0; // beta is discarded after the collective check
            if (samples.empty())
                return 0.0;
            cp = closest_point_on_interface(samples, Pq, kb);
            if (kb < 0)
                return 0.0; // no interface in the mesh -> harmless
        }
        return vel.eval(kb, cp, comp, op_id);
    };
    interpolate(Vh, beta.array(), fun_ext);

    if (!allow_surface_sample_fallback) {
        const unsigned long local_failures =
            static_cast<unsigned long>(failed_queries);
        unsigned long global_failures = 0;
        MPIcf::AllReduce(local_failures, global_failures, MPI_MAX);
        if (global_failures != 0)
            throw std::runtime_error(
                "HOCP closest-point extension failed at "
                + std::to_string(global_failures)
                + " local interpolation queries on at least one MPI rank; "
                  "no fallback applied");
    }
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
