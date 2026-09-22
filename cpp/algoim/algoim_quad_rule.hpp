#ifndef ALGOIM_QUAD_RULE_HPP
#define ALGOIM_QUAD_RULE_HPP

template <typeMesh Mesh>
struct AlgoimQuadratureRule {

    std::vector<typename Mesh::Rd> points;    // Physical coordinates
    std::vector<double> weights;              // Quadrature weights
    std::vector<typename Mesh::Rd> normals;   // Normals (for surfaces)
    // Correction diagnostics for the Mesh2 implicit-polynomial path.  Counts
    // refer to generator calls (one top-level call may recurse over children).
    int ibp_subdivisions = 0;
    int ibp_capped_failures = 0;
    bool ibp_correction_ok = true;
    double ibp_scaled_residual = 0.0;
    double ibp_acceptance_tolerance = 0.0;
    int ibp_svd_rank = 0;
    int ibp_svd_attempts = 0;
    int ibp_svd_selected = 0;
    double ibp_svd_discarded_rhs = 0.0;
    double ibp_svd_update_norm = 0.0;
    double ibp_relative_vector_change = 0.0;
    double ibp_relative_surface_mass_change = 0.0;
    double ibp_min_normal_alignment = 1.0;
    
    bool empty() const {return points.empty();}
    size_t size() const {return points.size();}

    void reserve(size_t n) {
        points.reserve(n);
        weights.reserve(n);
        normals.reserve(n);
    }
};

#endif
