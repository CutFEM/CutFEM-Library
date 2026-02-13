#ifndef ALGOIM_QUAD_RULE_HPP
#define ALGOIM_QUAD_RULE_HPP

template <typeMesh Mesh>
struct AlgoimQuadratureRule {

    std::vector<typename Mesh::Rd> points;    // Physical coordinates
    std::vector<double> weights;              // Quadrature weights
    std::vector<typename Mesh::Rd> normals;   // Normals (for surfaces)
    
    bool empty() const {return points.empty();}
    size_t size() const {return points.size();}

    void reserve(size_t n) {
        points.reserve(n);
        weights.reserve(n);
        normals.reserve(n);
    }
};

#endif