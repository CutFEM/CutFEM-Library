#ifdef CUTFEM_USE_EIGEN
#include <Eigen/Sparse>

namespace cutfem::eigen {

using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
using Trip  = Eigen::Triplet<double>;

template <class MapMat>
SpMat to_sparse(int n, const MapMat& Amap) {
    std::vector<Trip> triplets;
    triplets.reserve(static_cast<size_t>(Amap.size()));

    for (const auto& it : Amap) {
        const int i = it.first.first;
        const int j = it.first.second;
        const double aij = static_cast<double>(it.second);
        if (aij != 0.0) triplets.emplace_back(i, j, aij);
    }

    SpMat A(n, n);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    return A;
}

inline Eigen::VectorXd to_vec(const std::vector<double>& b) {
    Eigen::VectorXd v(b.size());
    for (int i = 0; i < (int)b.size(); ++i) v[i] = b[(size_t)i];
    return v;
}

inline void from_vec(const Eigen::VectorXd& x, std::vector<double>& out) {
    out.resize((size_t)x.size());
    for (int i = 0; i < x.size(); ++i) out[(size_t)i] = x[i];
}

} // namespace cutfem::eigen
#endif
