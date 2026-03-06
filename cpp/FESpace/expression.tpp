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
template <typename M>
FunFEM<M>::FunFEM(const FESpace &vh, const ExpressionVirtual &fh) : FunFEMVirtual(vh.NbDoF()), Vh(&vh) {

    // allocate memory for basis functions
    std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
    std::size_t n_chunk     = cutfem_get_max_threads();
    pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);

    assert(Vh->N == 1);

// #ifdef USE_MPI
//     double dataSend[Vh->nbDoF];
//     Rn_ fhSend(dataSend, Vh->nbDoF);
//     fhSend = 1e+50;
// #else

// #endif

    const int d   = Vh->N;
    const int nve = Vh->TFE(0)->NbPtforInterpolation;
    // KNM<R> Vpf(nve, d);                             // value of f at the interpolation points

// #pragma omp parallel
//     {
        std::vector<std::vector<double>> Vpf(d, std::vector<double>(nve)); // value of f at the interpolation points
        std::vector<double> ggf(Vh->MaxNbDFPerElement); // stock the values of the dof of the interpolate

// #pragma omp for
        for (int k = Vh->first_element(); k < Vh->last_element(); k += Vh->next_element()) {

            const FElement &FK((*Vh)[k]);
            const int nbdf   = FK.NbDoF(); // nof local
            const int domain = FK.get_domain();
            const int kb     = Vh->idxElementInBackMesh(k);

            for (int p = 0; p < FK.tfe->NbPtforInterpolation; p++) { // all interpolation points
                const Rd &P(FK.Pt(p));                               // the coordinate of P in K hat
                for (int i = 0; i < d; ++i) {
                    Vpf[i][p] = fh.evalOnBackMesh(kb, domain, P);
                }
            }
            FK.Pi_h(Vpf, ggf);

// #ifdef USE_MPI
// #pragma omp critical
//             for (int df = 0; df < nbdf; df++) {
//                 fhSend(FK.loc2glb(df)) = ggf[df];
//                 // fh[K(df)] =  ggf[df] ;
//             }
// #else
// #pragma omp critical
            for (int df = 0; df < nbdf; df++) {
                v[FK.loc2glb(df)] = ggf[df];
            }
// #endif
        }
    // }

// #ifdef USE_MPI
//     MPIcf::AllReduce(dataSend, v.data(), fhSend.size(), MPI_MIN);
// #endif
}

template <typename M>
FunFEM<M>::FunFEM(const FESpace &vh, const ExpressionVirtual &fh1, const ExpressionVirtual &fh2)
    : FunFEMVirtual(vh.NbDoF()), Vh(&vh) {

    // allocate memory for basis functions
    std::size_t databf_size = vh[0].NbDoF() * vh.N * 10;
    std::size_t n_chunk     = cutfem_get_max_threads();
    pool_databf             = std::make_unique<MemoryPool>(n_chunk, databf_size);

    assert(Vh->N == 2);
    double dataSend[Vh->nbDoF];
    Rn_ fhSend(dataSend, Vh->nbDoF);
    fhSend        = 1e+50;
    const int d   = Vh->N;
    const int nve = Vh->TFE(0)->NbPtforInterpolation;
    KNM<R> Vpf(nve, d);               // value of f at the interpolation points
    KN<R> ggf(Vh->MaxNbDFPerElement); // stock the values of the dof of the
                                      // interpolate

    for (int k = Vh->first_element(); k < Vh->last_element(); k += Vh->next_element()) {

        const FElement &FK((*Vh)[k]);
        const int nbdf   = FK.NbDoF(); // nof local
        const int domain = FK.whichDomain();
        const int kb     = Vh->idxElementInBackMesh(k);

        for (int p = 0; p < FK.tfe->NbPtforInterpolation; p++) { // all interpolation points
            const Rd &P(FK.Pt(p));                               // the coordinate of P in K hat
            for (int i = 0; i < d; ++i) {
                const ExpressionVirtual &fh = (d == 0) ? fh1 : fh2;
                Vpf(p, i)                   = fh.evalOnBackMesh(kb, domain, P);
            }
        }
        // std::cout << Vpf << std::endl;
        FK.Pi_h(Vpf, ggf);
        for (int df = 0; df < nbdf; df++) {
            fhSend(FK.loc2glb(df)) = ggf[df];
            // fh[K(df)] =  ggf[df] ;
        }
        // for(int i=0, j=FK.dfcbegin(ci);j<FK.dfcend(ci);++j,++i) {
        //   Rd mip = FK.Pt(i);
        //   fhSend(FK.loc2glb(j)) = fh.evalOnBackMesh(kb, domain, mip);
        // }
    }
#ifdef USE_MPI
    MPIcf::AllReduce(dataSend, v, fhSend.size(), MPI_MIN);
#else
    assert(0 && "need to fixe the output");
#endif
}

template <typename M> void FunFEM<M>::print() const { std::cout << v << std::endl; }

template <typename M> double FunFEM<M>::eval(const int k, const R *x, int cu, int op) const {
    const FElement &FK((*Vh)[k]);
    int ndf = FK.NbDoF();

    std::size_t i_thread         = omp_get_thread_num();
    std::span<double> mem_buffer = pool_databf->get_data(i_thread);
    RNMK_ w(mem_buffer.data(), ndf, Vh->N, op_dz + 1);
    FK.BF(FK.T.toKref(x), w);

    double val = 0.;

    for (int j = FK.dfcbegin(cu); j < FK.dfcend(cu); ++j) {
        val += v[FK(j)] * w(j, cu, op);
    }

    return val;
}

template <typename M> double FunFEM<M>::eval(const int k, const R *x, const R t, int cu, int op, int opt) const {

    if (!In)
        return eval(k, x, cu, op);

    const FElement &FK((*Vh)[k]);
    int ndf = FK.NbDoF();

    std::size_t i_thread         = omp_get_thread_num();
    std::span<double> mem_buffer = pool_databf->get_data(i_thread);

    RNMK_ w(mem_buffer.data(), ndf, Vh->N, op_dz + 1);
    KNMK<R> wt(In->NbDoF(), 1, op_dz);

    FK.BF(FK.T.toKref(x), w);
    In->BF(In->T.toKref(t), wt);

    double val = 0.;
    for (int jt = 0; jt < In->NbDoF(); ++jt) {
        for (int j = FK.dfcbegin(cu); j < FK.dfcend(cu); ++j) {
            val += v[FK(j) + jt * Vh->NbDoF()] * w(j, cu, op) * wt(jt, 0, opt);
        }
    }

    return val;
}

template <typename M> void FunFEM<M>::eval(R *u, const int k) const {
    assert(v.data() && u);
    const FElement &FK((*Vh)[k]);
    for (int ci = 0; ci < Vh->N; ++ci) { // loop over componant
        for (int j = FK.dfcbegin(ci); j < FK.dfcend(ci); ++j)
            u[j] = v[FK(j)];
    }
}

template <typename M> double FunFEM<M>::evalOnBackMesh(const int kb, int dom, const R *x, int cu, int op) const {
    int k = Vh->idxElementFromBackMesh(kb, dom);

    return eval(k, x, cu, op);
}

template <typename M>
double FunFEM<M>::evalOnBackMesh(const int kb, int dom, const R *x, const R t, int cu, int op, int opt) const {

    int k = Vh->idxElementFromBackMesh(kb, dom);

    return eval(k, x, t, cu, op, opt);
}

template <typename M> std::vector<std::shared_ptr<ExpressionFunFEM<M>>> FunFEM<M>::exprList(int n) const {
    if (n == -1)
        n = Vh->N;
    assert(n <= Vh->N);
    std::vector<std::shared_ptr<ExpressionFunFEM<Mesh>>> l;
    for (int i = 0; i < n; ++i) {
        l.push_back(std::make_shared<ExpressionFunFEM<Mesh>>(*this, i, op_id));
    }
    return l;
}

template <typename M> std::shared_ptr<ExpressionFunFEM<M>> FunFEM<M>::expr(int i0) const {
    assert(i0 < Vh->N);
    return std::make_shared<ExpressionFunFEM<Mesh>>(*this, i0, op_id);
}

template <typename M> std::vector<std::shared_ptr<ExpressionFunFEM<M>>> FunFEM<M>::exprList(int n, int i0) const {
    assert(n <= Vh->N);
    std::vector<std::shared_ptr<ExpressionFunFEM<Mesh>>> l;
    for (int i = 0; i < n; ++i) {
        l.push_back(std::make_shared<ExpressionFunFEM<Mesh>>(*this, i + i0, op_id));
    }
    return l;
}

template<typename M>
FunFEM<M>& FunFEM<M>::operator+=(const FunFEM<M>& other) {
    assert(v.size() == other.v.size() && "Size mismatch in operator+=");
    for (std::size_t i = 0; i < v.size(); ++i) {
        v[i] += other.v[i];
    }
    return *this;
}

template<typename M>
FunFEM<M>& FunFEM<M>::operator-=(const FunFEM<M>& other) {
    assert(v.size() == other.v.size() && "Size mismatch in operator+=");
    for (std::size_t i = 0; i < v.size(); ++i) {
        v[i] -= other.v[i];
    }
    return *this;
}

template<typename M>
FunFEM<M>& FunFEM<M>::operator*=(const double c) {
    for (std::size_t i = 0; i < v.size(); ++i) {
        v[i] *= c;
    }
    return *this;
}

template<typename M>
FunFEM<M>& FunFEM<M>::operator/=(const double c) {
    assert(c != 0.);
    for (std::size_t i = 0; i < v.size(); ++i) {
        v[i] /= c;
    }
    return *this;
}



// Algoim additions

// Compute barycentric coordinates of a physical point on a simplex.
// lam[0..nv-1] are filled with Interval values.
// Dl[i] = grad(lambda_i) from Gradlambda — pure doubles.
// A[j] = j-th coordinate of K[0] — pure double.
// Template params Itv, nv, d are all explicit to avoid deduction issues.
template<typename Itv, int nv, int d, typename RdType>
inline void computeBarycentricCoords(const algoim::uvector<Itv, d>& x, const RdType* Dl,
                                     const RdType& A, Itv* lam) {
    lam[0] = 1.0;
    for (int i = 1; i < nv; ++i) {
        lam[i] = 0.0;
        for (int j = 0; j < d; ++j)
            lam[i] += Dl[i][j] * (x(j) - A[j]);
        lam[0] -= lam[i];
    }
}

// operator()(Interval)
template<typename M>
algoim::Interval<M::Rd::d>
FunFEM<M>::operator()(const algoim::uvector<algoim::Interval<M::Rd::d>, M::Rd::d>& x) const {
    using Itv = algoim::Interval<M::Rd::d>;
    constexpr int d = M::Rd::d;

    assert(algoim_k_ >= 0);
    const FElement& FK = (*Vh)[algoim_k_];
    const Element&  K  = FK.T;

    switch (getBasisFctType()) {

    case BasisFctType::P1:
    case BasisFctType::P1dc: {
        if constexpr (std::is_same_v<Mesh, MeshQuad2>) {
            // Q1 bilinear on quad
            double x0 = K[0][0], y0 = K[0][1];
            double hx  = K[1][0] - x0;
            double hy  = K[3][1] - y0;
            Itv s = (x(0) - x0) / hx;
            Itv t = (x(1) - y0) / hy;
            Itv val = v[FK(0)]*(1.0-s)*(1.0-t) + v[FK(1)]*s*(1.0-t)
                    + v[FK(2)]*s*t             + v[FK(3)]*(1.0-s)*t;
            val.eps = 0.0;  // Polynomial evaluation is exact
            return val;
        } else {
            // P1 Lagrange on simplex (Mesh2 or Mesh3)
            Rd Dl[Element::nv];
            K.Gradlambda(Dl);
            const Rd& A = static_cast<const Rd&>(K[0]);
            Itv lam[Element::nv];
            computeBarycentricCoords<Itv, Element::nv, d>(x, Dl, A, lam);
            Itv val(0.0);
            for (int i = 0; i < Element::nv; ++i)
                val += v[FK(i)] * lam[i];
            val.eps = 0.0;  // Polynomial evaluation is exact
            return val;
        }
    }

    case BasisFctType::P2:
    case BasisFctType::P2dc: {
        if constexpr (std::is_same_v<Mesh, Mesh2>) {
            // P2 Lagrange on triangle: 3 vertex DOFs + 3 edge DOFs
            Rd Dl[3];
            K.Gradlambda(Dl);
            const Rd& A = static_cast<const Rd&>(K[0]);
            Itv l[3];
            computeBarycentricCoords<Itv, 3, d>(x, Dl, A, l);
            // vertex basis functions: phi_i = l[i] * (2*l[i] - 1)
            Itv val(0.0);
            for (int i = 0; i < 3; ++i)
                val += v[FK(i)] * l[i] * (2.0*l[i] - 1.0);
            // edge basis functions: phi = 4 * l[a] * l[b]
            for (int i = 0; i < 3; ++i) {
                int i0 = Element::nvedge[i][0], i1 = Element::nvedge[i][1];
                val += v[FK(3+i)] * 4.0 * l[i0] * l[i1];
            }
            val.eps = 0.0;  // Polynomial evaluation is exact
            return val;
        } else if constexpr (std::is_same_v<Mesh, MeshQuad2>) {
            // Q2 tensor-product on quad: 9 DOFs
            double x0 = K[0][0], y0 = K[0][1];
            double hx  = K[1][0] - x0;
            double hy  = K[3][1] - y0;
            Itv s = (x(0) - x0) / hx;
            Itv t = (x(1) - y0) / hy;
            Itv psi_x[3] = {1.0-3.0*s+2.0*s*s, s*(2.0*s-1.0), 4.0*s*(1.0-s)};
            Itv psi_y[3] = {1.0-3.0*t+2.0*t*t, t*(2.0*t-1.0), 4.0*t*(1.0-t)};
            const int ix[9] = {0,1,1,0,2,1,2,0,2};
            const int iy[9] = {0,0,1,1,0,2,1,2,2};
            Itv val(0.0);
            for (int i = 0; i < 9; ++i)
                val += v[FK(i)] * psi_x[ix[i]] * psi_y[iy[i]];
            // For polynomial evaluation on reference element, result is exact
            // Zero out accumulated eps from interval arithmetic
            val.eps = 0.0;
            return val;
        } else {
            assert(0 && "P2 Interval eval not yet implemented for this mesh type");
            return {};
        }
    }

    default:
        assert(0 && "Interval eval not implemented for this FE type");
        return {};
    }
}

// grad(Interval)
template<typename M>
algoim::uvector<algoim::Interval<M::Rd::d>, M::Rd::d>
FunFEM<M>::grad(const algoim::uvector<algoim::Interval<M::Rd::d>, M::Rd::d>& x) const {
    using Itv = algoim::Interval<M::Rd::d>;
    constexpr int d = M::Rd::d;

    assert(algoim_k_ >= 0);
    const FElement& FK = (*Vh)[algoim_k_];
    const Element&  K  = FK.T;
    algoim::uvector<Itv, d> g;  // zero-initialised by uvector default ctor

    switch (getBasisFctType()) {

    case BasisFctType::P1:
    case BasisFctType::P1dc: {
        if constexpr (std::is_same_v<Mesh, MeshQuad2>) {
            // Q1: gradient is bilinear — depends on x
            double x0 = K[0][0], y0 = K[0][1];
            double hx  = K[1][0] - x0;
            double hy  = K[3][1] - y0;
            Itv s = (x(0) - x0) / hx;
            Itv t = (x(1) - y0) / hy;
            double u0=v[FK(0)], u1=v[FK(1)], u2=v[FK(2)], u3=v[FK(3)];
            g(0) = ((u1-u0)*(1.0-t) + (u2-u3)*t) / hx;
            g(1) = ((u3-u0)*(1.0-s) + (u2-u1)*s) / hy;
            g(0).eps = 0.0;  // Polynomial gradient is exact
            g(1).eps = 0.0;
        } else {
            // P1 on simplex: gradient is element-constant (all Dl[i] are pure doubles)
            Rd Dl[Element::nv];
            K.Gradlambda(Dl);
            for (int j = 0; j < d; ++j) {
                g(j) = 0.0;
                for (int i = 0; i < Element::nv; ++i)
                    g(j) += v[FK(i)] * Dl[i][j];
                g(j).eps = 0.0;  // Polynomial gradient is exact
            }
        }
        break;
    }

    case BasisFctType::P2:
    case BasisFctType::P2dc: {
        if constexpr (std::is_same_v<Mesh, Mesh2>) {
            // P2 on triangle: gradient is linear in barycentric coords → depends on x
            Rd Dl[3];
            K.Gradlambda(Dl);
            const Rd& A = static_cast<const Rd&>(K[0]);
            Itv l[3];
            computeBarycentricCoords<Itv, 3, d>(x, Dl, A, l);
            // vertex DOFs: grad(phi_i) = (4*l[i] - 1) * Dl[i] * u_i
            for (int i = 0; i < 3; ++i) {
                Itv coeff = (4.0*l[i] - 1.0) * v[FK(i)];
                for (int j = 0; j < d; ++j)
                    g(j) += coeff * Dl[i][j];
            }
            // edge DOFs: grad(phi) = 4 * (l[b]*Dl[a] + l[a]*Dl[b]) * u_edge
            for (int i = 0; i < 3; ++i) {
                int i0 = Element::nvedge[i][0], i1 = Element::nvedge[i][1];
                double u_edge = v[FK(3+i)];
                for (int j = 0; j < d; ++j)
                    g(j) += 4.0 * u_edge * (l[i1]*Dl[i0][j] + l[i0]*Dl[i1][j]);
            }
            for (int j = 0; j < d; ++j)
                g(j).eps = 0.0;  // Polynomial gradient is exact
        } else if constexpr (std::is_same_v<Mesh, MeshQuad2>) {
            // Q2 tensor-product on quad: 9 DOFs
            double x0 = K[0][0], y0 = K[0][1];
            double hx  = K[1][0] - x0;
            double hy  = K[3][1] - y0;
            Itv s = (x(0) - x0) / hx;
            Itv t = (x(1) - y0) / hy;
            Itv psi_x[3]   = {1.0-3.0*s+2.0*s*s, s*(2.0*s-1.0), 4.0*s*(1.0-s)};
            Itv psi_y[3]   = {1.0-3.0*t+2.0*t*t, t*(2.0*t-1.0), 4.0*t*(1.0-t)};
            Itv dpsi_x[3]  = {(-3.0+4.0*s)/hx, (4.0*s-1.0)/hx, (4.0-8.0*s)/hx};
            Itv dpsi_y[3]  = {(-3.0+4.0*t)/hy, (4.0*t-1.0)/hy, (4.0-8.0*t)/hy};
            const int ix[9] = {0,1,1,0,2,1,2,0,2};
            const int iy[9] = {0,0,1,1,0,2,1,2,2};
            g(0) = 0.0; g(1) = 0.0;
            for (int i = 0; i < 9; ++i) {
                double ui = v[FK(i)];
                g(0) += ui * dpsi_x[ix[i]] * psi_y[iy[i]];
                g(1) += ui * psi_x[ix[i]]  * dpsi_y[iy[i]];
            }
            g(0).eps = 0.0;  // Polynomial gradient is exact
            g(1).eps = 0.0;
        } else {
            assert(0 && "P2 Interval grad not yet implemented for this mesh type");
        }
        break;
    }

    default:
        assert(0 && "Interval grad not implemented for this FE type");
    }

    return g;
}