
namespace algoim_cut_detail {

template <typename Phi>
class NegatedLevelSet {
public:
    explicit NegatedLevelSet(Phi &phi) : phi_(phi) {}

    template <typename X>
    auto operator()(const X &x) const {
        return -phi_(x);
    }

    template <typename X>
    auto grad(const X &x) const {
        auto g = phi_.grad(x);
        g *= -1.0;
        return g;
    }

    template <typename... Args>
    auto normal(Args &&...args) const {
        auto n = phi_.normal(std::forward<Args>(args)...);
        n *= -1.0;
        return n;
    }

private:
    Phi &phi_;
};

} // namespace algoim_cut_detail

// IBP-consistent quadrature correction (triangle multipoly path).
//
// Enforces the discrete divergence theorem on each cut element K:
//   int_{K cap {phi<0}} div F  =  int_{K cap {phi=0}} F.n  +  int_{edges cap {phi<0}} F.n
// exactly (to small dense-solver precision) for all polynomial fields F with
// deg F <= m := ProblemOption::algoim_ibp_degree_.  Rationale: the raw
// multipoly rules carry O(1e-3..1e-5) divergence-theorem defects on cut cells
// whose geometry is hard at the element scale; in pressure-robust Stokes these
// defects enter the momentum residual scaled by |p|/nu.  The element-edge
// (face) integrals are 1D integrals of the Bernstein level set's edge traces,
// computable to machine precision, and serve as the trusted reference.
//
// Stage 1 (correct_surface_rule): minimally correct the surface rule's VECTOR
// weights w_q n_q such that  sum_q F(y_q) . (w_q n_q) = - sum_faces int F.n
// for all divergence-free polynomial F = curl(psi), deg psi <= m+1.  This is
// the part of the divergence theorem that constrains the surface rule alone.
// The constraint set for the region {phi<0} and its complement coincide up to
// the sign of the weights (full-edge integrals of div-free fields vanish
// exactly), so the corrected rule is valid for both sides / either sign
// convention of the caller.
//
// Stage 2 (correct_volume_rule): after Stage 1 the volume moments
// T(p) = B(P) with div P = p are well-defined (independent of the chosen
// antiderivative P); minimally correct the volume weights to reproduce T(p)
// for all deg p <= m-1.
namespace algoim_ibp {

using real = algoim::real;

// Solve (M + ridge*I) y = r for symmetric PSD M (n x n, row-major, by value).
inline bool ridge_cholesky_solve(std::vector<double> M, const std::vector<double> &r, int n,
                                 std::vector<double> &y, double ridge_rel = 1e-12) {
    double dmax = 0;
    for (int i = 0; i < n; ++i) dmax = std::max(dmax, M[i * n + i]);
    const double ridge = std::max(dmax, 1.0) * ridge_rel;
    for (int i = 0; i < n; ++i) M[i * n + i] += ridge;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double s = M[i * n + j];
            for (int k = 0; k < j; ++k) s -= M[i * n + k] * M[j * n + k];
            if (i == j) {
                if (!(s > 0) || !std::isfinite(s)) return false;
                M[i * n + i] = std::sqrt(s);
            } else
                M[i * n + j] = s / M[j * n + j];
        }
    }
    y = r;
    for (int i = 0; i < n; ++i) {
        double s = y[i];
        for (int k = 0; k < i; ++k) s -= M[i * n + k] * y[k];
        y[i] = s / M[i * n + i];
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = y[i];
        for (int k = i + 1; k < n; ++k) s -= M[k * n + i] * y[k];
        y[i] = s / M[i * n + i];
    }
    return true;
}

// Symmetric Jacobi eigendecomposition: M (n x n, row-major, destroyed).
// Eigenvalues in lam; eigenvectors are the COLUMNS of Q (Q[i*n+j] = i-th
// component of eigenvector j).
inline void jacobi_eig(std::vector<double> &M, int n, std::vector<double> &lam,
                       std::vector<double> &Q) {
    Q.assign(n * n, 0.0);
    for (int i = 0; i < n; ++i) Q[i * n + i] = 1.0;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) off += M[p * n + q] * M[p * n + q];
        if (off < 1e-28) break;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) {
                const double apq = M[p * n + q];
                if (std::abs(apq) < 1e-300) continue;
                const double tau = (M[q * n + q] - M[p * n + p]) / (2 * apq);
                const double t = (tau >= 0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1 + tau * tau));
                const double c = 1 / std::sqrt(1 + t * t), s = t * c;
                for (int k = 0; k < n; ++k) {
                    const double mkp = M[k * n + p], mkq = M[k * n + q];
                    M[k * n + p] = c * mkp - s * mkq;
                    M[k * n + q] = s * mkp + c * mkq;
                }
                for (int k = 0; k < n; ++k) {
                    const double mpk = M[p * n + k], mqk = M[q * n + k];
                    M[p * n + k] = c * mpk - s * mqk;
                    M[q * n + k] = s * mpk + c * mqk;
                }
                for (int k = 0; k < n; ++k) {
                    const double qkp = Q[k * n + p], qkq = Q[k * n + q];
                    Q[k * n + p] = c * qkp - s * qkq;
                    Q[k * n + q] = s * qkp + c * qkq;
                }
            }
    }
    lam.resize(n);
    for (int i = 0; i < n; ++i) lam[i] = M[i * n + i];
}

// Truncated min-norm update: solve  min ||d||  s.t.  A d = r  through the
// eigendecomposition of A A^T, DROPPING eigendirections with lam < tol*lam_max
// instead of inverting them.  Near-singular constraint directions (high-degree
// moments on short arcs, tiny slivers) would otherwise amplify machine-level
// residuals into O(1e-6) weight noise; truncation leaves those directions
// uncorrected (their raw residual is retained, which is the lesser evil).
// A final trust-region cap ||d|| <= cap rejects pathological updates entirely;
// the caller then keeps the raw rule.
inline bool trust_region_update(const std::vector<double> &A, const std::vector<double> &r, int nc,
                                size_t nu, double cap, std::vector<double> &d) {
    std::vector<double> M(nc * nc, 0.0), lam, Q;
    for (int i = 0; i < nc; ++i)
        for (int j = 0; j <= i; ++j) {
            double s = 0;
            for (size_t k = 0; k < nu; ++k) s += A[i * nu + k] * A[j * nu + k];
            M[i * nc + j] = M[j * nc + i] = s;
        }
    jacobi_eig(M, nc, lam, Q);
    double lmax = 0;
    for (int i = 0; i < nc; ++i) lmax = std::max(lmax, lam[i]);
    if (!(lmax > 0) || !std::isfinite(lmax)) return false;

    std::vector<double> y(nc, 0.0);
    for (int j = 0; j < nc; ++j) {
        if (!(lam[j] > 1e-12 * lmax)) continue;
        double cj = 0;
        for (int i = 0; i < nc; ++i) cj += Q[i * nc + j] * r[i];
        const double coef = cj / lam[j];
        for (int i = 0; i < nc; ++i) y[i] += coef * Q[i * nc + j];
    }

    d.assign(nu, 0.0);
    double nrm2 = 0;
    for (size_t q = 0; q < nu; ++q) {
        double s = 0;
        for (int i = 0; i < nc; ++i) s += A[i * nu + q] * y[i];
        d[q] = s;
        nrm2 += s * s;
    }
    return std::isfinite(nrm2) && std::sqrt(nrm2) <= cap;
}

// 1D quadrature for the {phiB < 0} portions of the three (reference-)triangle
// edges, returned as physical points with OUTWARD vector weights (unit outward
// normal times arc measure).  phiB is the Bernstein level set on the unit
// square used by the volume/surface rules; its edge traces are polynomials
// (degree bern_deg on the two legs, 2*bern_deg on the diagonal), so the roots
// and the resulting portions are essentially exact.
inline void inside_face_rule(const Mesh2::Element &K, const algoim::xarray<real, 2> &phiB,
                             int bern_deg, int q1d, std::vector<std::array<double, 2>> &pts,
                             std::vector<std::array<double, 2>> &vw) {
    using R2 = typename Mesh2::Rd;
    pts.clear();
    vw.clear();

    const R2 v0(K.at(0)[0], K.at(0)[1]);
    const R2 v1(K.at(1)[0], K.at(1)[1]);
    const R2 v2(K.at(2)[0], K.at(2)[1]);

    for (int e = 0; e < 3; ++e) {
        // reference parametrization, physical endpoints A->B, opposite vertex C
        R2 A, B, C;
        if (e == 0) { A = v0; B = v1; C = v2; }        // ref (t, 0)
        else if (e == 1) { A = v0; B = v2; C = v1; }   // ref (0, t)
        else { A = v1; B = v2; C = v0; }               // ref (1-t, t)

        const int dtr = (e == 2) ? 2 * bern_deg : bern_deg;

        std::vector<real> gdata(dtr + 1);
        algoim::xarray<real, 1> g(gdata.data(), algoim::uvector<int, 1>(dtr + 1));
        algoim::bernstein::bernsteinInterpolate<1>(
            [&](const algoim::uvector<real, 1> &t) -> real {
                const real tt = t(0);
                algoim::uvector<real, 2> xi;
                if (e == 0)      xi = algoim::uvector<real, 2>(tt, real(0));
                else if (e == 1) xi = algoim::uvector<real, 2>(real(0), tt);
                else             xi = algoim::uvector<real, 2>(real(1) - tt, tt);
                return algoim::bernstein::evalBernsteinPoly(phiB, xi);
            },
            g);

        // Robust roots of the edge trace by monotone splitting: split [0,1] at
        // the roots of g' (monotone pieces), then bisect every sign change.
        // The fast Bernstein root-finder alone can MISS near-tangent double
        // crossings (tiny dips of the trace, common where the interface runs
        // nearly parallel to an edge), which would corrupt the face reference
        // by the dip length.
        auto geval = [&](double t) -> double {
            return double(algoim::bernstein::evalBernsteinPoly(g, algoim::uvector<real, 1>(t)));
        };
        std::vector<double> knots;
        knots.push_back(0.0);
        if (dtr >= 2) {
            std::vector<real> gd(dtr);
            for (int i = 0; i < dtr; ++i) gd[i] = dtr * (gdata[i + 1] - gdata[i]);
            std::vector<real> droots(std::max(dtr - 1, 1));
            const int ndr = algoim::bernstein::bernsteinUnitIntervalRealRoots(gd.data(), dtr, droots.data());
            for (int i = 0; i < ndr; ++i)
                if (droots[i] > 0.0 && droots[i] < 1.0) knots.push_back(double(droots[i]));
        }
        knots.push_back(1.0);
        std::sort(knots.begin(), knots.end());

        std::vector<double> brk;
        brk.push_back(0.0);
        for (size_t i = 0; i + 1 < knots.size(); ++i) {
            const double a = knots[i], b = knots[i + 1];
            if (b - a <= 1e-15) continue;
            double ga = geval(a), gb = geval(b);
            if (!(ga * gb < 0)) continue; // no sign change on this monotone piece
            double lo = a, hi = b;
            if (ga > 0) std::swap(lo, hi); // g(lo) < 0 < g(hi)
            for (int it = 0; it < 55; ++it) {
                const double mid = 0.5 * (lo + hi);
                if (geval(mid) < 0) lo = mid;
                else hi = mid;
            }
            brk.push_back(0.5 * (lo + hi));
        }
        brk.push_back(1.0);
        std::sort(brk.begin(), brk.end());

        // physical edge geometry: outward unit normal and length
        const double dxe = B[0] - A[0], dye = B[1] - A[1];
        const double len = std::sqrt(dxe * dxe + dye * dye);
        if (len <= 0) continue;
        double n0 = dye / len, n1 = -dxe / len;
        if (n0 * (C[0] - A[0]) + n1 * (C[1] - A[1]) > 0) { n0 = -n0; n1 = -n1; }

        for (size_t i = 0; i + 1 < brk.size(); ++i) {
            const double ta = brk[i], tb = brk[i + 1];
            if (tb - ta <= 1e-14) continue;
            const double tm = 0.5 * (ta + tb);
            if (algoim::bernstein::evalBernsteinPoly(g, algoim::uvector<real, 1>(tm)) >= 0)
                continue; // outside {phi<0}
            for (int j = 0; j < q1d; ++j) {
                const double t = ta + (tb - ta) * double(algoim::GaussQuad::x(q1d, j));
                const double w = (tb - ta) * len * double(algoim::GaussQuad::w(q1d, j));
                pts.push_back({A[0] + t * dxe, A[1] + t * dye});
                vw.push_back({w * n0, w * n1});
            }
        }
    }
}

// Stage 1: correct the surface rule's vector weights (see file-top comment).
// Convention: `rule` normals are +grad(phi)/|grad(phi)|, outward for the
// region {phi < 0} that phiB's negative side defines.
inline void correct_surface_rule(AlgoimQuadratureRule<Mesh2> &rule, const Mesh2::Element &K,
                                 const algoim::xarray<real, 2> &phiB, int bern_deg, int m) {
    using R2       = typename Mesh2::Rd;
    const size_t ns = rule.points.size();
    if (ns == 0) return;

    std::vector<std::array<double, 2>> fpts, fvw;
    inside_face_rule(K, phiB, bern_deg, 10, fpts, fvw);

    const double cx = (K.at(0)[0] + K.at(1)[0] + K.at(2)[0]) / 3.;
    const double cy = (K.at(0)[1] + K.at(1)[1] + K.at(2)[1]) / 3.;
    const double hs = K.get_h();

    // divergence-free basis F = curl(psi), psi = xi^a eta^b, 1 <= a+b <= m+1
    std::vector<std::pair<int, int>> basis;
    for (int d = 1; d <= m + 1; ++d)
        for (int a = 0; a <= d; ++a) basis.push_back({a, d - a});
    const int nc = (int)basis.size();

    auto Feval = [&](int i, double x, double y, double F[2]) {
        const auto [a, b] = basis[i];
        const double xi = (x - cx) / hs, eta = (y - cy) / hs;
        F[0] = (b == 0) ? 0.0 : b * std::pow(xi, a) * std::pow(eta, b - 1);
        F[1] = (a == 0) ? 0.0 : -a * std::pow(xi, a - 1) * std::pow(eta, b);
    };

    const size_t nu = 2 * ns; // unknowns: vector-weight components
    // constraint matrix (fixed: node positions don't change) and constant
    // face contributions; row-normalized for conditioning
    std::vector<double> A(nc * nu), face(nc), rownorm(nc);
    for (int i = 0; i < nc; ++i) {
        double F[2];
        for (size_t q = 0; q < ns; ++q) {
            Feval(i, rule.points[q][0], rule.points[q][1], F);
            A[i * nu + 2 * q]     = F[0];
            A[i * nu + 2 * q + 1] = F[1];
        }
        double fc = 0;
        for (size_t q = 0; q < fpts.size(); ++q) {
            Feval(i, fpts[q][0], fpts[q][1], F);
            fc += F[0] * fvw[q][0] + F[1] * fvw[q][1];
        }
        face[i] = fc;

        double rn = 0;
        for (size_t j = 0; j < nu; ++j) rn += A[i * nu + j] * A[i * nu + j];
        rownorm[i] = std::sqrt(rn);
        if (rownorm[i] > 0)
            for (size_t j = 0; j < nu; ++j) A[i * nu + j] /= rownorm[i];
    }

    // working state: vector weights
    std::vector<double> omega(nu);
    double wsum = 0;
    for (size_t q = 0; q < ns; ++q) {
        omega[2 * q]     = rule.weights[q] * rule.normals[q][0];
        omega[2 * q + 1] = rule.weights[q] * rule.normals[q][1];
        wsum += rule.weights[q];
    }

    // residual (in row-normalized units) of the constraints at a given omega
    auto residual = [&](const std::vector<double> &om, std::vector<double> &r) -> double {
        double rmax = 0;
        for (int i = 0; i < nc; ++i) {
            double acc = face[i];
            for (size_t j = 0; j < nu; ++j) acc += rownorm[i] * A[i * nu + j] * om[j];
            r[i] = (rownorm[i] > 0) ? -acc / rownorm[i] : 0.0;
            rmax = std::max(rmax, std::abs(r[i]));
        }
        return rmax;
    };

    // iterative refinement with monotone acceptance: each pass re-solves the
    // constraints from the freshly evaluated residual (classic refinement,
    // pushes the single-solve ~1e-12 relative accuracy toward machine -- for
    // problems with large pressure scales, e.g. coriolis w=1e4, that factor
    // shows up directly in the velocity error); a pass is committed only if
    // the residual actually decreases.
    const double cap = 0.25 * (wsum + hs * 1e-8);
    const double tol = 3e-16 * (1.0 + wsum);
    std::vector<double> r(nc), rtry(nc), d, otry(nu);
    double rmax = residual(omega, r);
    bool changed = false;
    for (int pass = 0; pass < 4 && rmax >= tol; ++pass) {
        if (!trust_region_update(A, r, nc, nu, cap, d)) break;
        for (size_t j = 0; j < nu; ++j) otry[j] = omega[j] + d[j];
        const double rmax_try = residual(otry, rtry);
        if (!(rmax_try < rmax)) break; // no improvement: keep previous state
        omega.swap(otry);
        r.swap(rtry);
        rmax    = rmax_try;
        changed = true;
    }
    if (!changed) return;

    for (size_t q = 0; q < ns; ++q) {
        const double o0 = omega[2 * q], o1 = omega[2 * q + 1];
        const double nn = std::sqrt(o0 * o0 + o1 * o1);
        if (!(nn > 1e-14 * hs) || !std::isfinite(nn)) continue; // keep original entry
        rule.weights[q] = nn;
        rule.normals[q] = R2(o0 / nn, o1 / nn);
    }
}

// Stage 2: moment-fit the volume weights to the divergence-theorem moments
// defined by the (Stage-1-corrected) surface rule plus the face rule.
// `surf` must carry outward normals for the region {phi < 0} (i.e. the rule
// produced by quadGenSurf for the same phi whose negative side `vol`
// integrates).
inline void correct_volume_rule(AlgoimQuadratureRule<Mesh2> &vol,
                                const AlgoimQuadratureRule<Mesh2> &surf, const Mesh2::Element &K,
                                const algoim::xarray<real, 2> &phiB, int bern_deg, int m) {
    const size_t nv = vol.points.size(), ns = surf.points.size();
    if (nv == 0 || ns == 0) return;

    std::vector<std::array<double, 2>> fpts, fvw;
    inside_face_rule(K, phiB, bern_deg, 10, fpts, fvw);

    const double cx = (K.at(0)[0] + K.at(1)[0] + K.at(2)[0]) / 3.;
    const double cy = (K.at(0)[1] + K.at(1)[1] + K.at(2)[1]) / 3.;
    const double hs = K.get_h();

    // moment basis p = xi^a eta^b, 0 <= a+b <= m-1; antiderivative
    // P = (hs * xi^{a+1} eta^b / (a+1), 0) satisfies div P = p.
    std::vector<std::pair<int, int>> basis;
    for (int d = 0; d <= m - 1; ++d)
        for (int a = 0; a <= d; ++a) basis.push_back({a, d - a});
    const int nc = (int)basis.size();

    auto peval = [&](int i, double x, double y) -> double {
        const auto [a, b] = basis[i];
        const double xi = (x - cx) / hs, eta = (y - cy) / hs;
        return std::pow(xi, a) * std::pow(eta, b);
    };
    auto P1eval = [&](int i, double x, double y) -> double {
        const auto [a, b] = basis[i];
        const double xi = (x - cx) / hs, eta = (y - cy) / hs;
        return hs * std::pow(xi, a + 1) * std::pow(eta, b) / (a + 1);
    };

    std::vector<double> V(nc * nv), r(nc);
    std::vector<double> T(nc), rownorm(nc);
    for (int i = 0; i < nc; ++i) {
        // target moment from the boundary rules (fixed across passes)
        double Ti = 0;
        for (size_t q = 0; q < ns; ++q)
            Ti += surf.weights[q] * P1eval(i, surf.points[q][0], surf.points[q][1]) * surf.normals[q][0];
        for (size_t q = 0; q < fpts.size(); ++q)
            Ti += P1eval(i, fpts[q][0], fpts[q][1]) * fvw[q][0];
        T[i] = Ti;

        for (size_t q = 0; q < nv; ++q)
            V[i * nv + q] = peval(i, vol.points[q][0], vol.points[q][1]);

        double rn = 0;
        for (size_t j = 0; j < nv; ++j) rn += V[i * nv + j] * V[i * nv + j];
        rownorm[i] = std::sqrt(rn);
        if (rownorm[i] > 0)
            for (size_t j = 0; j < nv; ++j) V[i * nv + j] /= rownorm[i];
    }

    std::vector<double> wq(vol.weights.begin(), vol.weights.end());
    double wsum = 0;
    for (size_t q = 0; q < nv; ++q) wsum += wq[q];

    auto residual = [&](const std::vector<double> &w, std::vector<double> &rr) -> double {
        double rmax = 0;
        for (int i = 0; i < nc; ++i) {
            double acc = 0;
            for (size_t q = 0; q < nv; ++q) acc += w[q] * V[i * nv + q];
            rr[i] = (rownorm[i] > 0) ? (T[i] / rownorm[i] - acc) : 0.0;
            rmax  = std::max(rmax, std::abs(rr[i]));
        }
        return rmax;
    };

    // iterative refinement with monotone acceptance (see correct_surface_rule)
    const double cap = 0.25 * (wsum + hs * hs * 1e-8);
    const double tol = 3e-16 * (1.0 + wsum);
    std::vector<double> rtry(nc), d, wtry(nv);
    double rmax  = residual(wq, r);
    bool changed = false;
    for (int pass = 0; pass < 4 && rmax >= tol; ++pass) {
        if (!trust_region_update(V, r, nc, nv, cap, d)) break;
        for (size_t q = 0; q < nv; ++q) wtry[q] = wq[q] + d[q];
        const double rmax_try = residual(wtry, rtry);
        if (!(rmax_try < rmax)) break;
        wq.swap(wtry);
        r.swap(rtry);
        rmax    = rmax_try;
        changed = true;
    }
    if (!changed) return;

    for (size_t q = 0; q < nv; ++q)
        if (std::isfinite(wq[q])) vol.weights[q] = wq[q];
}

} // namespace algoim_ibp

// Mesh2 specialization
template<typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenVol(const Mesh2::Element& K, Phi& phi, const ProblemOption& option) {

    // Quadrature generation strategy:
    // 1. Map the triangle K in physical coordinates to the reference triangle (0,0)-(1,0)-(0,1) called K_ref
    // 2. K_ref can be defined implicitly as the unit square intersected with {psi < 0} where psi = x+y-1 represents the hypothenuse of K_ref     
    // 3. Add the intersection with {phi < 0} where phi is the levelset function used to describe our domain
    // 4. Generate quadrature rules for volume and surface using Bernstein polynomials for phi (psi is always linear)
    // 5. Map quadrature nodes and weights back to the physical triangle using an affine map

    using real = algoim::real;
    using Vec2 = algoim::uvector<real, 2>;
    using R2   = typename Mesh2::Rd;

    std::array<Vec2,3> vertices = {
        Vec2(K.at(0)[0], K.at(0)[1]),
        Vec2(K.at(1)[0], K.at(1)[1]),
        Vec2(K.at(2)[0], K.at(2)[1])
    };

    AlgoimQuadratureRule<Mesh2> rule;
    
    real tol_phi = 1e-10;
    real tol_psi = 1e-12;

    // Create the affine map F(xi0, xi1) = v0 + xi0*e1 + xi1*e2 mapping from the reference triangle to the physical triangle
    Vec2 v0 = vertices[0];
    Vec2 e1 = vertices[1] - vertices[0];
    Vec2 e2 = vertices[2] - vertices[0];
    auto F = [&](const Vec2& xi) -> Vec2 { return v0 + e1*xi(0) + e2*xi(1); };  
    auto detJ = e1(0)*e2(1) - e1(1)*e2(0);

    // Define the reference level set psi_ref(xi) = xi0 + xi1 - 1
    auto psi_ref = [&](const Vec2& xi) -> real { return xi(0) + xi(1) - real(1); };

    // Interpolate phi(F(x)) on the unit square using Bernstein polynomials
    int bernstein_deg = option.algoim_bernstein_deg_;
    //int q1d = option.algoim_q1d_;
    int q1d = option.algoim_vol_quad_deg_;
    int n = bernstein_deg + 1;
    
    algoim::uvector<int,2> P(n,n);
    std::vector<real> phi_data(n*n), psi_data(n*n);
    
    algoim::xarray<real,2> phiB(phi_data.data(), P);
    algoim::xarray<real,2> psiB(psi_data.data(), P);
    
    algoim::bernstein::bernsteinInterpolate<2>(
        [&](const Vec2& xi) -> real { return phi(F(xi)); },
        phiB
    );
    algoim::bernstein::bernsteinInterpolate<2>(
        [&](const Vec2& xi) -> real { return psi_ref(xi); },
        psiB
    );

    // Compute quadrature using the implicit polynomial quadrature method
    algoim::ImplicitPolyQuadrature<2> ipquad(phiB, psiB);

    const auto strategy = static_cast<algoim::QuadStrategy>(option.algoim_quad_strategy_);
    ipquad.integrate(strategy, q1d, [&](const Vec2& xi, real w_ref)
    {
        // pick component: inside triangle AND phi<0
        if (psi_ref(xi) >= 0) return; // exact
        if (algoim::bernstein::evalBernsteinPoly(phiB, xi) >= 0) return;
        // if (phi(F(xi)) >= 0) return;

        Vec2 x = F(xi);
        rule.points.emplace_back(R2(x(0), x(1)));
        rule.weights.emplace_back(double(w_ref) * std::abs(detJ));
    });

    if (option.algoim_ibp_consistent_ && rule.points.size() > 0) {
        // Stage 2 of the IBP correction: fit the volume weights to the
        // divergence-theorem moments of the (Stage-1-corrected) boundary rule.
        // quadGenSurf below is called with the same phi, so its normals are
        // outward for the region {phi<0} integrated here, and (with the flag
        // set) it returns the Stage-1-corrected rule the assembly also uses.
        auto surf = quadGenSurf(K, phi, option);
        if (surf.points.size() > 0)
            algoim_ibp::correct_volume_rule(rule, surf, K, phiB, bernstein_deg,
                                            option.algoim_ibp_degree_);
    }

    return rule;
}

template<typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenVol(const Mesh2::Element& K, Phi& phi, const ProblemOption& option, int domain) {
    assert(domain == 0 || domain == 1);
    if (domain == 0) {
        algoim_cut_detail::NegatedLevelSet<Phi> positive_phi(phi);
        return quadGenVol(K, positive_phi, option);
    }
    return quadGenVol(K, phi, option);
}

template <typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenSurf(const Mesh2::Element& K, Phi& phi, const ProblemOption& option) {
    // Quadrature generation strategy:
    // 1. Map the triangle K in physical coordinates to the reference triangle (0,0)-(1,0)-(0,1) called K_ref
    // 2. K_ref can be defined implicitly as the unit square intersected with {psi < 0} where psi = x+y-1 represents the hypothenuse of K_ref     
    // 3. Add the intersection with {phi = 0} where phi is the levelset function used to describe our domain
    // 4. Generate quadrature rules for surface using Bernstein polynomials for phi (psi is always linear)
    // 5. Map quadrature nodes, weights and normals back to the physical triangle using an affine map

    using real = algoim::real;
    using Vec2 = algoim::uvector<real, 2>;
    using R2   = typename Mesh2::Rd;

    AlgoimQuadratureRule<Mesh2> rule;

    // --- physical triangle vertices
    const Vec2 v0(K.at(0)[0], K.at(0)[1]);
    const Vec2 v1(K.at(1)[0], K.at(1)[1]);
    const Vec2 v2(K.at(2)[0], K.at(2)[1]);

    // Affine map x = F(xi) = v0 + J*xi,   J = [e1 e2]
    const Vec2 e1 = v1 - v0;
    const Vec2 e2 = v2 - v0;

    const auto F = [&](const Vec2& xi) -> Vec2 {
        return v0 + e1*xi(0) + e2*xi(1);
    };

    const real detJ = e1(0)*e2(1) - e1(1)*e2(0);
    if (std::abs(detJ) < real(1e-30)) {
        // Degenerate element
        return rule;
    }

    // cofactor(J) = det(J) * J^{-T} = [[ e2y, -e1y ],
    //                                 [ -e2x,  e1x ]]
    // This maps (n ds)_ref -> (n ds)_phys for a curve under affine map.
    const auto cofJ_mul = [&](const Vec2& v) -> Vec2 {
        return Vec2(
            e2(1)*v(0) - e1(1)*v(1),
           -e2(0)*v(0) + e1(0)*v(1)
        );
    };

    // Reference triangle cut: psi(xi) = xi0 + xi1 - 1 < 0
    const auto psi_ref = [&](const Vec2& xi) -> real {
        return xi(0) + xi(1) - real(1);
    };

    // --- Bernstein interpolation on unit square [0,1]^2 (as algoim expects)
    const int bernstein_deg = option.algoim_bernstein_deg_;
    const int q1d           = option.algoim_surface_quad_deg_;

    const int n = bernstein_deg + 1;
    algoim::uvector<int,2> P(n,n);

    std::vector<real> phi_data(n*n), psi_data(n*n);
    algoim::xarray<real,2> phiB(phi_data.data(), P);
    algoim::xarray<real,2> psiB(psi_data.data(), P);

    algoim::bernstein::bernsteinInterpolate<2>(
        [&](const Vec2& xi) -> real { return phi(F(xi)); }, phiB
    );
    algoim::bernstein::bernsteinInterpolate<2>(
        [&](const Vec2& xi) -> real { return psi_ref(xi); }, psiB
    );

    algoim::ImplicitPolyQuadrature<2> ipquad(phiB, psiB);

    // Tolerances: only to filter out the psi=0 part (triangle edge)
    const real tol_inside = real(1e-14);
    const real tol_psi0   = real(1e-10);

    const auto strategy = static_cast<algoim::QuadStrategy>(option.algoim_quad_strategy_);
    ipquad.integrate_surf(strategy, q1d,
    [&](const Vec2& xi, real /*w_ref*/, const Vec2& wn_ref)
    {
        const real ph  = algoim::bernstein::evalBernsteinPoly(phiB, xi);
        // const real ph = phi(F(xi));
        const real ps  = psi_ref(xi); // exact
        if (ps >= 0) return;          // must be inside triangle

        // decide which constraint generated this boundary point
        if (std::abs(ps) < std::abs(ph)) return; // this is psi-boundary => skip

        // map point
        Vec2 x = F(xi);

        // map (n ds) using cof(J)
        Vec2 wn_phys = cofJ_mul(wn_ref);
        real w = std::sqrt(wn_phys(0)*wn_phys(0) + wn_phys(1)*wn_phys(1));
        if (w <= 0) return;

        Vec2 n = wn_phys / w;

        // orientation: make n outward for the phase "phi<0"
        // Vec2 g_ref = algoim::bernstein::evalBernsteinPolyGradient(phiB, xi);
        // Vec2 g_phys(
        //     ( e2(1)*g_ref(0) - e1(1)*g_ref(1)) / detJ,
        //     (-e2(0)*g_ref(0) + e1(0)*g_ref(1)) / detJ
        // );
        // if (n(0)*g_phys(0) + n(1)*g_phys(1) < 0) n = -n;
        
        
        // const R2 normal = phi.normal(R2(x(0), x(1)));
        rule.points.emplace_back(R2(x(0), x(1)));
        rule.weights.emplace_back(double(w));
        rule.normals.emplace_back(R2(n(0), n(1)));
        // rule.normals.emplace_back(normal);
    });

    if (option.algoim_ibp_consistent_ && rule.points.size() > 0) {
        // Stage 1 of the IBP correction: make the vector weights w*n satisfy
        // the divergence theorem for div-free polynomial fields against exact
        // 1D integrals over the element-edge portions.
        algoim_ibp::correct_surface_rule(rule, K, phiB, bernstein_deg,
                                         option.algoim_ibp_degree_);
    }

    return rule;
}

template<typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenFace(const Mesh2::Element& K, Phi& phi, const ProblemOption& option, int ifac) {
    // To be implemented
    assert(0);
    AlgoimQuadratureRule<Mesh2> rule;
    return rule;
}

template<typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenFace(const Mesh2::Element& K, Phi& phi, const ProblemOption& option, int ifac, int domain) {
    assert(domain == 0 || domain == 1);
    if (domain == 0) {
        algoim_cut_detail::NegatedLevelSet<Phi> positive_phi(phi);
        return quadGenFace(K, positive_phi, option, ifac);
    }
    return quadGenFace(K, phi, option, ifac);
}




// MeshQuad2 specialization
template<typename Phi>
AlgoimQuadratureRule<MeshQuad2> quadGenVol(const MeshQuad2::Element& K, Phi& phi, const ProblemOption& option) {
    // Get coordinates of current quadrilateral
    const auto &V0(K.at(0)); // vertex 0
    const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

    algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
    algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

    algoim::QuadratureRule<2> q =
        algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, option.algoim_vol_quad_deg_);

    AlgoimQuadratureRule<MeshQuad2> rule;
    for (size_t ipq = 0; ipq < q.nodes.size(); ++ipq) {
        rule.points.push_back(R2(q.nodes[ipq].x(0), q.nodes[ipq].x(1)));
        rule.weights.push_back(q.nodes[ipq].w);
    }

    // std::cout << "quad_rule.size() = " << rule.size() << " in quadGenVol for MeshQuad2\n";
    return rule;
}

template<typename Phi>
AlgoimQuadratureRule<MeshQuad2> quadGenVol(const MeshQuad2::Element& K, Phi& phi, const ProblemOption& option, int domain) {
    assert(domain == 0 || domain == 1);
    if (domain == 0) {
        algoim_cut_detail::NegatedLevelSet<Phi> positive_phi(phi);
        return quadGenVol(K, positive_phi, option);
    }
    return quadGenVol(K, phi, option);
}

template <typename Phi>
AlgoimQuadratureRule<MeshQuad2> quadGenSurf(const MeshQuad2::Element& K, Phi& phi, const ProblemOption& option) {
    // Get coordinates of current quadrilateral
    const auto &V0(K.at(0)); // vertex 0
    const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

    algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
    algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

    algoim::QuadratureRule<2> q =
        algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, option.algoim_surface_quad_deg_);
    
    AlgoimQuadratureRule<MeshQuad2> rule;
    for (size_t ipq = 0; ipq < q.nodes.size(); ++ipq) {
        rule.points.push_back(R2(q.nodes[ipq].x(0), q.nodes[ipq].x(1)));
        rule.weights.push_back(q.nodes[ipq].w);
        rule.normals.push_back(phi.normal(R2(q.nodes[ipq].x(0), q.nodes[ipq].x(1))));
    }

    return rule;
}
template<typename Phi>
AlgoimQuadratureRule<MeshQuad2> quadGenFace(const MeshQuad2::Element& K, Phi& phi, const ProblemOption& option, int ifac) {
    using R2 = typename MeshQuad2::Rd;
    
    const auto &V0(K[0]);
    const auto &V2(K[2]);

    algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
    algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

    int dim = -1, side = -1;
    switch (ifac) {
        case 0: dim = 1; side = 0; break; // bottom
        case 1: dim = 0; side = 1; break; // right
        case 2: dim = 1; side = 1; break; // top
        case 3: dim = 0; side = 0; break; // left
        default: assert(false && "Unexpected face index in isCutFace");
    }
    
    algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), dim, side, option.algoim_surface_quad_deg_);
    
    const R2 normal = (dim == 0) ? R2((side == 0) ? -1 : 1, 0) : R2(0, (side == 0) ? -1 : 1);

    AlgoimQuadratureRule<MeshQuad2> rule;
    for (size_t ipq = 0; ipq < q.nodes.size(); ++ipq) {
        rule.points.push_back(R2(q.nodes[ipq].x(0), q.nodes[ipq].x(1)));
        rule.weights.push_back(q.nodes[ipq].w);
        rule.normals.push_back(normal);
    }
    return rule;
}

template<typename Phi>
AlgoimQuadratureRule<MeshQuad2> quadGenFace(const MeshQuad2::Element& K, Phi& phi, const ProblemOption& option, int ifac, int domain) {
    assert(domain == 0 || domain == 1);
    if (domain == 0) {
        algoim_cut_detail::NegatedLevelSet<Phi> positive_phi(phi);
        return quadGenFace(K, positive_phi, option, ifac);
    }
    return quadGenFace(K, phi, option, ifac);
}


// Gustaf: MeshHexa specialization
template<typename Phi>
AlgoimQuadratureRule<MeshHexa> quadGenVol(const MeshHexa::Element& K, Phi& phi, const ProblemOption& option) {
    using R3 = typename MeshHexa::Rd;

    const auto& V0 = K.at(0); // vertex 0: min x, y, z
    const auto& V6 = K.at(6); // opposite vertex: max x, y, z

    algoim::uvector<double, 3> xyzmin{V0[0], V0[1], V0[2]};
    algoim::uvector<double, 3> xyzmax{V6[0], V6[1], V6[2]};

    algoim::QuadratureRule<3> q =
        algoim::quadGen<3>(
            phi,
            algoim::HyperRectangle<double, 3>(xyzmin, xyzmax),
            -1,
            -1,
            option.algoim_vol_quad_deg_
        );

    AlgoimQuadratureRule<MeshHexa> rule;
    for (size_t ipq = 0; ipq < q.nodes.size(); ++ipq) {
        rule.points.push_back(R3(q.nodes[ipq].x(0),
                                 q.nodes[ipq].x(1),
                                 q.nodes[ipq].x(2)));
        rule.weights.push_back(q.nodes[ipq].w);
    }

    return rule;
}

template<typename Phi>
AlgoimQuadratureRule<MeshHexa> quadGenVol(const MeshHexa::Element& K, Phi& phi, const ProblemOption& option, int domain) {
    assert(domain == 0 || domain == 1);
    if (domain == 0) {
        algoim_cut_detail::NegatedLevelSet<Phi> positive_phi(phi);
        return quadGenVol(K, positive_phi, option);
    }
    return quadGenVol(K, phi, option);
}

template<typename Phi>
AlgoimQuadratureRule<MeshHexa> quadGenSurf(const MeshHexa::Element& K, Phi& phi, const ProblemOption& option) {
    using R3 = typename MeshHexa::Rd;

    const auto& V0 = K.at(0); // vertex 0: min x, y, z
    const auto& V6 = K.at(6); // opposite vertex: max x, y, z

    algoim::uvector<double, 3> xyzmin{V0[0], V0[1], V0[2]};
    algoim::uvector<double, 3> xyzmax{V6[0], V6[1], V6[2]};

    algoim::QuadratureRule<3> q =
        algoim::quadGen<3>(
            phi,
            algoim::HyperRectangle<double, 3>(xyzmin, xyzmax),
            3,
            -1,
            option.algoim_surface_quad_deg_
        );

    AlgoimQuadratureRule<MeshHexa> rule;
    for (size_t ipq = 0; ipq < q.nodes.size(); ++ipq) {
        R3 x(q.nodes[ipq].x(0),
             q.nodes[ipq].x(1),
             q.nodes[ipq].x(2));

        rule.points.push_back(x);
        rule.weights.push_back(q.nodes[ipq].w);
        rule.normals.push_back(phi.normal(x));
    }

    return rule;
}

template<typename Phi>
AlgoimQuadratureRule<MeshHexa> quadGenFace(const MeshHexa::Element& K, Phi& phi, const ProblemOption& option, int ifac) {
    using R3 = typename MeshHexa::Rd;

    const auto& V0 = K.at(0);
    const auto& V6 = K.at(6);

    algoim::uvector<double, 3> xyzmin{V0[0], V0[1], V0[2]};
    algoim::uvector<double, 3> xyzmax{V6[0], V6[1], V6[2]};

    int dim  = -1;
    int side = -1;

    switch (ifac) {
        case 0: dim = 2; side = 0; break; // z-min: face {0,1,2,3}
        case 1: dim = 1; side = 0; break; // y-min: face {0,1,5,4}
        case 2: dim = 0; side = 1; break; // x-max: face {1,2,6,5}
        case 3: dim = 1; side = 1; break; // y-max: face {2,3,7,6}
        case 4: dim = 0; side = 0; break; // x-min: face {3,0,4,7}
        case 5: dim = 2; side = 1; break; // z-max: face {4,5,6,7}
        default:
            assert(false && "Unexpected face index in quadGenFace<MeshHexa>");
    }

    algoim::QuadratureRule<3> q =
        algoim::quadGen<3>(
            phi,
            algoim::HyperRectangle<double, 3>(xyzmin, xyzmax),
            dim,
            side,
            option.algoim_surface_quad_deg_
        );

    R3 normal(0., 0., 0.);
    normal[dim] = (side == 0) ? -1.0 : 1.0;

    AlgoimQuadratureRule<MeshHexa> rule;
    for (size_t ipq = 0; ipq < q.nodes.size(); ++ipq) {
        rule.points.push_back(R3(q.nodes[ipq].x(0),
                                 q.nodes[ipq].x(1),
                                 q.nodes[ipq].x(2)));
        rule.weights.push_back(q.nodes[ipq].w);
        rule.normals.push_back(normal);
    }

    return rule;
}

template<typename Phi>
AlgoimQuadratureRule<MeshHexa> quadGenFace(const MeshHexa::Element& K, Phi& phi, const ProblemOption& option, int ifac, int domain) {
    assert(domain == 0 || domain == 1);
    if (domain == 0) {
        algoim_cut_detail::NegatedLevelSet<Phi> positive_phi(phi);
        return quadGenFace(K, positive_phi, option, ifac);
    }
    return quadGenFace(K, phi, option, ifac);
}

// Gustaf: Mesh3 specialization
template<typename Phi>
AlgoimQuadratureRule<Mesh3> quadGenVol(const Mesh3::Element& K, Phi& phi, const ProblemOption& option) {
    // Quadrature generation strategy:
    // 1. Map the tetrahedron K in physical coordinates to the reference tetrahedron
    //    (0,0,0)-(1,0,0)-(0,1,0)-(0,0,1).
    // 2. K_ref is represented as the unit cube intersected with {psi < 0},
    //    where psi = x + y + z - 1.
    // 3. Add the intersection with {phi < 0}.
    // 4. Generate quadrature using Bernstein polynomials for phi and psi.
    // 5. Map quadrature nodes and weights back to the physical tetrahedron.

    using real = algoim::real;
    using Vec3 = algoim::uvector<real, 3>;
    using R3   = typename Mesh3::Rd;

    AlgoimQuadratureRule<Mesh3> rule;

    const R3 A = K.at(0);
    const R3 B = K.at(1);
    const R3 C = K.at(2);
    const R3 D = K.at(3);

    const R3 AB(A, B);
    const R3 AC(A, C);
    const R3 AD(A, D);

    const double detJ = det(AB, AC, AD);
    if (std::abs(detJ) < 1e-30) {
        return rule;
    }

    const auto F = [&](const Vec3& xi) -> R3 {
        return A + xi(0) * AB + xi(1) * AC + xi(2) * AD;
    };

    const auto psi_ref = [](const Vec3& xi) -> real {
        return xi(0) + xi(1) + xi(2) - real(1);
    };

    const int bernstein_deg = option.algoim_bernstein_deg_;
    const int q1d           = option.algoim_vol_quad_deg_;
    const int n             = bernstein_deg + 1;

    algoim::uvector<int, 3> P(n, n, n);

    std::vector<real> phi_data(n * n * n);
    std::vector<real> psi_data(n * n * n);

    algoim::xarray<real, 3> phiB(phi_data.data(), P);
    algoim::xarray<real, 3> psiB(psi_data.data(), P);

    algoim::bernstein::bernsteinInterpolate<3>(
        [&](const Vec3& xi) -> real {
            return phi(F(xi));
        },
        phiB
    );

    algoim::bernstein::bernsteinInterpolate<3>(
        [&](const Vec3& xi) -> real {
            return psi_ref(xi);
        },
        psiB
    );

    algoim::ImplicitPolyQuadrature<3> ipquad(phiB, psiB);

    ipquad.integrate(algoim::AutoMixed, q1d,
        [&](const Vec3& xi, real w_ref) {
            // Keep only points inside the reference tetrahedron.
            if (psi_ref(xi) >= real(0)) {
                return;
            }

            // Keep only the negative phase of phi.
            if (algoim::bernstein::evalBernsteinPoly(phiB, xi) >= real(0)) {
                return;
            }

            const R3 x = F(xi);

            rule.points.emplace_back(x);
            rule.weights.emplace_back(double(w_ref) * std::abs(detJ));
        }
    );

    return rule;
}

template<typename Phi>
AlgoimQuadratureRule<Mesh3> quadGenVol(const Mesh3::Element& K, Phi& phi, const ProblemOption& option, int domain) {
    assert(domain == 0 || domain == 1);
    if (domain == 0) {
        algoim_cut_detail::NegatedLevelSet<Phi> positive_phi(phi);
        return quadGenVol(K, positive_phi, option);
    }
    return quadGenVol(K, phi, option);
}

template<typename Phi>
AlgoimQuadratureRule<Mesh3> quadGenSurf(const Mesh3::Element& K, Phi& phi, const ProblemOption& option) {
    // Surface quadrature on {phi = 0} inside a tetrahedron.
    //
    // The tetrahedron is represented in the unit cube by psi_ref < 0,
    // where psi_ref = xi0 + xi1 + xi2 - 1.
    //
    // integrate_surf may return points on either phi = 0 or psi_ref = 0.
    // We keep only the phi = 0 part and discard the artificial tetrahedron boundary.

    using real = algoim::real;
    using Vec3 = algoim::uvector<real, 3>;
    using R3   = typename Mesh3::Rd;

    AlgoimQuadratureRule<Mesh3> rule;

    const R3 A = K.at(0);
    const R3 B = K.at(1);
    const R3 C = K.at(2);
    const R3 D = K.at(3);

    const R3 AB(A, B);
    const R3 AC(A, C);
    const R3 AD(A, D);

    const double detJ = det(AB, AC, AD);
    if (std::abs(detJ) < 1e-30) {
        return rule;
    }

    const auto F = [&](const Vec3& xi) -> R3 {
        return A + xi(0) * AB + xi(1) * AC + xi(2) * AD;
    };

    const auto psi_ref = [](const Vec3& xi) -> real {
        return xi(0) + xi(1) + xi(2) - real(1);
    };

    const int bernstein_deg = option.algoim_bernstein_deg_;
    const int q1d           = option.algoim_surface_quad_deg_;
    const int n             = bernstein_deg + 1;

    algoim::uvector<int, 3> P(n, n, n);

    std::vector<real> phi_data(n * n * n);
    std::vector<real> psi_data(n * n * n);

    algoim::xarray<real, 3> phiB(phi_data.data(), P);
    algoim::xarray<real, 3> psiB(psi_data.data(), P);

    algoim::bernstein::bernsteinInterpolate<3>(
        [&](const Vec3& xi) -> real {
            return phi(F(xi));
        },
        phiB
    );

    algoim::bernstein::bernsteinInterpolate<3>(
        [&](const Vec3& xi) -> real {
            return psi_ref(xi);
        },
        psiB
    );

    algoim::ImplicitPolyQuadrature<3> ipquad(phiB, psiB);

    ipquad.integrate_surf(algoim::AutoMixed, q1d,
        [&](const Vec3& xi, real /*w_ref*/, const Vec3& wn_ref) {
            const real ph = algoim::bernstein::evalBernsteinPoly(phiB, xi);
            const real ps = psi_ref(xi);

            // Must be inside the reference tetrahedron.
            if (ps >= real(0)) {
                return;
            }

            // Discard the artificial boundary psi_ref = 0.
            // This mirrors the Mesh2 logic.
            if (std::abs(ps) < std::abs(ph)) {
                return;
            }

            const R3 x = F(xi);

            // Cofactor columns for J = [AB AC AD]:
            //
            // cof(J) e0 = AC x AD
            // cof(J) e1 = AD x AB
            // cof(J) e2 = AB x AC
            //
            // This maps weighted reference normals to weighted physical normals:
            // (n dS)_phys = cof(J) (n dS)_ref.
            const R3 c0 = AC ^ AD;
            const R3 c1 = AD ^ AB;
            const R3 c2 = AB ^ AC;

            R3 wn_phys = wn_ref(0) * c0 + wn_ref(1) * c1 + wn_ref(2) * c2;

            const double w = wn_phys.norme();
            if (w <= 0.) {
                return;
            }

            R3 n_phys = wn_phys / w;

            rule.points.emplace_back(x);
            rule.weights.emplace_back(w);
            rule.normals.emplace_back(n_phys);
        }
    );

    return rule;
}

template<typename Phi>
AlgoimQuadratureRule<Mesh3> quadGenFace(const Mesh3::Element& K, Phi& phi, const ProblemOption& option, int ifac) {
    // Face quadrature on one triangular face of a tetrahedron.
    //
    // The face is parameterised by
    //
    //     F(u,v) = A + u AB + v AC
    //
    // over the reference triangle u >= 0, v >= 0, u + v <= 1.
    // The reference triangle is represented inside the unit square by
    //
    //     psi_face(u,v) = u + v - 1 < 0.
    //
    // We integrate over the part where phi(F(u,v)) < 0.

    using real = algoim::real;
    using Vec2 = algoim::uvector<real, 2>;
    using R3   = typename Mesh3::Rd;
    using Element = typename Mesh3::Element;

    assert(ifac >= 0 && ifac < 4);

    AlgoimQuadratureRule<Mesh3> rule;

    const int iv0 = Element::nvface[ifac][0];
    const int iv1 = Element::nvface[ifac][1];
    const int iv2 = Element::nvface[ifac][2];

    const R3 A = K.at(iv0);
    const R3 B = K.at(iv1);
    const R3 C = K.at(iv2);

    const R3 AB(A, B);
    const R3 AC(A, C);

    const double detJ_face = (AB ^ AC).norme();
    if (detJ_face < 1e-30) {
        return rule;
    }

    const auto F = [&](const Vec2& uv) -> R3 {
        return A + uv(0) * AB + uv(1) * AC;
    };

    const auto psi_face = [](const Vec2& uv) -> real {
        return uv(0) + uv(1) - real(1);
    };

    const int bernstein_deg = option.algoim_bernstein_deg_;
    const int q1d           = option.algoim_surface_quad_deg_;
    const int n             = bernstein_deg + 1;

    algoim::uvector<int, 2> P(n, n);

    std::vector<real> phi_data(n * n);
    std::vector<real> psi_data(n * n);

    algoim::xarray<real, 2> phiB(phi_data.data(), P);
    algoim::xarray<real, 2> psiB(psi_data.data(), P);

    algoim::bernstein::bernsteinInterpolate<2>(
        [&](const Vec2& uv) -> real {
            return phi(F(uv));
        },
        phiB
    );

    algoim::bernstein::bernsteinInterpolate<2>(
        [&](const Vec2& uv) -> real {
            return psi_face(uv);
        },
        psiB
    );

    algoim::ImplicitPolyQuadrature<2> ipquad(phiB, psiB);

    const R3 normal = K.n(ifac);

    ipquad.integrate(algoim::AutoMixed, q1d,
        [&](const Vec2& uv, real w_ref) {
            // Keep only points inside the triangular face.
            if (psi_face(uv) >= real(0)) {
                return;
            }

            // Keep only the negative phase of phi on this face.
            if (algoim::bernstein::evalBernsteinPoly(phiB, uv) >= real(0)) {
                return;
            }

            const R3 x = F(uv);

            rule.points.emplace_back(x);
            rule.weights.emplace_back(double(w_ref) * detJ_face);
            rule.normals.emplace_back(normal);
        }
    );

    return rule;
}

template<typename Phi>
AlgoimQuadratureRule<Mesh3> quadGenFace(const Mesh3::Element& K, Phi& phi, const ProblemOption& option, int ifac, int domain) {
    assert(domain == 0 || domain == 1);
    if (domain == 0) {
        algoim_cut_detail::NegatedLevelSet<Phi> positive_phi(phi);
        return quadGenFace(K, positive_phi, option, ifac);
    }
    return quadGenFace(K, phi, option, ifac);
}

// End of Gustaf


template <typeMesh Mesh, typename Phi>
void AlgoimCutFEMUnified<Mesh, Phi>::addElementContribution(const itemVFlist_t& VF, const int k, const TimeSlab* In, int itq, double cst_time) {
    const fespace_t& Vh(VF.get_spaceV(0));
    const auto& Th(Vh.get_mesh());
    const fe_element_t& FK(Vh[k]);
    const element_t&  K(FK.T);

    const int domain = FK.get_domain();
    const int kb     = Vh.idxElementInBackMesh(k);
    const int iam    = omp_get_thread_num();

    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;
    if constexpr (std::is_same_v<std::remove_cvref_t<Phi>, FunFEM<Mesh>>) {
        phi_.setTime(tid);
        phi_.setElementFromBackMesh(kb);
    } else {
        phi_.t = tid;
    }
    auto quad_rule = quadGenVol(K, phi_, options_, domain);

    if (quad_rule.points.size() == 0) {
        return;
    }

    // Data for the VirtualParameter coefficient lists (CutFEMParameter etc.).
    // These were previously NEVER evaluated in the algoim paths -- coefficients
    // like a per-domain viscosity eta were silently dropped on cut elements
    // (uncut elements go through the coefficient-aware standard BaseFEM path),
    // caught by the polynomial-reproduction MMS in
    // workfiles/src/active_surfaces/algoim_axisymmetric_stokes_mms.cpp.
    const double h_    = K.get_h();
    const double measK = K.measure();
    double meas_cut    = 0.0;
    for (const double w : quad_rule.weights) meas_cut += w;

    for (int l=0; l<VF.size(); ++l) {
        if (!VF[l].on(domain)) continue;

        const double coef_param = VF[l].computeCoefElement(h_, meas_cut, measK, meas_cut, domain);

        const fespace_t& Vhv(VF.get_spaceV(l));
        const fespace_t& Vhu(VF.get_spaceU(l));
        int kv = k, ku = k;
        if (&Vhv != &Vh) {
            kv = Vhv.idxElementFromBackMesh(kb, domain);
            if (kv < 0) continue;
        }
        if (&Vhu != &Vh) {
            ku = (&Vhu == &Vhv) ? kv : Vhu.idxElementFromBackMesh(kb, domain);
            if (ku < 0) continue;
        }
        const auto& FKv(Vhv[kv]);
        const auto& FKu(Vhu[ku]);
        this->initIndex(FKu, FKv);

        bool same  = (&Vhu == &Vhv);
        int lastop = getLastop(VF[l].du, VF[l].dv);

        long offset = iam * this->offset_bf_;
        RNMK_ fv(this->databf_ + offset, FKv.NbDoF(), FKv.N, lastop);
        RNMK_ fu(this->databf_ + offset + (same?0:FKv.NbDoF()*FKv.N*lastop),
                 FKu.NbDoF(), FKu.N, lastop);
        What_d Fop = Fwhatd(lastop);

        // Loop over quadrature points
        for (size_t ipq = 0; ipq < quad_rule.points.size(); ++ipq) {
            
            Rd mip = quad_rule.points[ipq];
            Rd cut_ip = K.mapToReferenceElement(mip);
            double Cint = quad_rule.weights[ipq]*cst_time;
            
            FKv.BF(Fop, cut_ip, fv);
            if (!same) FKu.BF(Fop, cut_ip, fu);

            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);
            Cint *= coef_param * VF[l].c;

            if (In) {
                if (VF.isRHS()) this->addToRHS(VF[l], *In, FKv, fv, Cint);
                else            this->addToMatrix(VF[l], *In, FKu, FKv, fu, fv, Cint);
            } else {
                if (VF.isRHS()) this->addToRHS(VF[l], FKv, fv, Cint);
                else            this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
            }
        }
    }
}


// template <typeMesh Mesh, typename Phi>
// template <typename Fct>
// void AlgoimCutFEMUnified<Mesh, Phi>::addElementContributionExact(const Fct &f, const itemVFlist_t& VF, const int k, const TimeSlab* In, int itq, double cst_time) {
//     const fespace_t& Vh(VF.get_spaceV(0));
//     const auto& Th(Vh.get_mesh());
//     const fe_element_t& FK(Vh[k]);
//     const element_t&  K(FK.T);

//     const int domain = FK.get_domain();
//     const int kb     = Vh.idxElementInBackMesh(k);
//     const int iam    = omp_get_thread_num();

//     auto tq    = this->get_quadrature_time(itq);
//     double tid = (In) ? (double)In->map(tq) : 0.;
//     phi_.t = tid; // if your L has time

//     // Get quadrature rule - template argument deduction from K and phi_
//     auto quad_rule = quadGenVol(K, phi_, options_);

//     if (quad_rule.points.size() == 0) {
//         std::cout << "Warning: no volume quadrature points for cut element element kb = " << kb << " in AlgoimCutFEMUnified::addElementContribution\n";
//         return;
//     }

//     // assemble (same as your current inner loops)
//     for (int l=0; l<VF.size(); ++l) {
//         if (!VF[l].on(domain)) continue;

//         const fespace_t& Vhv(VF.get_spaceV(l));
//         const fespace_t& Vhu(VF.get_spaceU(l));
//         const auto& FKv(Vhv[k]);
//         const auto& FKu(Vhu[k]);
//         this->initIndex(FKu, FKv);

//         bool same  = (&Vhu == &Vhv);
//         int lastop = getLastop(VF[l].du, VF[l].dv);

//         long offset = iam * this->offset_bf_;
//         RNMK_ fv(this->databf_ + offset, FKv.NbDoF(), FKv.N, lastop);
//         RNMK_ fu(this->databf_ + offset + (same?0:FKv.NbDoF()*FKv.N*lastop),
//                  FKu.NbDoF(), FKu.N, lastop);
//         What_d Fop = Fwhatd(lastop);

//         // Loop over quadrature points
//         for (size_t ipq = 0; ipq < quad_rule.points.size(); ++ipq) {
            
//             Rd mip = quad_rule.points[ipq];
//             Rd cut_ip = K.mapToReferenceElement(mip);
//             double Cint = quad_rule.weights[ipq]*cst_time;
            
//             FKv.BF(Fop, cut_ip, fv);
//             if (!same) FKu.BF(Fop, cut_ip, fu);

//             Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);
//                 Cint *= VF[l].c;
//                 Cint *= f(mip, VF[l].cv, tid); 

//                 // std::cout << "ALGOIMCUTFEM: k = " << k << ", Cint = " << Cint << std::endl;
//                 // getchar();
//                 if (In) {
//                     if (VF.isRHS()) this->addToRHS(VF[l], *In, FKv, fv, Cint);
//                     else            this->addToMatrix(VF[l], *In, FKu, FKv, fu, fv, Cint);
//                 } else {
//                     if (VF.isRHS()) this->addToRHS(VF[l], FKv, fv, Cint);
//                     else            this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
//                 }
//             }
//         }
//     }

template <typeMesh Mesh, typename Phi>
void AlgoimCutFEMUnified<Mesh, Phi>::addInterfaceContribution(const itemVFlist_t& VF, const Interface<mesh_t>& interface, int ifac, double tid,
                                  const TimeSlab* In, double cst_time, int itq) {

    //  GET IDX ELEMENT CONTAINING FACE ON backMes
    const int kb = interface.idxElementOfFace(ifac);
    const auto &K(interface.get_element(kb));

    // Select the surface cut-quadrature rule.
    AlgoimQuadratureRule<Mesh> generated_quad_rule;            // storage if regenerated
    const AlgoimQuadratureRule<Mesh> *quad_rule = nullptr;

    if (use_stored_interface_rule_) {
        // Read the rule stored on this (per-node) AlgoimInterface, built from the
        // node's own level set -- the standard-CutFEM TimeInterface behavior.  No
        // dependence on phi_ at all, so distinct time nodes integrate exactly
        // their own surfaces (needed for slab-to-slab mass conservation).
        if (const auto *algoim_interface =
                dynamic_cast<const AlgoimInterface<Mesh, std::remove_cvref_t<Phi>> *>(&interface)) {
            quad_rule = algoim_interface->get_cut_quadrature(kb);
        }
        if (quad_rule == nullptr || quad_rule->size() == 0) {
            // Not a stored-rule interface, or this cell is uncut there.
            return;
        }
    } else {
        if constexpr (std::is_same_v<std::remove_cvref_t<Phi>, FunFEM<Mesh>>) {
            phi_.setTime(tid);
            phi_.setElementFromBackMesh(kb, 0);
        } else {
            phi_.t = tid;
        }
        // Get quadrature rule - template argument deduction from K and phi_.
        // If regeneration misses a borderline cell that the AlgoimInterface already
        // accepted, fall back to the cached rule stored on the interface.
        generated_quad_rule = quadGenSurf(K, phi_, options_);
        quad_rule = algoim_surface_rule_or_cached<Mesh, Phi>(interface, kb, generated_quad_rule);

        if (quad_rule->size() == 0) {
            std::cout << "Warning: no surface quadrature points for cut element element kb = " << kb << " in AlgoimCutFEMUnified::addInterfaceContribution\n";
            return;
        }
    }

    // Data for the VirtualParameter coefficient lists (CutFEMParameter etc.);
    // previously never evaluated here, silently dropping e.g. the per-domain
    // viscosity eta in Nitsche consistency terms (caught by the MMS in
    // workfiles/src/active_surfaces/algoim_axisymmetric_stokes_mms.cpp).
    // NOTE: the per-side cut measures are not computed in this path (they
    // would need two extra volume rules); measK is passed instead, so
    // MEASURE-WEIGHTED kappa parameters are NOT supported here -- constant
    // kappas and per-domain parameters (the only ones used with algoim
    // drivers) are exact.
    const double h_        = K.get_h();
    const double measK     = K.measure();
    double meas_surf       = 0.0;
    for (const double w : quad_rule->weights) meas_surf += w;
    const std::array<double, 2> measCut{measK, measK};

    for (int l = 0; l < VF.size(); ++l) {

        // if(!VF[l].on(domain)) continue;

        // FINITE ELEMENT SPACES && ELEMENTS
        const fespace_t &Vhv(VF.get_spaceV(l));
        const fespace_t &Vhu(VF.get_spaceU(l));
        bool same = (VF.isRHS() || (&Vhu == &Vhv));

        std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
        std::vector<int> idxU = (same) ? idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());

        // Gustaf: Not having this casues stochastic segfaults, but I am not quite happy with the mathematical implications of adding this
        // std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());

        if (idxV.empty()) {
        #ifdef USE_MPI
            std::cerr << "Rank " << MPIcf::my_rank()
        #else
            std::cerr << "Rank serial"
        #endif
                    << ": empty idxV in addInterfaceContribution"
                    << ", kb = " << kb
                    << ", test domain = " << VF[l].get_domain_test_function()
                    << std::endl;
        #ifdef USE_MPI
            MPIcf::Abort(1);
        #else
            std::abort();
        #endif
        }

        // std::vector<int> idxU = same ? idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());

        if (idxU.empty()) {
        #ifdef USE_MPI
            std::cerr << "Rank " << MPIcf::my_rank()
        #else
            std::cerr << "Rank serial"
        #endif
                    << ": empty idxU in addInterfaceContribution"
                    << ", kb = " << kb
                    << ", trial domain = " << VF[l].get_domain_trial_function()
                    << std::endl;
        #ifdef USE_MPI
            MPIcf::Abort(1);
        #else
            std::abort();
        #endif
        }
        // end of gustaf block

        int kv = VF[l].onWhatElementIsTestFunction(idxV);
        int ku = VF[l].onWhatElementIsTrialFunction(idxU);
        // Gustaf
        if (kv < 0 || kv >= Vhv.get_nb_element()) {
            std::cerr << "Invalid kv = " << kv
                    << ", Vhv.get_nb_element() = " << Vhv.get_nb_element()
                    << ", kb = " << kb << std::endl;
        #ifdef USE_MPI
            MPIcf::Abort(1);
        #else
            std::abort();
        #endif
        }

        if (ku < 0 || ku >= Vhu.get_nb_element()) {
            std::cerr << "Invalid ku = " << ku
                    << ", Vhu.get_nb_element() = " << Vhu.get_nb_element()
                    << ", kb = " << kb << std::endl;
        #ifdef USE_MPI
            MPIcf::Abort(1);
        #else
            std::abort();
        #endif
        }

        const auto &FKu(Vhu[ku]);
        const auto &FKv(Vhv[kv]);
        int domu = FKu.get_domain();
        int domv = FKv.get_domain();
        this->initIndex(FKu, FKv);

        const double coef_param = VF[l].computeCoefInterface(h_, meas_surf, measK, measCut, {domu, domv});

        // BF MEMORY MANAGEMENT -
        int lastop = getLastop(VF[l].du, VF[l].dv);
        RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
        RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
        What_d Fop = Fwhatd(lastop);

        // Loop over quadrature points
        for (size_t ipq = 0; ipq < quad_rule->size(); ++ipq) {
            
            Rd mip = quad_rule->points[ipq];
            const double weight = quad_rule->weights[ipq];
            
            assert(weight > 0);
            
            const Rd face_ip = K.mapToReferenceElement(mip);
            double Cint      = weight * cst_time;

            // Algoim rules carry +grad(phi)/|grad(phi)| normals.  The standard-CutFEM
            // interface assembly (BaseFEM::addInterfaceContribution) uses
            // -interface.normal(ifac), i.e. the OPPOSITE orientation, and every
            // two-phase Nitsche form in the drivers is written for that convention.
            // Negate so both assembly paths agree; without this, all odd-in-n
            // interface terms (Nitsche consistency/adjoint, pressure-normal pairs)
            // get the wrong sign -- an O(1) consistency error at the interface
            // (caught by the polynomial-reproduction MMS in
            // workfiles/src/active_surfaces/algoim_axisymmetric_stokes_mms.cpp).
            const Rd normal = -quad_rule->normals[ipq];

            assert(std::fabs(normal.norm() - 1) < 1e-14);
            double coef = VF[l].computeCoefFromNormal(normal);

            FKv.BF(Fop, face_ip, fv);

            if (!same)
                FKu.BF(Fop, face_ip, fu);

            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domu, domv), mip, tid,
                                                           normal);
            Cint *= coef_param * coef * VF[l].c;

            // std::cout << "ALGOIMCUTFEM: kb = " << kb << ", Cint = " << Cint << std::endl;
            // getchar();
            if (In) {
                if (VF.isRHS())
                    this->addToRHS(VF[l], *In, FKv, fv, Cint);
                else
                    this->addToMatrix(VF[l], *In, FKu, FKv, fu, fv, Cint);
            } else {
                if (VF.isRHS()) {
                    this->addToRHS(VF[l], FKv, fv, Cint);
                } else
                    this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
            }
        }
    }
}

template <typeMesh Mesh, typename Phi>
void AlgoimCutFEMUnified<Mesh, Phi>::addFaceContribution(const itemVFlist_t &VF, const std::pair<int, int> &e1,
                             const std::pair<int, int> &e2, const TimeSlab *In, int itq, double cst_time) {

    
    // CHECK IF IT IS FOR RHS OR MATRIX
    // CONVENTION ki < kj
    bool to_rhs = VF.isRHS();
    int k = e1.first, ifac = e1.second;    
    int kn = e2.first, jfac = e2.second;

    // std::cout << "in addFaceContribution for element " << k << " ifac = " << ifac << std::endl;

    const fespace_t &Vh(VF.get_spaceV(0));
    const ActiveMesh<Mesh> &Th(Vh.get_mesh());
    const fe_element_t &FK(Vh[k]);
    const element_t &K(FK.T);
    int domain = FK.get_domain();

    int thread_id = omp_get_thread_num();
    
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    int kb = Vh.idxElementInBackMesh(k);

    if constexpr (std::is_same_v<std::remove_cvref_t<Phi>, FunFEM<Mesh>>) {
        phi_.setTime(tid);
        phi_.setElementFromBackMesh(kb);
    } else {
        phi_.t = tid;
    }

    auto quad_rule = quadGenFace(K, phi_, options_, ifac, domain);

    // VirtualParameter coefficient lists (e.g. eta-scaled ghost penalties);
    // previously never evaluated in the algoim face path.
    const double h_    = K.get_h();
    const double measK = K.measure();
    double meas_face   = 0.0;
    for (const double w : quad_rule.weights) meas_face += w;

    // Loop over the variational formulation items
    for (int l = 0; l < VF.size(); ++l) {
        if (!VF[l].on(domain))
            continue;

        const double coef_param = VF[l].computeCoefElement(h_, meas_face, measK, measK, domain);

        // FINITE ELEMENT SPACES && ELEMENTS
        const fespace_t &Vhv(VF.get_spaceV(l));
        const fespace_t &Vhu(VF.get_spaceU(l));
        assert(Vhv.get_nb_element() == Vhu.get_nb_element());
        const int kv = VF[l].onWhatElementIsTestFunction(k, kn);
        const int ku = VF[l].onWhatElementIsTrialFunction(k, kn);

        int kbv = Vhv.idxElementInBackMesh(kv);
        int kbu = Vhu.idxElementInBackMesh(ku);

        const fe_element_t &FKu(Vhu[ku]);
        const fe_element_t &FKv(Vhv[kv]);
        this->initIndex(FKu, FKv);

        // BF MEMORY MANAGEMENT -
        bool same  = (VF.isRHS() || (&Vhu == &Vhv && ku == kv));
        int lastop = getLastop(VF[l].du, VF[l].dv);

        // Calculate the offset for the current thread

        long offset = thread_id * this->offset_bf_;

        RNMK_ fv(this->databf_ + offset, FKv.NbDoF(), FKv.N, lastop);
        RNMK_ fu(this->databf_ + offset + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
        What_d Fop = Fwhatd(lastop);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < quad_rule.size(); ++ipq) {

            Rd mip = quad_rule.points[ipq];
            const double weight = quad_rule.weights[ipq];
            double Cint    = weight * cst_time;
            assert(weight > 0);

            const Rd normal(quad_rule.normals[ipq]);
            double coef = VF[l].computeCoefFromNormal(normal);

            assert(fabs(normal.norm() - 1) < 1e-14);

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, FKv.T.mapToReferenceElement(mip), fv);
            if (!same)
                FKu.BF(Fop, FKu.T.mapToReferenceElement(mip), fu);
            //   VF[l].applyFunNL(fu,fv);

            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kbu, kbv), std::make_pair(domain, domain),
                                                            mip, tid, normal);
            Cint *= coef_param * coef * VF[l].c;

            if (In) {
                if (VF.isRHS())
                    this->addToRHS(VF[l], *In, FKv, fv, Cint);
                else
                    this->addToMatrix(VF[l], *In, FKu, FKv, fu, fv, Cint);
            } else {
                if (VF.isRHS())
                    this->addToRHS(VF[l], FKv, fv, Cint);
                else
                    this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
            }
        }
        
    }
}


template <typeMesh Mesh, typename Phi>
void AlgoimCutFEMUnified<Mesh, Phi>::addBorderContribution(
    const itemVFlist_t &, const element_t &, const border_element_t &, int,
    const TimeSlab *, int, double) {
    throw std::logic_error(
        "AlgoimCutFEMUnified cut-boundary integration is not implemented; "
        "refusing the legacy linear get_cut_face/get_cut_part fallback");
}


template <typeMesh Mesh, typename Phi>
void AlgoimCutFEMUnified<Mesh, Phi>::addLagrangeContribution(const itemVFlist_t& VF, const int k, const TimeSlab* In, int itq, double cst_time) {
    const fespace_t& Vh(VF.get_spaceV(0));
    const fe_element_t& FK(Vh[k]);
    const element_t&  K(FK.T);

    const int domain = FK.get_domain();
    const int kb     = Vh.idxElementInBackMesh(k);
    const int iam    = omp_get_thread_num();

    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;
    if constexpr (std::is_same_v<std::remove_cvref_t<Phi>, FunFEM<Mesh>>) {
        phi_.setTime(tid);
        phi_.setElementFromBackMesh(kb);
    } else {
        phi_.t = tid;
    }
    auto quad_rule = quadGenVol(K, phi_, options_, domain);

    if (quad_rule.points.empty()) {
        return;
    }

    const double h_    = K.get_h();
    const double measK = K.measure();
    double meas_cut    = 0.0;
    for (const double w : quad_rule.weights) meas_cut += w;

    for (int l=0; l<VF.size(); ++l) {
        if (!VF[l].on(domain)) continue;

        const double coef_param = VF[l].computeCoefElement(h_, meas_cut, measK, meas_cut, domain);

        const fespace_t& Vhv(VF.get_spaceV(l));
        int kv = k;
        if (&Vhv != &Vh) {
            kv = Vhv.idxElementFromBackMesh(kb, domain);
            if (kv < 0) continue;
        }
        const auto& FKv(Vhv[kv]);
        this->initIndex(FKv, FKv);

        int lastop = getLastop(VF[l].du, VF[l].dv);

        long offset = iam * this->offset_bf_;
        RNMK_ fv(this->databf_ + offset, FKv.NbDoF(), FKv.N, lastop);
        What_d Fop = Fwhatd(lastop);

        for (size_t ipq = 0; ipq < quad_rule.points.size(); ++ipq) {
            const Rd mip = quad_rule.points[ipq];
            const Rd cut_ip = K.mapToReferenceElement(mip);
            double Cint = quad_rule.weights[ipq] * cst_time;

            FKv.BF(Fop, cut_ip, fv);
            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);
            Cint *= coef_param * VF[l].c;

            // A pressure-mean Lagrange item is represented as a linear form, but
            // its integral belongs in the multiplier row and column of the matrix.
            // The one-basis-function overload is the same path used by
            // BaseCutFEM::addLagrangeContribution; addToRHS/addToMatrix(FKu,FKv)
            // would assemble an entirely different problem.
            if (In) {
                this->addToMatrix(VF[l], *In, FKv, fv, Cint);
            } else {
                this->addToMatrix(VF[l], FKv, fv, Cint);
            }
        }
    }
}

// template <typeMesh Mesh, typename Phi>
// template <typename Fct>
// void AlgoimCutFEMUnified<Mesh, Phi>::addInterfaceContributionExact(const Fct &f, const itemVFlist_t& VF, const Interface<mesh_t>& interface, int ifac, double tid,
//                                   const TimeSlab* In, double cst_time, int itq) {


//     phi_.t = tid; // update time in level set function

//     //  GET IDX ELEMENT CONTAINING FACE ON backMes
//     const int kb = interface.idxElementOfFace(ifac);
//     const auto &K(interface.get_element(kb));

//     // Get quadrature rule - template argument deduction from K and phi_
//     auto quad_rule = quadGenSurf(K, phi_, options_);

//     if (quad_rule.points.size() == 0) {
//         std::cout << "Warning: no surface quadrature points for cut element element kb = " << kb << " in AlgoimCutFEMUnified::addInterfaceContribution\n";
//         return;
//     }

//     for (int l = 0; l < VF.size(); ++l) {

//         // if(!VF[l].on(domain)) continue;

//         // FINITE ELEMENT SPACES && ELEMENTS
//         const fespace_t &Vhv(VF.get_spaceV(l));
//         const fespace_t &Vhu(VF.get_spaceU(l));
//         bool same = (VF.isRHS() || (&Vhu == &Vhv));

//         std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
//         std::vector<int> idxU = (same) ? idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());

//         int kv = VF[l].onWhatElementIsTestFunction(idxV);
//         int ku = VF[l].onWhatElementIsTrialFunction(idxU);

//         const auto &FKu(Vhu[ku]);
//         const auto &FKv(Vhv[kv]);
//         int domu = FKu.get_domain();
//         int domv = FKv.get_domain();
//         this->initIndex(FKu, FKv);

//         // BF MEMORY MANAGEMENT -
//         int lastop = getLastop(VF[l].du, VF[l].dv);
//         RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
//         RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
//         What_d Fop = Fwhatd(lastop);

//         // Loop over quadrature points
//         for (size_t ipq = 0; ipq < quad_rule.points.size(); ++ipq) {
            
//             Rd mip = quad_rule.points[ipq];
//             const R weight = quad_rule.weights[ipq];
            
//             assert(weight > 0);
            
//             const Rd face_ip = K.mapToReferenceElement(mip);
//             double Cint      = weight * cst_time;

//             const Rd normal = quad_rule.normals[ipq];

//             assert(std::fabs(normal.norm() - 1) < 1e-14);
//             double coef = VF[l].computeCoefFromNormal(normal);

//             FKv.BF(Fop, face_ip, fv);

//             if (!same)
//                 FKu.BF(Fop, face_ip, fu);

//             Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domu, domv), mip, tid,
//                                                            normal);
//             Cint *= coef * VF[l].c;
//             Cint *= f(mip, VF[l].cv, tid);

//             // std::cout << "ALGOIMCUTFEM: kb = " << kb << ", Cint = " << Cint << std::endl;
//             // getchar();
//             if (In) {
//                 if (VF.isRHS())
//                     this->addToRHS(VF[l], *In, FKv, fv, Cint);
//                 else
//                     this->addToMatrix(VF[l], *In, FKu, FKv, fu, fv, Cint);
//             } else {
//                 if (VF.isRHS()) {
//                     this->addToRHS(VF[l], FKv, fv, Cint);
//                 } else
//                     this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
//             }
//         }
//     }
// }


// // Some integration functions
// template <typeMesh Mesh, typename Phi>
// void AlgoimCutFEMUnified<Mesh, Phi>::addBilinearAlgoim(const itemVFlist_t &VF, const ActiveMesh<mesh_t> &Th) {

// assert(!VF.isRHS());
//     progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
// #pragma omp parallel for num_threads(this->get_num_threads())
//     for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

//         bar += Th.next_element();

//         if (Th.isInactive(k, 0))
//             continue;

//         addElementContribution(VF, k, nullptr, 0, 1.);
        
//         this->addLocalContribution();
//     }
//     bar.end();
// }

// template <typeMesh Mesh, typename Phi>
// void AlgoimCutFEMUnified<Mesh, Phi>::addLinearAlgoim(const itemVFlist_t &VF, const ActiveMesh<mesh_t> &Th) {
// assert(VF.isRHS());
//     progress bar(" Add Linear CutMesh", Th.last_element(), globalVariable::verbose);
// #pragma omp parallel for num_threads(this->get_num_threads())
//     for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        
//         bar += Th.next_element();

//         if (Th.isInactive(k, 0))
//             continue;

//         addElementContribution(VF, k, nullptr, 0, 1.);
        
//     }
//     bar.end();
// }

// template <typeMesh M, typename Phi>
// template <typename Fct>
// void AlgoimCutFEMUnified<M, Phi>::addLinearExact(const Fct &f, const itemVFlist_t &VF, const ActiveMesh<M> &Th,
//                                             const TimeSlab &In) {
//     for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
//         assert(VF.isRHS());
//         auto tq    = this->get_quadrature_time(itq);
//         double tid = In.map(tq);

//         // KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
//         RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
//         In.BF(tq.x, bf_time); // compute time basic funtions
//         double cst_time = tq.a * In.get_measure();

//         for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

//             if (Th.isInactive(k, itq))
//                 continue;

//             addElementContributionExact(f, VF, k, &In, itq, cst_time);
//         }
//     }
// }

// // template <typeMesh M, typename Phi>
// // template <typename Fct>
// // void AlgoimCutFEMUnified<M, Phi>::addLinearExact(const Fct &f, const itemVFlist_t &VF, const ActiveMesh<M> &Th,
// //                                             const TimeSlab &In, const QuadratureFormular1d &qtime) {
// //     for (int itq = 0; itq < qtime.n; ++itq) {
// //         assert(VF.isRHS());
// //         auto tq    = qtime.at(itq);
// //         double tid = In.map(tq);

// //         //KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
// //         RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
// //         In.BF(tq.x, bf_time); // compute time basic funtions
// //         double cst_time = tq.a * In.get_measure();

// //         for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

// //             // if (Th.isInactive(k, itq))
// //             //     continue;

// //             addElementContributionExactSensitive(f, VF, k, &In, itq, cst_time);
// //         }
// //     }
// // }

// template <typeMesh M, typename Phi>
// template <typename Fct>
// void AlgoimCutFEMUnified<M, Phi>::addLinearExact(const Fct &f, const itemVFlist_t &VF, const ActiveMesh<M> &Th,
//                                             const int itq, const TimeSlab &In, const double scaling_time) {

//     assert(VF.isRHS());
//     auto tq    = this->get_quadrature_time(itq);
//     double tid = In.map(tq);
//     KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
//     RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
//     In.BF(tq.x, bf_time);
//     for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

//         if (Th.isInactive(k, itq))
//             continue;

//         addElementContributionExact(f, VF, k, &In, itq, 1.);
//     }
// }

// /**
//  * @brief Add bulk rhs integral in specific time quadrature point and scale with time.
//  *
//  * @tparam M
//  * @tparam L
//  * @tparam Fct
//  * @param f
//  * @param VF
//  * @param Th
//  * @param In
//  */
// template <typeMesh M, typename Phi>
// template <typename Fct>
// void AlgoimCutFEMUnified<M, Phi>::addLinearExact(const Fct &f, const itemVFlist_t &VF, const ActiveMesh<M> &Th,
//                                             const TimeSlab &In, const int itq) {

//     assert(VF.isRHS());
//     auto tq    = this->get_quadrature_time(itq);
//     double tid = In.map(tq);

//     KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
//     RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
//     In.BF(tq.x, bf_time); // compute time basic funtions
//     double cst_time = tq.a * In.get_measure();

//     for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

//         if (Th.isInactive(k, itq))
//             continue;

//         addElementContributionExact(f, VF, k, &In, itq, cst_time);
//     }
// }




// template <typeMesh M, typename Phi>
// template <typename Fct>
// void AlgoimCutFEMUnified<M, Phi>::addLinearExact(const Fct &f, const itemVFlist_t &VF, const TimeInterface<M> &gamma,
//                                             const TimeSlab &In) {
//     assert(VF.isRHS());

//     for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {

//         auto tq    = this->get_quadrature_time(itq);
//         double tid = In.map(tq);

//         KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
//         RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
//         In.BF(tq.x, bf_time); // compute time basic funtions
//         double cst_time = tq.a * In.get_measure();

//         for (int iface = gamma[itq]->first_element(); iface < gamma[itq]->last_element();
//              iface += gamma[itq]->next_element()) {
//             // const typename Interface<M>::Face &face = (*gamma[itq])[iface]; // the face

//             addInterfaceContributionExact(f, VF, *gamma[itq], iface, tid, &In, cst_time, itq);
//         }
//     }
// }

// /**
//  * @brief Add exact RHS in specific time quadrature point, without scaling in time.
//  *
//  * @tparam M
//  * @tparam L
//  * @tparam Fct
//  * @param f
//  * @param VF
//  * @param gamma
//  * @param In
//  * @param itq
//  */
// template <typeMesh M, typename Phi>
// template <typename Fct>
// void AlgoimCutFEMUnified<M, Phi>::addLinearExact(const Fct &f, const itemVFlist_t &VF, const Interface<M> &gamma,
//                                             const TimeSlab &In, const int itq) {
//     assert(VF.isRHS());

//     auto tq    = this->get_quadrature_time(itq);
//     double tid = In.map(tq);

//     KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
//     RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
//     In.BF(tq.x, bf_time); // compute time basic funtions

//     for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
//         const typename Interface<M>::Face &face = gamma[iface]; // the face

//         addInterfaceContributionExact(f, VF, gamma, iface, tid, &In, 1., itq);
//     }
// }

// /**
//  * @brief Add exact RHS in specific time quadrature point, WITH scaling in time.
//  *
//  * @tparam M
//  * @tparam L
//  * @tparam Fct
//  * @param f
//  * @param VF
//  * @param gamma
//  * @param In
//  * @param itq
//  */
// template <typeMesh M, typename Phi>
// template <typename Fct>
// void AlgoimCutFEMUnified<M, Phi>::addLinearExact(const Fct &f, const itemVFlist_t &VF, const TimeInterface<M> &gamma,
//                                             const TimeSlab &In, const int itq) {
//     assert(VF.isRHS());

//     auto tq    = this->get_quadrature_time(itq);
//     double tid = In.map(tq);

//     KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
//     RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
//     In.BF(tq.x, bf_time); // compute time basic funtions

//     double cst_time = tq.a * In.get_measure();

//     for (int iface = gamma[itq]->first_element(); iface < gamma[itq]->last_element();
//          iface += gamma[itq]->next_element()) {
//         const typename Interface<M>::Face &face = (*gamma[itq])[iface]; // the face

//         addInterfaceContributionExact(f, VF, *gamma[itq], iface, tid, &In, cst_time, itq);
//     }
// }
