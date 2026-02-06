#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "quadrature_multipoly.hpp"

namespace cuttri_multipoly {

using real = algoim::real;
using Vec2 = algoim::uvector<real,2>;

inline real norm2(const Vec2& v) { return std::sqrt(v(0)*v(0) + v(1)*v(1)); }
inline Vec2 perp(const Vec2& n)  { return Vec2(-n(1), n(0)); }

struct QuadScheme
{
    std::vector<std::array<real,2>> x_vol;
    std::vector<real>              w_vol;

    std::vector<std::array<real,2>> x_surf;
    std::vector<real>               w_surf;
};

// Affine map x = v0 + ξ e1 + η e2
struct AffineTriMap
{
    Vec2 v0, e1, e2;

    Vec2 map(const Vec2& xi) const { return v0 + e1*xi(0) + e2*xi(1); }

    real detJ() const { return e1(0)*e2(1) - e1(1)*e2(0); }

    Vec2 J_mul(const Vec2& t) const { return e1*t(0) + e2*t(1); }
};

// Phi functor requirement: real operator()(Vec2 x) const;
template <class Phi>
inline QuadScheme compute_cut_triangle_quadrature(
    const std::array<Vec2,3>& verts, // physical triangle vertices
    const Phi& phi,                  // physical level set
    int bernstein_deg,               // interpolation degree (>=1)
    int q1d,                         // algoim 1D quadrature order
    algoim::QuadStrategy strat = algoim::AutoMixed,
    real tol_phi = 1e-10,
    real tol_psi = 1e-12
)
{
    AffineTriMap F;
    F.v0 = verts[0];
    F.e1 = verts[1] - verts[0];
    F.e2 = verts[2] - verts[0];

    const real detJ = F.detJ();
    const real abs_detJ = std::abs(detJ);

    // Unit-square coords (ξ,η) ∈ [0,1]^2; triangle is ψ(ξ,η)=ξ+η-1 <= 0
    auto psi_ref = [&](const Vec2& xi) -> real { return xi(0) + xi(1) - real(1); };

    const int n = bernstein_deg + 1;
    algoim::uvector<int,2> P(n,n);

    std::vector<real> phi_data(n*n), psi_data(n*n);
    algoim::xarray<real,2> phiB(phi_data.data(), P);
    algoim::xarray<real,2> psiB(psi_data.data(), P);

    // Interpolate φ∘F on unit square
    algoim::bernstein::bernsteinInterpolate<2>(
        [&](const Vec2& xi) -> real { return phi(F.map(xi)); },
        phiB
    );

    // Interpolate ψ (linear)
    algoim::bernstein::bernsteinInterpolate<2>(
        [&](const Vec2& xi) -> real { return psi_ref(xi); },
        psiB
    );

    algoim::ImplicitPolyQuadrature<2> ipquad(phiB, psiB);

    QuadScheme out;

    // -------- Volume nodes: keep phi<0 and inside simplex (psi<0) --------
    ipquad.integrate(strat, q1d, [&](const Vec2& xi, real w_ref) {
        const real ps = psi_ref(xi);
        if (ps >= 0) return;

        const Vec2 x = F.map(xi);
        //const real ph = phi(x);   // before: this needs phi to be a polynomial
        const real phB = algoim::bernstein::evalBernsteinPoly(phiB, xi);    
        if (phB >= 0) return;

        const real w_phys = w_ref * abs_detJ;
        out.x_vol.push_back({x(0), x(1)});
        out.w_vol.push_back(w_phys);
    });

    // -------- Surface nodes: ipquad gives union of interfaces (phi=0 and psi=0) --------
    // We keep phi=0, reject psi=0, and require psi<0.
    ipquad.integrate_surf(strat, q1d, 
    [&](const Vec2& xi, real w_ref, const Vec2& wn_ref)
    {
        const real ps = psi_ref(xi);    // psi is linear, so exact
        if (ps >= 0) return;                // outside triangle

        //const real ph = phi(x);   // before: this needs phi to be a polynomial
        const real phB = algoim::bernstein::evalBernsteinPoly(phiB, xi);    // evaluate phi using its Bernstein poly, since phi may not always be a polynomial

        const bool on_phi = std::abs(phB) < tol_phi;  // robust: same polynomial as ipquad
        const bool on_psi = std::abs(ps)  < tol_psi;  // cheap & exact
        if (!on_phi || on_psi) return;

        const Vec2 x = F.map(xi);

        Vec2 n_ref = wn_ref;
        const real nrm = norm2(n_ref);
        if (nrm == 0) return;
        n_ref /= nrm;

        Vec2 t_ref = perp(n_ref);
        const real scale = norm2(F.J_mul(t_ref));

        const real w_phys = w_ref * scale;
        out.x_surf.push_back({x(0), x(1)});
        out.w_surf.push_back(w_phys);
    });

    return out;
}

} // namespace cuttri_multipoly
