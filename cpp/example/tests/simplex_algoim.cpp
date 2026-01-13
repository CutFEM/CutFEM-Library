
// main_cut_tri_multipoly.cpp
//
// Cut-triangle quadrature using algoim's multivariate polynomial routines (multipoly):
//  - volume integral over K ∩ {phi < 0}
//  - surface integral over K ∩ {phi = 0}
//  - exports quadrature points + weights (physical) for both volume and surface
//
// Geometry assumption (your simplification):
//  - right isosceles triangle aligned with axes:
//      V0 = (xmin, xmin), V1 = (xmax, xmin), V2 = (xmin, xmax)
//  - embed in the square [xmin,xmax]^2 and add linear cut
//      psi(x,y) = (x-xmin)/L + (y-xmin)/L - 1   (L = xmax-xmin)
//    so psi < 0 selects the triangle, psi = 0 is the hypotenuse.
//
// Build:
//   clang++ -O2 -std=c++17 -I/path/to/algoim main_cut_tri_multipoly.cpp -o cuttri
//
// Note: If you compile with C++20, you may need the xarraySlice constructor patch discussed earlier.

// #include <../algoim/quadrature_multipoly.hpp>
// #include <../algoim/bernstein.hpp>
// #include <../algoim/uvector.hpp>
// #include <../algoim/xarray.hpp>

#include <../algoim/cut_triangle_quadrature.hpp>

#include "../cutfem.hpp"
// #include "finiteElement.hpp"
// // #include "GenericMapping.hpp"
// #include "levelSet.hpp"
// #include "baseProblem.hpp"

// #include <array>
// #include <cmath>
// #include <fstream>
// #include <iomanip>
// #include <iostream>
// #include <string>
// #include <vector>

/*
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
*/


/*int main() {

    std::array<cuttri_multipoly::Vec2,3> tri_verts = {
        // cuttri_multipoly::Vec2(0.0, 0.0),
        // cuttri_multipoly::Vec2(1.0, 0.0),
        // cuttri_multipoly::Vec2(0.0, 1.0)
        cuttri_multipoly::Vec2(0.0, 0.0),
        cuttri_multipoly::Vec2(0.33, 0.33),
        cuttri_multipoly::Vec2(0.0, 1.0)
    };

    // Level set: circle of radius 0.6 centered at (0.5,0.5)
    struct PhiCircle {
        cuttri_multipoly::real operator()(const cuttri_multipoly::Vec2& x) const {
            const cuttri_multipoly::real xc = 0.5;
            const cuttri_multipoly::real yc = 0.5;
            const cuttri_multipoly::real r  = 0.6;
            const cuttri_multipoly::real dx = x(0) - xc;
            const cuttri_multipoly::real dy = x(1) - yc;
            return (dx*dx + dy*dy) - r*r;   //! needs to be a polynomial
        }
    } phi_circle;

    const int bernstein_deg = 8;
    const int q1d = 12;
    const auto quad_scheme = cuttri_multipoly::compute_cut_triangle_quadrature(
        tri_verts,
        phi_circle,
        bernstein_deg,
        q1d,
        algoim::AutoMixed
    );
    // Output volume quadrature points + weights
    std::ofstream fout_vol("cuttri_multipoly_vol_quad.txt");
    fout_vol << std::setprecision(12);
    for (size_t i=0; i<quad_scheme.x_vol.size(); ++i) {
        fout_vol << quad_scheme.x_vol[i][0] << " "
                 << quad_scheme.x_vol[i][1] << " "
                 << quad_scheme.w_vol[i] << "\n";
    }
    fout_vol.close();

    // Output surface quadrature points + weights
    std::ofstream fout_surf("cuttri_multipoly_surf_quad.txt");
    fout_surf << std::setprecision(12);
    for (size_t i=0; i<quad_scheme.x_surf.size(); ++i) {
        fout_surf << quad_scheme.x_surf[i][0] << " "
                  << quad_scheme.x_surf[i][1] << " "
                  << quad_scheme.w_surf[i] << "\n";
    }
    fout_surf.close();



    return 0;
}*/


using mesh_t       = Mesh2;
using barymesh_t   = BarycentricMesh2;
using funtest_t    = TestFunction<mesh_t>;
using fct_t        = FunFEM<mesh_t>;
using activemesh_t = ActiveMesh<mesh_t>;
using baryactivemesh_t = BarycentricActiveMesh2;
using space_t      = GFESpace<mesh_t>;
using cutspace_t   = CutFESpace<mesh_t>;
using fe_element_t = space_t::FElement;
using Rd           = fe_element_t::Rd;

template <typename L>
double integrate(const activemesh_t &Th, const fct_t &fh, L &phi) {
    
    double val = 0.0;
    int domain = 0;
    int cu = 0; // component index
    int op = 0; // operation index

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (Th.isInactive(k, 0))
            continue;

        int kb = Th.idxElementInBackMesh(k);

        const auto &T(Th.Th[kb]);   // get element in background mesh
        const auto &V0(T.at(0)); // vertex 0
        const auto &V1(T.at(1)); // vertex 1
        const auto &V2(T.at(2)); // vertex 2
        
        const std::array<algoim::uvector<double, 2>, 3> triangle{ 
            algoim::uvector<double, 2>{V0[0], V0[1]},
            algoim::uvector<double, 2>{V1[0], V1[1]},
            algoim::uvector<double, 2>{V2[0], V2[1]}
        };
        
        const auto q = cuttri_multipoly::compute_cut_triangle_quadrature(
            triangle,
            phi,
            8,   // bernstein degree
            12   // algoim 1D quadrature order
        );
        for (int ipq = 0; ipq < q.x_vol.size(); ++ipq) {

            const Rd mip(q.x_vol.at(ipq)[0], q.x_vol.at(ipq)[1]);
            const R weight = q.w_vol.at(ipq);

            val += weight * fh.evalOnBackMesh(kb, domain, mip, cu, op);
        }

        
    }
    return val;

}


struct PhiCircle {
    cuttri_multipoly::real operator()(const cuttri_multipoly::Vec2& x) const {
        const cuttri_multipoly::real xc = 0.5;
        const cuttri_multipoly::real yc = 0.5;
        const cuttri_multipoly::real r  = 0.6;
        const cuttri_multipoly::real dx = x(0) - xc;
        const cuttri_multipoly::real dy = x(1) - yc;
        //return (dx*dx + dy*dy) - r*r;   
        return std::sqrt(dx*dx + dy*dy) - r;
    }
};

double fun_phi(double *P, const int comp) {
    const double x = P[0];
    const double y = P[1];

    const double xc = 0.5;
    const double yc = 0.5;
    const double r  = 0.6;
    const double dx = x - xc;
    const double dy = y - yc;
    // return (dx*dx + dy*dy) - r*r; 
    return std::sqrt(dx*dx + dy*dy) - r;
}


double fun_test(double *P, const int comp) {
    const double x = P[0];
    const double y = P[1];

    return 1.;
}

int main() {

    int nx = 10;
    int ny = 10;
    double bottom_left_x = -1.0;
    double bottom_left_y = -1.0;
    double width_domain  = 2.5;
    double height_domain = 2.5;
    
    mesh_t Th(nx
            , ny
            , bottom_left_x
            , bottom_left_y
            , width_domain
            , height_domain);

    // barymesh_t Th(nx
    //         , ny
    //         , bottom_left_x
    //         , bottom_left_y
    //         , width_domain
    //         , height_domain);

            
    PhiCircle phi_circle;

    space_t Lh(Th, DataFE<mesh_t>::P1);
    fct_t phi_lin(Lh, fun_phi);
    InterfaceLevelSet<mesh_t> surface(Th, phi_lin);
    activemesh_t Th_active(Th);
    Th_active.truncate(surface, 1);

    space_t V_test(Th, DataFE<mesh_t>::P4);   
    fct_t test_fun(V_test, fun_test);


    double volV = integrate(Th_active, test_fun, phi_circle);

    // double surfA = integrate_surface(Th_active, test_fun, phi_lin, phi_circle);

    std::cout << "Integrated volume: " << volV << "\n";

    Paraview<mesh_t> paraview(Th, "simplex_algoim.vtk");
    paraview.add(phi_lin, "phi_lin", 0, 1);
    paraview.writeActiveMesh(Th_active, "active_mesh.vtk");

    
    return 0;
}