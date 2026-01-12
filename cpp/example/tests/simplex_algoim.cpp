// #include <iostream>
// #include <filesystem>
// #include "../cutfem.hpp"
// #include "../algoim/quadrature_general.hpp"
// #include "../algoim/quadrature_multipoly.hpp"


// using namespace algoim;

// using mesh_t       = Mesh2;
// using barymesh_t   = BarycentricMesh2;
// using funtest_t    = TestFunction<mesh_t>;
// using fct_t        = FunFEM<mesh_t>;
// using activemesh_t = ActiveMesh<mesh_t>;
// using baryactivemesh_t = BarycentricActiveMesh2;
// using space_t      = GFESpace<mesh_t>;
// using cutspace_t   = CutFESpace<mesh_t>;


// namespace Simplex {

//     std::string problem_name = "simplex";

//     size_t nx = 6; // Number of elements in x direction
//     size_t ny = 6; // Number of elements in y direction

//     double bottom_left_x = 0; // Bottom left corner x coordinate
//     double bottom_left_y = 0; // Bottom left corner y coordinate
//     double width_domain  = 1; // Width of the domain
//     double height_domain = 1; // Height of the domain

//     double Rad = std::sqrt(0.2); // Radius of the drop
//     R2 shift(0.5, 0.5);

//     double levelset(double* P, int i) {
        
//         return (P[0] - shift.x) * (P[0] - shift.x) +
//                 (P[1] - shift.y) * (P[1] - shift.y) - Rad*Rad;
//     }

//     double fun_normal(double* P, int i, int dom) {
//         if (i == 0)
//             return (P[0] - shift.x) / std::sqrt((P[0] - shift.x) * (P[0] - shift.x) + (P[1] - shift.y) * (P[1] - shift.y));
//         else
//             return (P[1] - shift.y) / std::sqrt((P[0] - shift.x) * (P[0] - shift.x) + (P[1] - shift.y) * (P[1] - shift.y));
//     }

//     template <int N> struct Levelset {
    
//         template <typename V> typename V::value_type operator()(const V &P) const {
//             R xc = shift.x, yc = shift.y;
//             return ((P[0] - xc) * (P[0] - xc) + (P[1] - yc) * (P[1] - yc) - Rad*Rad - 1e-14);
//         }

//         // gradient of level set function
//         template <typename T> algoim::uvector<T, N> grad(const algoim::uvector<T, N> &x) const {

//             return algoim::uvector<T, N>(2.0 * (x(0) - shift.x),
//                                         2.0 * (x(1) - shift.y));
//         }

//         // normal = grad(phi)/norm(grad(phi))
//         R2 normal(std::span<double> P) const {
//             R norm = sqrt(pow(2.0 * (P[0] - shift.x), 2) +
//                         pow(2.0 * (P[1] - shift.y), 2));
//             // R normalize = 1. / sqrt(4. * P[0] * P[0] + 4. * P[1] * P[1]);
//             return R2(2.0 * (P[0] - shift.x) / norm,
//                     2.0 * (P[1] - shift.y) / norm);
//         }
//     };

// }


// using namespace Simplex;


// double integral() {


//     auto ellipsoid = [](const uvector<real,3>& x)
//     {
//         return x(0)*x(0) + x(1)*x(1)*4 + x(2)*x(2)*9 - 1;
//     };


//     real data[27]; // or some other memory block elsewhere, on the stack, heap, etc.
//     xarray<real,3> phi(data, uvector<int,3>(3,3,3)); // (3,3,3) indicates the array extent, equal to the Bernstein degree plus one

//     bernstein::bernsteinInterpolate<N>([&](const uvector<real,N>& x) { return ellipsoid(xmin + x * (xmax - xmin)); }, phi);

//     ImplicitPolyQuadrature<3> ipquad(phi);
//     int q = 3; // q nodes per one-dimensional line integral
//     std::vector<uvector<real,N>> phase0, phase1; // stores quadrature nodes for the 'inside' and 'outside'
//     ipquad.integrate(AutoMixed, q, [&](const uvector<real,N>& x, real w)
//     {
//         if (bernstein::evalBernsteinPoly(phi, x) < 0)
//             phase0.push_back(x);
//         else
//             phase1.push_back(x);
//     });

//     std::cout << "Number of quadrature points inside: " << phase0.size() << "\n";
//     std::cout << "Number of quadrature points outside: " << phase1.size()

//                 << "\n";
//     return phase0.size();

// }


// int main() {


//     std::cout << "Integral inside the ellipsoid: " << integral() << "\n";

//     return 0;
// }


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

#include <../algoim/quadrature_multipoly.hpp>
#include <../algoim/bernstein.hpp>
#include <../algoim/uvector.hpp>
#include <../algoim/xarray.hpp>

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using real = algoim::real;
template <int N> using Vec = algoim::uvector<real, N>;

// ------------------------- Example level set -------------------------
// Circle: phi(x,y) = (x-cx)^2 + (y-cy)^2 - R^2
struct CirclePhi
{
    real cx, cy, R;
    real operator()(const Vec<2>& x) const
    {
        const real dx = x(0) - cx;
        const real dy = x(1) - cy;
        return dx * dx + dy * dy - R * R;
    }
};

// ------------------------- CSV helpers -------------------------
static void write_csv(const std::string& path,
                      const std::vector<std::array<real, 2>>& X,
                      const std::vector<real>& W)
{
    std::ofstream os(path);
    os << std::setprecision(17);
    os << "x,y,w\n";
    for (size_t i = 0; i < X.size(); ++i)
        os << X[i][0] << "," << X[i][1] << "," << W[i] << "\n";
}

// ------------------------- Main -------------------------
int main()
{
    // ----- physical element (right isosceles triangle) -----
    const real xmin = 0.0;
    const real xmax = 1.0;
    const real L    = xmax - xmin;

    // triangle vertices: (xmin,xmin), (xmax,xmin), (xmin,xmax)
    const std::array<std::array<real, 2>, 3> tri = {{
        {xmin, xmin},
        {xmax, xmin},
        {xmin, xmax},
    }};

    // ----- choose level set in physical coordinates -----
    CirclePhi phi_phys{/*cx=*/0.55, /*cy=*/0.40, /*R=*/0.35};

    // ----- linear cut selecting the triangle inside the square -----
    auto psi_phys = [&](const Vec<2>& x) -> real
    {
        // psi < 0   => inside triangle
        // psi = 0   => hypotenuse
        const real xi  = (x(0) - xmin) / L;
        const real eta = (x(1) - xmin) / L;
        return xi + eta - real(1);
    };

    // ----- algoim settings -----
    // Bernstein degree P: pick >= polynomial degree of phi composed with mapping.
    // For non-polynomial phi (e.g. circle), you are approximating it by interpolation.
    const int deg = 6;
    const int q1d = 6; // 1D Gauss order used internally

    const algoim::uvector<int, 2> P(deg + 1, deg + 1);

    // ----- build Bernstein polynomials on unit square [0,1]^2 -----
    // x_unit in [0,1]^2, map to physical square: x_phys = xmin + L * x_unit
    std::vector<real> phi_data((deg + 1) * (deg + 1));
    std::vector<real> psi_data((deg + 1) * (deg + 1));

    algoim::xarray<real, 2> phiB(phi_data.data(), P);
    algoim::xarray<real, 2> psiB(psi_data.data(), P);

    algoim::bernstein::bernsteinInterpolate<2>(
        [&](const Vec<2>& x_unit) -> real {
            Vec<2> x_phys(xmin + L * x_unit(0), xmin + L * x_unit(1));
            return phi_phys(x_phys);
        },
        phiB);

    // psi is linear, but we still interpolate it (degree deg is fine)
    algoim::bernstein::bernsteinInterpolate<2>(
        [&](const Vec<2>& x_unit) -> real {
            Vec<2> x_phys(xmin + L * x_unit(0), xmin + L * x_unit(1));
            return psi_phys(x_phys);
        },
        psiB);

    // ----- build multipoly quadrature hierarchy -----
    algoim::ImplicitPolyQuadrature<2> ipquad(phiB, psiB);

    // ----- results to output -----
    real I_vol  = 0;
    real I_surf = 0;

    std::vector<std::array<real, 2>> X_vol, X_surf;
    std::vector<real> W_vol, W_surf;

    auto map_to_phys = [&](const Vec<2>& x_unit) -> Vec<2> {
        return Vec<2>(xmin + L * x_unit(0), xmin + L * x_unit(1));
    };

    // ----- VOLUME: integrate over {phi<0} ∩ {psi<0} -----
    ipquad.integrate(algoim::AutoMixed, q1d, [&](const Vec<2>& x_unit, real w_unit) {
        Vec<2> x = map_to_phys(x_unit);

        const real ph = phi_phys(x);
        const real ps = psi_phys(x);

        if (ph < 0 && ps < 0) {
            const real w_phys = w_unit * (L * L); // area scaling (uniform)
            X_vol.push_back({x(0), x(1)});
            W_vol.push_back(w_phys);

            // volume integrand: change this as needed
            const real f = 1.0;
            I_vol += w_phys * f;
        }
    });

    // ----- SURFACE: integrate_surf returns union of interfaces (phi=0 and psi=0).
    // Keep only points on phi=0 (and inside triangle psi<0), reject psi=0.
    const real tol_phi = 1e-10;  // classification tolerance
    const real tol_psi = 1e-10;

    ipquad.integrate_surf(algoim::AutoMixed, q1d,
                          [&](const Vec<2>& x_unit, real w_unit, const Vec<2>& /*wn*/) {
        Vec<2> x = map_to_phys(x_unit);

        const real ph = phi_phys(x);
        const real ps = psi_phys(x);

        const bool on_phi = std::abs(ph) < tol_phi;
        const bool on_psi = std::abs(ps) < tol_psi;

        if (on_phi && !on_psi && ps < 0) {
            const real w_phys = w_unit * L; // arc-length scaling (uniform)
            X_surf.push_back({x(0), x(1)});
            W_surf.push_back(w_phys);

            // surface integrand: change this as needed
            const real g = 1.0;
            I_surf += w_phys * g;
        }
    });

    // ----- print summary -----
    std::cout << std::setprecision(17);
    std::cout << "Triangle vertices:\n";
    for (int i = 0; i < 3; ++i)
        std::cout << "  V" << i << " = (" << tri[i][0] << ", " << tri[i][1] << ")\n";

    std::cout << "\nResults:\n";
    std::cout << "  I_vol  = " << I_vol << "   (∫_{K∩{phi<0}} 1 dx)\n";
    std::cout << "  I_surf = " << I_surf << "   (∫_{K∩{phi=0}} 1 ds)\n";

    std::cout << "\nQuadrature sizes:\n";
    std::cout << "  volume nodes  = " << X_vol.size() << "\n";
    std::cout << "  surface nodes = " << X_surf.size() << "\n";

    // ----- dump points/weights -----
    write_csv("vol_quad.csv", X_vol, W_vol);
    write_csv("surf_quad.csv", X_surf, W_surf);

    std::cout << "\nWrote:\n  vol_quad.csv\n  surf_quad.csv\n";
    return 0;
}
