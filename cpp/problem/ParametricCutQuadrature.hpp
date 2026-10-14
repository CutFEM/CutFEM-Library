#ifndef PARAMETRIC_CUT_QUADRATURE_HPP
#define PARAMETRIC_CUT_QUADRATURE_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

// Cut quadrature on triangles from a polynomial reconstruction of the interface.
//
// On a triangle whose boundary the level set crosses exactly twice, {phi = 0}
// is replaced by a graph over the chord between the two edge roots,
//     gamma(s) = p0 + s (p1 - p0) + alpha(s) n_c,      alpha in P_r,
// where alpha interpolates the roots of phi along the chord normal n_c at the
// interior Chebyshev-Lobatto points.  Each side of the cut triangle is then
// bounded by straight pieces of the triangle boundary and by gamma.
//
//   boundary rule: Gauss-Legendre in the parameter of every piece, with vector
//                  weights W = w gamma'^perp;
//   volume rule:   every boundary node y is joined to a star point a, with
//                  Gauss nodes t on [0,1], points a + t (y - a) and weights
//                  w t (y - a).W  (scaled boundary cubature, Chin & Sukumar,
//                  CMAME 2021).
//
// Let D = ProblemOption::algoim_ibp_degree_.  For F in [P_D]^2 the integrand
// F(gamma).gamma'^perp is a polynomial in s, so the boundary rule is exact,
// and the volume rule applied to div F equals the boundary rule applied to a
// polynomial field with the same divergence.  Hence, on each side of every
// triangle and up to rounding,
//     sum_volume w div F = sum_boundary F.W        for all F in [P_D]^2,
// however the interface cuts the triangle.  A triangle that is not crossed
// exactly twice, whose arc is not a graph over its chord, or whose sides are
// not star-shaped about a star point is split into its four midpoint
// children; star-shapedness is what makes the volume weights positive.  Edge
// roots are never discarded and children inherit the roots of their parent's
// edges, so neighbouring triangles agree on the inside part of a shared edge
// up to the rounding of their two level-set evaluations.
namespace parametric_cut {

using R2 = Mesh2::Rd;

inline double dot(const R2 &a, const R2 &b) { return a[0] * b[0] + a[1] * b[1]; }
inline double cross(const R2 &a, const R2 &b) { return a[0] * b[1] - a[1] * b[0]; }
inline bool lex_less(const R2 &a, const R2 &b) { return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]); }

template <typename Phi> double evaluate(Phi &phi, const R2 &x) {
    return double(phi(algoim::uvector<algoim::real, 2>(x[0], x[1])));
}

// Roots of phi on the open segment (A, B), ordered from A to B.  The segment
// is traversed in a direction fixed by its end points alone, so that two
// triangles sharing it see the same trace.  The Bernstein interpolant of
// degree `degree` (exact for a level set in P_degree) only locates the
// stationary points of the trace; phi itself is bisected between them.
template <typename Phi>
std::vector<R2> roots_on_segment(Phi &phi, const R2 &A, const R2 &B, int degree) {
    using real = algoim::real;
    const bool flip = lex_less(B, A);
    const R2 P = flip ? B : A, Q = flip ? A : B;
    auto g = [&](double t) { return evaluate(phi, (1.0 - t) * P + t * Q); };

    std::vector<double> knots{0.0, 1.0};
    if (degree >= 2) {
        std::vector<real> c(degree + 1), dc(degree), stationary(degree);
        algoim::xarray<real, 1> trace(c.data(), algoim::uvector<int, 1>(degree + 1));
        algoim::bernstein::bernsteinInterpolate<1>(
            [&](const algoim::uvector<real, 1> &t) -> real { return g(double(t(0))); }, trace);
        for (int i = 0; i < degree; ++i) dc[i] = degree * (c[i + 1] - c[i]);
        const int n = algoim::bernstein::bernsteinUnitIntervalRealRoots(dc.data(), degree, stationary.data());
        for (int i = 0; i < n; ++i)
            if (stationary[i] > 0.0 && stationary[i] < 1.0) knots.push_back(double(stationary[i]));
        std::sort(knots.begin(), knots.end());
    }

    std::vector<R2> roots;
    double a = knots[0], ga = g(a);
    for (size_t i = 1; i < knots.size(); ++i) {
        const double b = knots[i], gb = g(b);
        if ((ga < 0) != (gb < 0)) {
            double lo = a, hi = b; // g(lo) has the sign of ga
            for (int it = 0; it < 60; ++it) {
                const double mid = 0.5 * (lo + hi);
                ((g(mid) < 0) == (ga < 0) ? lo : hi) = mid;
            }
            const double t = 0.5 * (lo + hi);
            roots.push_back((1.0 - t) * P + t * Q);
        }
        a  = b;
        ga = gb;
    }
    if (flip) std::reverse(roots.begin(), roots.end());
    return roots;
}

// Gauss nodes and vector weights of the degree-r arc from `from` to `to`.  The
// weights w gamma'^perp point to the right of the direction of travel.  The
// arc is built in an orientation fixed by its end points, so that both sides
// of the interface use the same nodes.  Returns false unless every chord-normal
// line through an interior interpolation point crosses the level set exactly
// once inside the triangle V.  Pass r = 1 for the chord itself.
template <typename Phi>
bool build_arc(Phi &phi, const std::array<R2, 3> &V, const R2 &from, const R2 &to, int r, int D,
               int degree, std::vector<R2> &points, std::vector<R2> &weights) {
    const bool flip = lex_less(to, from);
    const R2 p0 = flip ? to : from, chord = (flip ? from : to) - p0;
    const double L = std::sqrt(dot(chord, chord));
    const R2 n(chord[1] / L, -chord[0] / L);

    std::vector<double> s(r + 1), alpha(r + 1, 0.0);
    for (int i = 0; i <= r; ++i) s[i] = 0.5 * (1.0 - std::cos(M_PI * i / r));
    for (int i = 1; i < r; ++i) {
        const R2 c = p0 + s[i] * chord;
        // part [lo, hi] of the line c + a n inside the triangle
        double lo = 0.0, hi = 0.0;
        for (int e = 0; e < 3; ++e) {
            const R2 edge = V[(e + 1) % 3] - V[e], to_edge = V[e] - c;
            const double det = cross(n, edge);
            if (std::abs(det) < 1e-12 * std::sqrt(dot(edge, edge))) continue;
            const double beta = cross(to_edge, n) / det;
            if (beta < -1e-12 || beta > 1.0 + 1e-12) continue;
            const double a = cross(to_edge, edge) / det;
            lo             = std::min(lo, a);
            hi             = std::max(hi, a);
        }
        if (!(hi > lo)) return false;
        const auto roots = roots_on_segment(phi, c + lo * n, c + hi * n, degree);
        if (roots.size() != 1) return false;
        alpha[i] = dot(roots[0] - c, n);
    }

    const int ns = (r * (D + 1) + 1) / 2;
    for (int q = 0; q < ns; ++q) {
        const double x = algoim::GaussQuad::x(ns, q), w = algoim::GaussQuad::w(ns, q);
        double a = 0.0, da = 0.0; // alpha(x) and alpha'(x) in Lagrange form
        for (int i = 1; i < r; ++i) {
            double li = 1.0, dli = 0.0;
            for (int j = 0; j <= r; ++j) {
                if (j == i) continue;
                dli = (dli * (x - s[j]) + li) / (s[i] - s[j]);
                li *= (x - s[j]) / (s[i] - s[j]);
            }
            a  += alpha[i] * li;
            da += alpha[i] * dli;
        }
        const R2 tangent = chord + da * n;
        points.push_back(p0 + x * chord + a * n);
        weights.push_back((flip ? -w : w) * R2(tangent[1], -tangent[0]));
    }
    return true;
}

// One accepted triangle.  Index 0 is its part in {phi < 0}, index 1 in {phi > 0}.
struct Leaf {
    std::array<std::vector<R2>, 2> boundary; // counter-clockwise pieces of the triangle boundary
    std::array<R2, 2> star;                  // star points of the two parts
    std::vector<R2> arc_points, arc_weights; // arc closing boundary[0]; weights outward for {phi < 0}
};

// `roots[i]` holds the roots on the edge from V[i] to V[i+1], in that order.
template <typename Phi>
void collect(Phi &phi, const std::array<R2, 3> &V, const std::array<std::vector<R2>, 3> &roots,
             const ProblemOption &option, int depth, std::vector<Leaf> &leaves, int &subdivisions) {
    const int degree = option.algoim_bernstein_deg_;

    // vertices and edge roots in counter-clockwise order
    std::vector<R2> ring;
    std::vector<size_t> root_at;
    double h2 = 0.0;
    for (int i = 0; i < 3; ++i) {
        const R2 edge = V[(i + 1) % 3] - V[i];
        h2            = std::max(h2, dot(edge, edge));
        ring.push_back(V[i]);
        for (const R2 &x : roots[i]) {
            root_at.push_back(ring.size());
            ring.push_back(x);
        }
    }

    Leaf leaf;
    bool accepted   = false;
    auto whole_side = [&](int side) {
        leaf                = Leaf();
        leaf.boundary[side] = {V[0], V[1], V[2], V[0]};
        leaf.star[side]     = (V[0] + V[1] + V[2]) / 3.0;
        accepted            = true;
    };
    if (root_at.empty()) {
        whole_side(evaluate(phi, V[0]) < 0 ? 0 : 1);
    } else if (root_at.size() == 2) {
        // the two chains of the ring between the roots; each part owns one
        std::array<std::vector<R2>, 2> chain;
        chain[0].assign(ring.begin() + root_at[0], ring.begin() + root_at[1] + 1);
        chain[1].assign(ring.begin() + root_at[1], ring.end());
        chain[1].insert(chain[1].end(), ring.begin(), ring.begin() + root_at[0] + 1);
        const int with_vertex = chain[0].size() > 2 ? 0 : 1;
        const int negative    = evaluate(phi, chain[with_vertex][1]) < 0 ? with_vertex : 1 - with_vertex;
        leaf.boundary[0]      = chain[negative];
        leaf.boundary[1]      = chain[1 - negative];

        const R2 from = leaf.boundary[0].back(), to = leaf.boundary[0].front();
        const double chord2 = dot(to - from, to - from);
        // A chord below 1e-8 h deviates from the level set by less than
        // rounding, so the arc is the chord and both parts are convex.
        const int r = chord2 < 1e-16 * h2 ? 1 : option.parametric_cut_arc_degree_;
        if (chord2 == 0.0) {
            // the two roots coincide: the level set only touches the triangle
            whole_side(leaf.boundary[0].size() > leaf.boundary[1].size() ? 0 : 1);
        } else if (build_arc(phi, V, from, to, r, option.algoim_ibp_degree_, degree, leaf.arc_points,
                             leaf.arc_weights)) {
            // Star point of each part: the mean of its boundary points or, where
            // the arc bends so far into the part that it is star-shaped only
            // about its own boundary, one of those points.
            accepted = true;
            for (int side = 0; side < 2 && accepted; ++side) {
                R2 mean(0.0, 0.0);
                for (const R2 &x : leaf.boundary[side]) mean += x;
                std::vector<R2> candidates{mean / double(leaf.boundary[side].size())};
                candidates.insert(candidates.end(), leaf.boundary[side].begin(), leaf.boundary[side].end());

                accepted = false;
                for (const R2 &a : candidates) {
                    bool star_shaped = true; // a straight arc bounds a convex part
                    for (size_t q = 0; r > 1 && q < leaf.arc_points.size() && star_shaped; ++q) {
                        const double height = dot(leaf.arc_points[q] - a, leaf.arc_weights[q]);
                        star_shaped         = (side == 0 ? height : -height) > 0;
                    }
                    if (star_shaped) {
                        leaf.star[side] = a;
                        accepted        = true;
                        break;
                    }
                }
            }
        }
    }
    if (accepted) {
        leaves.push_back(std::move(leaf));
        return;
    }

    if (depth == 0) {
        std::ostringstream message;
        message << "parametric_cut: interface not resolved on the triangle (" << V[0] << "), (" << V[1]
                << "), (" << V[2] << "); raise algoim_subdivision_depth_ or refine the mesh";
        throw std::runtime_error(message.str());
    }
    ++subdivisions;
    // The roots on the two halves of an edge are those already found on the
    // edge; only the three edges between the midpoints are new.
    const std::array<R2, 3> M = {0.5 * (V[0] + V[1]), 0.5 * (V[1] + V[2]), 0.5 * (V[2] + V[0])};
    std::array<std::vector<R2>, 3> first, second, inner, reversed;
    for (int i = 0; i < 3; ++i) {
        const R2 edge = V[(i + 1) % 3] - V[i];
        for (const R2 &x : roots[i]) (2.0 * dot(x - V[i], edge) < dot(edge, edge) ? first[i] : second[i]).push_back(x);
        inner[i] = roots_on_segment(phi, M[i], M[(i + 1) % 3], degree);
        reversed[i].assign(inner[i].rbegin(), inner[i].rend());
    }
    collect(phi, {V[0], M[0], M[2]}, {first[0], reversed[2], second[2]}, option, depth - 1, leaves, subdivisions);
    collect(phi, {M[0], V[1], M[1]}, {second[0], first[1], reversed[0]}, option, depth - 1, leaves, subdivisions);
    collect(phi, {M[2], M[1], V[2]}, {reversed[1], second[1], first[2]}, option, depth - 1, leaves, subdivisions);
    collect(phi, {M[0], M[1], M[2]}, inner, option, depth - 1, leaves, subdivisions);
}

template <typename Phi>
std::vector<Leaf> leaves_of(const Mesh2::Element &K, Phi &phi, const ProblemOption &option, int &subdivisions) {
    std::array<R2, 3> V = {R2(K.at(0)), R2(K.at(1)), R2(K.at(2))};
    if (cross(V[1] - V[0], V[2] - V[0]) < 0) std::swap(V[1], V[2]);
    std::array<std::vector<R2>, 3> roots;
    for (int i = 0; i < 3; ++i) roots[i] = roots_on_segment(phi, V[i], V[(i + 1) % 3], option.algoim_bernstein_deg_);
    std::vector<Leaf> leaves;
    collect(phi, V, roots, option, option.algoim_subdivision_depth_, leaves, subdivisions);
    return leaves;
}

// Rule on K cap {phi = 0}; normals point out of {phi < 0}.
template <typename Phi>
AlgoimQuadratureRule<Mesh2> surface(const Mesh2::Element &K, Phi &phi, const ProblemOption &option) {
    AlgoimQuadratureRule<Mesh2> rule;
    for (const Leaf &leaf : leaves_of(K, phi, option, rule.ibp_subdivisions))
        for (size_t q = 0; q < leaf.arc_points.size(); ++q) {
            const R2 &W    = leaf.arc_weights[q];
            const double w = std::sqrt(dot(W, W));
            rule.points.push_back(leaf.arc_points[q]);
            rule.weights.push_back(w);
            rule.normals.push_back(W / w);
        }
    return rule;
}

// Rule on K cap {phi > 0} (domain 0) or K cap {phi < 0} (domain 1).
template <typename Phi>
AlgoimQuadratureRule<Mesh2> volume(const Mesh2::Element &K, Phi &phi, const ProblemOption &option, int domain) {
    const int side = domain == 0 ? 1 : 0;
    const int nt   = (option.algoim_ibp_degree_ + 2) / 2;

    AlgoimQuadratureRule<Mesh2> rule;
    for (const Leaf &leaf : leaves_of(K, phi, option, rule.ibp_subdivisions)) {
        const R2 &a = leaf.star[side];
        // Segment from a to the boundary node y; height = (y - a).W, which is
        // not positive only at rounding level or for a piece in line with a.
        auto cone = [&](const R2 &y, double height) {
            if (!(height > 0)) return;
            for (int j = 0; j < nt; ++j) {
                const double t = algoim::GaussQuad::x(nt, j);
                rule.points.push_back(a + t * (y - a));
                rule.weights.push_back(algoim::GaussQuad::w(nt, j) * t * height);
            }
        };
        const auto &boundary = leaf.boundary[side];
        for (size_t i = 0; i + 1 < boundary.size(); ++i) {
            const R2 d = boundary[i + 1] - boundary[i];
            // (y - a).W is constant along a straight piece
            const double height = dot(boundary[i] - a, R2(d[1], -d[0]));
            for (int q = 0; q < nt; ++q)
                cone(boundary[i] + algoim::GaussQuad::x(nt, q) * d, algoim::GaussQuad::w(nt, q) * height);
        }
        for (size_t q = 0; q < leaf.arc_points.size(); ++q) {
            const double height = dot(leaf.arc_points[q] - a, leaf.arc_weights[q]);
            cone(leaf.arc_points[q], side == 0 ? height : -height);
        }
    }
    return rule;
}

} // namespace parametric_cut

#endif
