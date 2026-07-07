
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

    ipquad.integrate(algoim::AutoMixed, q1d, [&](const Vec2& xi, real w_ref)
    {
        // pick component: inside triangle AND phi<0
        if (psi_ref(xi) >= 0) return; // exact
        if (algoim::bernstein::evalBernsteinPoly(phiB, xi) >= 0) return;
        // if (phi(F(xi)) >= 0) return;

        Vec2 x = F(xi);
        rule.points.emplace_back(R2(x(0), x(1)));
        rule.weights.emplace_back(double(w_ref) * std::abs(detJ));
    });

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

    ipquad.integrate_surf(algoim::AutoMixed, q1d,
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
        const auto& FKv(Vhv[k]);
        const auto& FKu(Vhu[k]);
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


// Gustaf
template <typeMesh Mesh, typename Phi>
void AlgoimCutFEMUnified<Mesh, Phi>::addBilinear(
    const itemVFlist_t& VF,
    const Interface<mesh_t>& interface
) {
    assert(!VF.isRHS());

    progress bar("Add Bilinear Algoim Interface",
                 interface.last_element(),
                 globalVariable::verbose);

#pragma omp parallel for num_threads(this->get_num_threads())
    for (int ifac = interface.first_element();
         ifac < interface.last_element();
         ifac += interface.next_element()) {

        bar += interface.next_element();

        this->addInterfaceContribution(VF, interface, ifac, 0., nullptr, 1., 0);
        this->addLocalContribution();
    }

    bar.end();
}

template <typeMesh Mesh, typename Phi>
void AlgoimCutFEMUnified<Mesh, Phi>::addLinear(
    const itemVFlist_t& VF,
    const Interface<mesh_t>& interface
) {
    assert(VF.isRHS());

    progress bar("Add Linear Algoim Interface",
                 interface.last_element(),
                 globalVariable::verbose);

#pragma omp parallel for num_threads(this->get_num_threads())
    for (int ifac = interface.first_element();
         ifac < interface.last_element();
         ifac += interface.next_element()) {

        bar += interface.next_element();

        this->addInterfaceContribution(VF, interface, ifac, 0., nullptr, 1., 0);
    }

    bar.end();
}
