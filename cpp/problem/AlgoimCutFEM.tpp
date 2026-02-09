
template<typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenVol(const Mesh2::Element& K, Phi& phi, const ProblemOption& option, double time) {

    // Quadrature generation strategy:
    // 1. Map the triangle K in physical coordinates to the reference triangle (0,0)-(1,0)-(0,1) called K_ref
    // 2. K_ref can be defined implicitly as the unit square intersected with {psi < 0} where psi = x+y-1 represents the hypothenuse of K_ref     
    // 3. Add the intersection with {phi < 0} where phi is the levelset function used to describe our domain
    // 4. Generate quadrature rules for volume and surface using Bernstein polynomials for phi (psi is always linear)
    // 5. Map quadrature nodes and weights back to the physical triangle using an affine map

    using real = algoim::real;
    using Vec2 = algoim::uvector<real, 2>;
    using R2   = typename Mesh2::Rd;

    phi.t = time; // update time in level set function

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
    int q1d = option.order_space_element_quadrature_;
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

    // Generate quadrature rule on K_ref \cap {phi < 0} \cap {psi < 0}
    ipquad.integrate(algoim::AutoMixed, q1d, [&](const Vec2& xi, real w_ref) {

        // Evaluate exact psi
        if (psi_ref(xi) >= 0) return;

        // Evaluate Bernstein interpolation of phi
        if (algoim::bernstein::evalBernsteinPoly(phiB, xi) >= 0) return;

        const Vec2 x = F(xi);   // physical coordinates of quadrature node
        const double w_phys = w_ref * std::abs(detJ);  // physical quadrature weight
        // const R2 normal = phi.normal(R2(x(0), x(1)));
        
        rule.points.push_back(R2(x(0), x(1)));
        rule.weights.push_back(w_phys);
        // rule.normals.push_back(normal);

    });

    return rule;
}


template <typename Phi>
AlgoimQuadratureRule<Mesh2> quadGenSurf(const Mesh2::Element& K, Phi& phi, const ProblemOption& option, double time) {
    // Quadrature generation strategy:
    // 1. Map the triangle K in physical coordinates to the reference triangle (0,0)-(1,0)-(0,1) called K_ref
    // 2. K_ref can be defined implicitly as the unit square intersected with {psi < 0} where psi = x+y-1 represents the hypothenuse of K_ref     
    // 3. Add the intersection with {phi < 0} where phi is the levelset function used to describe our domain
    // 4. Generate quadrature rules for surface using Bernstein polynomials for phi (psi is always linear)
    // 5. Map quadrature nodes, weights and normals back to the physical triangle using an affine map

    using real = algoim::real;
    using Vec2 = algoim::uvector<real, 2>;
    using R2   = typename Mesh2::Rd;

    phi.t = time; // update time in level set function

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
    auto J_mul = [&](const Vec2& t) -> Vec2 { return e1*t(0) + e2*t(1); };
    
    // Inverse transpose of Jacobian for transforming normals: n_phys = J^{-T} * n_ref
    auto J_inv_T_mul = [&](const Vec2& n) -> Vec2 {
        return Vec2(e2(1)*n(0) - e1(1)*n(1), -e2(0)*n(0) + e1(0)*n(1)) / detJ;
    };

    // Define the reference level set psi_ref(xi) = xi0 + xi1 - 1
    auto psi_ref = [&](const Vec2& xi) -> real { return xi(0) + xi(1) - real(1); };

    // Interpolate phi(F(x)) on the unit square using Bernstein polynomials
    int bernstein_deg = option.algoim_bernstein_deg_;
    int q1d = option.order_space_element_quadrature_;
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

    // Generate quadrature rule on K_ref \cap {phi = 0} \cap {psi < 0} (surface)
    ipquad.integrate_surf(algoim::AutoMixed, q1d, [&](const Vec2& xi, real w_ref, const Vec2& wn_ref) {

        // Evaluate exact psi
        real ps = psi_ref(xi);
        if (ps >= 0) return;  // outside triangle

        const real phB = algoim::bernstein::evalBernsteinPoly(phiB, xi);

        const bool on_phi = std::abs(phB) < tol_phi;
        const bool on_psi = std::abs(ps) < tol_psi;
        if (!on_phi || on_psi) return;

        // Transform normal from reference to physical coordinates
        Vec2 n_ref = wn_ref;
        const real nrm_ref = std::sqrt(n_ref(0)*n_ref(0) + n_ref(1)*n_ref(1));
        if (nrm_ref == 0) return;
        n_ref /= nrm_ref;  // normalize reference normal
        
        Vec2 n_phys = J_inv_T_mul(n_ref);  // transform to physical space
        const real nrm_phys = std::sqrt(n_phys(0)*n_phys(0) + n_phys(1)*n_phys(1));
        n_phys /= nrm_phys;  // normalize physical normal
        // const R2 n_phys(phi.normal(R2(F(xi)(0), F(xi)(1))));  // use normal from level set function in physical space
        
        // Compute tangent and surface Jacobian (same as before for weight computation)
        Vec2 t_ref = Vec2(-n_ref(1), n_ref(0));  // perp(n_ref)
        Vec2 t_phys = J_mul(t_ref);
        const real scale = std::sqrt(t_phys(0)*t_phys(0) + t_phys(1)*t_phys(1));

        const Vec2 x = F(xi);   // physical coordinates
        const double w_phys = w_ref * scale;  // physical quadrature weight
        
        rule.points.push_back(R2(x(0), x(1)));
        rule.weights.push_back(w_phys);
        rule.normals.push_back(R2(n_phys(0), n_phys(1)));  // use transformed normal from Bernstein poly
        // rule.normals.push_back(n_phys);  // use normal from level set function in physical space

    });

    return rule;

    
}



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
    phi_.t = tid; // if your L has time

    // Get quadrature rule - template argument deduction from K and phi_
    auto quad_rule = quadGenVol(K, phi_, options_, tid);

    if (quad_rule.points.size() == 0) {
        std::cout << "Warning: no volume quadrature points for cut element element kb = " << kb << ", returning\n";
        return;
    }

    // assemble (same as your current inner loops)
    for (int l=0; l<VF.size(); ++l) {
        if (!VF[l].on(domain)) continue;

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
            Cint *= VF[l].c;

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


template <typeMesh Mesh, typename Phi>
void AlgoimCutFEMUnified<Mesh, Phi>::addInterfaceContribution(const itemVFlist_t& VF, const Interface<mesh_t>& interface, int ifac, double tid,
                                  const TimeSlab* In, double cst_time, int itq) {


    phi_.t = tid; // update time in level set function

    //  GET IDX ELEMENT CONTAINING FACE ON backMes
    const int kb = interface.idxElementOfFace(ifac);
    const auto &K(interface.get_element(kb));

    // Get quadrature rule - template argument deduction from K and phi_
    auto quad_rule = quadGenSurf(K, phi_, options_, tid);

    if (quad_rule.points.size() == 0) {
        std::cout << "Warning: no surface quadrature points for cut element element kb = " << kb << ", returning\n";
        return;
    }

    for (int l = 0; l < VF.size(); ++l) {

        // if(!VF[l].on(domain)) continue;

        // FINITE ELEMENT SPACES && ELEMENTS
        const fespace_t &Vhv(VF.get_spaceV(l));
        const fespace_t &Vhu(VF.get_spaceU(l));
        bool same = (VF.isRHS() || (&Vhu == &Vhv));

        std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
        std::vector<int> idxU = (same) ? idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());

        int kv = VF[l].onWhatElementIsTestFunction(idxV);
        int ku = VF[l].onWhatElementIsTrialFunction(idxU);

        const auto &FKu(Vhu[ku]);
        const auto &FKv(Vhv[kv]);
        int domu = FKu.get_domain();
        int domv = FKv.get_domain();
        this->initIndex(FKu, FKv);

        // BF MEMORY MANAGEMENT -
        int lastop = getLastop(VF[l].du, VF[l].dv);
        RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
        RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
        What_d Fop = Fwhatd(lastop);

        // Loop over quadrature points
        for (size_t ipq = 0; ipq < quad_rule.points.size(); ++ipq) {
            
            Rd mip = quad_rule.points[ipq];
            const R weight = quad_rule.weights[ipq];
            
            assert(weight > 0);
            
            const Rd face_ip = K.mapToReferenceElement(mip);
            double Cint      = weight * cst_time;

            const Rd normal = quad_rule.normals[ipq];

            assert(std::fabs(normal.norm() - 1) < 1e-14);
            double coef = VF[l].computeCoefFromNormal(normal);

            FKv.BF(Fop, face_ip, fv);

            if (!same)
                FKu.BF(Fop, face_ip, fu);

            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domu, domv), mip, tid,
                                                           normal);
            Cint *= coef * VF[l].c;

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

// Stub implementation for addFaceContribution (to be implemented)
template <typeMesh Mesh, typename Phi>
void AlgoimCutFEMUnified<Mesh, Phi>::addFaceContribution(const itemVFlist_t& VF, const std::pair<int, int>& e1, 
                                                          const std::pair<int, int>& e2, const TimeSlab* In, 
                                                          int itq, double cst_time) {
    // TODO: Implement face contribution using Algoim quadrature
    // For now, this is a stub to resolve linker error
    std::cout << "Warning: addFaceContribution not yet implemented for AlgoimCutFEMUnified\n";
}

