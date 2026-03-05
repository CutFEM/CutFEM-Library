#include "baseProblem.hpp" 
#include "baseCutProblem.hpp" 

template <> void BaseCutFEM<Mesh2>::addFaceStabilization(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(!VF.isRHS());
    assert(Th.get_nb_domain() == 1);
    progress bar(" Add Face Stabilization CutMesh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        int domain = Th.get_domain_element(k);
        
        // const int active_micro_idx_d = k - ((domain == 0) ? 0 : Th.get_nb_element(0));  // the to domain local active mesh index
        // const int k_macro = Th.inverse_active_macro_map_d[domain].at(active_micro_idx_d);
        const int k_macro = Th.macro_of_micro(k);
        // std::cout << "Element " << k << " (local in active mesh = " << k - idx_shift << "), belongs to domain " << domain << " and macro element " << k_macro << "\n";
        if (!Th.is_macro_cut(k_macro, domain, 0))
            continue;

        
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn == -1) continue;

            //if (Th.get_domain_element(kn) != domain) continue; // <- critical fix

            const int kn_macro = Th.macro_of_micro(kn);

            if ((kn < k) && Th.is_macro_cut(kn_macro, domain, 0))
                continue;
            

            int kb  = Th.idxElementInBackMesh(k);
            int kbn = Th.idxElementInBackMesh(kn);

            std::pair<int, int> e1 = std::make_pair(kb, ifac);
            std::pair<int, int> e2 = std::make_pair(kbn, jfac);
            BaseFEM<Mesh2>::addFaceContribution(VF, e1, e2, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <> void BaseCutFEM<Mesh2>::addPatchStabilization(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Patch Stabilization BarycentricActiveMesh2", Th.last_element(), globalVariable::verbose);

    size_t num_stab_faces = 0;

    // Loop over micro element
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        int domain = Th.get_domain_element(k);
        
        // const int active_micro_idx_d = k - ((domain == 0) ? 0 : Th.get_nb_element(0));  // the to domain local active mesh index
        // const int k_macro = Th.inverse_active_macro_map_d[domain].at(active_micro_idx_d);
        const int k_macro = Th.macro_of_micro(k);
        // std::cout << "Element " << k << " (local in active mesh = " << k - idx_shift << "), belongs to domain " << domain << " and macro element " << k_macro << "\n";
        if (!Th.is_macro_cut(k_macro, domain, 0))
            continue;

        // if (!Th.isCut(k, 0) && !Th.isInactive(k, 0))
        //     continue;
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
        
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            // if (kn < k)
            //     continue;
            if (kn == -1) continue;

            // int kn_macro = Th.inverse_active_macro_map_d[domain][kn];
            // if (kn < k && Th.is_macro_cut(kn, domain, 0))

            if (Th.get_domain_element(kn) != domain) continue; // <- critical fix

            const int kn_macro = Th.macro_of_micro(kn);

            if ((kn < k) && Th.is_macro_cut(kn_macro, domain, 0))
                continue;
            
                // std::pair<int, int> e1 = std::make_pair(k, ifac);
            // std::pair<int, int> e2 = std::make_pair(kn, jfac);
            const int kb  = Th.idxElementInBackMesh(k, domain);
            const int kbn = Th.idxElementInBackMesh(kn, domain);
            
            // const int kbn_macro = Th.macro_idx_in_background_mesh_[domain][kn_macro];
            // const int kb_macro = Th.macro_idx_in_background_mesh_[domain][k_macro];

            // std::cout << "Domain " << domain << " stabilizing face between micro elements " << kb << " and " << kbn << " in the macros " << kb_macro << ", " << kbn_macro <<  "\n";
            // std::cout << "The element " << kb << " is (" << Th.Th[kb][0] << "), (" << Th.Th[kb][1] << "), (" << Th.Th[kb][2] << ")\n";
            BaseFEM<Mesh2>::addPatchContribution(VF, k, kn, nullptr, 0, 1.);
            num_stab_faces++;
        }
        this->addLocalContribution();
    }
    bar.end();
    //std::cout << "Number of stabilized faces: " << num_stab_faces << "\n";
}

template <> void BaseCutFEM<Mesh2>::addPatchStabilization(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th, const TimeSlab &In) {
    
    int number_of_quadrature_points = this->get_nb_quad_point_time();
    int num_stab_faces = 0;

    // Loop through time quadrature points
    for (int itq = 0; itq < number_of_quadrature_points; ++itq) {
        assert(!VF.isRHS());

        // Compute contribution from time basis functions
        auto tq    = this->get_quadrature_time(itq);
        double tid = In.map(tq);
        KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
        RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
        In.BF(tq.x, bf_time); // compute time basic funtions
        double cst_time = tq.a * In.get_measure();

        progress bar(" Add Patch Stabilization BarycentricActiveMesh2", Th.last_element(), globalVariable::verbose);

        // Loop over micro element
        for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
            bar += Th.next_element();

            int domain = Th.get_domain_element(k);
            
            // const int active_micro_idx_d = k - ((domain == 0) ? 0 : Th.get_nb_element(0));  // the to domain local active mesh index
            // const int k_macro = Th.inverse_active_macro_map_d[domain].at(active_micro_idx_d);
            const int k_macro = Th.macro_of_micro(k);
            // std::cout << "Element " << k << " (local in active mesh = " << k - idx_shift << "), belongs to domain " << domain << " and macro element " << k_macro << "\n";
            if (!Th.stabilize_macro(k_macro))
                continue;

            // std::cout << "Stabilizing macro element " << k_macro << " in domain " << domain << "\n";
      
            for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

                int jfac = ifac;
                int kn   = Th.ElementAdj(k, jfac);
            
                // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
                // if (kn < k)
                //     continue;
                if (kn == -1) continue;

                // int kn_macro = Th.inverse_active_macro_map_d[domain][kn];
                // if (kn < k && Th.is_macro_cut(kn, domain, 0))

                if (Th.get_domain_element(kn) != domain) continue; // <- critical fix

                const int kn_macro = Th.macro_of_micro(kn);

                if ((kn < k) && Th.stabilize_macro(kn_macro))
                    continue;
                
                // std::cout << " Stabilizing face between micro elements " << k << " and " << kn << " in the macros " << k_macro << ", " << kn_macro <<  "\n";
                
                    // std::pair<int, int> e1 = std::make_pair(k, ifac);
                // std::pair<int, int> e2 = std::make_pair(kn, jfac);
                const int kb  = Th.idxElementInBackMesh(k, domain);
                const int kbn = Th.idxElementInBackMesh(kn, domain);
                
                // const int kbn_macro = Th.macro_idx_in_background_mesh_[domain][kn_macro];
                // const int kb_macro = Th.macro_idx_in_background_mesh_[domain][k_macro];

                // std::cout << "Domain " << domain << " stabilizing face between micro elements " << kb << " and " << kbn << " in the macros " << kb_macro << ", " << kbn_macro <<  "\n";
                // std::cout << "The element " << kb << " is (" << Th.Th[kb][0] << "), (" << Th.Th[kb][1] << "), (" << Th.Th[kb][2] << ")\n";
                BaseFEM<Mesh2>::addPatchContribution(VF, k, kn, &In, itq, cst_time);
                num_stab_faces++;
            }
            this->addLocalContribution();
        }
        bar.end();

    }

    // std::cout << "Number of stabilized faces: " << num_stab_faces / number_of_quadrature_points << "\n";
}


template <> void BaseCutFEM<Mesh2>::addFaceStabilizationMixed(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(!VF.isRHS());
    assert(Th.get_nb_domain() == 1);
    progress bar(" Add Face Stabilization CutMesh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        int domain = Th.get_domain_element(k);
        
        // const int active_micro_idx_d = k - ((domain == 0) ? 0 : Th.get_nb_element(0));  // the to domain local active mesh index
        // const int k_macro = Th.inverse_active_macro_map_d[domain].at(active_micro_idx_d);
        const int k_macro = Th.macro_of_micro(k);
        // std::cout << "Element " << k << " (local in active mesh = " << k - idx_shift << "), belongs to domain " << domain << " and macro element " << k_macro << "\n";
        if (!Th.is_macro_cut(k_macro, domain, 0))
            continue;

        
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn == -1) continue;

            //if (Th.get_domain_element(kn) != domain) continue; // <- critical fix

            const int kn_macro = Th.macro_of_micro(kn);

            if ((kn < k) && Th.is_macro_cut(kn_macro, domain, 0))
                continue;
            

            int kb  = Th.idxElementInBackMesh(k);
            int kbn = Th.idxElementInBackMesh(kn);

            std::pair<int, int> e1 = std::make_pair(kb, ifac);
            std::pair<int, int> e2 = std::make_pair(kbn, jfac);
            BaseFEM<Mesh2>::addFaceContributionMixed(VF, e1, e2, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <> void BaseCutFEM<Mesh2>::addPatchStabilizationMixed(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Patch Stabilization BarycentricActiveMesh2", Th.last_element(), globalVariable::verbose);

    size_t num_stab_faces = 0;

    // Loop over micro element
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        int domain = Th.get_domain_element(k);
        
        // const int active_micro_idx_d = k - ((domain == 0) ? 0 : Th.get_nb_element(0));  // the to domain local active mesh index
        // const int k_macro = Th.inverse_active_macro_map_d[domain].at(active_micro_idx_d);
        const int k_macro = Th.macro_of_micro(k);
        // std::cout << "Element " << k << " (local in active mesh = " << k - idx_shift << "), belongs to domain " << domain << " and macro element " << k_macro << "\n";
        if (!Th.is_macro_cut(k_macro, domain, 0))
            continue;

        // if (!Th.isCut(k, 0) && !Th.isInactive(k, 0))
        //     continue;
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
        
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            // if (kn < k)
            //     continue;
            if (kn == -1) continue;

            // int kn_macro = Th.inverse_active_macro_map_d[domain][kn];
            // if (kn < k && Th.is_macro_cut(kn, domain, 0))

            // if (Th.get_domain_element(kn) != domain) continue; // <- critical fix

            const int kn_macro = Th.macro_of_micro(kn);

            if ((kn < k) && Th.is_macro_cut(kn_macro, domain, 0))
                continue;
            
                // std::pair<int, int> e1 = std::make_pair(k, ifac);
            // std::pair<int, int> e2 = std::make_pair(kn, jfac);
            const int kb  = Th.idxElementInBackMesh(k, domain);
            const int kbn = Th.idxElementInBackMesh(kn, domain);
            
            // const int kbn_macro = Th.macro_idx_in_background_mesh_[domain][kn_macro];
            // const int kb_macro = Th.macro_idx_in_background_mesh_[domain][k_macro];

            // std::cout << "Domain " << domain << " stabilizing face between micro elements " << kb << " and " << kbn << " in the macros " << kb_macro << ", " << kbn_macro <<  "\n";
            // std::cout << "The element " << kb << " is (" << Th.Th[kb][0] << "), (" << Th.Th[kb][1] << "), (" << Th.Th[kb][2] << ")\n";
            BaseFEM<Mesh2>::addPatchContributionMixed(VF, kb, kbn, nullptr, 0, 1.);
            num_stab_faces++;
        }
        this->addLocalContribution();
    }
    bar.end();
    //std::cout << "Number of stabilized faces: " << num_stab_faces << "\n";
}


template <> void BaseCutFEM<Mesh2>::addBilinearInner(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);

    const int domain = 0;
    assert(Th.get_nb_domain() == 1);
    for (int km = 0; km < Th.active_macro_elements_d[domain].size(); ++km) {
        
        if (Th.is_macro_cut(km))
            continue;

        assert(Th.is_macro_interior(km));
        
        const auto& micro_elements = Th.active_macro_elements_d[domain][km];

        // Loop over all micro elements in the macro element
        for (int k_micro : micro_elements) 
            BaseFEM<Mesh2>::addElementContribution(VF, k_micro, nullptr, 0, 1.);
        
        this->addLocalContribution();
    }
    bar.end();
}

template <> void BaseCutFEM<Mesh2>::addLinearInner(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
    
    const int domain = 0;
    assert(Th.get_nb_domain() == 1);
    for (int km = 0; km < Th.active_macro_elements_d[domain].size(); ++km) {
        
        if (Th.is_macro_cut(km))
            continue;

        assert(Th.is_macro_interior(km));
        
        const auto& micro_elements = Th.active_macro_elements_d[domain][km];

        // Loop over all micro elements in the macro element
        for (int k_micro : micro_elements) 
            BaseFEM<Mesh2>::addElementContribution(VF, k_micro, nullptr, 0, 1.);
        
        this->addLocalContribution();
    }
    bar.end();
}

template <> void BaseCutFEM<Mesh2>::addBilinearInnerBorder(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(!VF.isRHS());

    const int domain = 0;
    assert(Th.get_nb_domain() == 1);
    for (int km = 0; km < Th.active_macro_elements_d[domain].size(); ++km) {

        // Loop only over cut elements for efficiency
        if (!Th.is_macro_cut(km))
            continue;

        const auto& micro_elements = Th.active_macro_elements_d[domain][km];

        for (int k_micro : micro_elements) {
            
            for (int ifac = 0; ifac < Element::nea; ++ifac) {
                int jfac = ifac;
                int kn_micro = Th.ElementAdj(k_micro, jfac);

                if (kn_micro == -1) 
                    continue;
                int kn_macro = Th.inverse_active_macro_map_d[domain][kn_micro];
                if (Th.is_macro_cut(kn_macro))  // 2)
                    continue;
                
                // Skip if in same macro element
                if (kn_macro == km)
                    continue;
                
                // assert(!Th.isCut(kn_micro, 0));
                assert(Th.is_macro_interior(kn_macro));

                // std::cout << "Integrating on edge between micro elements " << k_micro << " (macro element " << Th.inverse_active_macro_map[k_micro] << ") and " << kn_micro << " (macro element " << Th.inverse_active_macro_map[kn_micro] << ")\n";
                BaseFEM<Mesh2>::addInnerBorderContribution(VF, k_micro, ifac, nullptr, 0, 1.);
            
            }
        
        }

    }

}

template <> void BaseCutFEM<Mesh2>::addLinearInnerBorder(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(VF.isRHS());

    const int domain = 0;
    assert(Th.get_nb_domain() == 1);
    for (int km = 0; km < Th.active_macro_elements_d[domain].size(); ++km) {

        // Loop only over cut elements for efficiency
        if (!Th.is_macro_cut(km))
            continue;

        const auto& micro_elements = Th.active_macro_elements_d[domain][km];

        for (int k_micro : micro_elements) {
            
            for (int ifac = 0; ifac < Element::nea; ++ifac) {
                int jfac = ifac;
                int kn_micro = Th.ElementAdj(k_micro, jfac);

            // a neighboring element can be 1) outside of the domain, 2) another cut element, 
            // 3) inside of the domain but in the same macro element, or 4) inside of the domain and in another macro element
                if (kn_micro == -1) 
                    continue;
                int kn_macro = Th.inverse_active_macro_map_d[domain][kn_micro];
                if (Th.is_macro_cut(kn_macro))  // 2)
                    continue;
                
                // Skip if in same macro element
                if (kn_macro == km)
                    continue;
                
                // assert(!Th.isCut(kn_micro, 0));
                assert(Th.is_macro_interior(kn_macro));

                // std::cout << "Integrating on edge between micro elements " << k_micro << " (macro element " << Th.inverse_active_macro_map[k_micro] << ") and " << kn_micro << " (macro element " << Th.inverse_active_macro_map[kn_micro] << ")\n";
                BaseFEM<Mesh2>::addInnerBorderContribution(VF, k_micro, ifac, nullptr, 0, 1.);
            
            }
        
        }

    }

}

template <> void BaseCutFEM<Mesh2>::addBilinearInnerBorderMixed(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(!VF.isRHS());

    const int domain = 0;
    assert(Th.get_nb_domain() == 1);
    for (int km = 0; km < Th.active_macro_elements_d[domain].size(); ++km) {

        // Loop only over cut elements for efficiency
        if (!Th.is_macro_cut(km))
            continue;

        const auto& micro_elements = Th.active_macro_elements_d[domain][km];

        for (int k_micro : micro_elements) {
            
            for (int ifac = 0; ifac < Element::nea; ++ifac) {
                int jfac = ifac;
                int kn_micro = Th.ElementAdj(k_micro, jfac);

                // a neighboring element can be 1) outside of the domain, 2) another cut macro element, 
                // 3) inside of the domain but in the same macro element, or 4) inside of the domain and in another macro element
                // we only want to integrate on the face to an element that satisfies 4)
                if (kn_micro == -1) // 1)
                    continue;
                
                int kn_macro = Th.inverse_active_macro_map_d[domain][kn_micro];
                if (Th.is_macro_cut(kn_macro))  // 2)
                    continue;
                
                // Skip if in same macro element
                if (kn_macro == km) // 3)
                    continue;
                
                const int kbn_macro = Th.macro_idx_in_background_mesh_[domain][kn_macro];
                const int kb_macro = Th.macro_idx_in_background_mesh_[domain][km];
                // std::cout << "Integrating on edge between micro elements " << k_micro << " (macro element " << Th.inverse_active_macro_map_d[domain][k_micro] << " (bg = " << kb_macro <<  ")) and " << kn_micro << " (macro element " << Th.inverse_active_macro_map_d[domain][kn_micro] << " (bg = " << kbn_macro <<  "))\n";
                

                // assert(!Th.isCut(kn_micro, 0));
                assert(Th.is_macro_interior(kn_macro));

                int kb_micro = Th.idxElementInBackMesh(k_micro);

                
                BaseFEM<Mesh2>::addInnerBorderContributionMixed(VF, kb_micro, ifac, nullptr, 0, 1.);
            
            }
        
        }

    }

}

template <> void BaseCutFEM<Mesh2>::addLinearInnerBorderMixed(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(VF.isRHS());

    const int domain = 0;
    assert(Th.get_nb_domain() == 1);
    for (int km = 0; km < Th.active_macro_elements_d[domain].size(); ++km) {

        // Loop only over cut elements for efficiency
        if (!Th.is_macro_cut(km))
            continue;

        const auto& micro_elements = Th.active_macro_elements_d[domain][km];

        for (int k_micro : micro_elements) {
            
            for (int ifac = 0; ifac < Element::nea; ++ifac) {
                int jfac = ifac;
                int kn_micro = Th.ElementAdj(k_micro, jfac);

                // a neighboring element can be 1) outside of the domain, 2) another cut element, 
                // 3) inside of the domain but in the same macro element, or 4) inside of the domain and in another macro element
                if (kn_micro == -1) 
                    continue;
                int kn_macro = Th.inverse_active_macro_map_d[domain][kn_micro];
                if (Th.is_macro_cut(kn_macro))  // 2)
                    continue;
                
                // Skip if in same macro element
                if (kn_macro == km)
                    continue;
                
                // assert(!Th.isCut(kn_micro, 0));
                assert(Th.is_macro_interior(kn_macro));

                int kb_micro = Th.idxElementInBackMesh(k_micro);

                // std::cout << "Integrating on edge between micro elements " << k_micro << " (macro element " << Th.inverse_active_macro_map[k_micro] << ") and " << kn_micro << " (macro element " << Th.inverse_active_macro_map[kn_micro] << ")\n";
                BaseFEM<Mesh2>::addInnerBorderContributionMixed(VF, kb_micro, ifac, nullptr, 0, 1.);
            
            }
        
        }

    }

}

// Outer border of the active mesh
template<>
void BaseCutFEM<Mesh2>::addBilinearOuterBorder(const itemVFlist_t& VF,
                                               const BarycentricActiveMesh2& Th)
{
    using Element = typename BarycentricActiveMesh2::Element;
    assert(!VF.isRHS());

    const int domain = 0;
    assert(Th.get_nb_domain() == 1);

    for (int km = 0; km < (int)Th.active_macro_elements_d[domain].size(); ++km) {
        if (!Th.is_macro_cut(km)) continue;

        const auto& micro = Th.active_macro_elements_d[domain][km];
        for (int k_micro : micro) {
            for (int ifac = 0; ifac < Element::nea; ++ifac) {
                int jfac = ifac;
                const int kn_micro = Th.ElementAdj(k_micro, jfac);
                if (kn_micro != -1) continue;         // only outer boundary

                BaseFEM<Mesh2>::addOuterBorderContribution(VF, k_micro, ifac, nullptr, 0, 1.0);
            }
        }
    }
}
template<>
void BaseCutFEM<Mesh2>::addBilinearOuterBorderMixed(const itemVFlist_t& VF,
                                               const BarycentricActiveMesh2& Th)
{
    using Element = typename BarycentricActiveMesh2::Element;
    assert(!VF.isRHS());

    const int domain = 0;
    assert(Th.get_nb_domain() == 1);

    for (int km = 0; km < (int)Th.active_macro_elements_d[domain].size(); ++km) {
        if (!Th.is_macro_cut(km)) continue;

        const auto& micro = Th.active_macro_elements_d[domain][km];
        for (int k_micro : micro) {
            for (int ifac = 0; ifac < Element::nea; ++ifac) {
                int jfac = ifac;
                const int kn_micro = Th.ElementAdj(k_micro, jfac);
                if (kn_micro != -1) continue;         // only outer boundary

                int kb_micro = Th.idxElementInBackMesh(k_micro);

                BaseFEM<Mesh2>::addOuterBorderContributionMixed(VF, kb_micro, ifac, nullptr, 0, 1.0);
            }
        }
    }
}


template <> void BaseCutFEM<Mesh2>::addLinearOuterBorder(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(VF.isRHS());

    const int domain = 0;
    assert(Th.get_nb_domain() == 1);

    for (int km = 0; km < Th.active_macro_elements_d[domain].size(); ++km) {

        // Loop only over cut elements for efficiency
        if (!Th.is_macro_cut(km))
            continue;

        const auto& micro_elements = Th.active_macro_elements_d[domain][km];

        for (int k_micro : micro_elements) {
            
            for (int ifac = 0; ifac < Element::nea; ++ifac) {
                int jfac = ifac;
                int kn_micro = Th.ElementAdj(k_micro, jfac);

                // if ((kn_micro != -1) ||(Th.inverse_active_macro_map[k_micro] == Th.inverse_active_macro_map[kn_micro]))
                //     continue;
                if (kn_micro != -1)
                    continue;
                
                // std::cout << "Integrating on edge between micro elements " << k_micro << " (macro element " << Th.inverse_active_macro_map[k_micro] << ") and " << kn_micro << " (macro element " << Th.inverse_active_macro_map[kn_micro] << ")\n";
                BaseFEM<Mesh2>::addOuterBorderContribution(VF, k_micro, ifac, nullptr, 0, 1.);
            
            }
        
        }

    }

}
template <> void BaseCutFEM<Mesh2>::addLinearOuterBorderMixed(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(VF.isRHS());

    const int domain = 0;
    assert(Th.get_nb_domain() == 1);

    for (int km = 0; km < Th.active_macro_elements_d[domain].size(); ++km) {

        // Loop only over cut elements for efficiency
        if (!Th.is_macro_cut(km))
            continue;

        const auto& micro_elements = Th.active_macro_elements_d[domain][km];

        for (int k_micro : micro_elements) {
            
            for (int ifac = 0; ifac < Element::nea; ++ifac) {
                int jfac = ifac;
                int kn_micro = Th.ElementAdj(k_micro, jfac);

                // if ((kn_micro != -1) ||(Th.inverse_active_macro_map[k_micro] == Th.inverse_active_macro_map[kn_micro]))
                //     continue;
                if (kn_micro != -1)
                    continue;
                
                int kb_micro = Th.idxElementInBackMesh(k_micro);

                // std::cout << "Integrating on edge between micro elements " << k_micro << " (macro element " << Th.inverse_active_macro_map[k_micro] << ") and " << kn_micro << " (macro element " << Th.inverse_active_macro_map[kn_micro] << ")\n";
                BaseFEM<Mesh2>::addOuterBorderContributionMixed(VF, kb_micro, ifac, nullptr, 0, 1.);
            
            }
        
        }

    }

}






// template <> void BaseCutFEM<Mesh2>::addElementStabilization(const itemVFlist_t &VF, const Interface<Mesh2> &interface, const CutMesh &Th) {
//     assert(!VF.isRHS());
//     progress bar("Add Bilinear Mesh", Th.last_element(), globalVariable::verbose);
// #pragma omp parallel for num_threads(this->get_num_threads())
//     for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
//         bar += Th.next_element();
        
//         // Compute parameter coonected to the mesh.
//         // on can take the one from the first test function
//         const FESpace &Vh(*VF[0].fespaceV);
//         const FElement &FK(Vh[k]);
//         const Element &K(FK.T);
//         double meas = K.measure();
//         double h    = K.get_h();
//         int domain  = FK.get_domain();
//         int kb      = Vh.idxElementInBackMesh(k);
//         int iam     = omp_get_thread_num();

//         // GET THE QUADRATURE RULE
//         const QF &qf(this->get_quadrature_formular_K());
    
//         double tid = 0.;

//         // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
//         for (int l = 0; l < VF.size(); ++l) {
//             if (!VF[l].on(domain))
//                 continue; 

//             // FINTE ELEMENT SPACES && ELEMENTS
//             const FESpace &Vhv(VF.get_spaceV(l));
//             const FESpace &Vhu(VF.get_spaceU(l));
//             const FElement &FKv(Vhv[k]);
//             const FElement &FKu(Vhu[k]);
//             this->initIndex(FKu, FKv); //, iam);
//             // BF MEMORY MANAGEMENT -
//             bool same   = (&Vhu == &Vhv);
//             int lastop  = getLastop(VF[l].du, VF[l].dv);
//             // RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop); //  the value for
//             // basic fonction RNMK_ fu(this->databf_+ (same
//             // ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value
//             // for basic fonction
//             long offset = iam * this->offset_bf_;
//             RNMK_ fv(this->databf_ + offset, FKv.NbDoF(), FKv.N,
//                     lastop); //  the value for basic fonction
//             RNMK_ fu(this->databf_ + offset + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
//             What_d Fop = Fwhatd(lastop);

//             // COMPUTE COEFFICIENT
//             double coef = VF[l].computeCoefElement(h, meas, meas, meas, domain);

//             // LOOP OVER QUADRATURE IN SPACE
//             for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
//                 typename QF::QuadraturePoint ip(qf[ipq]);
//                 const Rd mip = K.mapToPhysicalElement(ip);
//                 double Cint  = meas * ip.getWeight() * cst_time;

//                 // EVALUATE THE BASIS FUNCTIONS
//                 FKv.BF(Fop, ip, fv);
//                 if (!same)
//                     FKu.BF(Fop, ip, fu);
//                 //   VF[l].applyFunNL(fu,fv);

//                 // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
//                 Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);
//                 Cint *= coef * VF[l].c;

//                 if (In) {
//                     if (VF.isRHS())
//                         this->addToRHS(VF[l], *In, FKv, fv, Cint);
//                     else
//                         this->addToMatrix(VF[l], *In, FKu, FKv, fu, fv, Cint);
//                 } else {
//                     if (VF.isRHS())
//                         this->addToRHS(VF[l], FKv, fv, Cint);
//                     else
//                         this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint); //, iam);
//                 }
//             }
//         }

//         this->addLocalContribution();
//     }
//     bar.end();
// }