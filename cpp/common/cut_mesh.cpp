#include "time_interface.hpp"
#include "cut_mesh.hpp"



BarycentricActiveMesh2::BarycentricActiveMesh2(const BarycentricMesh2 &th)
    : ActiveMesh<Mesh2>(static_cast<const Mesh2 &>(th)) {}


BarycentricActiveMesh2::BarycentricActiveMesh2(const BarycentricMesh2 &th, const Interface<Mesh2> &interface) : ActiveMesh<Mesh2>(static_cast<const Mesh2 &>(th)) {
    // Since this calls the ActiveMesh(const Mesh &) constructor, we need to clear all the containers
    idx_in_background_mesh_.clear();
    idx_from_background_mesh_.clear();
    interface_id_.clear();
    idx_element_domain.clear();
    not_in_active_mesh_.clear();

    idx_in_background_mesh_.resize(2);
    idx_from_background_mesh_.resize(2);
    interface_id_.resize(1);
    // not_in_active_mesh_.resize(2);
    nb_quadrature_time_ = 1;
    

    macro_idx_in_background_mesh_.resize(2);
    macro_idx_from_background_mesh_.resize(2);
    active_macro_elements_d.resize(2);
    inverse_active_macro_map_d.resize(2);
    nb_active_macros_d.resize(2);

    not_in_active_mesh_.resize(21);
    for (int i = 0; i < 21; ++i)
        not_in_active_mesh_[i].resize(nb_quadrature_time_);

    init(interface);


}
    

// --- Stationary --------------------------------------------------------------
void BarycentricActiveMesh2::truncate(const Interface<Mesh2>& interface,
                                      int sign_domain_remove)
{
    // reset macro containers state
    macro_idx_in_background_mesh_.clear();
    macro_idx_from_background_mesh_.clear();
    inverse_active_macro_map_d.clear();
    active_macro_elements_d.clear();
    nb_active_macros_d.clear();

    // reset active mesh containers
    idx_in_background_mesh_.clear();
    idx_from_background_mesh_.clear();
    interface_id_.clear();
    idx_element_domain.clear();
    not_in_active_mesh_.clear();

    // initialize containers
    macro_idx_in_background_mesh_.resize(1);
    macro_idx_from_background_mesh_.resize(1);
    active_macro_elements_d.resize(1);
    inverse_active_macro_map_d.resize(1);
    nb_active_macros_d.resize(1);

    idx_in_background_mesh_.resize(1);
    idx_in_background_mesh_[0].reserve(Th.nt);
    idx_from_background_mesh_.resize(1);
    interface_id_.resize(1);
    nb_quadrature_time_ = 1;
    

    not_in_active_mesh_.resize(21);
    for (int i = 0; i < 21; ++i)
        not_in_active_mesh_[i].resize(nb_quadrature_time_);

    
    int domain = 0;

    const auto& Th_bary = static_cast<const BarycentricMesh2&>(this->Th);
    
    macro_idx_in_background_mesh_[domain].reserve(Th_bary.macro_elements.size());

    // macro bookkeeping
    int nt = 0;        // number of active elements in the domain
    int nt_macros = 0; // number of active macro elements 

    // --- decide kept macros ---
    for (int macro_k = 0; macro_k < (int)Th_bary.macro_elements.size(); ++macro_k) {
        const auto& sub = Th_bary.macro_elements[macro_k]; // 3 background subelements

        bool keep_macro = false;
        for (int i = 0; i < 3 && !keep_macro; ++i) {
            const int sub_k = sub[i];
            const auto signK = interface.get_SignElement(sub_k);
            if (interface.isCut(sub_k) || signK.sign() != sign_domain_remove)
                keep_macro = true;
        }
        if (!keep_macro) continue;

        macro_idx_in_background_mesh_[domain].push_back(macro_k);
        macro_idx_from_background_mesh_[domain][macro_k] = nt_macros;
        
        std::array<int, 3> sub_elems_active;

        // keep all three subelements (into every kept domain; adapt if you keep a subset)
        for (int i = 0; i < 3; ++i) {
            const int sub_k = sub[i];
            const bool is_cut = interface.isCut(sub_k);
            const int sign   = interface.get_SignElement(sub_k).sign();

            
            const int local_id = nt;     // index in active mesh
            idx_in_background_mesh_[domain].push_back(sub_k);
            idx_from_background_mesh_[domain][sub_k] = local_id;
            sub_elems_active[i] = local_id;
            inverse_active_macro_map_d[domain].push_back(nt_macros);

            if (is_cut) {
                interface_id_[0][{domain, local_id}].emplace_back(&interface, -sign_domain_remove);
            } else if (sign == sign_domain_remove) {
                // mark "outside" for stationary t=0
                not_in_active_mesh_[domain][0][local_id] = true;
            }
            ++nt;
            
        }
        active_macro_elements_d[domain].push_back(sub_elems_active);
        nt_macros++;
    }

    nb_active_macros_d[domain] = nt_macros;
    // finalize prefix sums
    idx_element_domain.push_back(0);
    idx_element_domain.push_back(nt);
    idx_in_background_mesh_[domain].shrink_to_fit();
    

}




// --- Time-dependent ----------------------------------------------------------
void BarycentricActiveMesh2::truncate(const TimeInterface<Mesh2>& interface,
                                      int sign_domain_remove)
{
    
    // reset macro containers state
    macro_idx_in_background_mesh_.clear();
    macro_idx_from_background_mesh_.clear();
    inverse_active_macro_map_d.clear();
    active_macro_elements_d.clear();
    nb_active_macros_d.clear();

    // reset active mesh containers
    idx_in_background_mesh_.clear();
    idx_from_background_mesh_.clear();
    interface_id_.clear();
    idx_element_domain.clear();
    not_in_active_mesh_.clear();

    // initialize containers
    macro_idx_in_background_mesh_.resize(1);
    macro_idx_from_background_mesh_.resize(1);
    active_macro_elements_d.resize(1);
    inverse_active_macro_map_d.resize(1);
    nb_active_macros_d.resize(1);

    int n_tid = interface.size();
    nb_quadrature_time_ = n_tid;

    interface_id_.resize(21);
    interface_id_.assign(nb_quadrature_time_, {}); // time slices
    not_in_active_mesh_.resize(10);
    for (int i = 0; i < 10; ++i)
        not_in_active_mesh_[i].resize(nb_quadrature_time_);

    idx_in_background_mesh_.resize(1);
    idx_in_background_mesh_[0].reserve(Th.nt);
    idx_from_background_mesh_.resize(1);

    const int domain = 0;
    int nt = 0;
    int nt_macros = 0; // number of active macro elements 

    const auto &Th_bary = static_cast<const BarycentricMesh2 &>(this->Th);

    macro_idx_in_background_mesh_[domain].reserve(Th_bary.macro_elements.size());


    // --- decide kept macros (exists t where inside OR cut OR sign flips) ---
    for (int macro_k = 0; macro_k < (int)Th_bary.macro_elements.size(); ++macro_k) {
        const auto& sub = Th_bary.macro_elements[macro_k];

        bool keep_macro = false;
        for (int i = 0; i < 3 && !keep_macro; ++i) {
            const int sub_k = sub[i];

            for (int it = 0; it + 1 < n_tid; ++it) {
                const int s_i  = interface(it)->get_SignElement(sub_k).sign();
                const int s_ip = interface(it+1)->get_SignElement(sub_k).sign();
                const bool cut_i  = interface(it)->isCut(sub_k);
                const bool cut_ip = interface(it+1)->isCut(sub_k);

                if (cut_i || cut_ip || s_i != sign_domain_remove || s_i * s_ip <= 0) {
                    keep_macro = true;
                    break;
                }
            }
        }
        if (!keep_macro) continue;

        macro_idx_in_background_mesh_[domain].push_back(macro_k);
        macro_idx_from_background_mesh_[domain][macro_k] = nt_macros;
        
        std::array<int, 3> sub_elems_active;

        // keep all three subelements
        for (int i = 0; i < 3; ++i) {
            const int sub_k = sub[i];

            
            const int local_id = nt;
            sub_elems_active[i] = local_id;
            idx_in_background_mesh_[domain].push_back(sub_k);
            idx_from_background_mesh_[domain][sub_k] = local_id;
            inverse_active_macro_map_d[domain].push_back(nt_macros);

            // per-time bookkeeping
            for (int it = 0; it < n_tid; ++it) {
                const bool is_cut = interface(it)->isCut(sub_k);
                const int sign    = interface(it)->get_SignElement(sub_k).sign();

                if (is_cut) {
                    interface_id_[it][{domain, local_id}].emplace_back(interface[it], -sign_domain_remove);
                } else if (sign == sign_domain_remove) {
                    not_in_active_mesh_[domain][it][local_id] = true; // outside at this t
                }
            }

            ++nt;
            
        }
        active_macro_elements_d[domain].push_back(sub_elems_active);
        nt_macros++;
    }

    nb_active_macros_d[domain] = nt_macros;

    // finalize prefix sums
    idx_element_domain.push_back(0);
    idx_element_domain.push_back(nt);
    idx_in_background_mesh_[domain].shrink_to_fit();

}


void BarycentricActiveMesh2::createSurfaceMesh(const Interface<Mesh2> &interface) {
    // reset macro containers state
    macro_idx_in_background_mesh_.clear();
    macro_idx_from_background_mesh_.clear();
    inverse_active_macro_map_d.clear();
    active_macro_elements_d.clear();
    nb_active_macros_d.clear();

    // reset active mesh containers
    idx_in_background_mesh_.clear();
    idx_from_background_mesh_.clear();
    interface_id_.clear();
    idx_element_domain.clear();
    not_in_active_mesh_.clear();

    // initialize containers
    macro_idx_in_background_mesh_.resize(1);
    macro_idx_from_background_mesh_.resize(1);
    active_macro_elements_d.resize(1);
    inverse_active_macro_map_d.resize(1);
    nb_active_macros_d.resize(1);

    idx_in_background_mesh_.resize(1);
    idx_in_background_mesh_[0].reserve(Th.nt);
    idx_from_background_mesh_.resize(1);
    interface_id_.resize(1);
    nb_quadrature_time_ = 1;
    

    not_in_active_mesh_.resize(21);
    for (int i = 0; i < 21; ++i)
        not_in_active_mesh_[i].resize(nb_quadrature_time_);


    int nt = 0;
    int nt_macros = 0; // number of active macro elements 

    const int domain = 0;

    const auto &Th_bary = static_cast<const BarycentricMesh2 &>(this->Th);

    macro_idx_in_background_mesh_[domain].reserve(Th_bary.macro_elements.size());

    for (int macro_k = 0; macro_k < Th_bary.macro_elements.size(); ++macro_k) {
        const auto &sub_elems = Th_bary.macro_elements[macro_k];

        bool keep_macro = false;
        for (int i = 0; i < 3; ++i) {
            int sub_k = sub_elems[i];
            const auto signK = interface.get_SignElement(sub_k);

            if (interface.isCut(sub_k)) {
                keep_macro = true;
                break;
            }
        }

        if (!keep_macro)
            continue;

        macro_idx_in_background_mesh_[domain].push_back(macro_k);
        macro_idx_from_background_mesh_[domain][macro_k] = nt_macros;
        
        std::array<int, 3> sub_elems_active;

        // Keep all three subelements
        for (int i = 0; i < 3; ++i) {
            int sub_k = sub_elems[i];
            const auto signK = interface.get_SignElement(sub_k);

            const int local_id = nt;
            sub_elems_active[i] = local_id;
            idx_in_background_mesh_[domain].push_back(sub_k);
            idx_from_background_mesh_[domain][sub_k] = local_id;
            inverse_active_macro_map_d[domain].push_back(nt_macros);


            if (interface.isCut(sub_k)) {
                interface_id_[0][{domain, local_id}].emplace_back(&interface, 0);
            } else {
                // std::cout << "sub_k = " << sub_k << "\n";
                not_in_active_mesh_[domain][0][local_id] = true; //! to exclude elements completely outside of domain from isInactive
            }

            ++nt;
            
        }
        active_macro_elements_d[domain].push_back(sub_elems_active);
        nt_macros++;
    }

    nb_active_macros_d[domain] = nt_macros;

    idx_element_domain.push_back(0);
    idx_element_domain.push_back(nt);
    idx_in_background_mesh_[domain].shrink_to_fit();

}

void BarycentricActiveMesh2::createSurfaceMesh(const TimeInterface<Mesh2> &interface) {
    
    // reset macro containers state
    macro_idx_in_background_mesh_.clear();
    macro_idx_from_background_mesh_.clear();
    inverse_active_macro_map_d.clear();
    active_macro_elements_d.clear();
    nb_active_macros_d.clear();

    // reset active mesh containers
    idx_in_background_mesh_.clear();
    idx_from_background_mesh_.clear();
    interface_id_.clear();
    idx_element_domain.clear();
    not_in_active_mesh_.clear();

    // initialize containers
    macro_idx_in_background_mesh_.resize(1);
    macro_idx_from_background_mesh_.resize(1);
    active_macro_elements_d.resize(1);
    inverse_active_macro_map_d.resize(1);
    nb_active_macros_d.resize(1);

    int n_tid = interface.size();
    nb_quadrature_time_ = n_tid;

    interface_id_.resize(21);
    interface_id_.assign(nb_quadrature_time_, {}); // time slices
    not_in_active_mesh_.resize(10);
    for (int i = 0; i < 10; ++i)
        not_in_active_mesh_[i].resize(nb_quadrature_time_);

    idx_in_background_mesh_.resize(1);
    idx_in_background_mesh_[0].reserve(Th.nt);
    idx_from_background_mesh_.resize(1);

    const int domain = 0;
    int nt = 0;
    int nt_macros = 0; // number of active macro elements 

    const auto &Th_bary = static_cast<const BarycentricMesh2 &>(this->Th);

    macro_idx_in_background_mesh_[domain].reserve(Th_bary.macro_elements.size());

    for (int macro_k = 0; macro_k < Th_bary.macro_elements.size(); ++macro_k) {
        const auto &sub_elems = Th_bary.macro_elements[macro_k];

        bool macro_is_active = false;

        // Check all subelements and time quadrature points  
        for (int i = 0; i < 3 && !macro_is_active; ++i) {
            int sub_k = sub_elems[i];
            for (int it = 0; it < n_tid - 1; ++it) {
                const auto signKi  = interface(it)->get_SignElement(sub_k);
                const auto signKii = interface(it + 1)->get_SignElement(sub_k);
                bool Ki_cut  = interface(it)->isCut(sub_k);
                bool Kii_cut = interface(it + 1)->isCut(sub_k);

                if (Ki_cut || Kii_cut || signKi.sign() * signKii.sign() <= 0) {
                    macro_is_active = true;
                    break;
                }
            }
        }

        if (!macro_is_active)
            continue;
        
        macro_idx_in_background_mesh_[domain].push_back(macro_k);
        macro_idx_from_background_mesh_[domain][macro_k] = nt_macros;
        
        std::array<int, 3> sub_elems_active;

        // Keep all subelements if macro is active
        for (int i = 0; i < 3; ++i) {
            int sub_k = sub_elems[i];

            int local_id = nt;
            sub_elems_active[i] = local_id;

            idx_in_background_mesh_[domain].push_back(sub_k);
            idx_from_background_mesh_[domain][sub_k] = local_id;
            inverse_active_macro_map_d[domain].push_back(nt_macros);


            for (int t = 0; t < n_tid; ++t) {
                bool is_cut = interface(t)->isCut(sub_k);

                if (is_cut) {
                    interface_id_[t][{domain, local_id}].emplace_back(interface[t], 0);
                } else {
                    not_in_active_mesh_[domain][t][local_id] = true;
                }
            }

            ++nt;
            
        }
        active_macro_elements_d[domain].push_back(sub_elems_active);
        nt_macros++;
    }

    nb_active_macros_d[domain] = nt_macros;
    
    idx_element_domain.push_back(0);
    idx_element_domain.push_back(nt);
    idx_in_background_mesh_[domain].shrink_to_fit();
        
    
}


void BarycentricActiveMesh2::init(const Interface<Mesh2>& interface)
{
    // Reserve memory for the indices in the background mesh for positive and negative domains
    idx_in_background_mesh_[0].reserve(Th.nt);
    idx_in_background_mesh_[1].reserve(Th.nt);

    // Push 0 as the first index for the element domain
    idx_element_domain.push_back(0);

    // Initialize the counters for positive and negative domains
    int nt0 = 0, nt1 = 0;
    int nt0_macro = 0, nt1_macro = 0;

    const auto &Th_bary = static_cast<const BarycentricMesh2 &>(this->Th);
    macro_idx_in_background_mesh_[0].reserve(Th_bary.macro_elements.size());
    macro_idx_in_background_mesh_[1].reserve(Th_bary.macro_elements.size());

    // Loop through all the macro elements in the background mesh
    // for (int k = 0; k < Th.nt; ++k) {
    for (int macro_k = 0; macro_k < Th_bary.macro_elements.size(); ++macro_k) {
        
        // Micro elements
        std::array<int, 3> micro_elems = Th_bary.macro_elements[macro_k];

        // Find out if element is cut
        bool is_cut = false;
        for (int i = 0; i < 3; ++i) {
            int micro_k = micro_elems[i];

            // macro is cut if any of its micros are cut
            if (interface.isCut(micro_k)) {
                is_cut = true;
                break;
            }
        }

        // Find out if which sign the element has of the level set function if not cut
        int sign = 0;      // get the sign of the interface     
        if (!is_cut) {
            // since the macro element is not cut, none of its micro elements are cut a
            // should all have the same sign 
            for (int i = 0; i < 3; ++i) {
                int micro_k = micro_elems[i];

                sign += interface.get_SignElement(micro_k).sign();  
            }
            sign /= 3;  // divide with number of micro-elements

            assert((sign == -1) || (sign == 1));
        }

        // If the element is cut, add it to both domains
        if (is_cut) {
            assert(sign == 0);
            // std::cout << "Macroelement " << macro_k << " is cut, adding index " << nt0_macro << " to domain 0 and index " << nt1_macro << " to domain 1\n";

            std::array<int, 3> active_micros_0;   // index of micro elements in active mesh 0
            std::array<int, 3> active_micros_1;   // index of micro elements in active mesh 1

            // Push the indices of the micro elements into both positive and negative domain index arrays
            macro_idx_in_background_mesh_[0].push_back(macro_k);
            macro_idx_from_background_mesh_[0][macro_k] = nt0_macro;
            macro_idx_in_background_mesh_[1].push_back(macro_k);
            macro_idx_from_background_mesh_[1][macro_k] = nt1_macro;

            for (int i = 0; i < 3; ++i) {

                int micro_k = micro_elems[i];
                idx_in_background_mesh_[0].push_back(micro_k);
                idx_from_background_mesh_[0][micro_k] = nt0;
                idx_in_background_mesh_[1].push_back(micro_k);
                idx_from_background_mesh_[1][micro_k] = nt1;
                active_micros_0[i] = nt0;
                active_micros_1[i] = nt1;
                inverse_active_macro_map_d[0].push_back(nt0_macro);
                inverse_active_macro_map_d[1].push_back(nt1_macro);

                // Update the interface_id array for both positive and negative domains
                if (interface.isCut(micro_k)) {
                    interface_id_[0][std::make_pair(0, nt0)].push_back(std::make_pair(&interface, 1));
                    interface_id_[0][std::make_pair(1, nt1)].push_back(std::make_pair(&interface, -1));
                } else if (interface.get_SignElement(micro_k).sign() < 0) {
                    // opposite sign -> inactive in domain 0
                    not_in_active_mesh_[0][0][nt0] = true;
                } else if (interface.get_SignElement(micro_k).sign() > 0) {
                    // inactive in domain 1
                    not_in_active_mesh_[1][0][nt1] = true;
                } else {
                    std::cout << "ERROR! QUTTING \n";
                    return;
                }

                // std::cout << "Microelement " << micro_k << " is cut, adding index " << nt0_macro << " to domain 0 and index " << nt1_macro << " to domain 1\n";
                
                // Increment the counters for both positive and negative domain elements
                nt0++;
                nt1++;
            }
            active_macro_elements_d[0].push_back(active_micros_0);
            active_macro_elements_d[1].push_back(active_micros_1);

            nt0_macro++;
            nt1_macro++;
        }
        // If the element is not cut, add it to the appropriate domain
        else {
        
            std::array<int, 3> active_micros; 

            // Set the counters to the appropriate domain counters
            int &nnt = (sign > 0) ? nt0 : nt1;
            int &nnt_macro = (sign > 0) ? nt0_macro : nt1_macro;
        
            // std::cout << "Macro element " << macro_k << " has sign " << sign << ", adding index " << nnt_macro << " to domain " << (sign < 0) << "\n";

            macro_idx_in_background_mesh_[(sign < 0)].push_back(macro_k);
            macro_idx_from_background_mesh_[(sign < 0)][macro_k] = nnt_macro;

            for (int i = 0; i < 3; ++i) {
                int micro_k = micro_elems[i];
                
                // Push the indices of the corresponding background mesh micro elements into the appropriate domain index array
                idx_in_background_mesh_[(sign < 0)].push_back(micro_k);
                // Push the index of the corresponding active mesh micro elements into the appropriate domain index array
                idx_from_background_mesh_[(sign < 0)][micro_k] = nnt;

                active_micros[i] = nnt;
                inverse_active_macro_map_d[(sign < 0)].push_back(nnt_macro);

                // std::cout << "Micro element " << micro_k << " has sign " << sign << ", pushing back index " << nnt_macro << " to domain " << (sign < 0) << "\n";
                nnt++;
            }

            active_macro_elements_d[(sign < 0)].push_back(active_micros);            
            nnt_macro++;
        
        }
    }

    nb_active_macros_d[0] = nt0_macro;
    nb_active_macros_d[1] = nt1_macro;

    // Shrink the size of the index arrays to the actual size
    idx_in_background_mesh_[0].resize(nt0);
    idx_in_background_mesh_[1].resize(nt1);
    idx_in_background_mesh_[0].shrink_to_fit();
    idx_in_background_mesh_[1].shrink_to_fit();

    // Push the total number of elements in the positive and negative domains as the last index of the element domain
    // array
    assert(nt0 == 3*active_macro_elements_d[0].size());
    assert(nt1 == 3*active_macro_elements_d[1].size());
    idx_element_domain.push_back(nt0);
    idx_element_domain.push_back(nt0 + nt1);
    
}


int BarycentricActiveMesh2::macro_of_micro(int k) const {
    const int d     = this->get_domain_element(k);            // 0 or 1, uses idx_element_domain
    const int k_loc = k - this->idx_element_domain[d];
    assert(0 <= k_loc && k_loc < (int)inverse_active_macro_map_d[d].size());
    return inverse_active_macro_map_d[d][k_loc];
}


/*
k_active = active macro element index
*/
int BarycentricActiveMesh2::get_macro_in_background_mesh(int k_active) const {
    const int domain = this->get_domain_element(k_active);
    assert(0 <= k_active && k_active < (int)active_macro_elements_d[domain].size());

    return macro_idx_in_background_mesh_[domain][k_active];
}

int BarycentricActiveMesh2::get_macro_in_active_mesh(int k_bg, int domain) const {
    if (this->get_nb_domain() == 1) {
        domain = 0;
    }
    auto it = macro_idx_from_background_mesh_[domain].find(k_bg);
    if (it == macro_idx_from_background_mesh_[domain].end()) {
        return -1;
    }
    return this->idx_element_domain[domain] + it->second; // 
}


bool BarycentricActiveMesh2::is_macro_cut(int macro_k, int domain, int t) const {
    assert(0 <= macro_k && macro_k < (int)active_macro_elements_d[domain].size());
    std::array<int,3> sub = active_macro_elements_d[domain][macro_k]; // local micro ids in domain d
    const int shift = idx_element_domain[domain];              // 0 for d=0, nt0 for d=1
    // Loop through the micro elements
    for (int k_loc : sub) {
        const int k_global = k_loc + shift; // map to global active index from per-domain index
        if (this->isCut(k_global, t)) return true;
    }
    return false;
}




bool BarycentricActiveMesh2::is_macro_cut(int macro_k, int t /*=0*/) const {
    const int domain = 0;   // assumes single domain
    assert(0 <= macro_k && macro_k < (int)active_macro_elements_d[domain].size());
    std::array<int,3> sub = active_macro_elements_d[domain][macro_k]; // 3 active micro ids

    for (int k_micro : sub) {
        assert(this->get_domain_element(k_micro) == domain);
        if (this->isCut(k_micro, t)) return true;
    }
        
    return false;
}

// "interior micro" := !cut && !inactive
bool BarycentricActiveMesh2::is_macro_interior(int macro_k, int t /*=0*/) const {
    const int domain = 0;
    assert(0 <= macro_k && macro_k < (int)active_macro_elements_d[domain].size());
    std::array<int,3> sub = active_macro_elements_d[domain][macro_k];

    for (int k_micro : sub) {
        assert(this->get_domain_element(k_micro) == domain);
        if (this->isCut(k_micro, t) || this->isInactive(k_micro, t)) return false;
    }
        
    return true;
}

// "macro has no active micros" at time t
bool BarycentricActiveMesh2::is_macro_inactive(int macro_k, int t /*=0*/) const {
    const int domain = 0;
    assert(0 <= macro_k && macro_k < (int)active_macro_elements_d[domain].size());
    std::array<int,3> sub = active_macro_elements_d[domain][macro_k];

    for (int k_micro : sub) {
        assert(this->get_domain_element(k_micro) == domain);
        if (!this->isInactive(k_micro, t)) return false;
    }
        
    return true;
}

// Stabilize if macro is cut OR inactive at ANY time instance
bool BarycentricActiveMesh2::stabilize_macro(int macro_k) const {
    const int domain = 0;
    if (!(0 <= macro_k && macro_k < (int)active_macro_elements_d[domain].size())) {
        std::cout << "macro_k = " << macro_k << " index out of range";
        assert(0);
    }
    // assert(0 <= macro_k && macro_k < (int)active_macro_elements.size());
    for (int t = 0; t < nb_quadrature_time_; ++t)
        if (is_macro_cut(macro_k, t) || is_macro_inactive(macro_k, t)) return true;
    return false;
}

// a macro element is exterior if it's inactive in each time instance (this can only happen using the global_truncate function for active meshes)
bool BarycentricActiveMesh2::is_macro_exterior(int macro_k) const {
    const int domain = 0;
    if (!(0 <= macro_k && macro_k < (int)active_macro_elements_d[domain].size())) {
        std::cout << "macro_k = " << macro_k << " index out of range";
        assert(0);
    }
    // assert(0 <= macro_k && macro_k < (int)active_macro_elements.size());
    for (int t = 0; t < nb_quadrature_time_; ++t)
        if (!is_macro_inactive(macro_k, t)) return false;
    return true;
}


int BarycentricActiveMesh2::macro_adjacent(const int macro_k, const int iface_adj) const {
    const int domain = 0;
    if (!(0 <= macro_k && macro_k < (int)active_macro_elements_d[domain].size())) {
        std::cout << "macro_k = " << macro_k << " index out of range";
        assert(0);
    }
    int macro_kb = get_macro_in_background_mesh(macro_k);
 
    const auto& Th_bary = static_cast<const BarycentricMesh2&>(this->Th);

    int macro_kbn = Th_bary.element_adj(macro_kb, iface_adj);

    if (macro_kbn == -1)    
        return -1;
    
    return get_macro_in_active_mesh(macro_kbn, domain);

}

