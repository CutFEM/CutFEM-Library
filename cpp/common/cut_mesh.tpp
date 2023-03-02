/**
 * @file cpp/common/cut_mesh.tpp
 *
 * @brief Contains implementation of methods
 *
 */
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

#ifndef COMMON_CUT_MESH_TPP
#define COMMON_CUT_MESH_TPP

//* --- Cut_Part class --- *//

template <typename E> Cut_Part<E>::Cut_Part(const Partition<E> p, int s) : ip(p), sign_cut_(s), pp(p.T) {
    partition_ = &ip;
}

template <typename E> Cut_Part<E>::Cut_Part(const Physical_Partition<E> p, int s) : pp(p), sign_cut_(s), ip(p.T) {
    partition_ = &pp;
}

template <typename E> Cut_Part<E>::Cut_Part(const Cut_Part<E> &p) : pp(p.pp), ip(p.ip) {
    if (p.partition_ == &p.pp)
        partition_ = &pp;
    else
        partition_ = &ip;
}

template <typename E> Cut_Part<E> &Cut_Part<E>::operator=(const Cut_Part<E> &p) {
    // Copy the partition and the index from the given Cut_Part.
    pp = p.pp;
    ip = p.ip;
    // If the given Cut_Part is using the
    // partition, use the copy of the partition.
    if (p.partition_ == &p.pp)
        partition_ = &pp;
    else
        // Otherwise, use the copy of the index.
        partition_ = &ip;
    return *this;
}

// GETTERS
template <typename E> int Cut_Part<E>::get_sign() const { return sign_cut_; }

template <typename E> int Cut_Part<E>::get_sign_node(int i) const { return partition_->get_sign_node(i); }

template <typename E> void Cut_Part<E>::get_list_node(std::vector<typename E::Rd> &node) const {
    partition_->get_list_node(node, sign_cut_);
}

template <typename E> CutElement<E> Cut_Part<E>::get_element(int k) const { return partition_->get_element(k); }

template <typename E> typename Cut_Part<E>::Rd Cut_Part<E>::get_vertex(const_element_iterator it, const int i) const {
    return partition_->get_vertex(it, i);
}

template <typename E> int Cut_Part<E>::get_nb_element() const { return partition_->nb_element(sign_cut_); }

template <typename E> int Cut_Part<E>::get_local_domain_id() const {
    if (sign_cut_ == 0)
        return -1; // not cout
    else
        return (sign_cut_ == -1);
}

// OTHER METHODS
// GIVE THE MEASURE OF THE CUT PART IN Rd

template <typename E> double Cut_Part<E>::measure() const { return partition_->measure(sign_cut_); }

template <typename E> double Cut_Part<E>::measure(const_element_iterator it) const { return partition_->measure(it); }

// //GIVE THE MEASURE OF CUT PART OF A FACE IN RdBord
// double measureBord(int ifac) const {return
// partition_->measureBord(sign_cut_, ifac);}

template <typename E> bool Cut_Part<E>::multi_interface() const { return partition_ == &pp; }

template <typename E>
typename Cut_Part<E>::Rd Cut_Part<E>::mapToPhysicalElement(const_element_iterator it, const RdHat Phat) const {
    return partition_->mapToPhysicalElement(it, Phat);
}

// ITERATORS

template <typename E> typename Cut_Part<E>::const_element_iterator Cut_Part<E>::element_begin() const {
    return partition_->element_begin(sign_cut_);
}

template <typename E> typename Cut_Part<E>::const_element_iterator Cut_Part<E>::element_end() const {
    return partition_->element_end(sign_cut_);
}

template <typename E> typename Cut_Part<E>::const_element_iterator Cut_Part<E>::other_side_element_begin() const {
    assert(sign_cut_ != 0);
    int other_side_cut = (sign_cut_ == 1);
    return partition_->element_begin(other_side_cut);
}

template <typename E> typename Cut_Part<E>::const_element_iterator Cut_Part<E>::other_side_element_end() const {
    assert(sign_cut_ != 0);
    int other_side_cut = (sign_cut_ == 1);
    return partition_->element_end(other_side_cut);
}

//* --- ActiveMesh class --- *//

//* Public Members *//

template <typename Mesh> ActiveMesh<Mesh>::ActiveMesh(const Mesh &th) : Th(th) {
    idx_in_background_mesh_.reserve(10);
    idx_from_background_mesh_.reserve(10);
    idx_in_background_mesh_.resize(1);
    idx_from_background_mesh_.resize(1);
    idx_in_background_mesh_[0].resize(Th.nt);

    interface_id_.resize(10);

    nb_quadrature_time_ = 1;        // by default, the active mesh of the background mesh is stationary
    
    // set the active mesh indexing to the same as the background element indexing
    for (int k = 0; k < Th.nt; ++k) {
        idx_in_background_mesh_[0][k]   = k;
        idx_from_background_mesh_[0][k] = k;
    }

    idx_element_domain.push_back(0);
    idx_element_domain.push_back(Th.nt);
    in_active_mesh_.resize(10);
    for (int i = 0; i < 10; ++i)
        in_active_mesh_[i].resize(nb_quadrature_time_);
}

// Give the background mesh and a sign Function defined on the mesh nodes
// Will create 2 subdomains
template <typename Mesh> ActiveMesh<Mesh>::ActiveMesh(const Mesh &th, const Interface<Mesh> &interface) : Th(th) {
    idx_in_background_mesh_.reserve(10);
    idx_from_background_mesh_.reserve(10);
    idx_in_background_mesh_.resize(2);
    idx_from_background_mesh_.resize(2);
    interface_id_.resize(1);
    nb_quadrature_time_ = 1;
    this->init(interface);
}

template <typename Mesh> ActiveMesh<Mesh>::ActiveMesh(const Mesh &th, const TimeInterface<Mesh> &interface) : Th(th) {
    idx_in_background_mesh_.reserve(10);
    idx_from_background_mesh_.reserve(10);
    idx_in_background_mesh_.resize(2);
    idx_from_background_mesh_.resize(2);
    nb_quadrature_time_ = interface.size();
    interface_id_.resize(nb_quadrature_time_);
    in_active_mesh_.resize(10);
    for (int i = 0; i < 10; ++i)
        in_active_mesh_[i].resize(nb_quadrature_time_);
    this->init(interface);
}



template <typename Mesh> void ActiveMesh<Mesh>::truncate(const Interface<Mesh> &interface, int sign_domain_remove) {

    // Get number of subdomains of resulting mesh //?
    int dom_size = this->get_nb_domain();
    idx_element_domain.resize(0);

    {
        // Iterate through number of remaining subdomains
        for (int d = 0; d < dom_size; ++d) {
            idx_in_background_mesh_[d].resize(0);
            // Compute number of elements in subdomain d
            int nt_max = idx_from_background_mesh_[d].size();
            // Reserve memory for these elements
            idx_in_background_mesh_[d].reserve(nt_max);
        }
    }

    // Vector to hold number of elements in each subdomain
    std::vector<int> nt(dom_size, 0.);
    for (int d = 0; d < dom_size; ++d) {
        for (auto it_k = idx_from_background_mesh_[d].begin(); it_k != idx_from_background_mesh_[d].end();) {

            int kb = it_k->first;   // background mesh element index
            int k  = it_k->second;  // active mesh element index

            // Get interface segment
            auto it_gamma                                               = interface_id_[0].find(std::make_pair(d, k));
            // Get the sign of the level set function in the element kb
            const SignElement<typename ActiveMesh<Mesh>::Element> signK = interface.get_SignElement(kb);

            // Remove the Mesh::Element in the domain corresopnding to sign_domain_remove
            if (signK.sign() == sign_domain_remove) {

                it_k = idx_from_background_mesh_[d].erase(it_k);

                continue;
            }

            // Save and erase old interfaces 
            int nb_interface = (it_gamma == interface_id_[0].end()) ? 0 : it_gamma->second.size();
            std::vector<const Interface<Mesh> *> old_interface(nb_interface);
            std::vector<int> ss(nb_interface);
            for (int i = 0; i < nb_interface; ++i)
                old_interface[i] = it_gamma->second[i].first;
            for (int i = 0; i < nb_interface; ++i)
                ss[i] = it_gamma->second[i].second;
            if (it_gamma != interface_id_[0].end()) {
                auto ittt = interface_id_[0].erase(it_gamma);
            }

            // SET NEW INDICES AND PUT BACK INTERFACES
            idx_in_background_mesh_[d].push_back(kb);
            it_k->second = nt[d];
            for (int i = 0; i < nb_interface; ++i) {
                interface_id_[0][std::make_pair(d, nt[d])].push_back(std::make_pair(old_interface[i], ss[i]));
            }
            // IS CUT SO NEED TO ADD INTERFACE AND SIGN
            if (signK.cut()) {
                interface_id_[0][std::make_pair(d, nt[d])].push_back(std::make_pair(&interface, -sign_domain_remove));
            }
            nt[d]++;
            it_k++;
        }
    }

    idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        idx_in_background_mesh_[d].resize(nt[d]);
        idx_in_background_mesh_[d].shrink_to_fit();
        int sum_nt = idx_element_domain[d] + nt[d];
        idx_element_domain.push_back(sum_nt);
    }
}

template <typename Mesh> void ActiveMesh<Mesh>::truncate(const TimeInterface<Mesh> &interface, int sign_domain_remove) {

    int n_tid = interface.size();
    assert(n_tid < interface_id_.size());
    nb_quadrature_time_ = n_tid;
    in_active_mesh_.resize(10);
    for (int i = 0; i < 10; ++i)
        in_active_mesh_[i].resize(nb_quadrature_time_);

    int dom_size = this->get_nb_domain();
    idx_element_domain.resize(0);

    {
        for (int d = 0; d < dom_size; ++d) {
            idx_in_background_mesh_[d].resize(0);
            int nt_max = idx_from_background_mesh_[d].size();
            idx_in_background_mesh_[d].reserve(nt_max);
        }
    }

    std::vector<int> nt(dom_size, 0.); //! nt is never changed after this line? why are its components just 0.?

    // Loop over all subdomains
    for (int d = 0; d < dom_size; ++d) {
        // Loop over all elements in the active mesh of subdomain d
        for (auto it_k = idx_from_background_mesh_[d].begin(); it_k != idx_from_background_mesh_[d].end();) {

            // Index of the element in the background mesh
            int kb = it_k->first;
            // Index of the element in the active mesh
            int k  = it_k->second;

            // Variable to check status if the element is active or inactive
            bool active_element = false;
            // Temporary variable to hold the sign of the element
            int s;

            // Loop over all time quadrature points
            for (int t = 0; t < interface.size() - 1; ++t) {
                // Get SignElement at time t
                const SignElement<typename ActiveMesh<Mesh>::Element> signKi  = interface(t)->get_SignElement(kb);
                // Get SignElement at time t+1
                const SignElement<typename ActiveMesh<Mesh>::Element> signKii = interface(t + 1)->get_SignElement(kb);
                // Save the sign of element at time t
                s                                                             = signKi.sign();
                // Check if the element is cut in either t or t+1 or if the sign of the element changes
                // -> that means it's active
                if (signKi.cut() || signKii.cut() || signKi.sign() * signKii.sign() <= 0) {
                    active_element = true;
                    break;
                }
            }

            //! CONTINUE HERE
            // REMOVE THE Mesh::Element IN THE INPUT DOMAIN
            if (s == sign_domain_remove && !active_element) {
                it_k = idx_from_background_mesh_[d].erase(it_k);
                continue;
            }

            // SAVE AND ERASE OLD INTERFACES
            for (int it = 0; it < n_tid; ++it) {
                auto it_gamma = interface_id_[it].find(std::make_pair(d, k));

                int nb_interface = (it_gamma == interface_id_[it].end()) ? 0 : it_gamma->second.size();
                std::vector<const Interface<Mesh> *> old_interface(nb_interface);
                std::vector<int> ss(nb_interface);
                for (int i = 0; i < nb_interface; ++i)
                    old_interface[i] = it_gamma->second[i].first;
                for (int i = 0; i < nb_interface; ++i)
                    ss[i] = it_gamma->second[i].second;
                if (it_gamma != interface_id_[it].end()) {
                    auto ittt = interface_id_[it].erase(it_gamma);
                }
                // PUT BACK INTERFACES
                for (int i = 0; i < nb_interface; ++i) {
                    interface_id_[it][std::make_pair(d, nt[d])].push_back(std::make_pair(old_interface[i], ss[i]));
                }
                // IS CUT SO NEED TO ADD INTERFACE AND SIGN
                const SignElement<typename ActiveMesh<Mesh>::Element> signK = interface(it)->get_SignElement(kb);
                if (signK.cut()) {
                    interface_id_[it][std::make_pair(d, nt[d])].push_back(
                        std::make_pair(interface[it], -sign_domain_remove));
                } else if (signK.sign() == sign_domain_remove) {
                    in_active_mesh_[d][it][nt[d]] = false;
                }
            }
            // // SET NEW INDICES AND PUT BACK INTERFACES
            idx_in_background_mesh_[d].push_back(kb);
            it_k->second = nt[d];
            nt[d]++;
            it_k++;
        }
    }
    return;
    idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        idx_in_background_mesh_[d].resize(nt[d]);
        idx_in_background_mesh_[d].shrink_to_fit();
        int sum_nt = idx_element_domain[d] + nt[d];
        idx_element_domain.push_back(sum_nt);
    }
}

//  TODO Comment: construct a subdomain corresponding to the positive sign
template <typename Mesh> void ActiveMesh<Mesh>::add(const Interface<Mesh> &interface, int sign_domain) {

    int dom_size = this->get_nb_domain();
    idx_element_domain.resize(0);
    // int sign_domain = -1;
    // Initialize the first new subdomain domain

    // and clear old array with Mesh::Element indices.
    {
        idx_in_background_mesh_.resize(dom_size + 1);
        idx_from_background_mesh_.resize(dom_size + 1);
        for (int d = 0; d < dom_size + 1; ++d) {
            idx_in_background_mesh_[d].resize(0);
        }
        int nt_max = idx_from_background_mesh_[0].size();
        idx_in_background_mesh_[dom_size].reserve(nt_max);
    }
    std::vector<int> nt(2 * dom_size, 0.);
    int new_dom_id = dom_size;
    for (int d = 0; d < dom_size; ++d) {
        bool sub_is_cut = false;
        for (auto it_k = idx_from_background_mesh_[d].begin(); it_k != idx_from_background_mesh_[d].end();) {

            int kb = it_k->first;
            int k  = it_k->second;

            auto it_gamma                                               = interface_id_[0].find(std::make_pair(d, k));
            const SignElement<typename ActiveMesh<Mesh>::Element> signK = interface.get_SignElement(kb);

            int nb_interface = (it_gamma == interface_id_[0].end()) ? 0 : it_gamma->second.size();
            std::vector<const Interface<Mesh> *> old_interface(nb_interface);
            std::vector<int> ss(nb_interface);
            for (int i = 0; i < nb_interface; ++i)
                old_interface[i] = it_gamma->second[i].first;
            for (int i = 0; i < nb_interface; ++i)
                ss[i] = it_gamma->second[i].second;

            if (it_gamma != interface_id_[0].end()) {
                auto ittt = interface_id_[0].erase(it_gamma);
            }

            if (signK.sign() == sign_domain || signK.cut()) {

                // Initialize first time we find a cut Mesh::Element
                if (!sub_is_cut && new_dom_id != dom_size) {
                    new_dom_id++;
                    idx_in_background_mesh_.resize(new_dom_id + 1);
                    idx_from_background_mesh_.resize(new_dom_id + 1);
                    idx_in_background_mesh_[new_dom_id].resize(0);
                    int nt_max = idx_from_background_mesh_[d + 1].size();
                    idx_in_background_mesh_[new_dom_id].reserve(nt_max);
                }

                sub_is_cut = true;
                idx_in_background_mesh_[new_dom_id].push_back(kb);
                idx_from_background_mesh_[new_dom_id][kb] = nt[new_dom_id];

                for (int i = 0; i < nb_interface; ++i) {
                    interface_id_[0][std::make_pair(new_dom_id, nt[new_dom_id])].push_back(
                        std::make_pair(old_interface[i], ss[i]));
                }

                if (!signK.cut()) {
                    it_k = idx_from_background_mesh_[d].erase(it_k);
                } else {
                    idx_in_background_mesh_[d].push_back(kb);
                    it_k->second = nt[d];

                    // need to change k and add new interface
                    // attach all interface to new Mesh::Element
                    for (int i = 0; i < nb_interface; ++i) {
                        interface_id_[0][std::make_pair(d, nt[d])].push_back(std::make_pair(old_interface[i], ss[i]));
                    }
                    interface_id_[0][std::make_pair(new_dom_id, nt[new_dom_id])].push_back(
                        std::make_pair(&interface, sign_domain));
                    interface_id_[0][std::make_pair(d, nt[d])].push_back(std::make_pair(&interface, -sign_domain));
                    nt[d]++;
                    it_k++;
                }

                nt[new_dom_id]++;
            } else {
                // std::cout << " in old domain " << std::endl;
                idx_in_background_mesh_[d].push_back(kb);
                it_k->second = nt[d];

                // change the key of the interface_id map
                for (int i = 0; i < nb_interface; ++i) {
                    interface_id_[0][std::make_pair(d, nt[d])].push_back(std::make_pair(old_interface[i], ss[i]));
                }
                it_k++;
                nt[d]++;
            }
        }

        // if(sub_is_cut && d+1 !=dom_size){
        //    new_dom_id++;
        //    idx_in_background_mesh_.resize(new_dom_id+1);
        //    idx_from_background_mesh_.resize(new_dom_id+1);
        //    idx_in_background_mesh_[new_dom_id].resize(0);
        //    int nt_max = idx_from_background_mesh_[d+1].size();
        //    idx_in_background_mesh_[new_dom_id].reserve(nt_max);
        //  }
    }

    idx_element_domain.push_back(0);
    for (int d = 0; d < new_dom_id + 1; ++d) {
        idx_in_background_mesh_[d].resize(nt[d]);
        idx_in_background_mesh_[d].shrink_to_fit();
        int sum_nt = idx_element_domain[d] + nt[d];
        idx_element_domain.push_back(sum_nt);
    }
}

template <typename Mesh> void ActiveMesh<Mesh>::createSurfaceMesh(const Interface<Mesh> &interface) {

    int dom_size = this->get_nb_domain();
    idx_element_domain.resize(0);

    {
        for (int d = 0; d < dom_size; ++d) {
            idx_in_background_mesh_[d].resize(0);
            int nt_max = idx_from_background_mesh_[d].size();
            idx_in_background_mesh_[d].reserve(nt_max);
        }
    }

    std::vector<int> nt(dom_size, 0.);
    for (int d = 0; d < dom_size; ++d) {
        for (auto it_k = idx_from_background_mesh_[d].begin(); it_k != idx_from_background_mesh_[d].end();) {

            int kb = it_k->first;
            int k  = it_k->second;

            // std::cout << "domain \t" << d << " Mesh::Element back " << kb << "\t =>
            // loc id " << k << std::endl;
            auto it_gamma                                               = interface_id_[0].find(std::make_pair(d, k));
            const SignElement<typename ActiveMesh<Mesh>::Element> signK = interface.get_SignElement(kb);

            // REMOVE THE Mesh::Element IN THE INPUT DOMAIN
            if (!signK.cut()) {

                it_k = idx_from_background_mesh_[d].erase(it_k);
                continue;
            }

            // SAVE AND ERASE OLD INTERFACES
            int nb_interface = (it_gamma == interface_id_[0].end()) ? 0 : it_gamma->second.size();
            std::vector<const Interface<Mesh> *> old_interface(nb_interface);
            std::vector<int> ss(nb_interface);
            for (int i = 0; i < nb_interface; ++i)
                old_interface[i] = it_gamma->second[i].first;
            for (int i = 0; i < nb_interface; ++i)
                ss[i] = it_gamma->second[i].second;
            if (it_gamma != interface_id_[0].end()) {
                auto ittt = interface_id_[0].erase(it_gamma);
            }

            // SET NEW INDICES AND PUT BACK INTERFACES
            idx_in_background_mesh_[d].push_back(kb);
            it_k->second = nt[d];
            for (int i = 0; i < nb_interface; ++i) {
                interface_id_[0][std::make_pair(d, nt[d])].push_back(std::make_pair(old_interface[i], ss[i]));
            }
            // IS CUT SO NEED TO ADD INTERFACE AND SIGN
            interface_id_[0][std::make_pair(d, nt[d])].push_back(std::make_pair(&interface, 0));
            nt[d]++;
            it_k++;
        }
    }

    idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        idx_in_background_mesh_[d].resize(nt[d]);
        idx_in_background_mesh_[d].shrink_to_fit();
        int sum_nt = idx_element_domain[d] + nt[d];
        idx_element_domain.push_back(sum_nt);
    }
}

template <typename Mesh> void ActiveMesh<Mesh>::createSurfaceMesh(const TimeInterface<Mesh> &interface) {

    int n_tid           = interface.size();
    nb_quadrature_time_ = n_tid;
    in_active_mesh_.resize(10);
    for (int i = 0; i < 10; ++i)
        in_active_mesh_[i].resize(nb_quadrature_time_);

    int dom_size = this->get_nb_domain();
    idx_element_domain.resize(0);

    {
        for (int d = 0; d < dom_size; ++d) {
            idx_in_background_mesh_[d].resize(0);
            int nt_max = idx_from_background_mesh_[d].size();
            idx_in_background_mesh_[d].reserve(nt_max);
        }
    }

    std::vector<int> nt(dom_size, 0.);
    for (int d = 0; d < dom_size; ++d) {
        for (auto it_k = idx_from_background_mesh_[d].begin(); it_k != idx_from_background_mesh_[d].end();) {

            int kb = it_k->first;
            int k  = it_k->second;

            bool active_element = false;
            for (int t = 0; t < interface.size() - 1; ++t) {
                const SignElement<typename ActiveMesh<Mesh>::Element> signKi  = interface(t)->get_SignElement(kb);
                const SignElement<typename ActiveMesh<Mesh>::Element> signKii = interface(t + 1)->get_SignElement(kb);

                if (signKi.cut() || signKii.cut() || signKi.sign() * signKii.sign() <= 0) {
                    active_element = true;
                    break;
                }
            }

            if (!active_element) {
                it_k = idx_from_background_mesh_[d].erase(it_k);
                continue;
            }

            // SAVE AND ERASE OLD INTERFACES EACH TIME STEP
            for (int it = 0; it < n_tid; ++it) {
                auto it_gamma    = interface_id_[it].find(std::make_pair(d, k));
                int nb_interface = (it_gamma == interface_id_[it].end()) ? 0 : it_gamma->second.size();
                std::vector<const Interface<Mesh> *> old_interface(nb_interface);
                std::vector<int> ss(nb_interface);
                for (int i = 0; i < nb_interface; ++i)
                    old_interface[i] = it_gamma->second[i].first;
                for (int i = 0; i < nb_interface; ++i)
                    ss[i] = it_gamma->second[i].second;
                if (it_gamma != interface_id_[it].end()) {
                    auto ittt = interface_id_[it].erase(it_gamma);
                }
                // PUT BACK INTERFACES
                for (int i = 0; i < nb_interface; ++i) {
                    interface_id_[it][std::make_pair(d, nt[d])].push_back(std::make_pair(old_interface[i], ss[i]));
                }
                // IS CUT SO NEED TO ADD INTERFACE AND SIGN
                const SignElement<typename ActiveMesh<Mesh>::Element> signK = interface(it)->get_SignElement(kb);

                if (signK.cut()) {
                    interface_id_[it][std::make_pair(d, nt[d])].push_back(std::make_pair(interface[it], 0));
                } else {
                    in_active_mesh_[d][it][nt[d]] = false;
                }
            }
            it_k->second = nt[d];
            idx_in_background_mesh_[d].push_back(kb);
            nt[d]++;
            it_k++;
        }
    }

    idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        idx_in_background_mesh_[d].resize(nt[d]);
        idx_in_background_mesh_[d].shrink_to_fit();
        int sum_nt = idx_element_domain[d] + nt[d];
        idx_element_domain.push_back(sum_nt);
    }
}



//* Private Members *//


//  constructor for basic 2 subdomains problem {1, -1}
template <typename Mesh> void ActiveMesh<Mesh>::init(const TimeInterface<Mesh> &interface) {

    int n_tid           = interface.size();
    nb_quadrature_time_ = n_tid;
    in_active_mesh_.resize(10);
    for (int i = 0; i < 10; ++i)
        in_active_mesh_[i].resize(nb_quadrature_time_);

    idx_in_background_mesh_[0].reserve(Th.nt);
    idx_in_background_mesh_[1].reserve(Th.nt);
    idx_element_domain.push_back(0);
    int nt0 = 0, nt1 = 0;
    for (int k = 0; k < Th.nt; ++k) {

        bool active_element = false;
        int s;
        for (int t = 0; t < interface.size() - 1; ++t) {
            const SignElement<typename ActiveMesh<Mesh>::Element> signKi  = interface(t)->get_SignElement(k);
            const SignElement<typename ActiveMesh<Mesh>::Element> signKii = interface(t + 1)->get_SignElement(k);
            s                                                             = signKi.sign();
            if (signKi.cut() || signKii.cut() || signKi.sign() * signKii.sign() <= 0) {
                active_element = true;
                break;
            }
        }

        if (active_element) {
            for (int it = 0; it < n_tid; ++it) {
                const SignElement<typename ActiveMesh<Mesh>::Element> signK = interface(it)->get_SignElement(k);
                int st                                                      = signK.sign();
                if (!signK.cut() && st == -1) {
                    in_active_mesh_[0][it][nt0] = false;
                }
                if (!signK.cut() && st == 1) {
                    in_active_mesh_[1][it][nt1] = false;
                }

                if (signK.cut()) {
                    interface_id_[it][std::make_pair(0, nt0)].push_back(std::make_pair(interface[it], 1));
                    interface_id_[it][std::make_pair(1, nt1)].push_back(std::make_pair(interface[it], -1));
                }
            }
            idx_in_background_mesh_[0].push_back(k);
            idx_from_background_mesh_[0][k] = nt0;
            idx_in_background_mesh_[1].push_back(k);
            idx_from_background_mesh_[1][k] = nt1;
            nt0++;
            nt1++;
        } else {
            int dom_add = (s < 0);
            int dom_rm  = (s > 0);
            int &nnt    = (s > 0) ? nt0 : nt1;
            idx_in_background_mesh_[dom_add].push_back(k);
            idx_from_background_mesh_[dom_add][k] = nnt;
            nnt++;
        }
    }
    idx_in_background_mesh_[0].resize(nt0);
    idx_in_background_mesh_[1].resize(nt1);
    idx_in_background_mesh_[0].shrink_to_fit();
    idx_in_background_mesh_[1].shrink_to_fit();

    idx_element_domain.push_back(nt0);
    idx_element_domain.push_back(nt0 + nt1);
}




// Check if a given cell index exists in the active mesh
template <typename Mesh> bool ActiveMesh<Mesh>::check_exist(int k, int dom) const {
    // Get the iterator corresponding to the given key
    const auto it = idx_from_background_mesh_[dom].find(k);

    // If the iterator is equal to the end of the map, then the key does not exist
    if (it == idx_from_background_mesh_[dom].end())
        return false;
    else
        return true;
}

/**
 * @brief Get the Element of a given index in the active mesh
 *
 * @tparam Mesh
 * @param i index in the active mesh
 * @return const ActiveMesh<Mesh>::Element&
 */
template <typename Mesh> const typename ActiveMesh<Mesh>::Element &ActiveMesh<Mesh>::operator[](int i) const {

    // Get the index of the Mesh::Element in the back mesh
    int k = idxElementInBackMesh(i);

    // Return the Mesh::Element
    return Th[k];
}

template <typename Mesh> int ActiveMesh<Mesh>::idxK_begin(int i) const { return this->idx_element_domain[i]; }

template <typename Mesh> int ActiveMesh<Mesh>::idxK_in_domain(int k, int i) const {
    return k - this->idx_element_domain[i];
}

template <typename Mesh> int ActiveMesh<Mesh>::nbElmts() const { return get_nb_element(); }

template <typename Mesh> int ActiveMesh<Mesh>::NbElement() const { return get_nb_element(); }

template <typename Mesh> int ActiveMesh<Mesh>::nbBrdElmts() const { return Th.nbe; }

template <typename Mesh> int ActiveMesh<Mesh>::nbVertices() const { return Th.nv; }

template <typename Mesh> const typename ActiveMesh<Mesh>::BorderElement &ActiveMesh<Mesh>::be(int i) const {
    return Th.be(i);
}

template <typename Mesh> int ActiveMesh<Mesh>::get_nb_domain() const { return idx_in_background_mesh_.size(); }

template <typename Mesh> int ActiveMesh<Mesh>::get_nb_element(int i) const {
    assert(i >= 0 && i < this->get_nb_domain());
    return idx_in_background_mesh_[i].size();
}

template <typename Mesh> int ActiveMesh<Mesh>::get_nb_element() const {
    int s = 0;
    for (int i = 0; i < this->get_nb_domain(); ++i) {
        s += get_nb_element(i);
    }
    return s;
}

/**
 * @brief Get the domain of a given element in the active mesh //! What if element belongs to several domains????
 *
 * @tparam Mesh Mesh
 * @param k element index in active mesh
 * @return int
 */
template <typename Mesh> int ActiveMesh<Mesh>::get_domain_element(const int k) const {
    for (int i = 0; i < this->get_nb_domain(); ++i) {
        if (k < idx_element_domain[i + 1]) {
            return i;
        }
    }
    std::cout << " Mesh::Element " << k << std::endl;
    assert(0);
}

template <typename Mesh> bool ActiveMesh<Mesh>::isCut(int k, int itq) const {
    int domain = get_domain_element(k);
    int kloc   = idxK_in_domain(k, domain);
    auto it    = interface_id_[itq].find(std::make_pair(domain, kloc));
    return (it != interface_id_[itq].end());
}

template <typename Mesh> bool ActiveMesh<Mesh>::isCutFace(int k, int ifac, int itq) const {
    // check if the element is cut
    if (!isCut(k, itq))
        return false;
    // get the domain id and the local id of the element in the domain
    int domain = get_domain_element(k);
    int kloc   = idxK_in_domain(k, domain);
    // find the element in the interfaces
    auto it    = interface_id_[itq].find(std::make_pair(domain, kloc));
    assert(it != interface_id_[itq].end());
    // get the interface id and the local element id in the interface
    int s = it->second.at(0).second; // 0 because no multi cut
    if (s == 0)
        return false; // means surface mesh
    // get the element id in the back mesh
    int kb = this->idxElementInBackMesh(k);
    // check if the face is cut
    return it->second.at(0).first->isCutFace(kb, ifac);
}

template <typename Mesh> bool ActiveMesh<Mesh>::isStabilizeElement(int k) const {
    for (int i = 0; i < nb_quadrature_time_; ++i) {
        // element is cut or not always active
        if (this->isCut(k, i) || this->isInactive(k, i)) {
            return true;
        }
    }
    return false;
}

template <typename Mesh> bool ActiveMesh<Mesh>::isInactive(int k, int t) const {
    // Get the domain index for the element
    int domain = get_domain_element(k);
    // Get the local element index within that domain
    int kloc   = idxK_in_domain(k, domain);
    // Search for the element in the active mesh at time t
    auto it    = in_active_mesh_[domain][t].find(kloc);
    // If the element was found, return true; otherwise, return false
    if (it == in_active_mesh_[domain][t].end())
        return false;
    return true;
}

template <typename Mesh> const Interface<Mesh> &ActiveMesh<Mesh>::get_interface(int k, int t) const {
    // Find the domain number of the element k
    int domain = get_domain_element(k);
    // Find the local index of the element k in the domain
    int kloc   = idxK_in_domain(k, domain);
    // Find the iterator to the interface_id_ map for the time instance t
    auto it    = interface_id_[t].find(std::make_pair(domain, kloc));
    // Make sure we found the interface_id_ map for the time instance t
    assert(it != interface_id_[t].end());
    // Return the first interface in the list associated with the iterator it
    return *(it->second.at(0).first);
}

/**
 * @brief Get the sub-partition of a given element in the active mesh?
 *
 * @tparam Mesh Mesh
 * @param k Element index in active mesh.
 * @param t The time step.
 * @return Partition<typename ActiveMesh<Mesh>::Element>
 */
template <typename Mesh>
Partition<typename ActiveMesh<Mesh>::Element> ActiveMesh<Mesh>::get_partition(int k, int t) const {
    // get the domain of the given element
    int domain = get_domain_element(k);
    // get the local index of the given element in the domain
    int kloc   = idxK_in_domain(k, domain);
    // find the interface id of the given element in the domain
    auto it    = interface_id_[t].find(std::make_pair(domain, kloc));
    // check that the interface id is found
    assert(it != interface_id_[t].end());
    // get the global index of the given element in the back mesh
    int kb = this->idxElementInBackMesh(k);
    // return the partition of the given element in the back mesh
    return it->second.at(0).first->get_partition(kb);
}

/**
 * @brief Returns the sign of the cut at the given element in the given time step.
 * The sign is +1 if the cut is in the positive side of the element,
 * and -1 if the cut is in the negative side of the element.
 * @tparam Mesh The mesh type.
 * @param k The element index.
 * @param t The time step.
 * @return The sign of the cut (-1 or +1).
 */
template <typename Mesh> int ActiveMesh<Mesh>::get_sign_cut(int k, int t) const {
    // get the domain of the given element
    int domain = get_domain_element(k);
    // get the local index of the given element in the domain
    int kloc   = idxK_in_domain(k, domain);
    // find the element in the given time step
    auto it    = interface_id_[t].find(std::make_pair(domain, kloc));
    // check that the element was found
    assert(it != interface_id_[t].end());
    // return the sign of the cut
    return it->second.at(0).second;
}

/**
 * @brief Get the cut part of the element k in the time t
 *
 * @tparam Mesh The mesh type
 * @param k The element index in the active mesh
 * @param t The time step
 * @return Cut_Part<typename ActiveMesh<Mesh>::Element>
 */
template <typename Mesh>
Cut_Part<typename ActiveMesh<Mesh>::Element> ActiveMesh<Mesh>::get_cut_part(int k, int t) const {
    // find the domain of the element k
    int domain = get_domain_element(k);
    // find the local index of the element k in the domain
    int kloc   = idxK_in_domain(k, domain);
    // find the cut of the element k in the time t
    auto it    = interface_id_[t].find(std::make_pair(domain, kloc));
    // if not cut build a partition that consider full Mesh::Element
    if (it == interface_id_[t].end()) {
        return Cut_Part<typename ActiveMesh<Mesh>::Element>(Partition<typename ActiveMesh<Mesh>::Element>((*this)[k]),
                                                            -1);
    }
    // find the index of the element k in the back mesh
    int kb = this->idxElementInBackMesh(k);
    // if the element k is cut in one interface
    if (it->second.size() == 1)
        // return the partition of the element k in the time t
        // and the local index of the interface
        return Cut_Part<typename ActiveMesh<Mesh>::Element>(it->second.at(0).first->get_partition(kb),
                                                            it->second.at(0).second);
    else
        // return the partition of the element k in the time t
        // and the local index of the interface
        return Cut_Part<typename ActiveMesh<Mesh>::Element>(this->build_local_partition(k), 0);
}

/**
 * @brief Get the face of the element k that is cut in time t
 *
 * @tparam Mesh Mesh
 * @param face Face
 * @param k Element index in active mesh
 * @param ifac Face index
 * @param t Time instance
 * @return Cut_Part<typename ActiveMesh<Mesh>::Element::Face>
 */
template <typename Mesh>
Cut_Part<typename ActiveMesh<Mesh>::Element::Face> ActiveMesh<Mesh>::get_cut_face(Face &face, int k, int ifac,
                                                                                  int t) const {

    // BUILD THE FACE
    // In the class mesh the inner faces are not built
    int kb = this->idxElementInBackMesh(k);
    int iv[Face::nv];
    for (int i = 0; i < Face::nv; ++i)
        iv[i] = Th(kb, Mesh::Element::nvhyperFace[ifac][i]);
    face.set(Th.vertices, iv, 0);

    // GET THE INTERFACE
    int domain = get_domain_element(k);
    int kloc   = idxK_in_domain(k, domain);
    auto it    = interface_id_[t].find(std::make_pair(domain, kloc));
    assert(it != interface_id_[t].end());

    if (it->second.size() == 1)
        return Cut_Part<Face>(it->second.at(0).first->get_partition_face(face, kb, ifac), it->second.at(0).second);
    else
        return Cut_Part<Face>(this->build_local_partition(face, k, ifac), 0);
}

/**
 *
 * @brief Returns the index of all elements in the back mesh corresponding to the given element in the given domain.//!
 * ????
 * @tparam Mesh Type of the mesh.
 * @param k Index of the element in the mesh.
 * @param d Domain of the element. If -1, all domains are considered.
 * @return Vector of indices of the corresponding elements in the back mesh.
 * If the domain is specified, the method returns the index of the element in the back mesh corresponding to the
 * specified domain. If the domain is not specified, the method returns the indices of the elements in the back mesh
 * corresponding to each domain.
 */
template <typename Mesh> std::vector<int> ActiveMesh<Mesh>::idxAllElementFromBackMesh(int k, int d) const {
    std::vector<int> idx(0);
    if (d != -1) {
        // If the domain is specified, get the index of the element
        //    in the back mesh corresponding to the specified domain.
        int ret = idxElementFromBackMesh(k, d);
        assert(ret != -1);
        idx.push_back(ret);
        return idx;
    }
    for (int i = 0; i < get_nb_domain(); ++i) {
        // Otherwise, get the index of the element in the back mesh
        //    corresponding to each domain.
        int ret = idxElementFromBackMesh(k, i);
        if (ret != -1)
            idx.push_back(ret);
    }

    // Assert that the number of domains is less than 3.
    assert(idx.size() > 0 && idx.size() < 3);
    return idx;
}

/**
 * @brief //! ????
 *
 * @tparam Mesh
 * @param k
 * @return std::vector<int>
 */
template <typename Mesh> std::vector<int> ActiveMesh<Mesh>::getAllDomainId(int k) const {
    // create a vector of size 0
    std::vector<int> idx(0);
    // loop over all domains
    for (int i = 0; i < get_nb_domain(); ++i) {
        // get the index of the element in the back mesh
        int ret = idxElementFromBackMesh(k, i);
        // if the element exists in the back mesh, add the domain id to the vector
        if (ret != -1)
            idx.push_back(i);
    }
    // check that the element is in at most 2 domains
    assert(idx.size() > 0 && idx.size() < 3);
    return idx;
}

/**
 * @brief //! ????
 *
 * @tparam Mesh
 * @param k
 * @return int
 */
template <typename Mesh> int ActiveMesh<Mesh>::idxElementFromBackMesh(int k) const {
    assert(0);
    return -1;
}

/**
 * @brief Get the index of the element k in the active mesh from element index in the background mesh //! ????
 *
 * @tparam Mesh
 * @param k index of element in the background mesh
 * @param i domain //! ????
 * @return int index of element in the active mesh
 */
template <typename Mesh> int ActiveMesh<Mesh>::idxElementFromBackMesh(int k, int i) const {
    if (i == -1)
        assert(0);
    if (get_nb_domain() == 1) {
        i = 0;
    } // fix bug but not satisfying
    auto it = idx_from_background_mesh_[i].find(k);
    if (it == idx_from_background_mesh_[i].end())
        return -1;
    return idxK_begin(i) + it->second;
}

/**
 * @brief Get index of element in the background mesh
 *
 * @tparam Mesh
 * @param k index of element in the active mesh
 * @return int index of element in the background mesh
 */
template <typename Mesh> int ActiveMesh<Mesh>::idxElementInBackMesh(const int k) const {
    int i = this->get_domain_element(k);
    int l = idxK_in_domain(k, i);
    return idx_in_background_mesh_[i][l];
}

/**
 * @brief Get index of element in the background mesh
 *
 * @tparam Mesh
 * @param k index of element in the active mesh
 * @param i domain index
 * @return int index of element in the background mesh
 */
template <typename Mesh> int ActiveMesh<Mesh>::idxElementInBackMesh(const int k, int i) const {
    int l = idxK_in_domain(k, i);
    return idx_in_background_mesh_[i][l];
}

/**
 * @brief Get the element adjacent to a given element in the active mesh across a given face
 *
 * @tparam Mesh
 * @param k element in the active mesh
 * @param j face index across which the adjacent element is searched
 * @return int index of adjacent element in the active mesh
 */
template <typename Mesh> int ActiveMesh<Mesh>::ElementAdj(const int k, int &j) const {
    int domain = get_domain_element(k);
    int kb     = this->idxElementInBackMesh(k);
    int kbn    = this->Th.ElementAdj(kb, j);
    if (kbn == -1)
        return -1;

    return this->idxElementFromBackMesh(kbn, domain);
}

/**
 * @brief Print information about the active mesh
 *
 * @tparam Mesh
 * @return prints number of subdomains, number of elements in each subdomain, and total number of elements
 */
template <typename Mesh> void ActiveMesh<Mesh>::info() const {
    std::cout << " ---  INFO CUT MESH  --- " << std::endl;
    std::cout << " Number of subdomains          :\t" << get_nb_domain() << std::endl;
    for (int i = 0; i < get_nb_domain(); ++i) {
        std::cout << " Number of elements in Omega_" << i << " :\t" << this->get_nb_element(i) << std::endl;
    }
    std::cout << " Total number of elements      :\t" << this->get_nb_element() << std::endl;
}

/**
 * @brief Constructs active mesh for problem with two subdomains {1, -1}
 *
 * @tparam Mesh
 * @param interface interface between the two subdomains
 */
template <typename Mesh> void ActiveMesh<Mesh>::init(const Interface<Mesh> &interface) {
    // Reserve memory for the indices in the background mesh for positive and negative domains
    idx_in_background_mesh_[0].reserve(Th.nt);
    idx_in_background_mesh_[1].reserve(Th.nt);

    // Push 0 as the first index for the element domain
    idx_element_domain.push_back(0);

    // Initialize the counters for positive and negative domains
    int nt0 = 0, nt1 = 0;

    // Loop through all the elements in the background mesh
    for (int k = 0; k < Th.nt; ++k) {

        // Get the sign of the current element
        const SignElement<typename ActiveMesh<Mesh>::Element> signK = interface.get_SignElement(k);

        // If the element is cut, add it to both domains
        if (signK.cut()) {
            // Push the index of the element into both positive and negative domain index arrays
            idx_in_background_mesh_[0].push_back(k);
            idx_from_background_mesh_[0][k] = nt0;
            idx_in_background_mesh_[1].push_back(k);
            idx_from_background_mesh_[1][k] = nt1;

            // Update the interface_id array for both positive and negative domains
            interface_id_[0][std::make_pair(0, nt0)].push_back(std::make_pair(&interface, 1));
            interface_id_[0][std::make_pair(1, nt1)].push_back(std::make_pair(&interface, -1));

            // Increment the counters for both positive and negative domain elements
            nt0++;
            nt1++;

        }
        // If the element is not cut, add it to the appropriate domain
        else {
            // Get the sign of the element
            int s    = signK.sign();
            // Get the index of the element in the appropriate domain
            int &nnt = (s > 0) ? nt0 : nt1;
            // Push the index of the corresponding background mesh element into the appropriate domain index array
            idx_in_background_mesh_[(s < 0)].push_back(k);
            // Push the index of the corresponding active mesh element into the appropriate domain index array
            idx_from_background_mesh_[(s < 0)][k] = nnt;

            nnt++;
        }
    }
    // Shrink the size of the index arrays to the actual size
    idx_in_background_mesh_[0].resize(nt0);
    idx_in_background_mesh_[1].resize(nt1);
    idx_in_background_mesh_[0].shrink_to_fit();
    idx_in_background_mesh_[1].shrink_to_fit();

    // Push the total number of elements in the positive and negative domains as the last index of the element domain
    // array
    idx_element_domain.push_back(nt0);
    idx_element_domain.push_back(nt0 + nt1);

    in_active_mesh_.resize(10);
    for (int i = 0; i < 10; ++i)
        in_active_mesh_[i].resize(nb_quadrature_time_);
}

template <typename Mesh> void ActiveMesh<Mesh>::addArtificialInterface(const Interface<Mesh> &interface) {

    int dom_size = this->get_nb_domain();
    assert(dom_size == 1);
    int sign_domain = 0;

    bool sub_is_cut = false;
    for (auto it_k = idx_from_background_mesh_[0].begin(); it_k != idx_from_background_mesh_[0].end(); ++it_k) {

        int kb = it_k->first;
        int k  = it_k->second;
        assert(kb == k);

        // auto it_gamma = interface_id_[0].find(std::make_pair(0, k));
        const SignElement<typename ActiveMesh<Mesh>::Element> signK = interface.get_SignElement(kb);

        if (signK.cut()) {
            interface_id_[0][std::make_pair(0, k)].push_back(std::make_pair(&interface, sign_domain));
        }
    }
}

template <typename Mesh>
Physical_Partition<typename ActiveMesh<Mesh>::Element> ActiveMesh<Mesh>::build_local_partition(const int k,
                                                                                               int t) const {

    typedef typename ActiveMesh<Mesh>::Element Element;

    int nvc = Element::Rd::d + 1;
    typedef SortArray<Ubyte, Element::Rd::d + 1> ElementIdx;
    // GET THE Mesh::Element TO BE CUT
    const Element &K((*this)[k]);
    Physical_Partition<Element> partition(K);
    std::vector<ElementIdx> elements_idx;

    Ubyte iv[Element::nb_ntcut][nvc];

    // INITIALIZE THE LOCAL_PARTITION WITH K
    for (int e = 0; e < Element::nb_ntcut; e++) {
        for (int i = 0; i < nvc; ++i) {
            iv[e][i] = (i + 2 * e) % Element::nv;
        }
    }
    for (int i = 0; i < Element::nv; ++i) {
        partition.add_node(K[i]);
    }
    for (int e = 0; e < Element::nb_ntcut; e++) {
        // partition.add_element(ElementIdx(iv[e]));
        elements_idx.push_back(ElementIdx(iv[e]));
    }

    if (!isCut(k, t))
        return partition;

    // START CUTTING PROCEDURE
    std::list<int> erased_element;
    std::vector<ElementIdx> new_element_idx;

    int domain = get_domain_element(k);
    int kloc   = idxK_in_domain(k, domain);
    auto it    = interface_id_[t].find(std::make_pair(domain, kloc));

    // std::cout << " Mesh::Element \t" << k << "in domain " << domain << std::endl;
    // std::cout << it->second.size() << " interfaces " << std::endl;
    for (int i = 0; i < it->second.size(); ++i) {

        // FRACTURE i DOES NOT CUT THIS Mesh::Element => NEXT
        //   if(!cut_lines_[i]->is_cut_element(k)) continue;
        //   std::cout << " is cut by fracture \t" << i << std::endl;
        //   // WE PERFORM THE CUT
        int s = it->second[i].second;
        const Interface<Mesh> &interface(*(it->second[i].first));
        interface.cut_partition(partition, new_element_idx, erased_element, s);

        // SORT THE LIST FROM HIGHER TO LOWER INDEX
        // TO NOT DESTROYED INDICES WHEN ERASING
        erased_element.sort(std::greater<int>());
        // std::cout << " nb K to erase " << erased_element.size() << std::endl;
        // ERASING THE CUT ELEMENTS
        for (auto it = erased_element.begin(); it != erased_element.end(); ++it) {
            // std::cout << " erased Mesh::Element " << *it << std::endl;
            // partition.erase_element(*it);
            elements_idx.erase(elements_idx.begin() + (*it));
        }
        // CREATING THE NEW ELEMENTS
        for (auto it = new_element_idx.begin(); it != new_element_idx.end(); ++it) {
            // partition.add_element(*it);
            elements_idx.push_back(ElementIdx(*it));
        }
        // std::cout << " iteration " << i << "nb K " << partition.nb_element()
        // << std::endl;
    }
    if (elements_idx.size() >= partition.max_size()) {
        std::cout << " Need to increase maximum size of array element_idx in "
                     "Physical partition"
                  << std::endl;
        assert(0);
    }
    for (int i = 0; i < elements_idx.size(); ++i) {
        partition.set_element(i, elements_idx[i]);
    }
    return partition;
}

template <typename Mesh>
Physical_Partition<typename ActiveMesh<Mesh>::Element::Face>
ActiveMesh<Mesh>::build_local_partition(Face &face, const int k, int ifac, int t) const {

    Physical_Partition<typename ActiveMesh<Mesh>::Element::Face> partition(face);

    return partition;
}

#endif
