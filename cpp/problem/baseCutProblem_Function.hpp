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

// INTEGRATION ON FULL ELEMENT

/**
 * @brief This function adds a bilinear form to the local matrix integrated over the cut mesh.
 * @tparam M The type of the matrix
 *
 * @tparam M Mesh
 * @param VF: Inner products as a list of vector of basis functions
 * @param Th: Active mesh
 * @note
 * The function allows for OpenMP parallelization, and loops over all elements in the active mesh.
 * For each element, if the element is a cut element, the function addElementContribution is called
 * from the BaseCutFEM class, otherwise the function addElementContribution from the BaseFEM class is called.
 * Finally, the function addLocalContribution is called to add the local contribution to the matrix.
 */

template <typename M> void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        bar += Th.next_element();

        if (Th.isInactive(k, 0))
            continue;

        if (Th.isCut(k, 0)) {
            addElementContribution(VF, k, nullptr, 0, 1.);
        } else {
            BaseFEM<M>::addElementContribution(VF, k, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M> 
const std::map<std::pair<int,int>, int>& BaseCutFEM<M>::get_dof_index_data(const FESpace &Vh, const CutMesh &Th) {
    dof_index_data.clear();

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (Th.isInactive(k, 0))
            continue;

        const int domain = Th.get_domain_element(k); // 0 or 1
        int cut_index    = -1;

        // Determine the current element's cut_index based on the 0/1/2/3 rule
        const bool current_is_cut_el = Th.isCut(k, 0);

        // actually i want to consider a dof if it belongs to a cut element or its neighbor element is cut
        bool neighbor_is_cut_el = false;
        for (int ifac = 0; ifac < Element::nea; ++ifac) { 
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            if (Th.isCut(kn, 0)) {
                neighbor_is_cut_el = true;
                break;
            }
        }

        const bool is_cut = current_is_cut_el || neighbor_is_cut_el;
        

        if (domain == 0) {
            cut_index = is_cut ? 1 : 0; // 0: Not Cut, Omega0; 1: Cut, Omega0
        } else if (domain == 1) {
            cut_index = is_cut ? 3 : 2; // 2: Not Cut, Omega1; 3: Cut, Omega1
        } else {
             // Handle unexpected domain index (optional, but safe)
             assert(false && "Unexpected Domain Index"); 
        }


        assert(cut_index >= 0);

        const bool current_is_cut = (cut_index == 1 || cut_index == 3);

        const FElement &FK(Vh[k]);

        // Loop over the DoF pairs
        for (int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
            for (int j = FK.dfcbegin(0); j < FK.dfcend(0); ++j) {

                auto key = std::make_pair(FK.loc2glb(i), FK.loc2glb(j));
                auto it  = dof_index_data.find(key);

                if (it == dof_index_data.end()) {
                    // 1. First time: always insert
                    dof_index_data.emplace(key, cut_index);
                } else {
                    // Already have a value: Check overwrite rules
                    // std::cout << "here1\n"; // Debug print

                    const int existing_cut_index = it->second;
                    const bool existing_is_cut   = (existing_cut_index == 1 || existing_cut_index == 3);

                    if (current_is_cut && !existing_is_cut) {
                        // Rule 1: Upgrade from interior/exterior (0 or 2) to cut (1 or 3)
                        it->second = cut_index;
                    } 
                    // Rule 2: If existing is cut (1 or 3), we do nothing. The cut status prevails.
                    // Rule 3: If both are non-cut (0 or 2), we keep the existing domain index.

                    // We do NOT need to check domains here because the cut status (1 or 3) already implies the domain (0 or 1).
                    // If we overwrite 0->1, the domain is preserved (Omega0).
                    // If we overwrite 2->3, the domain is preserved (Omega1).
                }
            }
        }
    }

    return dof_index_data;
}

template <typename M> 
std::map<int, std::pair<typename BaseCutFEM<M>::Rd, int>>& BaseCutFEM<M>::get_dof_vertex_data(const FESpace &Vh, const CutMesh &Th) {
    dof_vertex_data.clear();

    // dof_vertex_data[i] = (global_idx, (global_dof, cut_index))

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (Th.isInactive(k, 0))
            continue;

        const int domain = Th.get_domain_element(k); // 0 or 1
        int cut_index    = -1;

        // Determine the current element's cut_index based on the 0/1/2/3 rule
        const bool current_is_cut_el = Th.isCut(k, 0);

        // actually i want to consider a dof if it belongs to a cut element or its neighbor element is cut
        bool neighbor_is_cut_el = false;
        for (int ifac = 0; ifac < Element::nea; ++ifac) { 
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            if (Th.isCut(kn, 0)) {
                neighbor_is_cut_el = true;
                break;
            }
        }

        const bool is_cut = current_is_cut_el || neighbor_is_cut_el;
        

        if (domain == 0) {
            cut_index = is_cut ? 1 : 0; // 0: Not Cut, Omega0; 1: Cut, Omega0
        } else if (domain == 1) {
            cut_index = is_cut ? 3 : 2; // 2: Not Cut, Omega1; 3: Cut, Omega1
        } else {
             // Handle unexpected domain index (optional, but safe)
             assert(false && "Unexpected Domain Index"); 
        }

        assert(cut_index >= 0);
        
        const bool current_is_cut = (cut_index == 1 || cut_index == 3);

        // std::cout << "k = " << k << ", cut_index = " << cut_index << ", current_is_cut = " << current_is_cut << "\n";

        const FElement &FK(Vh[k]);

        // Loop over the DoFs (nodes) associated with vertices
        for (int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
            
            const Rd P = FK.Pt(i); // Get the physical coordinate of the DoF
            // auto key = std::make_pair(FK.loc2glb(i), std::make_pair(P, cut_index));
            
            auto it  = dof_vertex_data.find(FK.loc2glb(i));
            
            if (it == dof_vertex_data.end()) {
                // 1. First time: always insert
                dof_vertex_data.emplace(FK.loc2glb(i), std::make_pair(P, cut_index));
            } else {
                // Already have a value: Check overwrite rules
                // std::cout << "here\n"; // Debug print

                
                const int existing_cut_index = it->second.second;
                const bool existing_is_cut   = (existing_cut_index == 1 || existing_cut_index == 3);

                if (current_is_cut && !existing_is_cut) {
                    // Rule 1: Upgrade from interior/exterior (0 or 2) to cut (1 or 3)
                    it->second.second = cut_index;
                }
                // Rule 2: If existing is cut (1 or 3), we do nothing (preserve cut status).
                // Rule 3: If both are non-cut (0 or 2), we preserve the existing domain index (0 or 2).
            }
        }
    }

    return dof_vertex_data;
}



/**
 * @brief Bilinear form integrated over cut mesh and over time slab
 *
 * @tparam M Mesh
 * @param VF Bilinear form
 * @param Th Active mesh
 * @param In Time slab
 */
template <typename M> void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In) {

    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {

        addBilinear(VF, Th, In, itq);

        // MPIcf::Barrier();
    }
}

/**
 * @brief Bilinear form integrated over cut mesh in specific time quadrature point in the time slab.
 *
 * @tparam M Mesh
 * @param VF Bilinear form
 * @param Th Active mesh
 * @param In Time slab
 * @note The function is scaled with the time quadrature weight, and should therefore not be called directly.
 */
template <typename M>
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In, int itq) {
    // Check if the input VF is not a RHS (right-hand side)
    assert(!VF.isRHS());

    // Get the time quadrature
    auto tq = this->get_quadrature_time(itq);

    // Map the time quadrature to the time slab
    double tid = In.map(tq);

    // Start parallel region
#pragma omp parallel default(shared) num_threads(this->get_num_threads())
    {
        // Get the thread ID
        // Check if the number of threads matches the number of threads set in this object
        // assert(this->get_nb_thread() == omp_get_num_threads());
        int thread_id = omp_get_thread_num();

        // Calculate the offset for the current thread
        long offset = thread_id * this->offset_bf_time;

        // Create a matrix to store the time basis functions for this thread
        RNMK_ bf_time(this->databf_time_ + offset, In.NbDoF(), 1, op_dz);

        // Compute the time basis functions for this time quadrature point
        In.BF(tq.x, bf_time);

        // Calculate the time integration constant
        double cst_time = tq.a * In.get_measure();

        // Create a progress bar for this thread
        std::string title = " Add Bilinear CutMesh, In(" + std::to_string(itq) + ")";
        progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

        // Loop over all elements in the active mesh
#pragma omp for
        for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
            bar += Th.next_element();

            // Skip inactive elements for this time quadrature (if the element is inactive, it will be active for
            // another time quadrature)
            if (Th.isInactive(k, itq))
                continue;

            // Check if the element is cut
            if (Th.isCut(k, itq))
                // If the element is cut, add contribution to local matrix using addElementContribution function from
                // BaseCutFEM class
                addElementContribution(VF, k, &In, itq, cst_time);
            else
                // If the element is not cut, add contribution to local matrix using addElementContribution function
                // from BaseFEM class
                BaseFEM<M>::addElementContribution(VF, k, &In, itq, cst_time);

            // Add the contribution to the local matrix
            this->addLocalContribution();
        }
        bar.end();
    }
    // End parallel region
}

/**
 * @brief Bilinear form integrated over cut mesh in specific time quadrature point in the time slab.
 *
 * @tparam M Mesh
 * @param VF Bilinear form
 * @param Th Active mesh
 * @param In Time slab
 * @note The function is *NOT* scaled with the time quadrature weight, and may therefore be called directly.
 */
template <typename M>
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th, int itq, const TimeSlab &In) {
    // Assert that the input is not a RHS
    assert(!VF.isRHS());
    // Get the quadrature time for iteration "itq"
    auto tq = this->get_quadrature_time(itq);

    // Calculate the time using the map function of the TimeSlab "In"
    double tid = In.map(tq);

    // Allocate memory for the time-dependent basis functions
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);

    // Compute the time basic functions
    In.BF(tq.x, bf_time);

    // Set the title for the progress bar
    std::string title = " Add Bilinear Kh, In(" + std::to_string(itq) + ")";

    // Initialize the progress bar
    progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

    // Loop over each element of the active mesh
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        // Increment the progress bar
        bar += Th.next_element();

        // Skip the element if it is inactive at iteration "itq"
        if (Th.isInactive(k, itq))
            continue;

        // If the element is cut, add its contribution using BaseCutFEM
        if (Th.isCut(k, itq))
            addElementContribution(VF, k, &In, itq, 1.);
        // Else, add its contribution using BaseFEM
        else
            BaseFEM<M>::addElementContribution(VF, k, &In, itq, 1.);

        // Add the local contribution
        this->addLocalContribution();
    }

    // End the progress bar
    bar.end();
}

/**
 * @brief Bilinear form integrated over given time quadrature points in a subset of the time slab.
 *
 * @tparam M Mesh
 * @param VF Bilinear form
 * @param Th Active mesh
 * @param In Time slab, needed here because the basis functions in time live here
 * @param quad_pts quadrature points Jn subset In
 * @param phi vector of level-set functions evaluated in quad_pts 
 * @note the interfaces that are associated with the active mesh Th do not coincide with the interfaces
 * in phi, but the time quadrature points quad_pts must be a subset in the time interval In associated with Th
 */
// template <typename M>
// void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In, const std::vector<double> &quad_pts, const TimeInterface<Mesh> &interface) {
//     // Assert that the input is not a RHS
//     assert(!VF.isRHS());

//     assert(quad_pts.size() != 0);

//     // get end points of In on the reference interval [0,1]
//     double tq0 = In.map(this->get_quadrature_time(0));
//     double tqN = In.map(this->get_quadrature_time(this->get_nb_quad_point_time()-1));

//     assert(tq0 <= quad_pts[0]);
//     assert(quad_pts.back() <= tqN);

//     for (int itq = 0; itq < quad_pts.size(); ++itq) {

//         double tid = quad_pts[itq];

//         double t_ref = (tid-tq0) / (tqN-tq0);   // This maps tid from *In* to [0,1] for evaluation of the time basis functions
        
//         // Allocate memory for the time-dependent basis functions
//         RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);

//         // Compute the time basic functions
//         In.BF(t_ref, bf_time);

//         // Set the title for the progress bar
//         std::string title = " Add Bilinear Kh, In(" + std::to_string(itq) + ")";

//         // Initialize the progress bar
//         progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

//         // Loop over each element of the active mesh
//         for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
//             // Increment the progress bar
//             bar += Th.next_element();

//             // Skip the element if it is inactive at iteration "itq"
//             if (Th.isInactive(k, itq))
//                 continue;

//             // If the element is cut, add its contribution using BaseCutFEM
//             if (Th.isCut(k, itq))
//                 addElementContribution(VF, k, &In, itq, cst_time);
//             // Else, add its contribution using BaseFEM
//             else
//                 BaseFEM<M>::addElementContribution(VF, k, &In, itq, cst_time);

//             // Add the local contribution
//             this->addLocalContribution();
//         }

//         // End the progress bar
//         bar.end();
//     }
// }


/**
 * @brief Add bilinear in specific time
 *
 * @tparam M
 * @param VF
 * @param Th
 * @param t
 */
// template <typename M>
// void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const double t, const TimeSlab &In) {
//     // Assert that the input is not a RHS
//     assert(!VF.isRHS());

//     // Allocate memory for the time-dependent basis functions
//     RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);

//     // Compute the time basic functions
//     In.BF(tq.x, bf_time);

//     // Set the title for the progress bar
//     std::string title = " Add Bilinear Kh, In(" + std::to_string(itq) + ")";

//     // Initialize the progress bar
//     progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

//     // Loop over each element of the active mesh
//     for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
//         // Increment the progress bar
//         bar += Th.next_element();

//         // Skip the element if it is inactive at iteration "itq"
//         if (Th.isInactive(k, itq))
//             continue;

//         // If the element is cut, add its contribution using BaseCutFEM
//         if (Th.isCut(k, itq))
//             addElementContribution(VF, k, &In, itq, 1.);
//         // Else, add its contribution using BaseFEM
//         else
//             BaseFEM<M>::addElementContribution(VF, k, &In, itq, 1.);

//         // Add the local contribution
//         this->addLocalContribution();
//     }

//     // End the progress bar
//     bar.end();
// }

template <typename M> void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(VF.isRHS());
    progress bar(" Add Linear CutMesh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        if (Th.isInactive(k, 0))
            continue;

        if (Th.isCut(k, 0)) {
            addElementContribution(VF, k, nullptr, 0, 1.);
        } else {
            BaseFEM<M>::addElementContribution(VF, k, nullptr, 0, 1.);
        }
        // if(Th.isCut(k, 0))  BaseCutFEM<M>::addElementContribution(VF,
        // k,nullptr, 0, 1.); else BaseFEM<M>::addElementContribution(VF,
        // k,nullptr, 0, 1.);
    }
    bar.end();
}

template <typename M>
template <typename Fct>
void BaseCutFEM<M>::addLinear(const Fct &f, const itemVFlist_t &VF, const CutMesh &Th) {
    assert(VF.isRHS());

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (Th.isInactive(k, 0))
            continue;

        if (Th.isCut(k, 0)) {
            addElementContribution(f, VF, k, nullptr, 0, 1.);
        } else {
            BaseFEM<M>::addElementContribution(f, VF, k, nullptr, 0, 1.);
        }
        // if(Th.isCut(k, 0))  BaseCutFEM<M>::addElementContribution(VF,
        // k,nullptr, 0, 1.); else BaseFEM<M>::addElementContribution(VF,
        // k,nullptr, 0, 1.);
    }
}

template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &Th, int itq, const TimeSlab &In) {
    assert(VF.isRHS());
    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);
    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time);
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (Th.isInactive(k, itq))
            continue;

        if (Th.isCut(k, itq))
            addElementContribution(VF, k, &In, itq, 1.);
        else
            BaseFEM<M>::addElementContribution(VF, k, &In, itq, 1.);

        this->addLocalContribution();
    }
}


template <typename M>
template <typename Fct>
void BaseCutFEM<M>::addLinear(const Fct &f, const itemVFlist_t &VF, const CutMesh &Th, int itq, const TimeSlab &In) {
    //! DOES NOT WORK SAFELY
    // assert(0);
    assert(VF.isRHS());
    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);
    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time);
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (Th.isInactive(k, itq))
            continue;

        if (Th.isCut(k, itq))
            addElementContribution(f, VF, k, &In, itq, 1.);
        else
            BaseFEM<M>::addElementContribution(f, VF, k, &In, itq, 1.);

        this->addLocalContribution();
    }
}

template <typename M>
template <typename Fct>
void BaseCutFEM<M>::addLinear(const Fct &f, const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In) {
    //! DOES NOT WORK SAFELY, CHECK E.G. time_nitsche.cpp
    assert(0);
    assert(VF.isRHS());
    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        auto tq    = this->get_quadrature_time(itq);
        double tid = In.map(tq);
        KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
        RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
        In.BF(tq.x, bf_time);
        double cst_time = tq.a * In.get_measure();

        for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

            if (Th.isInactive(k, itq))
                continue;

            if (Th.isCut(k, itq))
                addElementContribution(f, VF, k, &In, itq, cst_time);
            else
                BaseFEM<M>::addElementContribution(f, VF, k, &In, itq, cst_time);

            this->addLocalContribution();
        }
    }
}

template <typename M> void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In) {
    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addLinear(VF, Th, In, itq);
    }
}

template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In, int itq) {
    assert(VF.isRHS());
    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time = tq.a * In.get_measure();

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (Th.isInactive(k, itq))
            continue;

        if (Th.isCut(k, itq))
            addElementContribution(VF, k, &In, itq, cst_time);
        else
            BaseFEM<M>::addElementContribution(VF, k, &In, itq, cst_time);
    }
}

template <typename M>
void BaseCutFEM<M>::addElementContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq,
                                           double cst_time) {

    // GET CUT AND COMPUTE PARAMETERS
    const FESpace &Vh(VF.get_spaceV(0));
    const CutMesh &Th(Vh.get_mesh());
    const Cut_Part<Element> cutK(Th.get_cut_part(k, itq));
    const FElement &FK(Vh[k]);
    const Element &K(FK.T);

    if (cutK.multi_interface()) {
        assert(0);
    } // not handled yet

    double meas = K.measure();
    double h    = K.get_h();
    int domain  = FK.get_domain();
    int kb      = Vh.idxElementInBackMesh(k);
    int iam     = omp_get_thread_num();
    // GET THE QUADRATURE RULE
    const QF &qf(this->get_quadrature_formular_cutK());
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    // LOOP OVER ELEMENTS IN THE CUT
    for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

        double meas_cut = cutK.measure(it);
        // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
        for (int l = 0; l < VF.size(); ++l) {
            if (!VF[l].on(domain))
                continue;

            // FINTE ELEMENT SPACES && ELEMENTS
            const FESpace &Vhv(VF.get_spaceV(l));
            const FESpace &Vhu(VF.get_spaceU(l));
            const FElement &FKv(Vhv[k]);
            const FElement &FKu(Vhu[k]);
            this->initIndex(FKu, FKv);

            // BF MEMORY MANAGEMENT -
            bool same   = (&Vhu == &Vhv);
            int lastop  = getLastop(VF[l].du, VF[l].dv);
            // RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop); //  the value for
            // basic fonction RNMK_ fu(this->databf_+ (same
            // ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the
            // value for basic fonction
            long offset = iam * this->offset_bf_;
            RNMK_ fv(this->databf_ + offset, FKv.NbDoF(), FKv.N,
                     lastop); //  the value for basic function
            RNMK_ fu(this->databf_ + offset + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N,
                     lastop); //  the value for basic function
            What_d Fop = Fwhatd(lastop);

            // COMPUTE COEFFICIENT
            double coef = VF[l].computeCoefElement(h, meas_cut, meas, meas_cut, domain);

            // LOOP OVER QUADRATURE IN SPACE
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

                typename QF::QuadraturePoint ip(qf[ipq]);
                Rd mip      = cutK.mapToPhysicalElement(it, ip); // to the physical cut part
                Rd cut_ip   = K.mapToReferenceElement(mip);      // back to the cut part in reference element
                double Cint = meas_cut * ip.getWeight() * cst_time;

                // EVALUATE THE BASIS FUNCTIONS
                FKv.BF(Fop, cut_ip, fv);
                if (!same)
                    FKu.BF(Fop, cut_ip, fu);
                //   VF[l].applyFunNL(fu,fv);

                // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
                Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);
                Cint *= coef * VF[l].c;

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
}

template <typename M>
template <typename Fct>
void BaseCutFEM<M>::addElementContribution(const Fct &f, const itemVFlist_t &VF, const int k, const TimeSlab *In,
                                           int itq, double cst_time) {

    // GET CUT AND COMPUTE PARAMETERS
    const FESpace &Vh(VF.get_spaceV(0));
    const CutMesh &Th(Vh.get_mesh());
    const Cut_Part<Element> cutK(Th.get_cut_part(k, itq));
    const FElement &FK(Vh[k]);
    const Element &K(FK.T);

    if (cutK.multi_interface()) {
        assert(0);
    } // not handled yet

    double meas = K.measure();
    double h    = K.get_h();
    int domain  = FK.get_domain();
    int kb      = Vh.idxElementInBackMesh(k);
    int iam     = omp_get_thread_num();
    // GET THE QUADRATURE RULE
    const QF &qf(this->get_quadrature_formular_cutK());
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    // LOOP OVER ELEMENTS IN THE CUT
    for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

        double meas_cut = cutK.measure(it);
        // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
        for (int l = 0; l < VF.size(); ++l) {
            if (!VF[l].on(domain))
                continue;

            // FINTE ELEMENT SPACES && ELEMENTS
            const FESpace &Vhv(VF.get_spaceV(l));
            const FESpace &Vhu(VF.get_spaceU(l));
            const FElement &FKv(Vhv[k]);
            const FElement &FKu(Vhu[k]);
            this->initIndex(FKu, FKv);

            // BF MEMORY MANAGEMENT -
            bool same   = (&Vhu == &Vhv);
            int lastop  = getLastop(VF[l].du, VF[l].dv);
            // RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop); //  the value for
            // basic fonction RNMK_ fu(this->databf_+ (same
            // ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the
            // value for basic fonction
            long offset = iam * this->offset_bf_;
            RNMK_ fv(this->databf_ + offset, FKv.NbDoF(), FKv.N,
                     lastop); //  the value for basic function
            RNMK_ fu(this->databf_ + offset + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N,
                     lastop); //  the value for basic function
            What_d Fop = Fwhatd(lastop);

            // COMPUTE COEFFICIENT
            double coef = VF[l].computeCoefElement(h, meas_cut, meas, meas_cut, domain);

            // LOOP OVER QUADRATURE IN SPACE
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

                typename QF::QuadraturePoint ip(qf[ipq]);
                Rd mip      = cutK.mapToPhysicalElement(it, ip); // to the physical cut part
                Rd cut_ip   = K.mapToReferenceElement(mip);      // back to the cut part in reference element
                double Cint = meas_cut * ip.getWeight() * cst_time;

                // EVALUATE THE BASIS FUNCTIONS
                FKv.BF(Fop, cut_ip, fv);
                if (!same)
                    FKu.BF(Fop, cut_ip, fu);
                //   VF[l].applyFunNL(fu,fv);

                // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
                Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);
                Cint *= coef * VF[l].c;
                Cint *= f(mip, VF[l].cv, domain);

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
}

template <typename M>
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const CExtension &ext, const int espE) {
    assert(!VF.isRHS());
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (Th.isInactive(k, 0))
            continue;

        if (Th.isCut(k, 0)) {
            addElementContribution(VF, k, nullptr, 0, 1.);
            addElementContributionOtherSide(VF, k, nullptr, 0, espE);
        } else
            BaseFEM<M>::addElementContribution(VF, k, nullptr, 0, 1.);

        this->addLocalContribution();
    }
}
template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &Th, const CExtension &ext, const int espE) {
    assert(VF.isRHS());
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        if (Th.isInactive(k, 0))
            continue;

        if (Th.isCut(k, 0)) {
            addElementContribution(VF, k, nullptr, 0, 1.);
            addElementContributionOtherSide(VF, k, nullptr, 0, espE);
        } else
            BaseFEM<M>::addElementContribution(VF, k, nullptr, 0, 1.);
    }
}
template <typename M>
void BaseCutFEM<M>::addElementContributionOtherSide(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq,
                                                    double cst_time) {

    // GET CUT AND COMPUTE PARAMETERS
    const FESpace &Vh(VF.get_spaceV(0));
    const CutMesh &Th(Vh.get_mesh());
    const Cut_Part<Element> cutK(Th.get_cut_part(k, itq));
    const FElement &FK(Vh[k]);
    const Element &K(FK.T);

    if (cutK.multi_interface()) {
        assert(0);
    } // not handled yet

    double meas = K.measure();
    double h    = K.get_h();
    int domain  = FK.get_domain();
    int kb      = Vh.idxElementInBackMesh(k);

    // GET THE QUADRATURE RULE
    const QF &qf(this->get_quadrature_formular_cutK());
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    // LOOP OVER ELEMENTS IN THE CUT
    for (auto it = cutK.other_side_element_begin(); it != cutK.other_side_element_end(); ++it) {

        double meas_cut = cutK.measure(it);

        // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
        for (int l = 0; l < VF.size(); ++l) {
            if (!VF[l].on(domain))
                continue;

            // FINTE ELEMENT SPACES && ELEMENTS
            const FESpace &Vhv(VF.get_spaceV(l));
            const FESpace &Vhu(VF.get_spaceU(l));
            const FElement &FKv(Vhv[k]);
            const FElement &FKu(Vhu[k]);
            this->initIndex(FKu, FKv);

            // BF MEMORY MANAGEMENT -
            bool same  = (&Vhu == &Vhv);
            int lastop = getLastop(VF[l].du, VF[l].dv);
            RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N,
                     lastop); //  the value for basic fonction
            RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N,
                     lastop); //  the value for basic fonction
            What_d Fop = Fwhatd(lastop);

            // COMPUTE COEFFICIENT
            double coef = VF[l].computeCoefElement(h, meas_cut, meas, meas_cut, domain);

            // LOOP OVER QUADRATURE IN SPACE
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

                typename QF::QuadraturePoint ip(qf[ipq]);
                Rd mip      = cutK.mapToPhysicalElement(it, ip); // to the physical cut part
                Rd cut_ip   = K.mapToReferenceElement(mip);      // back to the cut part in reference element
                double Cint = meas_cut * ip.getWeight() * cst_time;

                // EVALUATE THE BASIS FUNCTIONS
                FKv.BF(Fop, cut_ip, fv);
                if (!same)
                    FKu.BF(Fop, cut_ip, fu);
                //   VF[l].applyFunNL(fu,fv);

                // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
                Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);
                Cint *= coef * VF[l].c;
                // Cint = 1.;/coef;
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
}

// INTEGRATION ON INNER FACE
template <typename M> void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const CFacet &b) {
    assert(!VF.isRHS());
    progress bar(" Add Bilinear Face", Th.last_element(), globalVariable::verbose);

#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        // bar += Th.next_element();
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn == -1 || kn < k)
                continue;

            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);
            // CHECK IF IT IS A CUT EDGE
            if (Th.isCutFace(k, ifac, 0))
                addFaceContribution(VF, e1, e2, nullptr, 0, 1.);
            else
                BaseFEM<M>::addFaceContribution(VF, e1, e2, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const CFacet &b, const TimeSlab &In) {
    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addBilinear(VF, Th, b, In, itq);
    }
}

template <typename M>
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const CFacet &b, const TimeSlab &In,
                                int itq) {
    assert(!VF.isRHS());

    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time = tq.a * In.get_measure();

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        if (Th.isInactive(k, itq))
            continue;

        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn == -1 || kn < k)
                continue;

            if (Th.isInactive(kn, itq))
                continue;

            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);
            // CHECK IF IT IS A CUT EDGE
            if (Th.isCutFace(k, ifac, itq)) {
                // std::cout << k << "\t" << ifac << std::endl;
                addFaceContribution(VF, e1, e2, &In, itq, cst_time);

            } else
                BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
        }
        this->addLocalContribution();
    }
}

template <typename M> void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &Th, const CFacet &b) {
    assert(VF.isRHS());
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn == -1 || kn < k)
                continue;

            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);
            // CHECK IF IT IS A CUT EDGE
            if (Th.isCutFace(k, ifac, 0))
                addFaceContribution(VF, e1, e2, nullptr, 0, 1.);
            else
                BaseFEM<M>::addFaceContribution(VF, e1, e2, nullptr, 0, 1.);
        }
    }
}

template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &Th, const CFacet &b, const TimeSlab &In) {
    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addLinear(VF, Th, b, In, itq);
    }
}

template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &Th, const CFacet &b, const TimeSlab &In, int itq) {
    assert(VF.isRHS());

    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time = tq.a * In.get_measure();

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        if (Th.isInactive(k, itq))
            continue;

        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn == -1 || kn < k)
                continue;
            if (Th.isInactive(kn, itq))
                continue;

            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);
            // CHECK IF IT IS A CUT EDGE
            if (Th.isCutFace(k, ifac, itq))
                addFaceContribution(VF, e1, e2, &In, itq, cst_time);
            else
                BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
        }
    }
}

template <typename M>
void BaseCutFEM<M>::addFaceContribution(const itemVFlist_t &VF, const std::pair<int, int> &e1,
                                        const std::pair<int, int> &e2, const TimeSlab *In, int itq, double cst_time) {

    typedef typename FElement::RdHatBord RdHatBord;

    // CHECK IF IT IS FOR RHS OR MATRIX
    // CONVENTION ki < kj
    bool to_rhs = VF.isRHS();
    int ki = e1.first, ifac = e1.second;
    int kj = e2.first, jfac = e2.second;
    // Compute parameter coonected to the mesh.
    // on can take the one from the first test function
    const FESpace &Vh(*VF[0].fespaceV);
    const CutMesh &Th(Vh.get_mesh());
    const FElement &FKi(Vh[ki]);
    const FElement &FKj(Vh[kj]);
    const Element &Ki(FKi.T);
    const Element &Kj(FKj.T);

    const Cut_Part<Element> cutKi(Th.get_cut_part(ki, itq));
    const Cut_Part<Element> cutKj(Th.get_cut_part(kj, itq));
    typename Element::Face face;
    const Cut_Part<typename Element::Face> cutFace(Th.get_cut_face(face, ki, ifac, itq));

    double measK = cutKi.measure() + cutKj.measure();
    double measF = Ki.mesureBord(ifac);
    double h     = 0.5 * (Ki.get_h() + Kj.get_h());
    int domain   = FKi.get_domain();
    Rd normal    = Ki.N(ifac);

    int thread_id = omp_get_thread_num();
    // GET THE QUADRATURE RULE
    const QFB &qfb(this->get_quadrature_formular_cutFace());

    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    // LOOP OVER ELEMENTS IN THE CUT
    for (auto it = cutFace.element_begin(); it != cutFace.element_end(); ++it) {

        double meas = cutFace.measure(it);
        // std::cout << meas << std::endl;
        for (int l = 0; l < VF.size(); ++l) {
            if (!VF[l].on(domain))
                continue;

            // FINITE ELEMENT SPACES && ELEMENTS
            const FESpace &Vhv(VF.get_spaceV(l));
            const FESpace &Vhu(VF.get_spaceU(l));
            assert(Vhv.get_nb_element() == Vhu.get_nb_element());
            const int kv = VF[l].onWhatElementIsTestFunction(ki, kj);
            const int ku = VF[l].onWhatElementIsTrialFunction(ki, kj);

            int kbv = Vhv.idxElementInBackMesh(kv);
            int kbu = Vhu.idxElementInBackMesh(ku);

            const FElement &FKu(Vhu[ku]);
            const FElement &FKv(Vhv[kv]);
            this->initIndex(FKu, FKv);

            // BF MEMORY MANAGEMENT -
            bool same  = (VF.isRHS() || (&Vhu == &Vhv && ku == kv));
            int lastop = getLastop(VF[l].du, VF[l].dv);

            // Calculate the offset for the current thread

            long offset = thread_id * this->offset_bf_;

            RNMK_ fv(this->databf_ + offset, FKv.NbDoF(), FKv.N, lastop);
            RNMK_ fu(this->databf_ + offset + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
            What_d Fop = Fwhatd(lastop);

            // COMPUTE COEFFICIENT && NORMAL
            double coef = VF[l].computeCoefElement(h, meas, measK, measK, domain);
            coef *= VF[l].computeCoefFromNormal(normal);

            // LOOP OVER QUADRATURE IN SPACE
            for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
                typename QFB::QuadraturePoint ip(qfb[ipq]);
                const Rd mip = cutFace.mapToPhysicalElement(it, (RdHatBord)ip);
                double Cint  = meas * ip.getWeight() * cst_time;

                // EVALUATE THE BASIS FUNCTIONS
                FKv.BF(Fop, FKv.T.mapToReferenceElement(mip), fv);
                if (!same)
                    FKu.BF(Fop, FKu.T.mapToReferenceElement(mip), fu);
                //   VF[l].applyFunNL(fu,fv);

                Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kbu, kbv), std::make_pair(domain, domain),
                                                               mip, tid, normal);
                Cint *= coef * VF[l].c;
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
}

// INTEGRATION ON BOUNDARY
template <typename M>
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &cutTh, const CBorder &b, std::list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);

    for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element();
         idx_be += cutTh.next_boundary_element()) {

        int ifac;
        const int kb          = cutTh.Th.BoundaryElement(idx_be, ifac);
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);
        if (idxK.size() == 0)
            continue;

        const Element &K(cutTh.Th[kb]);
        const BorderElement &BE(cutTh.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {

            // for(int i=0;i<idxK.size();++i){
            //   std::cout << idxK[i] << "\t";
            // }
            // std::cout << std::endl;

            // CHECK IF IT IS A CUT EDGE
            if (cutTh.isCutFace(idxK[0], ifac, 0))
                addBorderContribution(VF, K, BE, ifac, nullptr, 0, 1.);
            else {
                if (idxK.size() == 1) {
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac, nullptr, 0, 1.);
                } else {
                    assert(cutTh.get_nb_domain() < 3);

                    // need to find out in what domain is the BOUNDARY
                    int k0          = idxK[0];
                    const auto cutK = cutTh.get_cut_part(k0, 0);
                    int ss          = cutK.get_sign();
                    int nv          = Element::nvface[ifac][0];
                    int id_sub      = (cutK.get_sign_node(nv) * ss == 1) ? 0 : 1;
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac + id_sub * Element::nea, nullptr, 0, 1.);
                }
            }
            this->addLocalContribution();
        }
    }
}

template <typename M>
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const CBorder &b, const TimeSlab &In,
                                std::list<int> label) {
    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addBilinear(VF, Th, b, In, itq, label);
    }
}

template <typename M>
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &cutTh, const CBorder &b, const TimeSlab &In,
                                int itq, std::list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);

    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time   = tq.a * In.get_measure();
    std::string title = " Add Bilinear Boundary, In(" + std::to_string(itq) + ")";
    progress bar(title.c_str(), cutTh.last_boundary_element(), globalVariable::verbose);

    for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element();
         idx_be += cutTh.next_boundary_element()) {
        bar += cutTh.next_boundary_element();
        int ifac;
        const int kb          = cutTh.Th.BoundaryElement(idx_be, ifac);
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);

        const Element &K(cutTh.Th[kb]);
        const BorderElement &BE(cutTh.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {

            // CHECK IF IT IS A CUT EDGE
            if (cutTh.isCutFace(idxK[0], ifac, itq)) {
                addBorderContribution(VF, K, BE, ifac, &In, itq, cst_time);
            } else {
                if (idxK.size() == 1) {
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac, &In, itq, cst_time);
                } else if (!cutTh.isCut(idxK[0], itq)) {
                    int id_sub = (cutTh.isInactive(idxK[0], itq)) ? 1 : 0;
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac + id_sub * Element::nea, &In, itq, cst_time);
                } else {
                    assert(cutTh.get_nb_domain() < 3);

                    // need to find out in what domain is the BOUNDARY
                    int k0          = idxK[0];
                    const auto cutK = cutTh.get_cut_part(k0, itq);
                    int ss          = cutK.get_sign();
                    int nv          = Element::nvface[ifac][0];
                    int id_sub      = (cutK.get_sign_node(nv) * ss == 1) ? 0 : 1;
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac + id_sub * Element::nea, &In, itq, cst_time);
                }
            }
            this->addLocalContribution();
        }
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &cutTh, const CBorder &b, std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);

    for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element();
         idx_be += cutTh.next_boundary_element()) {

        int ifac;
        const int kb          = cutTh.Th.BoundaryElement(idx_be, ifac);
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);
        if (idxK.size() == 0)
            continue;
        const Element &K(cutTh.Th[kb]);
        const BorderElement &BE(cutTh.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {

            // CHECK IF IT IS A CUT EDGE
            if (cutTh.isCutFace(idxK[0], ifac, 0)) {
                BaseCutFEM<M>::addBorderContribution(VF, K, BE, ifac, nullptr, 0, 1.);

            } else {
                if (idxK.size() == 1) {
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac, nullptr, 0, 1.);
                } else {
                    assert(cutTh.get_nb_domain() < 3);
                    // need to find out in what domain is the BOUNDARY
                    int k0          = idxK[0];
                    const auto cutK = cutTh.get_cut_part(k0, 0);
                    int ss          = cutK.get_sign();
                    int nv          = Element::nvface[ifac][0];
                    int id_sub      = (cutK.get_sign_node(nv) * ss == 1) ? 0 : 1;
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac + id_sub * Element::nea, nullptr, 0, 1.);
                }
            }
        }
    }
}

template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &Th, const CBorder &b, const TimeSlab &In,
                              std::list<int> label) {
    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addLinear(VF, Th, b, In, itq, label);
    }
}

template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &cutTh, const CBorder &b, const TimeSlab &In,
                              int itq, std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);

    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);

    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time = tq.a * In.get_measure();
    for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element();
         idx_be += cutTh.next_boundary_element()) {

        int ifac;
        const int kb          = cutTh.Th.BoundaryElement(idx_be, ifac);
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);

        const Element &K(cutTh.Th[kb]);
        const BorderElement &BE(cutTh.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {

            // CHECK IF IT IS A CUT EDGE
            if (cutTh.isCutFace(idxK[0], ifac, itq))
                addBorderContribution(VF, K, BE, ifac, &In, itq, cst_time);
            else {
                if (idxK.size() == 1) {
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac, &In, itq, cst_time);
                } else if (!cutTh.isCut(idxK[0], itq)) {
                    int id_sub = (cutTh.isInactive(idxK[0], itq)) ? 1 : 0;
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac + id_sub * Element::nea, &In, itq, cst_time);
                } else {
                    assert(cutTh.get_nb_domain() < 3);

                    // need to find out in what domain is the BOUNDARY
                    int k0          = idxK[0];
                    const auto cutK = cutTh.get_cut_part(k0, itq);
                    int ss          = cutK.get_sign();
                    int nv          = Element::nvface[ifac][0];
                    int id_sub      = (cutK.get_sign_node(nv) * ss == 1) ? 0 : 1;
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac + id_sub * Element::nea, &In, itq, cst_time);
                }
            }
        }
    }
}

template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &cutTh, const CBorder &b, int itq,
                              const TimeSlab &In, std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);

    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);

    // KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    // RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);

    // In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time = 1.; // do NOT scale with time

    for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element();
         idx_be += cutTh.next_boundary_element()) {

        int ifac;
        const int kb          = cutTh.Th.BoundaryElement(idx_be, ifac);
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);

        const Element &K(cutTh.Th[kb]);
        const BorderElement &BE(cutTh.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {

            // CHECK IF IT IS A CUT EDGE
            if (cutTh.isCutFace(idxK[0], ifac, itq))
                addBorderContribution(VF, K, BE, ifac, &In, itq, cst_time);
            else {
                if (idxK.size() == 1) {
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac, &In, itq, cst_time);
                } else if (!cutTh.isCut(idxK[0], itq)) {
                    int id_sub = (cutTh.isInactive(idxK[0], itq)) ? 1 : 0;
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac + id_sub * Element::nea, &In, itq, cst_time);
                } else {
                    assert(cutTh.get_nb_domain() < 3);

                    // need to find out in what domain is the BOUNDARY
                    int k0          = idxK[0];
                    const auto cutK = cutTh.get_cut_part(k0, itq);
                    int ss          = cutK.get_sign();
                    int nv          = Element::nvface[ifac][0];
                    int id_sub      = (cutK.get_sign_node(nv) * ss == 1) ? 0 : 1;
                    BaseFEM<M>::addBorderContribution(VF, K, BE, ifac + id_sub * Element::nea, &In, itq, cst_time);
                }
            }
        }
    }
}

template <typename M>
void BaseCutFEM<M>::addBorderContribution(const itemVFlist_t &VF, const Element &K, const BorderElement &BE, int ifac,
                                          const TimeSlab *In, int itq, double cst_time) {

    typedef typename FElement::RdHatBord RdHatBord;

    // Compute parameter connected to the mesh.
    double measK = K.measure();
    double h     = K.get_h();
    Rd normal    = K.N(ifac);

    // U and V HAS TO BE ON THE SAME MESH
    const FESpace &Vh(VF.get_spaceV(0));
    const CutMesh &Th(Vh.get_mesh());
    int kb                = Vh.Th(K);
    std::vector<int> idxK = Vh.idxAllElementFromBackMesh(kb, -1);
    // assert(idxK.size() == 2);

    // GET THE QUADRATURE RULE
    const QFB &qfb(this->get_quadrature_formular_cutFace());
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    int Ne = idxK.size();
    for (int e = 0; e < Ne; ++e) {

        int k = idxK[e];
        typename Element::Face face;
        const Cut_Part<typename Element::Face> cutFace(Th.get_cut_face(face, k, ifac, itq));
        const Cut_Part<Element> cutK(Th.get_cut_part(k, itq));

        double meas_cut = cutK.measure();

        // LOOP OVER ELEMENTS IN THE CUT
        for (auto it = cutFace.element_begin(); it != cutFace.element_end(); ++it) {

            double meas = cutFace.measure(it);

            for (int l = 0; l < VF.size(); ++l) {
                // if(!VF[l].on(domain)) continue;

                // FINITE ELEMENT SPACES && ELEMENTS
                const FESpace &Vhv(VF.get_spaceV(l));
                const FESpace &Vhu(VF.get_spaceU(l));
                assert(Vhv.get_nb_element() == Vhu.get_nb_element());
                bool same = (VF.isRHS() || (&Vhu == &Vhv));
                const FElement &FKu(Vhu[k]);
                const FElement &FKv(Vhv[k]);
                int domain = FKv.get_domain();
                this->initIndex(FKu, FKv);

                // BF MEMORY MANAGEMENT -
                int lastop = getLastop(VF[l].du, VF[l].dv);
                RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
                RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
                What_d Fop = Fwhatd(lastop);

                // COMPUTE COEFFICIENT && NORMAL
                double coef = VF[l].computeCoefElement(h, meas, measK, meas_cut, domain);
                coef *= VF[l].computeCoefFromNormal(normal);

                // LOOP OVER QUADRATURE IN SPACE
                for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
                    typename QFB::QuadraturePoint ip(qfb[ipq]);
                    const Rd mip    = cutFace.mapToPhysicalElement(it, (RdHatBord)ip);
                    const Rd cut_ip = K.mapToReferenceElement(mip);
                    double Cint     = meas * ip.getWeight() * cst_time;

                    // EVALUATE THE BASIS FUNCTIONS
                    FKv.BF(Fop, cut_ip, fv);
                    if (!same)
                        FKu.BF(Fop, cut_ip, fu);
                    //   VF[l].applyFunNL(fu,fv);

                    Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid, normal);
                    Cint *= coef * VF[l].c;

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
    }
}

template <typename Mesh>
void BaseCutFEM<Mesh>::setDirichlet(const FunFEM<Mesh> &gh, const CutMesh &cutTh, std::list<int> label) {

    bool all_label = (label.size() == 0);
    std::map<int, double> dof2set;
    const FESpace &Vh(gh.getSpace());

    // for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element();
    //      idx_be += cutTh.next_boundary_element()) {
    for (int idx_be = 0; idx_be < cutTh.last_boundary_element(); idx_be++) {

        int ifac;
        const int kb          = cutTh.Th.BoundaryElement(idx_be, ifac);
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);
        int k                 = idxK[0];

        assert(idxK.size() == 1);

        const Element &K(cutTh.Th[kb]);
        const BorderElement &BE(cutTh.be(idx_be));
        const FElement &FK(Vh[k]);

        if (util::contain(label, BE.lab) || all_label) {

            // for( int ic=0; ic<Vh.N;++ic) {
            for (int ic = 0; ic < 1; ++ic) {
                for (int df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df) {

                    int id_item = FK.DFOnWhat(df);

                    if (id_item < K.nv) {
                        // DOF is node (0,1,2)
                        assert(0);
                    }
                    // DOF is edge (3,4,5)
                    else if (id_item < K.nv + K.ne) {

                        // std::cout << " on edge  " <<FK.DFOnWhat(df) << std::endl;
                        int id_face = id_item - K.nv; // id of boundary face
                        if (id_face == ifac) {
                            // edge is boundary edge

                            int df_glob = FK.loc2glb(df);
                            // dof2set.insert({df_glob, gh(df_glob)});
                            dof2set.insert({df_glob, gh(df_glob)});

                            // std::cout << df_glob << std::endl;
                        }
                    } else {
                        // std::cout << " on face  " << FK.DFOnWhat(df) << std::endl;
                    }
                }
            }
        }
        // getchar();
    }

    assert(this->pmat_.size() == 1);
    eraseAndSetRow(this->get_nb_dof(), *(this->pmat_[0]), this->rhs_, dof2set);
}

/**
 * Strong boundary conditions for BDM1 in time
 */
template <typename Mesh>
void BaseCutFEM<Mesh>::setDirichlet(const FunFEM<Mesh> &gh, const CutMesh &cutTh, const TimeSlab &In,
                                    std::list<int> label) {

    bool all_label = (label.size() == 0);
    std::map<int, double> dof2set;
    const FESpace &Vh(gh.getSpace());

    // Loop through all fitted boundary elements
    for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element();
         idx_be += cutTh.next_boundary_element()) {

        int ifac;
        const int kb          = cutTh.Th.BoundaryElement(idx_be, ifac);
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);
        int k                 = idxK[0];

        assert(idxK.size() == 1);

        const Element &K(cutTh.Th[kb]);
        const BorderElement &BE(cutTh.be(idx_be));
        const FElement &FK(Vh[k]);

        if (util::contain(label, BE.lab) || all_label) {

            // for( int ic=0; ic<Vh.N;++ic) {

            // DOF for BDM is only one component (v*n)
            for (int ic = 0; ic < 1; ++ic) {

                // Loop over degrees of freedom of current element
                for (int df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df) {

                    int id_item = FK.DFOnWhat(df);

                    if (id_item < K.nv) {
                        assert(0);
                    } else if (id_item < K.nv + K.ne) {
                        // std::cout << " on edge  " <<FK.DFOnWhat(df) << std::endl;
                        int id_face = id_item - K.nv;
                        if (id_face == ifac) {

                            for (int dof_time; dof_time < In.NbDoF(); ++dof_time) {
                                int df_glob = FK.loc2glb(df, dof_time);
                                // dof2set.insert({df_glob, gh(df_glob)});
                                dof2set.insert({df_glob, gh(df_glob)});

                                // std::cout << df_glob << std::endl;
                            }
                        }
                    } else {
                        // std::cout << " on face  " << FK.DFOnWhat(df) << std::endl;
                    }
                }
            }
        }
        // getchar();
    }

    assert(this->pmat_.size() == 1);
    eraseAndSetRow(this->get_nb_dof(), *(this->pmat_[0]), this->rhs_, dof2set);
}

// set Dirichlet BC strongly for H(curl), ie edges (3D)
template <typename Mesh>
void BaseCutFEM<Mesh>::setDirichletHcurl(const FunFEM<Mesh> &gh, const CutMesh &cutTh, std::list<int> label) {
    std::cout << "WARNING (setDirichletHcurl): This sets H(curl) DOFs only (in 3D!) This space needs to be added first to the CutFEM object." << std::endl;

    bool all_label = (label.size() == 0);    // Check if we need to apply to all labels
    std::map<int, double> dof2set;    // Map to store degrees of freedom (DOFs) and their corresponding values to be set
    const FESpace &Vh(gh.getSpace());    // Get the finite element space from the provided FunFEM object
    int dof_init = this->mapIdx0_[&Vh]; // Get the FIRST index of the finite element space in the list of finite element spaces

    // Iterate over each boundary element in the cut mesh
    for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element(); idx_be += cutTh.next_boundary_element()) {
        int idx_bdry_face; // Index of the boundary face in the boundary element
        const int kb = cutTh.Th.BoundaryElement(idx_be, idx_bdry_face); // Get the boundary element index in the original mesh
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1); // Get all the indices of the active elements corresponding to the boundary element in the cut mesh

        // Ensure there's exactly one active element corresponding to the boundary element
        assert(idxK.size() == 1);
        int k = idxK[0];
        const FElement &FK(Vh[k]);  // Get the finite element for the current element in the finite element space

        // Get the current element and its boundary element
        const Element &K(cutTh.Th[kb]);
        const BorderElement &BE(cutTh.be(idx_be));

        // Check if the current boundary element label is in the list of labels to process, or if all labels are being processed
        if (util::contain(label, BE.lab) || all_label) {
            // Skip processing if the element is cut
            if (cutTh.isCut(k, 0)) {
                continue;
            }
            for (int df = FK.dfcbegin(0); df < FK.dfcend(0); ++df) { // Iterate over each DOF, there are 1 per edge, so FK.dfcend(0)=6
                // std::cout << "df = " << df << "\n";
                int id_item = FK.DFOnWhat(df); // Get the item (vertex/edge/face) that the DOF is associated with
                // std::cout << "id_item = " << id_item << "\n";
                // Make sure the DOF isn't associated with a vertex or a face
                if (id_item < K.nv)
                    continue;
                // assert(!(id_item < K.nv || id_item >= K.nv + K.ne));
                int id_edge = id_item - K.nv; // [4,9] -> [0,5]
                // assert(df == id_edge); // int id_edge = df;
                auto edges = K.edgeOfFace[idx_bdry_face]; // Check that edge is on boundary
                if (id_edge == edges[0] || id_edge == edges[1] || id_edge == edges[2]) {
                    // Get the global index of the DOF
                    int df_glob = dof_init + FK.loc2glb(df);
                    // Insert the DOF and its corresponding value into the map
                    dof2set.insert({df_glob, gh(df_glob)});
                }
            }
        }
    }

    // Ensure there's exactly one matrix in the problem
    // assert(this->pmat_.size() == 1);

    // Modify the matrix and the right-hand side to enforce the Dirichlet boundary conditions
    
    eraseAndSetRowCol(this->get_nb_dof(), this->mat_, this->rhs_, dof2set);
    
}


// set Dirichlet BC strongly for H¹, ie nodes (3D)
template <typename Mesh>
void BaseCutFEM<Mesh>::setDirichletHone(const FunFEM<Mesh> &gh, const CutMesh &cutTh, std::list<int> label) {
    std::cout << "WARNING (setDirichletH1): This sets H1 DOFs only (in 3D!)" << std::endl;

    bool all_label = (label.size() == 0);    // Check if we need to apply to all labels
    std::map<int, double> dof2set;    // Map to store degrees of freedom (DOFs) and their corresponding values to be set
    const FESpace &Vh(gh.getSpace());    // Get the finite element space from the provided FunFEM object

    int dof_init = this->mapIdx0_[&Vh]; // Get the FIRST index of the finite element space in the list of finite element spaces

    // Iterate over each boundary element in the cut mesh
    for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element(); idx_be += cutTh.next_boundary_element()) {
        int idx_bdry_face; // Index of the boundary face in the boundary element
        const int kb = cutTh.Th.BoundaryElement(idx_be, idx_bdry_face); // Get the boundary element index in the original mesh
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1); // Get all the indices of the active elements corresponding to the boundary element in the cut mesh

        // Ensure there's exactly one active element corresponding to the boundary element
        assert(idxK.size() == 1);
        int k = idxK[0];
        const FElement &FK(Vh[k]);  // Get the finite element for the current element in the finite element space

        // Get the current element and its boundary element
        const Element &K(cutTh.Th[kb]);
        const BorderElement &BE(cutTh.be(idx_be));

        // Check if the current boundary element label is in the list of labels to process, or if all labels are being processed
        if (util::contain(label, BE.lab) || all_label) {
            // Skip processing if the element is cut
            if (cutTh.isCut(k, 0)) {
                continue;
            }

            for (int df = FK.dfcbegin(0); df < FK.dfcend(0); ++df) { // Iterate over each degree of freedom (DOF)
                int id_item = FK.DFOnWhat(df); // Get the item (vertex/edge/face) that the DOF is associated with

                // Ensure the DOF is NOT associated with an edge or a face
                assert(!(id_item >= K.nv));
                // std::cout << K.edgeOfFace[idx_face] << id_edge << std::endl;
                auto nodes = K.nvface[idx_bdry_face]; // Check that edge is on boundary
                if (id_item == nodes[0] || id_item == nodes[1] || id_item == nodes[2]) {
                    // Get the global index of the DOF
                    // std::cout << (dof_init-1) + FK.loc2glb(df) << " : " << this->get_nb_dof() << std::endl;
                    int dof = FK.loc2glb(df);
                    int df_glob = dof_init + dof; // correct place in Matrix
                    // Insert the DOF and its corresponding value into the map
                    dof2set.insert({df_glob, gh(dof)});
                }
            }
        }
    }

    // Ensure there's exactly one matrix in the problem
    // assert(this->pmat_.size() == 1);

    // Modify the matrix and the right-hand side to enforce the Dirichlet boundary conditions
    // eraseAndSetRow(this->get_nb_dof(), *(this->pmat_[0]), this->rhs_, dof2set);

    eraseAndSetRowCol(this->get_nb_dof(), *(this->pmat_), this->rhs_, dof2set);
}

template <typename Mesh> void BaseCutFEM<Mesh>::removeDofForHansbo(const FESpace &Vh) {

    assert(Vh.basisFctType == BasisFctType::P0 || Vh.basisFctType == BasisFctType::P1dc);

    std::set<int> dof2rm;
    const CutMesh &Th(Vh.get_mesh());
    int idx0 = this->mapIdx0_[&Vh];
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (Th.isCut(k, 0)) {
        } else {
            if (Vh.basisFctType == BasisFctType::P0) {
                int df_glob0 = 2 * k + idx0;
                int df_glob1 = 2 * k + 1 + idx0;

                dof2rm.insert(df_glob0);
                dof2rm.insert(df_glob1);
            } else {
                for (int i = 0; i < 6; ++i) {
                    int df_glob = 6 * k + i + idx0;
                    dof2rm.insert(df_glob);
                }
            }
        }
    }
    int N = this->get_nb_dof();
    eraseRow(N, *this->pmat_, this->rhs_, dof2rm);
}

// On Ridges
template <typename M>
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const Interface<M> &gamma, const CRidge &innerRidge,
                                std::list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceRidgeContribution(VF, gamma, iface, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
}
template <typename M>
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const TimeInterface<M> &gamma, const CRidge &innerRidge,
                                const TimeSlab &In, std::list<int> label) {
    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addBilinear(VF, gamma, innerRidge, In, itq, label);
    }
}
template <typename M>
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const TimeInterface<M> &interface, const CRidge &innerRidge,
                                const TimeSlab &In, int itq, std::list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);
    auto tq        = this->get_quadrature_time(itq);
    double tid     = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time = tq.a * In.get_measure();

    const Interface<M> &gamma(*interface[itq]);

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceRidgeContribution(VF, gamma, iface, &In, itq, cst_time);
        }
        this->addLocalContribution();
    }
}

template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const Interface<M> &gamma, const CRidge &innerRidge,
                              std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceRidgeContribution(VF, gamma, iface, nullptr, 0, 1.);
        }
    }
}
template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const TimeInterface<M> &gamma, const CRidge &innerRidge,
                              const TimeSlab &In, std::list<int> label) {
    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addLinear(VF, gamma, innerRidge, In, itq, label);
    }
}
template <typename M>
void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const TimeInterface<M> &interface, const CRidge &innerRidge,
                              const TimeSlab &In, int itq, std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);
    auto tq        = this->get_quadrature_time(itq);
    double tid     = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time = tq.a * In.get_measure();
    const Interface<M> &gamma(*interface[itq]);
    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceRidgeContribution(VF, gamma, iface, &In, itq, cst_time);
        }
    }
}

// FACE STABILIZATION
template <typename M> void BaseCutFEM<M>::addFaceStabilization(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Face Stabilization CutMesh", Th.last_element(), globalVariable::verbose);

#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        if (!Th.isCut(k, 0) && !Th.isInactive(k, 0))
            continue;

        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if(kn == -1) continue;
            if (kn < k && Th.isCut(kn,0)) 
                continue;

            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);
            BaseFEM<M>::addFaceContribution(VF, e1, e2, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

// // For mixed finite element spaces on surface
// template <typename M> void BaseCutFEM<M>::addFaceStabilizationSurface(const itemVFlist_t &VF, const CutMesh &Th) {
//     assert(!VF.isRHS());
//     progress bar(" Add Face Stabilization CutMesh", Th.last_element(), globalVariable::verbose);

//     for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
//         bar += Th.next_element();

//         if (!Th.isCut(k, 0) && !Th.isInactive(k, 0))
//             continue;
//         for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

//             int jfac = ifac;
//             int kn   = Th.ElementAdj(k, jfac);      // -1 if outside active mesh
            
//             // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
//             if (kn < k)
//                 continue;

//             int kb  = Th.idxElementInBackMesh(k);
//             int kbn = Th.idxElementInBackMesh(kn);

//             std::pair<int, int> e1 = std::make_pair(kb, ifac);
//             std::pair<int, int> e2 = std::make_pair(kbn, jfac);
//             BaseFEM<M>::addFaceContributionSurface(VF, e1, e2, nullptr, 0, 1.);
//         }
//         this->addLocalContribution();
//     }
//     bar.end();
// }

template <typename M> void BaseCutFEM<M>::addFaceStabilizationMixed(const itemVFlist_t &VF, const CutMesh &Th) {
    //assert(0);  //! DO NOT USE, INNER EDGE SPECIAL CASE NOT FIXED
    assert(!VF.isRHS());
    assert(Th.get_nb_domain() == 1);
    progress bar(" Add Face Stabilization CutMesh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        if (!Th.isCut(k, 0) && !Th.isInactive(k, 0))
            continue;
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn == -1) continue;
            if ((kn < k) && Th.isCut(kn, 0))
                continue;

            int kb  = Th.idxElementInBackMesh(k);
            int kbn = Th.idxElementInBackMesh(kn);

            std::pair<int, int> e1 = std::make_pair(kb, ifac);
            std::pair<int, int> e2 = std::make_pair(kbn, jfac);
            BaseFEM<M>::addFaceContributionMixed(VF, e1, e2, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addFaceStabilization(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In) {

    number_of_stabilized_edges      = 0;
    int number_of_quadrature_points = this->get_nb_quad_point_time();
    for (int itq = 0; itq < number_of_quadrature_points; ++itq) {
        addFaceStabilization(VF, Th, In, itq);
    }
    number_of_stabilized_edges /= number_of_quadrature_points;
}

template <typename M>
void BaseCutFEM<M>::addFaceStabilization(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In, int itq) {
    assert(!VF.isRHS());
    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time   = tq.a * In.get_measure();
    std::string title = " Add Face Stab, In(" + std::to_string(itq) + ")";
    progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();
        if (!Th.isStabilizeElement(k))
            continue;
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn == -1)
                continue;
            if ((kn < k) && Th.isStabilizeElement(kn))
                continue;

            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);
            number_of_stabilized_edges += 1;
            BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M> void BaseCutFEM<M>::addPatchStabilization(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Patch Stabilization CutMesh", Th.last_element(), globalVariable::verbose);

    size_t num_stab_faces = 0;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        if (!Th.isCut(k, 0) && !Th.isInactive(k, 0))
            continue;
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            // if (kn < k)
            //     continue;
            if (kn == -1) continue;
            if (kn < k && Th.isCut(kn,0)) 
                continue;
        

            // std::pair<int, int> e1 = std::make_pair(k, ifac);
            // std::pair<int, int> e2 = std::make_pair(kn, jfac);
            BaseFEM<M>::addPatchContribution(VF, k, kn, nullptr, 0, 1.);
            num_stab_faces++;
        }
        this->addLocalContribution();
    }
    bar.end();
    //std::cout << "Number of stabilized faces: " << num_stab_faces << "\n";
}

template <typename M> void BaseCutFEM<M>::addPatchStabilizationMixed(const itemVFlist_t &VF, const CutMesh &Th) {
    //assert(0);  //! DO NOT USE, INNER EDGE SPECIAL CASE NOT FIXED
    assert(!VF.isRHS());
    progress bar(" Add Mixed Patch Stabilization CutMesh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        if (!Th.isCut(k, 0) && !Th.isInactive(k, 0))
            continue;
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn < k)
                continue;
            
            int kb  = Th.idxElementInBackMesh(k);
            int kbn = Th.idxElementInBackMesh(kn);

            std::pair<int, int> e1 = std::make_pair(kb, ifac);
            std::pair<int, int> e2 = std::make_pair(kbn, jfac);
            BaseFEM<M>::addPatchContributionMixed(VF, kb, kbn, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addPatchStabilization(const itemVFlist_t &VF, const CutMesh &Th, const MacroElement<M> &macro) {
    assert(!VF.isRHS());

    progress bar(" Add Macro Stabilization CutMesh", macro.macro_element.size(), globalVariable::verbose);

    for (auto me = macro.macro_element.begin(); me != macro.macro_element.end(); ++me) {
        bar += 1;
        for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {

            int k    = it->first;
            int ifac = it->second;
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            BaseFEM<M>::addPatchContribution(VF, k, kn, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addPatchStabilization(const itemVFlist_t &VF, const CutMesh &Th, const MacroElementSurface<M> &macro) {
    assert(!VF.isRHS());

    progress bar(" Add Macro Stabilization CutMesh", macro.macro_element.size(), globalVariable::verbose);

    for (auto me = macro.macro_element.begin(); me != macro.macro_element.end(); ++me) {
        bar += 1;
        for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {

            int k    = it->first;
            int ifac = it->second;
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            BaseFEM<M>::addPatchContribution(VF, k, kn, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addPatchStabilization(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In) {

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

        std::string title = " Add Patch Stab, In(" + std::to_string(itq) + ")";
        progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

        // Loop through active mesh elements
        for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
            bar += Th.next_element();

            // Exclude elements whose edges do not need stabilization
            if (!Th.isStabilizeElement(k))
                continue;

            // std::cout << "Stabilized element " << k << "\n";

            // Loop through the element's edges
            for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

                int jfac = ifac;
                int kn   = Th.ElementAdj(k, jfac); // get neighbor element's index

                // By skipping neighbors with smaller indices, we avoid adding contribution to the same edge twice
                // However, we don't want this condition to apply if kn is an interior element, since then we won't reach
                // this edge from the other side of the edge again, since the interior element won't be looped over
                //if ((kn < k) && (Th.isCut(kn,itq) || Th.isInactive(kn,itq)))

                if (kn == -1) continue;

                if ((kn < k) && (Th.isStabilizeElement(kn)))
                    continue;

                std::pair<int, int> e1 = std::make_pair(k, ifac);  // (element index, edge index) current element
                std::pair<int, int> e2 = std::make_pair(kn, jfac); // (element index, edge index) neighbor element
                num_stab_faces++;

                // Add patch contribution
                // BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
                BaseFEM<M>::addPatchContribution(VF, k, kn, &In, itq, cst_time);
            }
            this->addLocalContribution();
        }
        bar.end();
    }
    // std::cout << "Number of stabilized faces: " << num_stab_faces / number_of_quadrature_points << "\n";
}


template <typename M>
void BaseCutFEM<M>::addPatchStabilizationExterior(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In) {

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

        std::string title = " Add Patch Stab, In(" + std::to_string(itq) + ")";
        progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

        // Loop through active mesh elements
        for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
            bar += Th.next_element();

            // Exclude elements who are not exterior
            if (!Th.isExterior(k))
                continue;

            // Loop through the element's edges
            for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

                int jfac = ifac;
                int kn   = Th.ElementAdj(k, jfac); // get neighbor element's index

                // By skipping neighbors with smaller indices, we avoid adding contribution to the same edge twice
                // However, we don't want this condition to apply if kn is an interior element, since then we won't reach
                // this edge from the other side of the edge again, since the interior element won't be looped over
                //if ((kn < k) && (Th.isCut(kn,itq) || Th.isInactive(kn,itq)))

                if (kn == -1) continue;

                if ((kn < k) && (!Th.isExterior(kn)))
                    continue;

                std::pair<int, int> e1 = std::make_pair(k, ifac);  // (element index, edge index) current element
                std::pair<int, int> e2 = std::make_pair(kn, jfac); // (element index, edge index) neighbor element
                num_stab_faces++;

                // Add patch contribution
                // BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
                BaseFEM<M>::addPatchContribution(VF, k, kn, &In, itq, cst_time);
            }
            this->addLocalContribution();
        }
        bar.end();
    }
    // std::cout << "Number of stabilized faces: " << num_stab_faces / number_of_quadrature_points << "\n";
}




/**
 * @brief Patch stabilization in specific time quadrature point
 *
 * @tparam M Mesh
 * @param VF stabilization integrand
 * @param Th Active mesh
 * @param In Time slab
 * @param itq Time quadrature point
 */
template <typename M>
void BaseCutFEM<M>::addPatchStabilization(const itemVFlist_t &VF, const ActiveMesh<M> &Th, const TimeSlab &In,
                                          const int itq) {

    assert(!VF.isRHS());

    // Compute contribution from time basis functions
    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);
    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    // double cst_time = tq.a * In.get_measure();

    std::string title = " Add Patch Stab, In(" + std::to_string(itq) + ")";
    progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

    // Loop through active mesh elements
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        // Exclude elements whose edges do not need stabilization
        // if (!Th.isStabilizeElement(k))
        //     continue;

        if (!(Th.isCut(k, itq) || Th.isInactive(k, itq)))
            continue;

        // Loop through the element's edges
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac); // get neighbor element's index

            // By skipping neighbors with smaller indices, we avoid adding contribution to the same edge twice
            if (kn == -1) 
                continue;
            if ((kn < k) && (!(Th.isCut(kn, itq) || Th.isInactive(kn, itq))))
                continue;

            std::pair<int, int> e1 = std::make_pair(k, ifac);  // (element index, edge index) current element
            std::pair<int, int> e2 = std::make_pair(kn, jfac); // (element index, edge index) neighbor element

            // Add patch contribution
            // BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
            BaseFEM<M>::addPatchContribution(VF, k, kn, &In, itq, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M>
template <typename L>
void BaseCutFEM<M>::addPatchStabilization(const itemVFlist_t &VF, const ActiveMesh<M> &Th, const TimeSlab &In,
                                          const AlgoimMacro<M, L> &macro) {

    // number_of_stabilized_edges      = 0;
    int number_of_quadrature_points = this->get_nb_quad_point_time();
    for (int itq = 0; itq < number_of_quadrature_points; ++itq) {

        assert(!VF.isRHS());
        auto tq    = this->get_quadrature_time(itq);
        double tid = In.map(tq);

        KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
        RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
        In.BF(tq.x, bf_time); // compute time basic funtions
        double cst_time = tq.a * In.get_measure();

        for (auto me = macro.macro_element.begin(); me != macro.macro_element.end(); ++me) {

            for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {
                int k    = it->first;
                int ifac = it->second;
                int jfac = ifac;
                int kn   = Th.ElementAdj(k, jfac);

                // std::pair<int, int> e1 = std::make_pair(k, ifac);
                // std::pair<int, int> e2 = std::make_pair(kn, jfac);

                // number_of_stabilized_edges += 1;
                BaseFEM<M>::addPatchContribution(VF, k, kn, &In, itq, cst_time);
            }
            this->addLocalContribution();
        }
    }
    // number_of_stabilized_edges /= number_of_quadrature_points;
}


template <typename M>
void BaseCutFEM<M>::addPatchStabilization(const itemVFlist_t &VF, const ActiveMesh<M> &Th, const TimeSlab &In,
                                          const MacroElementPartition<M> &macro) {

    // number_of_stabilized_edges      = 0;
    int number_of_quadrature_points = this->get_nb_quad_point_time();
    for (int itq = 0; itq < number_of_quadrature_points; ++itq) {

        assert(!VF.isRHS());
        auto tq    = this->get_quadrature_time(itq);
        double tid = In.map(tq);

        KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
        RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
        In.BF(tq.x, bf_time); // compute time basic funtions
        double cst_time = tq.a * In.get_measure();

        for (auto me = macro.macro_element.begin(); me != macro.macro_element.end(); ++me) {

            for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {
                int k    = it->first;
                int ifac = it->second;
                int jfac = ifac;
                int kn   = Th.ElementAdj(k, jfac);

                // std::pair<int, int> e1 = std::make_pair(k, ifac);
                // std::pair<int, int> e2 = std::make_pair(kn, jfac);

                // number_of_stabilized_edges += 1;
                BaseFEM<M>::addPatchContribution(VF, k, kn, &In, itq, cst_time);
            }
            this->addLocalContribution();
        }
    }
    // number_of_stabilized_edges /= number_of_quadrature_points;
}

/**
 * @brief This only stabilize in the faces corresponding to the active mesh in quadrature point itq.
 *
 * @tparam M
 * @param VF
 * @param Th
 * @param itq
 * @param In
 */
template <typename M>
void BaseCutFEM<M>::addFaceStabilization(const itemVFlist_t &VF, const CutMesh &Th, int itq, const TimeSlab &In) {
    assert(!VF.isRHS());
    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    // double cst_time   = tq.a * In.get_measure();
    std::string title = " Add Face Stab, In(" + std::to_string(itq) + ")";
    progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();
        if (!(Th.isCut(k, itq)) && !Th.isInactive(k, itq))
            continue;
        // if (!Th.isStabilizeElement(k))
        //     continue;
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn == -1)
                continue;
            if ((kn < k) && (!(Th.isCut(kn, itq)) && !Th.isInactive(kn, itq)))
                continue;

            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);
            number_of_stabilized_edges += 1;
            BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addFaceStabilizationSpecial(const itemVFlist_t &VF, const CutMesh &Th, int itq,
                                                const TimeSlab &In) {

    assert(!VF.isRHS());
    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    // double cst_time   = tq.a * In.get_measure();
    std::string title = " Add Face Stab, In(" + std::to_string(itq) + ")";
    progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();
        if (!(Th.isCut(k, itq)) && !Th.isInactive(k, itq))
            continue;
        // if (!Th.isStabilizeElement(k))
        //     continue;
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn == -1)
                continue;
            if ((kn < k) && Th.isStabilizeElement(kn))
                continue;
            
            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);
            number_of_stabilized_edges += 1;
            BaseFEM<M>::addFaceContributionSpecial(VF, e1, e2, &In, itq, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addFaceStabilization(const itemVFlist_t &VF, const CutMesh &Th, const MacroElement<M> &macro) {

    progress bar(" Add Macro Stabilization CutMesh", macro.macro_element.size(), globalVariable::verbose);

    for (auto me = macro.macro_element.begin(); me != macro.macro_element.end(); ++me) {
        bar += 1;
        for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {

            int k    = it->first;
            int ifac = it->second;
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);

            BaseFEM<M>::addFaceContribution(VF, e1, e2, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addFaceStabilization(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In,
                                         const TimeMacroElement<M> &macro) {

    number_of_stabilized_edges      = 0;
    int number_of_quadrature_points = this->get_nb_quad_point_time();
    for (int itq = 0; itq < number_of_quadrature_points; ++itq) {

        assert(!VF.isRHS());
        auto tq    = this->get_quadrature_time(itq);
        double tid = In.map(tq);

        KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
        RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
        In.BF(tq.x, bf_time); // compute time basic funtions
        double cst_time = tq.a * In.get_measure();

        for (auto me = macro.macro_element.begin(); me != macro.macro_element.end(); ++me) {

            for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {
                int k    = it->first;
                int ifac = it->second;
                int jfac = ifac;
                int kn   = Th.ElementAdj(k, jfac);

                std::pair<int, int> e1 = std::make_pair(k, ifac);
                std::pair<int, int> e2 = std::make_pair(kn, jfac);

                number_of_stabilized_edges += 1;
                BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
            }
            this->addLocalContribution();
        }
    }
    number_of_stabilized_edges /= number_of_quadrature_points;
}

template <typename M>
void BaseCutFEM<M>::addFaceStabilization(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In,
                                         const MacroElementPartition<M> &macro) {

    number_of_stabilized_edges      = 0;
    int number_of_quadrature_points = this->get_nb_quad_point_time();
    for (int itq = 0; itq < number_of_quadrature_points; ++itq) {

        assert(!VF.isRHS());
        auto tq    = this->get_quadrature_time(itq);
        double tid = In.map(tq);

        KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
        RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
        In.BF(tq.x, bf_time); // compute time basic funtions
        double cst_time = tq.a * In.get_measure();

        for (auto me = macro.macro_element.begin(); me != macro.macro_element.end(); ++me) {

            for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {
                int k    = it->first;
                int ifac = it->second;
                int jfac = ifac;
                int kn   = Th.ElementAdj(k, jfac);

                std::pair<int, int> e1 = std::make_pair(k, ifac);
                std::pair<int, int> e2 = std::make_pair(kn, jfac);

                number_of_stabilized_edges += 1;
                BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
            }
            this->addLocalContribution();
        }
    }
    number_of_stabilized_edges /= number_of_quadrature_points;
}

template <typename M>
template <typename L>
void BaseCutFEM<M>::addFaceStabilization(const itemVFlist_t &VF, const ActiveMesh<M> &Th, const TimeSlab &In,
                                         const AlgoimMacro<M, L> &macro) {

    // number_of_stabilized_edges      = 0;
    int number_of_quadrature_points = this->get_nb_quad_point_time();
    for (int itq = 0; itq < number_of_quadrature_points; ++itq) {

        assert(!VF.isRHS());
        auto tq    = this->get_quadrature_time(itq);
        double tid = In.map(tq);

        KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
        RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
        In.BF(tq.x, bf_time); // compute time basic funtions
        double cst_time = tq.a * In.get_measure();

        for (auto me = macro.macro_element.begin(); me != macro.macro_element.end(); ++me) {

            for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {
                int k    = it->first;
                int ifac = it->second;
                int jfac = ifac;
                int kn   = Th.ElementAdj(k, jfac);

                std::pair<int, int> e1 = std::make_pair(k, ifac);
                std::pair<int, int> e2 = std::make_pair(kn, jfac);

                // number_of_stabilized_edges += 1;
                BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
            }
            this->addLocalContribution();
        }
    }
    // number_of_stabilized_edges /= number_of_quadrature_points;
}

template <typename M>
void BaseCutFEM<M>::addFaceStabilization(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In,
                                         const TimeMacroElementSurface<M> &macro) {
    number_of_stabilized_edges      = 0;
    int number_of_quadrature_points = this->get_nb_quad_point_time();
    for (int itq = 0; itq < number_of_quadrature_points; ++itq) {

        assert(!VF.isRHS());
        auto tq    = this->get_quadrature_time(itq);
        double tid = In.map(tq);

        KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
        RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
        In.BF(tq.x, bf_time); // compute time basic funtions
        double cst_time = tq.a * In.get_measure();

        for (auto me = macro.macro_element.begin(); me != macro.macro_element.end(); ++me) {

            for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {
                int k    = it->first;
                int ifac = it->second;
                int jfac = ifac;
                int kn   = Th.ElementAdj(k, jfac);

                std::pair<int, int> e1 = std::make_pair(k, ifac);
                std::pair<int, int> e2 = std::make_pair(kn, jfac);
                number_of_stabilized_edges += 1;
                BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
            }

            this->addLocalContribution();
        }
    }
    number_of_stabilized_edges /= number_of_quadrature_points;
}

template <typename M>
void BaseCutFEM<M>::addFaceStabilizationRHS(const itemVFlist_t &VF, const CutMesh &Th, const MacroElement<M> &macro) {

    progress bar(" Add Maro Stabilization RHS CutMesh", macro.macro_element.size(), globalVariable::verbose);

    assert(VF.isRHS());
    for (auto me = macro.macro_element.begin(); me != macro.macro_element.end(); ++me) {
        bar += 1;
        for (auto it = me->second.inner_edge.begin(); it != me->second.inner_edge.end(); ++it) {

            int k    = it->first;
            int ifac = it->second;
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);

            BaseFEM<M>::addFaceContribution(VF, e1, e2, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

// LAGRANGE MULTIPLIER
template <typename M> void BaseCutFEM<M>::addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &Th) {
    assert(VF.isRHS());
    int ndf = this->rhs_.size();
    this->rhs_.resize(ndf + 1);
    this->rhs_.at(ndf) = val;
    this->nb_dof_ = ndf + 1;  // Update the number of DOFs
    progress bar(" Add Lagrange Multiplier Kh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        if (Th.isInactive(k, 0))
            continue;

        if (Th.isCut(k, 0))
            addLagrangeContribution(VF, k, nullptr, 0, 1);
        else
            BaseFEM<M>::addLagrangeContribution(VF, k, nullptr, 0, 1);

        this->addLocalContributionLagrange(ndf);
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &Th, const int k) {
    assert(VF.isRHS());
    int ndf = this->rhs_.size();
    this->rhs_.resize(ndf + 1);
    this->rhs_[ndf] = val;
    this->nb_dof_ = ndf + 1;  // Update the number of DOFs

    if (Th.isCut(k, 0))
        BaseCutFEM<M>::addLagrangeContribution(VF, k, nullptr, 0, 1);
    else
        BaseFEM<M>::addLagrangeContribution(VF, k, nullptr, 0, 1);

    this->addLocalContributionLagrange(ndf);
}

template <typename M>
void BaseCutFEM<M>::addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &Th, const TimeSlab &In) {

    int ndf = this->rhs_.size();
    this->rhs_.resize(ndf + 1);
    this->rhs_[ndf] = val;
    this->nb_dof_ = ndf + 1;  // Update the number of DOFs
    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {

        addLagrangeMultiplier(VF, val, Th, In, itq, false);
    }
}

template <typename M>
void BaseCutFEM<M>::addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &Th, const TimeSlab &In,
                                          int itq, bool init) {
    assert(VF.isRHS());
    int ndf = this->rhs_.size() - 1;
    if (init) {
        ndf++;
        this->rhs_.resize(ndf + 1);
        this->rhs_[ndf] = val;
        this->nb_dof_ = ndf + 1;  // Update the number of DOFs
    }

    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time   = tq.a * In.get_measure();
    std::string title = " Add Lagrange Multiplier Kh, In(" + std::to_string(itq) + ")";
    progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();
        if (Th.isInactive(k, itq))
            continue;

        if (Th.isCut(k, itq))
            addLagrangeContribution(VF, k, &In, itq, cst_time);
        else
            BaseFEM<M>::addLagrangeContribution(VF, k, &In, itq, cst_time);

        this->addLocalContributionLagrange(ndf);
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &Th, int itq,
                                          const TimeSlab &In, bool init) {
    assert(VF.isRHS());
    int ndf = this->rhs_.size() - 1;
    if (init) {
        ndf++;
        this->rhs_.resize(ndf + 1);
        this->rhs_(ndf) = val;
        this->nb_dof_ = ndf + 1;  // Update the number of DOFs
    }

    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time   = tq.a * In.get_measure();
    std::string title = " Add Lagrange Multiplier Kh, In(" + std::to_string(itq) + ")";
    progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();
        if (Th.isInactive(k, itq))
            continue;

        if (Th.isCut(k, itq))
            addLagrangeContribution(VF, k, &In, itq, cst_time);
        else
            BaseFEM<M>::addLagrangeContribution(VF, k, &In, itq, 1.);

        this->addLocalContributionLagrange(ndf);
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addLagrangeContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq,
                                            double cst_time) {

    // GET CUT AND COMPUTE PARAMETERS
    const FESpace &Vh(VF.get_spaceV(0));
    const CutMesh &Th(Vh.get_mesh());
    const Cut_Part<Element> cutK(Th.get_cut_part(k, itq));
    const FElement &FK(Vh[k]);
    const Element &K(FK.T);

    if (cutK.multi_interface()) {
        assert(0);
    } // not handled yet

    double meas = K.measure();
    double h    = K.get_h();
    int domain  = FK.get_domain();
    int kb      = Vh.idxElementInBackMesh(k);

    // GET THE QUADRATURE RULE
    const QF &qf(this->get_quadrature_formular_cutK());
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    // LOOP OVER ELEMENTS IN THE CUT
    for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

        double meas_cut = cutK.measure(it);

        // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
        for (int l = 0; l < VF.size(); ++l) {
            if (!VF[l].on(domain))
                continue;
            // FINTE ELEMENT SPACES && ELEMENTS
            const FESpace &Vhv(VF.get_spaceV(l));
            const FElement &FKv(Vhv[k]);
            this->initIndex(FKv, FKv);

            // BF MEMORY MANAGEMENT -
            int lastop = getLastop(VF[l].du, VF[l].dv);
            RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N,
                     lastop); //  the value for basic fonction
            What_d Fop = Fwhatd(lastop);

            // COMPUTE COEFFICIENT
            double coef = VF[l].computeCoefElement(h, meas_cut, meas, meas_cut, domain);

            // LOOP OVER QUADRATURE IN SPACE
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

                typename QF::QuadraturePoint ip(qf[ipq]);
                Rd mip      = cutK.mapToPhysicalElement(it, ip); // to the physical cut part
                Rd cut_ip   = K.mapToReferenceElement(mip);      // back to the cut part in reference element
                double Cint = meas_cut * ip.getWeight() * cst_time;

                // EVALUATE THE BASIS FUNCTIONS
                FKv.BF(Fop, cut_ip, fv);
                //   VF[l].applyFunNL(fu,fv);

                // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
                // Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip);
                Cint *= coef * VF[l].c;
                // Cint = coef;
                if (In) {

                    this->addToMatrix(VF[l], *In, FKv, fv, Cint);
                } else {
                    this->addToMatrix(VF[l], FKv, fv, Cint);
                }
            }
        }
    }
}

template <typename M>
void BaseCutFEM<M>::addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &cutTh, const CBorder &b,
                                          std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);

    int ndf = this->rhs_.size();
    this->rhs_.resize(ndf + 1);
    this->rhs_[ndf] = val;
    this->nb_dof_ = ndf + 1;  // Update the number of DOFs

    for (int idx_be = cutTh.first_boundary_element(); idx_be < cutTh.last_boundary_element();
         idx_be += cutTh.next_boundary_element()) {

        int ifac;
        const int kb          = cutTh.Th.BoundaryElement(idx_be, ifac);
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);

        const Element &K(cutTh.Th[kb]);
        const BorderElement &BE(cutTh.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {

            // CHECK IF IT IS A CUT EDGE
            if (cutTh.isCutFace(idxK[0], ifac, 0))
                addLagrangeBorderContribution(VF, K, BE, ifac, nullptr, 0, 1.);
            else
                BaseFEM<M>::addLagrangeBorderContribution(VF, K, BE, ifac, nullptr, 0, 1.);
            this->addLocalContributionLagrange(ndf);
        }
    }
}

template <typename M>
void BaseCutFEM<M>::addLagrangeBorderContribution(const itemVFlist_t &VF, const Element &K, const BorderElement &BE,
                                                  int ifac, const TimeSlab *In, int itq, double cst_time) {

    // typedef typename FElement::RdHatBord RdHatBord;
    //
    // // Compute parameter connected to the mesh.
    // double measK = K.measure();
    // double h     = K.get_h();
    // Rd normal    = K.N(ifac);
    //
    // // U and V HAS TO BE ON THE SAME MESH
    // const FESpace& Vh(VF.get_spaceV(0));
    // const CutMesh& Th(Vh.get_mesh());
    // int kb = Vh.Th(K);
    // std::vector<int> idxK = Vh.idxAllElementFromBackMesh(kb,-1);
    // assert(idxK.size() == 2);
    //
    // // GET THE QUADRATURE RULE
    // const QFB& qfb(this->get_quadrature_formular_cutFace());
    // auto tq = this->get_quadrature_time(itq);
    // double tid = (In)? (double)In->map(tq) : 0.;
    //
    // for(int e=0;e<2;++e) {
    //
    //   int k = idxK[e];
    //   typename Element::Face face;
    //   const Cut_Part<typename Element::Face> cutFace(Th.get_cut_face(face, k,
    //   ifac)); const Cut_Part<Element> cutK(Th.get_cut_part(k, itq));
    //
    //   double meas_cut  = cutK.measure();
    //
    //   // LOOP OVER ELEMENTS IN THE CUT
    //   for(auto it = cutFace.element_begin();it != cutFace.element_end();
    //   ++it){
    //
    //     double meas  = cutFace.measure(it);
    //
    //     for(int l=0; l<VF.size();++l) {
    //       // if(!VF[l].on(domain)) continue;
    //
    //       // FINITE ELEMENT SPACES && ELEMENTS
    //       const FESpace& Vhv(VF.get_spaceV(l));
    //       const FESpace& Vhu(VF.get_spaceU(l));
    //       assert(Vhv.get_nb_element() == Vhu.get_nb_element());
    //       bool same = (VF.isRHS() || (&Vhu == &Vhv));
    //       const FElement& FKu(Vhu[k]);
    //       const FElement& FKv(Vhv[k]);
    //       int domain = FKv.get_domain();
    //       this->initIndex(FKu, FKv);
    //
    //
    //       // BF MEMORY MANAGEMENT -
    //       int lastop = getLastop(VF[l].du, VF[l].dv);
    //       RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop);
    //       RNMK_ fu(this->databf_+ (same ?0:FKv.NbDoF()*FKv.N*lastop)
    //       ,FKu.NbDoF(),FKu.N,lastop); What_d Fop = Fwhatd(lastop);
    //
    //       // COMPUTE COEFFICIENT && NORMAL
    //       double coef = VF[l].computeCoefElement(h,meas,measK,meas_cut,domain)
    //       ; coef *= VF[l].computeCoefFromNormal(normal);
    //
    //
    //       // LOOP OVER QUADRATURE IN SPACE
    //       for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
    //         typename QFB::QuadraturePoint ip(qfb[ipq]);
    //         const Rd mip = cutFace.mapToPhysicalElement(it, (RdHatBord)ip);
    //         const Rd cut_ip = K.mapToReferenceElement(mip);
    //         double Cint = meas * ip.getWeight() * cst_time;
    //
    //         // EVALUATE THE BASIS FUNCTIONS
    //         FKv.BF(Fop,cut_ip, fv);
    //         if(!same) FKu.BF(Fop, cut_ip, fu);
    //         //   VF[l].applyFunNL(fu,fv);
    //
    //         Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip,
    //         tid, normal); Cint *= coef * VF[l].c;
    //
    //         if( In ){
    //           if(VF.isRHS()) this->addToRHS(   VF[l], *In, FKv, fv, Cint);
    //           else           this->addToMatrix(VF[l], *In, FKu, FKv, fu, fv,
    //           Cint);
    //         }
    //         else {
    //           if(VF.isRHS()) this->addToRHS(   VF[l], FKv, fv, Cint);
    //           else           this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
    //         }
    //       }
    //     }
    //   }
    // }
}

template <typename M>
void BaseCutFEM<M>::addLagrangeMultiplier(const itemVFlist_t &VF, double val, const CutMesh &Th, const CExtension &ext,
                                          const int epsE) {
    assert(VF.isRHS());
    int ndf = this->rhs_.size();
    this->rhs_.resize(ndf + 1);
    this->rhs_[ndf] = val;
    this->nb_dof_ = ndf + 1;  // Update the number of DOFs

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (Th.isCut(k, 0)) {
            addLagrangeContribution(VF, k);
            addLagrangeContributionOtherSide(VF, k, epsE);
        } else
            BaseFEM<M>::addLagrangeContribution(VF, k);

        this->addLocalContributionLagrange(ndf);
    }
    // (*this)(ndf,ndf) = 1;
}
template <typename M>
void BaseCutFEM<M>::addLagrangeContributionOtherSide(const itemVFlist_t &VF, const int k, const int epsE) {

    // GET CUT AND COMPUTE PARAMETERS
    const FESpace &Vh(VF.get_spaceV(0));
    const CutMesh &Th(Vh.get_mesh());
    const Cut_Part<Element> cutK(Th.get_cut_part(k, 0));
    const FElement &FK(Vh[k]);
    const Element &K(FK.T);

    if (cutK.multi_interface()) {
        assert(0);
    } // not handled yet

    double meas = K.measure();
    double h    = K.get_h();
    int domain  = FK.get_domain();
    int kb      = Vh.idxElementInBackMesh(k);

    // GET THE QUADRATURE RULE
    const QF &qf(this->get_quadrature_formular_cutK());

    // LOOP OVER ELEMENTS IN THE CUT
    for (auto it = cutK.other_side_element_begin(); it != cutK.other_side_element_end(); ++it) {

        double meas_cut = cutK.measure(it);

        // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
        for (int l = 0; l < VF.size(); ++l) {
            if (!VF[l].on(domain))
                continue;
            // FINTE ELEMENT SPACES && ELEMENTS
            const FESpace &Vhv(VF.get_spaceV(l));
            const FElement &FKv(Vhv[k]);
            this->initIndex(FKv, FKv);

            // BF MEMORY MANAGEMENT -
            int lastop = getLastop(VF[l].du, VF[l].dv);
            RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N,
                     lastop); //  the value for basic fonction
            What_d Fop = Fwhatd(lastop);

            // COMPUTE COEFFICIENT
            double coef = VF[l].computeCoefElement(h, meas_cut, meas, meas_cut, domain);

            // LOOP OVER QUADRATURE IN SPACE
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

                typename QF::QuadraturePoint ip(qf[ipq]);
                Rd mip      = cutK.mapToPhysicalElement(it, ip); // to the physical cut part
                Rd cut_ip   = K.mapToReferenceElement(mip);      // back to the cut part in reference element
                double Cint = meas_cut * ip.getWeight();

                // EVALUATE THE BASIS FUNCTIONS
                FKv.BF(Fop, cut_ip, fv);
                //   VF[l].applyFunNL(fu,fv);

                // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
                // Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip);
                Cint *= coef * VF[l].c;
                // Cint = coef;
                this->addToMatrix(VF[l], FKv, fv, Cint);
            }
        }
    }
}

template <typename M>
void BaseCutFEM<M>::addLagrangeVecToRowAndCol(const std::span<double> vecRow, const std::span<double> vecCol,
                                              const R val_rhs) {
    int ndf = this->rhs_.size();
    this->rhs_.resize(ndf + 1);
    this->rhs_[ndf] = val_rhs;

    this->index_j0_[0] = 0;
    this->index_i0_[0] = 0;

    for (int idx = 0; idx < vecRow.size(); idx++) {
        // for (int idx = 0; idx < ndf; idx++) {
        //  this->mat_[0][std::make_pair(idx, ndf)] = vecCol[idx];
        //  this->mat_[0][std::make_pair(ndf, idx)] = vecRow[idx];
        (*this)(idx, ndf) += vecCol[idx];
        (*this)(ndf, idx) += vecRow[idx];
    }
}




template <typename M>
void BaseCutFEM<M>::addBilinearIntersection(const itemVFlist_t &VF, const CutMesh &funMesh, const CutMesh &activeMesh, const CBorder &b, std::list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);

    for (int idx_be = funMesh.first_boundary_element(); idx_be < funMesh.last_boundary_element();
         idx_be += funMesh.next_boundary_element()) {

        int ifac;
        const int kb          = funMesh.Th.BoundaryElement(idx_be, ifac);
        std::vector<int> idxK = activeMesh.idxAllElementFromBackMesh(kb, -1);
        if (idxK.size() == 0)
            continue;

        const Element &K(activeMesh.Th[kb]);
        const BorderElement &BE(activeMesh.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {

            // CHECK IF IT IS A CUT EDGE
            if (activeMesh.isCutFace(idxK[0], ifac, 0)) {
                addIntersectedBorderContribution(VF, activeMesh, K, BE, ifac, nullptr, 0, 1.);
                // std::cout << idxK[0] << " : " << ifac << std::endl;
            }
            else continue;
            this->addLocalContribution();
        }
    }
}

template <typename M>
void BaseCutFEM<M>::addLinearIntersection(const itemVFlist_t &VF, const CutMesh &funMesh, const CutMesh &activeMesh, const CBorder &b, std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);

    for (int idx_be = funMesh.first_boundary_element(); idx_be < funMesh.last_boundary_element();
         idx_be += funMesh.next_boundary_element()) {

        int ifac;
        const int kb          = funMesh.Th.BoundaryElement(idx_be, ifac);
        std::vector<int> idxK = activeMesh.idxAllElementFromBackMesh(kb, -1);
        if (idxK.size() == 0)
            continue;
        const Element &K(activeMesh.Th[kb]);
        const BorderElement &BE(activeMesh.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {

            // CHECK IF IT IS A CUT EDGE
            if (activeMesh.isCutFace(idxK[0], ifac, 0)) {
                BaseCutFEM<M>::addIntersectedBorderContribution(VF, activeMesh, K, BE, ifac, nullptr, 0, 1.);

            } else continue;
        }
    }
}


template <typename M>
void BaseCutFEM<M>::addIntersectedBorderContribution(const itemVFlist_t &VF, const CutMesh &activeMesh, const Element &K, const BorderElement &BE, int ifac,
                                          const TimeSlab *In, int itq, double cst_time) {

    typedef typename FElement::RdHatBord RdHatBord;

    // Compute parameter connected to the mesh.
    double measK = K.measure();
    double h     = K.get_h();
    Rd normal    = K.N(ifac);

    // std::cout << "K[0] = " << K[0] << ", K[1] = " << K[1] << ", K[2] = " << K[2] << "\n";
    // std::cout << "normal = " << normal << "\n";

    // U and V HAS TO BE ON THE SAME MESH
    const FESpace &Vh(VF.get_spaceV(0));
    const CutMesh &Th(Vh.get_mesh());
    int kb                = Vh.Th(K);
    // std::vector<int> idxK = Vh.idxAllElementFromBackMesh(kb, -1);
    std::vector<int> idxK_active = activeMesh.idxAllElementFromBackMesh(kb, -1);

    // assert(idxK.size() == 2);

    // GET THE QUADRATURE RULE
    const QFB &qfb(this->get_quadrature_formular_cutFace());
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    int Ne = idxK_active.size();
    assert(Ne == 1);
    for (int e = 0; e < Ne; ++e) {
        // int k_surf = idxK[e];
        int k_active = idxK_active[e];
        typename Element::Face face;
        const Cut_Part<typename Element::Face> cutFace(activeMesh.get_cut_face(face, k_active, ifac, itq));
        const Cut_Part<Element> cutK(activeMesh.get_cut_part(k_active, itq));
        // std::cout << "k_active = " << k_active << std::endl;

        double meas_cut = cutK.measure();
        // LOOP OVER ELEMENTS IN THE CUT
        for (auto it = cutFace.element_begin(); it != cutFace.element_end(); ++it) {

            double meas = cutFace.measure(it);


            for (int l = 0; l < VF.size(); ++l) {
                // if(!VF[l].on(domain)) continue;

                // FINITE ELEMENT SPACES && ELEMENTS
                const FESpace &Vhv(VF.get_spaceV(l));
                const FESpace &Vhu(VF.get_spaceU(l));
                std::vector<int> idxKu = Vhu.idxAllElementFromBackMesh(kb, -1);
                std::vector<int> idxKv = Vhv.idxAllElementFromBackMesh(kb, -1);
                int ku_surf = idxKu[e];
                int kv_surf = idxKv[e];
                // assert(Vhv.get_nb_element() == Vhu.get_nb_element());
                bool same = (VF.isRHS() || (&Vhu == &Vhv));
                const FElement &FKu(Vhu[ku_surf]);

                const FElement &FKv(Vhv[kv_surf]);

                int domain = FKv.get_domain();

                this->initIndex(FKu, FKv);
                // BF MEMORY MANAGEMENT -
                int lastop = getLastop(VF[l].du, VF[l].dv);
                RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
                RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
                What_d Fop = Fwhatd(lastop);

                // COMPUTE COEFFICIENT && NORMAL
                double coef = VF[l].computeCoefElement(h, meas, measK, meas_cut, domain);
                coef *= VF[l].computeCoefFromNormal(normal);
                // LOOP OVER QUADRATURE IN SPACE
                for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
                    typename QFB::QuadraturePoint ip(qfb[ipq]);
                    const Rd mip    = cutFace.mapToPhysicalElement(it, (RdHatBord)ip);
                    const Rd cut_ip = K.mapToReferenceElement(mip);
                    double Cint     = meas * ip.getWeight() * cst_time;

                    // EVALUATE THE BASIS FUNCTIONS
                    FKv.BF(Fop, cut_ip, fv);
                    if (!same)
                        FKu.BF(Fop, cut_ip, fu);
                    //   VF[l].applyFunNL(fu,fv);

                    Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid, normal);
                    Cint *= coef * VF[l].c;

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
    }
}




template <typename M> void BaseCutFEM<M>::addBilinearInner(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        bar += Th.next_element();

        if (Th.isCut(k, 0)) {
            continue;
        } else {
            BaseFEM<M>::addElementContribution(VF, k, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M> void BaseCutFEM<M>::addLinearInner(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        bar += Th.next_element();

        if (Th.isCut(k, 0)) {
            continue;
        } else {
            BaseFEM<M>::addElementContribution(VF, k, nullptr, 0, 1.);
        }
    }
    bar.end();
}
template <typename M> void BaseCutFEM<M>::addBilinearInnerBorder(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        bar += Th.next_element();

        // Loop only over cut elements for efficiency
        if (!Th.isCut(k, 0)) {
            continue;
        } 
        
        // Find the face that neighbors an interior element
        for (int ifac = 0; ifac < Element::nea; ++ifac) { 
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            // an neighboring element can be 1) outside of the domain, 2) another cut element or 3) inside of the domain 
            if (kn == -1) continue;
            if (Th.isCut(kn, 0)) continue;

            BaseFEM<M>::addInnerBorderContribution(VF, k, ifac, nullptr, 0, 1.);
            
        }
        
        this->addLocalContribution();
    }
    bar.end();
}
template <typename M> void BaseCutFEM<M>::addLinearInnerBorder(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        bar += Th.next_element();

        // Loop only over cut elements for efficiency
        if (!Th.isCut(k, 0)) {
            continue;
        } 
        
        // Find the face that neighbors an interior element
        for (int ifac = 0; ifac < Element::nea; ++ifac) { 
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            // an neighboring element can be 1) outside of the domain, 2) another cut element or 3) inside of the domain 
            if (kn == -1) continue;
            if (Th.isCut(kn, 0)) continue;

            BaseFEM<M>::addInnerBorderContribution(VF, k, ifac, nullptr, 0, 1.);
            
        }
        
        this->addLocalContribution();
    }
    bar.end();
}
template <typename M> void BaseCutFEM<M>::addBilinearInnerBorderMixed(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        bar += Th.next_element();

        // Loop only over cut elements for efficiency
        if (!Th.isCut(k, 0)) {
            continue;
        } 
        
        // Find the face that neighbors an interior element
        for (int ifac = 0; ifac < Element::nea; ++ifac) { 
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            // an neighboring element can be 1) outside of the domain, 2) another cut element or 3) inside of the domain 
            if (kn == -1) continue;
            if (Th.isCut(kn, 0)) continue;

            BaseFEM<M>::addInnerBorderContributionMixed(VF, k, ifac, nullptr, 0, 1.);
            
        }
        
        this->addLocalContribution();
    }
    bar.end();
}
template <typename M> void BaseCutFEM<M>::addLinearInnerBorderMixed(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        bar += Th.next_element();

        // Loop only over cut elements for efficiency
        if (!Th.isCut(k, 0)) {
            continue;
        } 
        
        // Find the face that neighbors an interior element
        for (int ifac = 0; ifac < Element::nea; ++ifac) { 
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            // an neighboring element can be 1) outside of the domain, 2) another cut element or 3) inside of the domain 
            if (kn == -1) continue;
            if (Th.isCut(kn, 0)) continue;

            BaseFEM<M>::addInnerBorderContributionMixed(VF, k, ifac, nullptr, 0, 1.);
            
        }
        
        this->addLocalContribution();
    }
    bar.end();
}


template <typename M> void BaseCutFEM<M>::addBilinearOuterBorder(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        bar += Th.next_element();

        // Loop only over cut elements for efficiency
        if (!Th.isCut(k, 0)) {
            continue;
        } 
        
        // Find the face that neighbors an exterior element
        for (int ifac = 0; ifac < Element::nea; ++ifac) { 
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            // an neighboring element can be 1) outside of the domain, 2) another cut element or 3) inside of the domain 
            if (kn != -1) continue;
            
            BaseFEM<M>::addOuterBorderContribution(VF, k, ifac, nullptr, 0, 1.);
            
        }
        
        this->addLocalContribution();
    }
    bar.end();
}
template <typename M> void BaseCutFEM<M>::addLinearOuterBorder(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        bar += Th.next_element();

        // Loop only over cut elements for efficiency
        if (!Th.isCut(k, 0)) {
            continue;
        } 
        
        // Find the face that neighbors an interior element
        for (int ifac = 0; ifac < Element::nea; ++ifac) { 
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            // an neighboring element can be 1) outside of the domain, 2) another cut element or 3) inside of the domain 
            if (kn != -1) continue;
            
            BaseFEM<M>::addOuterBorderContribution(VF, k, ifac, nullptr, 0, 1.);
            
        }
        
        this->addLocalContribution();
    }
    bar.end();
}
template <typename M> void BaseCutFEM<M>::addBilinearOuterBorderMixed(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        bar += Th.next_element();

        // Loop only over cut elements for efficiency
        if (!Th.isCut(k, 0)) {
            continue;
        } 
        
        // Find the face that neighbors an interior element
        for (int ifac = 0; ifac < Element::nea; ++ifac) { 
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            // an neighboring element can be 1) outside of the domain, 2) another cut element or 3) inside of the domain 
            if (kn != -1) continue;
            int kb = Th.idxElementInBackMesh(k);
            BaseFEM<M>::addOuterBorderContributionMixed(VF, kb, ifac, nullptr, 0, 1.);
            
        }
        
        this->addLocalContribution();
    }
    bar.end();
}
template <typename M> void BaseCutFEM<M>::addLinearOuterBorderMixed(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        bar += Th.next_element();

        // Loop only over cut elements for efficiency
        if (!Th.isCut(k, 0)) {
            continue;
        } 
        
        // Find the face that neighbors an interior element
        for (int ifac = 0; ifac < Element::nea; ++ifac) { 
            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);

            // an neighboring element can be 1) outside of the domain, 2) another cut element or 3) inside of the domain 
            if (kn != -1) continue;

            int kb = Th.idxElementInBackMesh(k);
            
            BaseFEM<M>::addOuterBorderContributionMixed(VF, kb, ifac, nullptr, 0, 1.);
            
        }
        
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M> void BaseCutFEM<M>::addElementStabilization(const itemVFlist_t &VF, const Interface<M> &gamma, const CutMesh &Th) {
    assert(Th.get_nb_domain()==1);
    assert(!VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
#pragma omp parallel for num_threads(this->get_num_threads())
    
    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        
        const int kb  = gamma.idxElementOfFace(iface);
        const int dom = 0;
        const int k   = Th.idxElementFromBackMesh(kb, dom);   // ASSUMES nb_domain = 1

        this->addElementInterfaceContribution(VF, gamma, iface, k, nullptr, 1., 0);

        this->addLocalContribution();
    }
    bar.end();
}



//! CHECK DOCUMENTATION ON THE TWO BELOW METHODS
/**
 * @brief Initializes the solution vector `u0` based on the `mapU0_` data.
 *
 * This function initializes the solution vector `u0` based on the `mapU0_` data,
 * which is a map of initial values for each degree of freedom. The map is cleared
 * at the end of the function.
 *
 * @tparam M The mesh type.
 * @param u0 The solution vector to be initialized.
 */

template <typename M> void BaseCutFEM<M>::initialSolution(std::span<double> u0) {
    // Note: this method changes the input vector u0

    // Get the number of degrees of freedom in time
    int nbTime = this->get_nb_dof_time();

    // Initialize u0 with the number of degrees of freedom
    // u0.init(this->get_nb_dof());

    assert(u0.size() == this->get_nb_dof());

    // If the mapU0_ is empty, return without performing any further operations
    if (this->mapU0_.size() == 0) {
        return;
    }

    // Initialize the id of the domain to 0
    int id_domain_0 = 0;

    // Loop through the solutions corresponding to different subdomains (and thus different FE spaces)
    for (auto q = this->mapIdx0_.begin(); q != this->mapIdx0_.end(); ++q) {

        // Get the FESpace object from the map
        const FESpace &Wh = *q->first;

        const int n0 = q->second;

        // Get the active mesh object from the FESpace object
        const ActiveMesh<M> &Th(Wh.get_mesh());

        // Get the back space from the FESpace object
        const FESpace &backVh = Wh.get_back_space();

        // Create a temporary variable u0S as a subarray of u0
        KN_<double> u0S = u0.subspan(n0, Wh.NbDoF() * nbTime);

        // Loop through all the elements of the active mesh
        for (int k = 0; k < Th.get_nb_element(); ++k) {

            // Only interested in elements that have intersection with the domain
            // in the first time quadrature point
            if (Th.isInactive(k, 0))
                continue; // has no intersection with domain at itq = 0

            // Get the FElement object for the current element
            const FElement &FK(Wh[k]);

            // Get the domain of the current element
            int domain = Th.get_domain_element(k);

            // If the domain is -1, set the id of the domain to id_domain_0
            // otherwise, set it to id_domain_0 + domain
            int id_domain = (domain == -1) ? id_domain_0 : id_domain_0 + domain;

            // Get the index of the current element in the back mesh
            int kb = Th.idxElementInBackMesh(k);

            // Get the FElement object for the corresponding element in the back mesh
            const FElement &FKback(backVh[kb]);

            // Loop over the components of the FE space (1 for scalar problems)
            for (int ic = 0; ic < Wh.N; ++ic) { // ESSAYER VH->N

                // Loop through all the degrees of freedom of the element
                for (int i = FK.dfcbegin(ic); i < FK.dfcend(ic); ++i) {
                    // Get the value from the map of initial conditions, with the key of the current domain and node in
                    // the back mesh.
                    //
                    u0S[FK(i)] = this->mapU0_[std::make_pair(id_domain, FKback(i))];
                }
            }
        }

        // Set id_domain_0 to the number of subdomains
        id_domain_0 += Th.get_nb_domain();
    }

    // Clear the map of initial conditions.
    this->mapU0_.clear();
}


template <typename M> void BaseCutFEM<M>::saveSolutionBackMesh(std::span<double> sol, FunFEM<M>& f_back) {


    std::span<double> u_back = f_back.array();
    std::ranges::fill(u_back, 0.);

    this->mapU0_.clear(); 

    int id_domain_0 = 0;                       
    int nbTime      = this->get_nb_dof_time(); 

    // Iterate over the finite element spaces in the finite element system
    for (typename std::map<const FESpace *, int>::const_iterator q = this->mapIdx0_.begin(); q != this->mapIdx0_.end();
         ++q) {
        const FESpace &Wh = *q->first;               // Get the finite element space.
        const int n0      = q->second;               // Get the starting index of the finite element space.
        const ActiveMesh<M> &Th(Wh.get_mesh());      // Get the mesh associated with the finite element space.
        const FESpace &backVh = Wh.get_back_space(); // Get the back space of the finite element space.

        assert(Th.get_nb_domain() == 1);
        // Pointer to the vector of coefficients of size nbTime*N_{h,i}^n
        // std::cout << "Wh.get_nb_dof() = " << Wh.get_nb_dof() << ", n0 = " << n0 << "\n";
        const KN_<double> solS = sol.subspan(n0, Wh.get_nb_dof() * nbTime);

        // Iterate over the elements in the current active mesh.
        for (int k = 0; k < Th.get_nb_element(); ++k) {

            const FElement &FK(Wh[k]); // Get the finite element associated with the current element in the mesh.
            const int domain = Th.get_domain_element(k); // Get the domain of the current element.
            int id_domain    = (domain == -1) ? id_domain_0 : id_domain_0 + domain; // Compute the domain ID.

            int kb = Th.idxElementInBackMesh(k); // Get the index of the current element in the back mesh.
            const FElement &FKback(backVh[kb]);  // Get the corresponding element in the back space.

            // Loop over the components of the FE space (1 for scalar problems)
            for (int ic = 0; ic < Wh.N; ++ic) {
                // Iterate over the degrees of freedom in the finite element (e.g. nodes).
                for (int i = FK.dfcbegin(ic); i < FK.dfcend(ic); ++i) {
                    R val = 0.; // Initialize the coefficient value to 0.

                    // Sum up the coefficients corresponding to the same space DOF in different time DOFs.
                    for (int it = 0; it < nbTime; ++it) {
                        val += solS[FK.loc2glb(i, it)];
                    }

                    // Store the coefficient value in the DOF in the background FE space, to be able to retreive it in
                    // the next time slab
                    this->mapU0_[std::make_pair(id_domain, FKback(i))] = val;

                    u_back[FKback(i)] = val;
                }
            }
        }

        id_domain_0 += Th.get_nb_domain();
    }

}


template <typename M> void BaseCutFEM<M>::saveSolution(std::span<double> sol) {
    
    this->mapU0_.clear(); 

    int id_domain_0 = 0;                       
    int nbTime      = this->get_nb_dof_time(); 

    // Iterate over the finite element spaces in the finite element system
    for (typename std::map<const FESpace *, int>::const_iterator q = this->mapIdx0_.begin(); q != this->mapIdx0_.end();
         ++q) {
        const FESpace &Wh = *q->first;               // Get the finite element space.
        const int n0      = q->second;               // Get the starting index of the finite element space.
        const ActiveMesh<M> &Th(Wh.get_mesh());      // Get the mesh associated with the finite element space.
        const FESpace &backVh = Wh.get_back_space(); // Get the back space of the finite element space.

        // Pointer to the vector of coefficients of size nbTime*N_{h,i}^n
        // std::cout << "Wh.get_nb_dof() = " << Wh.get_nb_dof() << ", n0 = " << n0 << "\n";
        const KN_<double> solS = sol.subspan(n0, Wh.get_nb_dof() * nbTime);

        // Iterate over the elements in the current active mesh.
        for (int k = 0; k < Th.get_nb_element(); ++k) {

            const FElement &FK(Wh[k]); // Get the finite element associated with the current element in the mesh.
            const int domain = Th.get_domain_element(k); // Get the domain of the current element.
            int id_domain    = (domain == -1) ? id_domain_0 : id_domain_0 + domain; // Compute the domain ID.

            int kb = Th.idxElementInBackMesh(k); // Get the index of the current element in the back mesh.
            const FElement &FKback(backVh[kb]);  // Get the corresponding element in the back space.

            // Loop over the components of the FE space (1 for scalar problems)
            for (int ic = 0; ic < Wh.N; ++ic) {
                // Iterate over the degrees of freedom in the finite element (e.g. nodes).
                for (int i = FK.dfcbegin(ic); i < FK.dfcend(ic); ++i) {
                    R val = 0.; // Initialize the coefficient value to 0.

                    // Sum up the coefficients corresponding to the same space DOF in different time DOFs.
                    for (int it = 0; it < nbTime; ++it) {
                        val += solS[FK.loc2glb(i, it)];
                    }

                    // Store the coefficient value in the DOF in the background FE space, to be able to retreive it in
                    // the next time slab
                    this->mapU0_[std::make_pair(id_domain, FKback(i))] = val;
                }
            }
        }

        id_domain_0 += Th.get_nb_domain();
    }
}

// dd
