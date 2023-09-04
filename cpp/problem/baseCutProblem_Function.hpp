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

    //  double t0 = MPIcf::Wtime();

#pragma omp parallel default(shared)
    {
#ifdef USE_OMP
        assert(this->get_nb_thread() == omp_get_num_threads());
        int thread_id = omp_get_thread_num();
#else
        int thread_id = 0;
#endif
        int verbose = (thread_id == 0) * globalVariable::verbose;
        progress bar(" Add Bilinear CutMesh", Th.last_element(), verbose, this->get_nb_thread());
#pragma omp for
        for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
            bar += Th.next_element();
            if (Th.isCut(k, 0)) {
                addElementContribution(VF, k, nullptr, 0, 1.);
            } else {
                BaseFEM<M>::addElementContribution(VF, k, nullptr, 0, 1.);
            }
            this->addLocalContribution();
        }
        bar.end();
    }
    // std::cout << " real time " << MPIcf::Wtime() - t0 << std::endl;
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
#pragma omp parallel default(shared)
    {
        // Get the thread ID
#ifdef USE_OMP
        // Check if the number of threads matches the number of threads set in this object
        assert(this->get_nb_thread() == omp_get_num_threads());
        int thread_id = omp_get_thread_num();
#else
        int thread_id = 0;
#endif

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
        int verbose       = (thread_id == 0) * globalVariable::verbose;
        progress bar(title.c_str(), Th.last_element(), verbose, this->get_nb_thread());

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

template <typename M> void BaseCutFEM<M>::addLinear(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(VF.isRHS());
    progress bar(" Add Linear CutMesh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

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
#ifdef USE_OMP
    int iam = omp_get_thread_num();
#else
    int iam = 0;
#endif
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
void BaseCutFEM<M>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th, const CExtension &ext, const int espE) {
    assert(!VF.isRHS());
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

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

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();
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
            RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
            RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
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
            for (int ic = 0; ic < 1; ++ic) {
                for (int df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df) {

                    int id_item = FK.DFOnWhat(df);

                    if (id_item < K.nv) {
                        assert(0);
                    } else if (id_item < K.nv + K.ne) {
                        // std::cout << " on edge  " <<FK.DFOnWhat(df) << std::endl;
                        int id_face = id_item - K.nv;
                        if (id_face == ifac) {
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
    progress bar("Add Face Stabilization CutMesh", Th.last_element(), globalVariable::verbose);

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

            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);
            BaseFEM<M>::addFaceContribution(VF, e1, e2, nullptr, 0, 1.);
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
            if (kn < k)
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


template <typename M>
void BaseCutFEM<M>::addPatchStabilization(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(!VF.isRHS());
    progress bar("Add Patch Stabilization CutMesh", Th.last_element(), globalVariable::verbose);

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

            std::pair<int, int> e1 = std::make_pair(k, ifac);
            std::pair<int, int> e2 = std::make_pair(kn, jfac);
            BaseFEM<M>::addPatchContribution(VF, k, kn, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M>
void BaseCutFEM<M>::addPatchStabilization(const itemVFlist_t &VF, const CutMesh &Th, const TimeSlab &In) {

    int number_of_quadrature_points = this->get_nb_quad_point_time();
    
    // Loop through time quadrature points
    for (int itq = 0; itq < number_of_quadrature_points; ++itq) {
        assert(!VF.isRHS());

        // Compute contribution from time basis functions
        auto tq    = this->get_quadrature_time(itq);
        double tid = In.map(tq);
        KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
        RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
        In.BF(tq.x, bf_time); // compute time basic funtions
        double cst_time   = tq.a * In.get_measure();
        
        std::string title = " Add Patch Stab, In(" + std::to_string(itq) + ")";
        progress bar(title.c_str(), Th.last_element(), globalVariable::verbose);

        // Loop through active mesh elements
        for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
            bar += Th.next_element();
            
            // Exclude elements whose edges do not need stabilization
            if (!Th.isStabilizeElement(k))
                continue;

            // Loop through the element's edges
            for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

                int jfac = ifac;
                int kn   = Th.ElementAdj(k, jfac);      // get neighbor element's index

                // By skipping neighbors with smaller indices, we avoid adding contribution to the same edge twice
                if (kn < k)
                    continue;

                std::pair<int, int> e1 = std::make_pair(k, ifac);   // (element index, edge index) current element
                std::pair<int, int> e2 = std::make_pair(kn, jfac);  // (element index, edge index) neighbor element

                // Add patch contribution
                //BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
                BaseFEM<M>::addPatchContribution(VF, k, kn, &In, itq, cst_time);
            }
            this->addLocalContribution();
        }
        bar.end();
    }
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
            if (kn < k)
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
            if (kn < k)
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

                std::pair<int, int> e1 = std::make_pair(k, ifac);
                std::pair<int, int> e2 = std::make_pair(kn, jfac);

                // number_of_stabilized_edges += 1;
                BaseFEM<M>::addPatchContribution(VF, k, kn, &In, itq, cst_time);
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
    this->rhs_(ndf) = val;
    progress bar(" Add Lagrange Multiplier Kh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();
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
    this->rhs_(ndf) = val;

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
    this->rhs_(ndf) = val;
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
        this->rhs_(ndf) = val;
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
    this->rhs_(ndf) = val;

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
    this->rhs_(ndf) = val;

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
    this->rhs_(ndf) = val_rhs;

    for (int idx = 0; idx < ndf; idx++) {
        this->mat_[0][std::make_pair(idx, ndf)] = vecCol[idx];
        this->mat_[0][std::make_pair(ndf, idx)] = vecRow[idx];
    }
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

template <typename M>
template <typename V>
    requires NonAllocVector<V> || std::is_same_v<V, KN<typename V::element_type>>
void BaseCutFEM<M>::initialSolution(V &u0) {
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

        // Get the second value from the map, which is n0
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

/**
 * @brief Save the coefficients of the solution.
 * @note The coefficients are stored in a map with keys representing the domain and the global index of the degrees of
 * freedom (DOFs) in the back space of the finite element space, and values representing the coefficients of the
 * solution.
 *
 * @param sol The vector of coefficients representing the numerical solution.
 * @tparam M The type of the mesh.
 * @tparam V The type of the vector for storing the coefficients.
 * @requires V to be a NonAllocVector or a vector with element type KN<typename V::element_type>.
 */
template <typename M>
template <typename V>
    requires NonAllocVector<V> || std::is_same_v<V, KN<typename V::element_type>>
void BaseCutFEM<M>::saveSolution(const V sol) {
    // Note: this method doesn't change the input sol

    this->mapU0_.clear(); // Clear the map of coefficients.

    int id_domain_0 = 0;                       // Initialize the domain ID to 0.
    int nbTime      = this->get_nb_dof_time(); // Get the number of degrees of freedom in time.

    // Iterate over the finite element spaces in the map of indices.
    for (typename std::map<const FESpace *, int>::const_iterator q = this->mapIdx0_.begin(); q != this->mapIdx0_.end();
         ++q) {
        const FESpace &Wh = *q->first;               // Get the finite element space.
        const int n0      = q->second;               // Get the starting index of the finite element space.
        const ActiveMesh<M> &Th(Wh.get_mesh());      // Get the mesh associated with the finite element space.
        const FESpace &backVh = Wh.get_back_space(); // Get the back space of the finite element space.

        // Pointer to the vector of coefficients of size nbTime*N_{h,i}^n
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
