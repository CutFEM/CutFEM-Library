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

/**
 * @brief This function adds the contribution to the local matrix from the given test and trial functions.
 * This function is called for each integration point, and is used for both the element and the boundary element case.
 */
template <typename M>
void BaseFEM<M>::addToMatrix(const itemVF_t &VFi, const FElement &FKu, const FElement &FKv, const RNMK_ &fu,
                             const RNMK_ &fv, double Cint) {
    // Get the current thread ID
#ifdef USE_OMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif

    // For each local DOF in the test function, i
    for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
        // For each local DOF in the trial function, j
        for (int j = FKu.dfcbegin(VFi.cu); j < FKu.dfcend(VFi.cu); ++j) {
            // Compute the contribution to the local matrix
            this->addToLocalContribution(FKv.loc2glb(i), FKu.loc2glb(j), thread_id) +=
                Cint * fv(i, VFi.cv, VFi.dv) * fu(j, VFi.cu, VFi.du);
        }
    }
}

/**
 * @brief Adds a contribution to the local matrix, using the time basis function.
 *
 * @tparam M The mesh
 * @param VFi A vector of basis functions.
 * @param In The time slab.
 * @param FKu The FElement of the trial function.
 * @param FKv The FElement of the test function.
 * @param fu The trial function.
 * @param fv The test function.
 * @param Cint The integration constant.
 *
 * @note
 * - The function loops over the time integration points, then over the spatial basis functions,
 *   and calls the function `addToLocalContribution` to add the contribution to the local matrix.
 */

template <typename M>
void BaseFEM<M>::addToMatrix(const itemVF_t &VFi, const TimeSlab &In, const FElement &FKu, const FElement &FKv,
                             const RNMK_ &fu, const RNMK_ &fv, double Cint) {
    // Get the thread id
#ifdef USE_OMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif

    // Get the offset
    long offset = thread_id * this->offset_bf_time;

    // Get the time basis function
    RNMK_ bf_time(this->databf_time_ + offset, In.NbDoF(), 1, op_dz);

    // Loop over the time integration points
    for (int it = In.dfcbegin(0); it < In.dfcend(0); ++it) {
        for (int jt = In.dfcbegin(0); jt < In.dfcend(0); ++jt) {

            // Get the value of the time basis function at the time integration points
            const R tval = bf_time(it, 0, VFi.dtv) * bf_time(jt, 0, VFi.dtu);
            
            // Loop over the spatial basis functions
            for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
                for (int j = FKu.dfcbegin(VFi.cu); j < FKu.dfcend(VFi.cu); ++j) {

                    // Add the contribution to the local matrix
                    this->addToLocalContribution(FKv.loc2glb(i, it), FKu.loc2glb(j, jt), thread_id) +=
                        Cint * tval * fv(i, VFi.cv, VFi.dv) * fu(j, VFi.cu, VFi.du);
                }
            }
        }
    }
}



template <typename M>
void BaseFEM<M>::addToMatrixSpecial(const itemVF_t &VFi, const TimeSlab &In, const FElement &FKu, const FElement &FKv,
                             const RNMK_ &fu, const RNMK_ &fv, double Cint, int from_time_dof_u, int from_time_dof_v) {
    // Get the thread id
#ifdef USE_OMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif

    // Get the offset
    long offset = thread_id * this->offset_bf_time;

    // Get the time basis function
    RNMK_ bf_time(this->databf_time_ + offset, In.NbDoF(), 1, op_dz);

    // Loop over the time integration points
    for (int it = from_time_dof_v; it < In.dfcend(0); ++it) {
        for (int jt = from_time_dof_u; jt < In.dfcend(0); ++jt) {

            // Get the value of the time basis function at the time integration points
            const R tval = bf_time(it, 0, VFi.dtv) * bf_time(jt, 0, VFi.dtu);

            // Loop over the spatial basis functions
            for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
                for (int j = FKu.dfcbegin(VFi.cu); j < FKu.dfcend(VFi.cu); ++j) {

                    // Add the contribution to the local matrix
                    this->addToLocalContribution(FKv.loc2glb(i, it), FKu.loc2glb(j, jt), thread_id) +=
                        Cint * tval * fv(i, VFi.cv, VFi.dv) * fu(j, VFi.cu, VFi.du);
                }
            }
        }
    }
}





/**
 * @brief Adds the contribution to the global matrix from the current cell, using the test function and
 * integration constant.
 *
 * @tparam M The mesh
 * @param VFi The object that contains information about the vector of basis functions and cell information
 * @param FKv The FElement of the test function
 * @param fv The test function values
 * @param Cint The integration constant
 *
 * The function loops over the local degrees of freedom (DOF) associated with the current cell and adds
 * the contribution of the current local DOF to the total local matrix using the function `addToLocalContribution`.
 */
template <typename M>
void BaseFEM<M>::addToMatrix(const itemVF_t &VFi, const FElement &FKv, const RNMK_ &fv, double Cint) {
    // Loop over the local degrees of freedom (DOF) associated with the
    // current cell
    for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
        // Add the contribution of the current local DOF to the global
        // matrix
        this->addToLocalContribution(FKv.loc2glb(i), 0) += Cint * fv(i, VFi.cv, VFi.dv);
    }
}

/**
 * @brief Adds the contribution to the global matrix from the current cell, using the test function,
 * integration constant, and time basis function.
 *
 * @tparam M The type of matrix
 * @param VFi The object that contains information about the vector of basis functions and cell information
 * @param In The time slab information
 * @param FKv The FElement of the test function
 * @param fv The test function values
 * @param Cint The integration constant
 *
 * The function loops over the local degrees of freedom (DOF) associated with the current cell and the time
 * integration points, and adds the contribution of the current local DOF to the global matrix, taking into
 * account the time basis function, using the function `addToLocalContribution`.
 */
template <typename M>
void BaseFEM<M>::addToMatrix(const itemVF_t &VFi, const TimeSlab &In, const FElement &FKv, const RNMK_ &fv,
                             double Cint) {

    // Step 1: Create a reference to the local time basis functions (bf_time).
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);

    // Step 2: Loop over the time degrees of freedom.
    for (int it = In.dfcbegin(0); it < In.dfcend(0); ++it) {
        
        // Step 3: Loop over the spatial degrees of freedom.
        for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {

            // Step 4: Add the local contribution to the global matrix.
            this->addToLocalContribution(FKv.loc2glb(i, it), 0) +=
                Cint * fv(i, VFi.cv, VFi.dv) * bf_time(it, 0, VFi.dtv);
        }
    }
}

template <typename M>
void BaseFEM<M>::addToMatrix_Opt(const itemVF_t &VFi, const FElement &FK, const RNMK_ &fv, double Cint) {

    for (int i = FK.dfcbegin(VFi.cv); i < FK.dfcend(VFi.cv); ++i) {
        for (int j = FK.dfcbegin(VFi.cu); j < FK.dfcend(VFi.cu); ++j) {
            this->addToLocalContribution_Opt(i, j) += Cint * fv(i, VFi.cv, VFi.dv) * fv(j, VFi.cu, VFi.du);
        }
    }
}

template <typename M>
void BaseFEM<M>::addToRHS(const itemVF_t &VFi, const FElement &FKv, const RNMK_ &fv, double Cint) {
    for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {

        (*this)(FKv.loc2glb(i)) += Cint * fv(i, VFi.cv, VFi.dv);
    }
}

template <typename M>
void BaseFEM<M>::addToRHS(const itemVF_t &VFi, const TimeSlab &In, const FElement &FKv, const RNMK_ &fv, double Cint) {
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);

    for (int it = In.dfcbegin(0); it < In.dfcend(0); ++it) {
        for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
            (*this)(FKv.loc2glb(i, it)) += Cint * fv(i, VFi.cv, VFi.dv) * bf_time(it, 0, VFi.dtv);
        }
    }
}

// Fun_h are always evaluated on backMesh
template <typename M> void BaseFEM<M>::addDiagonal(const FESpace &Qh, double val) {
    int idxBegin = this->mapIdx0_[&Qh];
    int idxEnd   = Qh.nbDoF;
    for (int i = idxBegin; i < idxEnd; ++i) {
        (*this->pmat_)[std::make_pair(i, i)] += val;
    }
}
template <typename M> void BaseFEM<M>::setDiagonal(const FESpace &Qh, double val) {
    int idxBegin = this->mapIdx0_[&Qh];
    int idxEnd   = Qh.nbDoF;
    for (int i = idxBegin; i < idxEnd; ++i) {
        (*this->pmat_)[std::make_pair(i, i)] = val;
    }
}

// INTEGRATION ON FULL ELEMENT
template <typename Mesh> void BaseFEM<Mesh>::addBilinear(const itemVFlist_t &VF, const Mesh &Th) {
    assert(!VF.isRHS());
    progress bar("Add Bilinear Mesh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();
        BaseFEM<Mesh>::addElementContribution(VF, k, nullptr, 0, 1.);

        this->addLocalContribution();
    }
    bar.end();
}

template <typename Mesh> void BaseFEM<Mesh>::addBilinear(const itemVFlist_t &VF, const CutMesh &Th) {
    assert(!VF.isRHS());
    progress bar("Add Bilinear Mesh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();
        BaseFEM<Mesh>::addElementContribution(VF, k, nullptr, 0, 1.);

        this->addLocalContribution();
    }
    bar.end();
}

template <typename Mesh> void BaseFEM<Mesh>::addLinear(const itemVFlist_t &VF, const Mesh &Th) {
    assert(VF.isRHS());
    progress bar("Add Linear Mesh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        BaseFEM<Mesh>::addElementContribution(VF, k, nullptr, 0, 1.);
    }
    bar.end();
}
template <typename M>
void BaseFEM<M>::addElementContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq,
                                        double cst_time) {

    // CHECK IF IT IS FOR RHS OR MATRIX
    bool to_rhs = VF.isRHS();

    // Compute parameter coonected to the mesh.
    // on can take the one from the first test function
    const FESpace &Vh(*VF[0].fespaceV);
    const FElement &FK(Vh[k]);
    const Element &K(FK.T);
    double meas = K.measure();
    double h    = K.get_h();
    int domain  = FK.get_domain();
    int kb      = Vh.idxElementInBackMesh(k);
#ifdef USE_OMP
    int iam = omp_get_thread_num();
#else
    int iam       = 0;
#endif

    // GET THE QUADRATURE RULE
    const QF &qf(this->get_quadrature_formular_K());
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
    for (int l = 0; l < VF.size(); ++l) {
        // if(!VF[l].on(domain)) continue;

        // FINTE ELEMENT SPACES && ELEMENTS
        const FESpace &Vhv(VF.get_spaceV(l));
        const FESpace &Vhu(VF.get_spaceU(l));
        const FElement &FKv(Vhv[k]);
        const FElement &FKu(Vhu[k]);
        this->initIndex(FKu, FKv); //, iam);

        // BF MEMORY MANAGEMENT -
        bool same   = (&Vhu == &Vhv);
        int lastop  = getLastop(VF[l].du, VF[l].dv);
        // RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop); //  the value for
        // basic fonction RNMK_ fu(this->databf_+ (same
        // ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value
        // for basic fonction
        long offset = iam * this->offset_bf_;
        RNMK_ fv(this->databf_ + offset, FKv.NbDoF(), FKv.N,
                 lastop); //  the value for basic fonction
        RNMK_ fu(this->databf_ + offset + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
        What_d Fop = Fwhatd(lastop);

        // COMPUTE COEFFICIENT
        double coef = VF[l].computeCoefElement(h, meas, meas, meas, domain);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
            typename QF::QuadraturePoint ip(qf[ipq]);
            const Rd mip = K.map(ip);
            double Cint  = meas * ip.getWeight() * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, ip, fv);
            if (!same)
                FKu.BF(Fop, ip, fu);
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
                    this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint); //, iam);
            }
        }
    }
}

// INTEGRATION ON INNER FACE
template <typename Mesh> void BaseFEM<Mesh>::addBilinear(const itemVFlist_t &VF, const Mesh &Th, const CFacet &b) {
    assert(!VF.isRHS());
    progress bar("Add Bilinear InnerEdge", Th.last_element(), globalVariable::verbose);
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
            BaseFEM<Mesh>::addFaceContribution(VF, e1, e2, nullptr, 0, 1.);
        }
        this->addLocalContribution();
    }
    bar.end();
}
template <typename Mesh> void BaseFEM<Mesh>::addLinear(const itemVFlist_t &VF, const Mesh &Th, const CFacet &b) {
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

            BaseFEM<Mesh>::addFaceContribution(VF, e1, e2, nullptr, 0, 1.);
        }
    }
}

template <typename M>
void BaseFEM<M>::addFaceContribution(const itemVFlist_t &VF, const std::pair<int, int> &e1,
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
    const FElement &FKi(Vh[ki]);
    const FElement &FKj(Vh[kj]);
    const Element &Ki(FKi.T);
    const Element &Kj(FKj.T);
    double measK = Ki.measure() + Kj.measure();
    double meas  = Ki.mesureBord(ifac);
    double h     = 0.5 * (Ki.get_h() + Kj.get_h());
    int domain   = FKi.get_domain();
    Rd normal    = Ki.N(ifac);

    // GET THE QUADRATURE RULE
    const QFB &qfb(this->get_quadrature_formular_dK());

    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
    for (int l = 0; l < VF.size(); ++l) {
        // if(!VF[l].on(domain)) continue;

        // FINTE ELEMENT SPACES && ELEMENTS
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

        // COMPUTE COEFFICIENT
        double coef = VF[l].computeCoefElement(h, meas, measK, measK, domain);
        coef *= VF[l].computeCoefFromNormal(normal);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
            typename QFB::QuadraturePoint ip(qfb[ipq]);
            const Rd ip_edge = Ki.mapToReferenceElement((RdHatBord)ip, ifac);
            const Rd mip     = Ki.mapToPhysicalElement(ip_edge);
            double Cint      = meas * ip.getWeight() * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, FKv.T.mapToReferenceElement(mip), fv);
            if (!same)
                FKu.BF(Fop, FKu.T.mapToReferenceElement(mip), fu);

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

template <typename M>
void BaseFEM<M>::addFaceContributionSpecial(const itemVFlist_t &VF, const std::pair<int, int> &e1,
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
    const FElement &FKi(Vh[ki]);
    const FElement &FKj(Vh[kj]);
    const Element &Ki(FKi.T);
    const Element &Kj(FKj.T);
    double measK = Ki.measure() + Kj.measure();
    double meas  = Ki.mesureBord(ifac);
    double h     = 0.5 * (Ki.get_h() + Kj.get_h());
    int domain   = FKi.get_domain();
    Rd normal    = Ki.N(ifac);

    // GET THE QUADRATURE RULE
    const QFB &qfb(this->get_quadrature_formular_dK());

    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
    for (int l = 0; l < VF.size(); ++l) {
        // if(!VF[l].on(domain)) continue;

        // FINTE ELEMENT SPACES && ELEMENTS
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

        // COMPUTE COEFFICIENT
        double coef = VF[l].computeCoefElement(h, meas, measK, measK, domain);
        coef *= VF[l].computeCoefFromNormal(normal);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
            typename QFB::QuadraturePoint ip(qfb[ipq]);
            const Rd ip_edge = Ki.mapToReferenceElement((RdHatBord)ip, ifac);
            const Rd mip     = Ki.mapToPhysicalElement(ip_edge);
            double Cint      = meas * ip.getWeight() * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, FKv.T.mapToReferenceElement(mip), fv);
            if (!same)
                FKu.BF(Fop, FKu.T.mapToReferenceElement(mip), fu);

            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kbu, kbv), std::make_pair(domain, domain),
                                                           mip, tid, normal);
            Cint *= coef * VF[l].c;

            if (In) {
                if (VF.isRHS())
                    assert(0);
                else
                    this->addToMatrixSpecial(VF[l], *In, FKu, FKv, fu, fv, Cint, 1, 0);
            } else {
                assert(0);
            }
        }
    }
}


// INTEGRATION ON BOUNDARY
template <typename Mesh>
void BaseFEM<Mesh>::addBilinear(const itemVFlist_t &VF, const Mesh &Th, const CBorder &b, std::list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);
    progress bar("Add Bilinear border", Th.last_boundary_element(), globalVariable::verbose);

    for (int idx_be = Th.first_boundary_element(); idx_be < Th.last_boundary_element();
         idx_be += Th.next_boundary_element()) {
        bar += Th.next_boundary_element();
        int ifac;
        const int kb = Th.BoundaryElement(idx_be, ifac);
        const Element &K(Th[kb]);
        const BorderElement &BE(Th.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {
            BaseFEM<Mesh>::addBorderContribution(VF, K, BE, ifac, nullptr, 0, 1.);
            this->addLocalContribution();
        }
    }
    bar.end();
}
template <typename Mesh>
void BaseFEM<Mesh>::addLinear(const itemVFlist_t &VF, const Mesh &Th, const CBorder &b, std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);
    progress bar("Add Bilinear border", Th.last_boundary_element(), globalVariable::verbose);

    for (int idx_be = Th.first_boundary_element(); idx_be < Th.last_boundary_element();
         idx_be += Th.next_boundary_element()) {
        bar += Th.next_boundary_element();

        int ifac;
        const int kb = Th.BoundaryElement(idx_be, ifac);

        const Element &K(Th[kb]);
        const BorderElement &BE(Th.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {

            // CHECK IF IT IS A CUT EDGE
            BaseFEM<Mesh>::addBorderContribution(VF, K, BE, ifac, nullptr, 0, 1.);
        }
    }
    bar.end();
}
template <typename M>
void BaseFEM<M>::addBorderContribution(const itemVFlist_t &VF, const Element &K, const BorderElement &BE, int iiifac,
                                       const TimeSlab *In, int itq, double cst_time) {

    int subDomId = iiifac / Element::nea;
    int ifac     = iiifac % Element::nea;
    typedef typename FElement::RdHatBord RdHatBord;

    // Compute parameter connected to the mesh.
    double measK = K.measure();
    double meas  = K.mesureBord(ifac);
    double h     = K.get_h();
    Rd normal    = K.N(ifac);

    // U and V HAS TO BE ON THE SAME MESH
    const FESpace &Vh(VF.get_spaceV(0));
    int kb                = Vh.Th(K);
    std::vector<int> idxK = Vh.idxAllElementFromBackMesh(kb, -1);
    // assert(idxK.size() == 1);
    // if(idxK.size() != 1) return;
    // when many subdomain are involve. Cut element but not cut boundary
    int k                 = idxK[subDomId]; // VF[0].onWhatElementIsTestFunction (idxK);

    // GET THE QUADRATURE RULE
    const QFB &qfb(this->get_quadrature_formular_dK());
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
    for (int l = 0; l < VF.size(); ++l) {
        // if(!VF[l].on(domain)) continue;

        // FINTE ELEMENT SPACES && ELEMENTS
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

        // COMPUTE COEFFICIENT
        double coef = VF[l].computeCoefElement(h, meas, measK, measK, domain);
        coef *= VF[l].computeCoefFromNormal(normal);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
            typename QFB::QuadraturePoint ip(qfb[ipq]);
            const Rd mip     = BE.mapToPhysicalElement((RdHatBord)ip);
            const Rd ip_edge = K.mapToReferenceElement(mip);
            double Cint      = meas * ip.getWeight() * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, ip_edge, fv);
            if (!same)
                FKu.BF(Fop, ip_edge, fu);

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

template <typename Mesh>
void BaseFEM<Mesh>::setDirichlet(const FunFEM<Mesh> &gh, const Mesh &Th, std::list<int> label) {

    bool all_label = (label.size() == 0);
    std::map<int, double> dof2set;
    const FESpace &Vh(gh.getSpace());

    for (int idx_be = Th.first_boundary_element(); idx_be < Th.last_boundary_element();
         idx_be += Th.next_boundary_element()) {

        int ifac;
        const int k = Th.BoundaryElement(idx_be, ifac);

        // assert(idxK.size() == 1);
        const Element &K(Th[k]);
        const BorderElement &BE(Th.be(idx_be));
        const FElement &FK(Vh[k]);

        if (util::contain(label, BE.lab) || all_label) {

            // for( int ic=0; ic<Vh.N;++ic) {
            for (int ic = 0; ic < 1; ++ic) {
                for (int df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df) {

                    int id_item = FK.DFOnWhat(df);

                    if (id_item < K.nv) {
                        bool is_on_border = false;
                        for (int i = 0; i < Element::nva; ++i) {

                            int i_e = Element::nvedge.at(ifac).at(i);
                            if (i_e == id_item) {
                                is_on_border = true;
                                break;
                            }
                        }
                        if (is_on_border) {
                            int df_glob = FK.loc2glb(df);
                            dof2set.insert({df_glob, gh(df_glob)});
                        }

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

    int N = this->get_nb_dof();
    eraseAndSetRow(this->get_nb_dof(), *this->pmat_[0], this->rhs_, dof2set);
}

// INTEGRATION ON INTERFACE
// On Faces
template <typename M>
void BaseFEM<M>::addBilinear(const itemVFlist_t &VF, const Interface<M> &gamma, std::list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);
    progress bar(" Add Bilinear Interface", gamma.last_element(), globalVariable::verbose);

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        bar += gamma.next_element();
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, gamma, iface, 0., nullptr, 1., 0);
            this->addLocalContribution();
        }
    }
    bar.end();
}

template <typename M>
void BaseFEM<M>::addBilinear(const itemVFlist_t &VF, const Interface<M> &gamma, const TimeSlab &In, int itq,
                             std::list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);
    auto tq        = this->get_quadrature_time(itq);
    double tid     = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions

    std::string title = " Add Bilinear Interface, In(" + std::to_string(itq) + ")";
    progress bar(title.c_str(), gamma.last_element(), globalVariable::verbose);

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        bar += gamma.next_element();
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, gamma, iface, tid, &In, 1., itq);
            this->addLocalContribution();
        }
    }

    bar.end();
}

template <typename M>
void BaseFEM<M>::addBilinear(const itemVFlist_t &VF, const TimeInterface<M> &gamma, const TimeSlab &In,
                             std::list<int> label) {

    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addBilinear(VF, gamma, In, itq, label);
    }
}

template <typename M>
void BaseFEM<M>::addBilinear(const itemVFlist_t &VF, const TimeInterface<M> &gamma, const TimeSlab &In, int itq,
                             std::list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);
    auto tq        = this->get_quadrature_time(itq);
    double tid     = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time   = tq.a * In.get_measure();
    std::string title = " Add Bilinear Interface, In(" + std::to_string(itq) + ")";
    progress bar(title.c_str(), gamma[itq]->last_element(), globalVariable::verbose);

    for (int iface = gamma[itq]->first_element(); iface < gamma[itq]->last_element();
         iface += gamma[itq]->next_element()) {
        bar += gamma[itq]->next_element();
        const typename Interface<M>::Face &face = (*gamma[itq])[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, *gamma[itq], iface, tid, &In, cst_time, itq);
            this->addLocalContribution();
        }
    }
    bar.end();
}

template <typename M>
void BaseFEM<M>::addLinear(const itemVFlist_t &VF, const Interface<M> &gamma, std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);
    progress bar(" Add Linear Interface", gamma.last_element(), globalVariable::verbose);

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        bar += gamma.next_element();

        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, gamma, iface, 0., nullptr, 1., 0);
        }
    }

    bar.end();
}

template <typename M>
void BaseFEM<M>::addLinear(const itemVFlist_t &VF, const Interface<M> &gamma, const TimeSlab &In, int itq,
                           std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);
    auto tq        = this->get_quadrature_time(itq);
    double tid     = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, gamma, iface, tid, &In, 1., itq);
        }
        this->addLocalContribution();
    }
}

template <typename M>
template <typename Fct>
void BaseFEM<M>::addLinear(const Fct &f, const itemVFlist_t &VF, const Interface<M> &gamma, const TimeSlab &In, int itq,
                           std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);
    auto tq        = this->get_quadrature_time(itq);
    double tid     = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(f, VF, gamma, iface, tid, &In, 1., itq);
        }
        this->addLocalContribution();
    }
}

template <typename M>
void BaseFEM<M>::addLinear(const itemVFlist_t &VF, const TimeInterface<M> &gamma, const TimeSlab &In,
                           std::list<int> label) {
    assert(VF.isRHS());

    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addLinear(VF, gamma, In, itq, label);
    }
}

template <typename M>
void BaseFEM<M>::addLinear(const itemVFlist_t &VF, const TimeInterface<M> &gamma, const TimeSlab &In, int itq,
                           std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);
    auto tq        = this->get_quadrature_time(itq);
    double tid     = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time = tq.a * In.get_measure();

    for (int iface = gamma[itq]->first_element(); iface < gamma[itq]->last_element();
         iface += gamma[itq]->next_element()) {
        const typename Interface<M>::Face &face = (*gamma[itq])[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, *gamma[itq], iface, tid, &In, cst_time, itq);
        }
    }
}

template <typename M>
template <typename Fct>
void BaseFEM<M>::addLinear(const Fct &f, const itemVFlist_t &VF, const TimeInterface<M> &gamma, const TimeSlab &In,
                           std::list<int> label) {
    assert(VF.isRHS());

    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addLinear(f, VF, gamma, In, itq, label);
    }
}

template <typename M>
template <typename Fct>
void BaseFEM<M>::addLinear(const Fct &f, const itemVFlist_t &VF, const TimeInterface<M> &gamma, const TimeSlab &In,
                           int itq, std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);
    auto tq        = this->get_quadrature_time(itq);
    double tid     = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions
    double cst_time = tq.a * In.get_measure();

    for (int iface = gamma[itq]->first_element(); iface < gamma[itq]->last_element();
         iface += gamma[itq]->next_element()) {
        const typename Interface<M>::Face &face = (*gamma[itq])[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(f, VF, *gamma[itq], iface, tid, &In, cst_time, itq);
        }
    }
}

template <typename M>
void BaseFEM<M>::addInterfaceContribution(const itemVFlist_t &VF, const Interface<M> &interface, int ifac, double tid,
                                          const TimeSlab *In, double cst_time, int itq) {
    typedef typename FElement::RdHatBord RdHatBord;

    // GET IDX ELEMENT CONTAINING FACE ON backMes
    const int kb = interface.idxElementOfFace(ifac);
    const Element &K(interface.get_element(kb));
    double measK = K.measure();
    double h     = K.get_h();
    double meas  = interface.measure(ifac);
    std::array<double, 2> measCut;

    { // compute each cut part
        const FESpace &Vh(VF.get_spaceV(0));
        const ActiveMesh<M> &Th(Vh.get_mesh());
        std::vector<int> idxV = Vh.idxAllElementFromBackMesh(kb, -1);

        const Cut_Part<Element> cutK(Th.get_cut_part(idxV[0], itq));
        measCut[0] = cutK.measure();
        measCut[1] = measCut[0];
        if (idxV.size() == 2) {
            const Cut_Part<Element> cutK(Th.get_cut_part(idxV[1], itq));
            measCut[1] = cutK.measure();
        }
    }

    const Rd normal(-interface.normal(ifac));

    // GET THE QUADRATURE RULE
    const QFB &qfb(this->get_quadrature_formular_cutFace());

    for (int l = 0; l < VF.size(); ++l) {
        // if(!VF[l].on(domain)) continue;

        // FINITE ELEMENT SPACES && ELEMENTS
        const FESpace &Vhv(VF.get_spaceV(l));
        const FESpace &Vhu(VF.get_spaceU(l));
        bool same = (VF.isRHS() || (&Vhu == &Vhv));

        std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
        std::vector<int> idxU = (same) ? idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());

        int kv = VF[l].onWhatElementIsTestFunction(idxV);
        int ku = VF[l].onWhatElementIsTrialFunction(idxU);

        const FElement &FKu(Vhu[ku]);
        const FElement &FKv(Vhv[kv]);
        int domu = FKu.get_domain();
        int domv = FKv.get_domain();
        this->initIndex(FKu, FKv);

        // BF MEMORY MANAGEMENT -
        int lastop = getLastop(VF[l].du, VF[l].dv);
        RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
        RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
        What_d Fop = Fwhatd(lastop);

        // COMPUTE COEFFICIENT && NORMAL
        double coef = VF[l].computeCoefInterface(h, meas, measK, measCut, {domu, domv});
        coef *= VF[l].computeCoefFromNormal(normal);

        // // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
            
            typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
            const Rd mip     = interface.mapToPhysicalFace(ifac, (RdHatBord)ip);
            const Rd face_ip = K.mapToReferenceElement(mip);
            double Cint      = meas * ip.getWeight() * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
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
                    // std::cout << Cint << std::endl;
                    // std::cout << fv << std::endl;
                    this->addToRHS(VF[l], FKv, fv, Cint);
                } else
                    this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
            }
        }
    }
}

// ! Add addInterfaceContribution for function evaluation (rhs)
template <typename M>
template <typename Fct>
void BaseFEM<M>::addInterfaceContribution(const Fct &f, const itemVFlist_t &VF, const Interface<M> &interface, int ifac,
                                          double tid, const TimeSlab *In, double cst_time, int itq) {
    typedef typename FElement::RdHatBord RdHatBord;

    // GET IDX ELEMENT CONTAINING FACE ON backMes
    const int kb = interface.idxElementOfFace(ifac);
    const Element &K(interface.get_element(kb));
    double measK = K.measure();
    double h     = K.get_h();
    double meas  = interface.measure(ifac);
    std::array<double, 2> measCut;

    { // compute each cut part
        const FESpace &Vh(VF.get_spaceV(0));
        const ActiveMesh<M> &Th(Vh.get_mesh());
        std::vector<int> idxV = Vh.idxAllElementFromBackMesh(kb, -1);

        const Cut_Part<Element> cutK(Th.get_cut_part(idxV[0], itq));
        measCut[0] = cutK.measure();
        measCut[1] = measCut[0];
        if (idxV.size() == 2) {
            const Cut_Part<Element> cutK(Th.get_cut_part(idxV[1], itq));
            measCut[1] = cutK.measure();
        }
    }

    const Rd normal(-interface.normal(ifac));

    // GET THE QUADRATURE RULE
    const QFB &qfb(this->get_quadrature_formular_cutFace());

    for (int l = 0; l < VF.size(); ++l) {
        // if(!VF[l].on(domain)) continue;

        // FINITE ELEMENT SPACES && ELEMENTS
        const FESpace &Vhv(VF.get_spaceV(l));
        const FESpace &Vhu(VF.get_spaceU(l));
        bool same = (VF.isRHS() || (&Vhu == &Vhv));

        std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
        std::vector<int> idxU = (same) ? idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());

        int kv = VF[l].onWhatElementIsTestFunction(idxV);
        int ku = VF[l].onWhatElementIsTrialFunction(idxU);

        const FElement &FKu(Vhu[ku]);
        const FElement &FKv(Vhv[kv]);
        int domu = FKu.get_domain();
        int domv = FKv.get_domain();
        this->initIndex(FKu, FKv);

        // BF MEMORY MANAGEMENT -
        int lastop = getLastop(VF[l].du, VF[l].dv);
        RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
        RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
        What_d Fop = Fwhatd(lastop);

        // COMPUTE COEFFICIENT && NORMAL
        double coef = VF[l].computeCoefInterface(h, meas, measK, measCut, {domu, domv});
        coef *= VF[l].computeCoefFromNormal(normal);

        // // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {

            typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
            Rd mip           = interface.mapToPhysicalFace(ifac, (RdHatBord)ip);
            const Rd face_ip = K.mapToReferenceElement(mip);
            double Cint      = meas * ip.getWeight() * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, face_ip, fv);
            if (!same)
                FKu.BF(Fop, face_ip, fu);

            Cint *= f(mip, VF[l].cv, tid);
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
                    // std::cout << Cint << std::endl;
                    // std::cout << fv << std::endl;
                    this->addToRHS(VF[l], FKv, fv, Cint);
                } else
                    this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
            }
        }
    }
}

// WITH MAPPING
template <typename M>
void BaseFEM<M>::addBilinear(const itemVFlist_t &VF, const Interface<M> &gamma, const Mapping<M> &mapping,
                             std::list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);
    progress bar(" Add Bilinear Interface", gamma.last_element(), globalVariable::verbose);

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        bar += gamma.next_element();
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, gamma, iface, 0., nullptr, 1., 0, mapping);
            this->addLocalContribution();
        }
    }
    bar.end();
}
template <typename M>
void BaseFEM<M>::addLinear(const itemVFlist_t &VF, const Interface<M> &gamma, const Mapping<M> &mapping,
                           std::list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);
    progress bar(" Add Linear Interface", gamma.last_element(), globalVariable::verbose);

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        bar += gamma.next_element();

        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, gamma, iface, 0., nullptr, 1., 0, mapping);
        }
    }

    bar.end();
}

static R2 operator*(const KNM_<double> &mat, const R2 v) {
    R2 val;
    val.x = mat(0, 0) * v.x + mat(0, 1) * v.y;
    val.y = mat(1, 0) * v.x + mat(1, 1) * v.y;
    return val;
}

static R determinant(const KNM_<double> &a) {
    R sum = 0.;
    if (a.N() == 2 && a.M() == 2) {
        sum = a(0, 0) * a(1, 1) - a(1, 0) * a(0, 1);
    } else if (a.N() == 3 && a.M() == 3) {
        sum = a(0, 0) * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) - a(0, 1) * (a(1, 0) * a(2, 2) - a(2, 0) * a(1, 2)) +
              a(0, 2) * (a(1, 0) * a(2, 1) - a(2, 0) * a(1, 1));
    } else
        assert(0);

    assert(std::fabs(sum) > 1e-16);
    return std::fabs(sum);
}

template <typename M>
void BaseFEM<M>::addInterfaceContribution(const itemVFlist_t &VF, const Interface<M> &interface, int ifac, double tid,
                                          const TimeSlab *In, double cst_time, int itq, const Mapping<M> &mapping) {
    typedef typename FElement::RdHatBord RdHatBord;

    // GET IDX ELEMENT CONTAINING FACE ON backMes
    const int kb = interface.idxElementOfFace(ifac);
    const Element &K(interface.get_element(kb));
    double measK = K.measure();
    double h     = K.get_h();
    double meas  = interface.measure(ifac);
    std::array<double, 2> measCut;
    KNM<double> invJ(Rd::d, Rd::d);

    { // compute each cut part
        const FESpace &Vh(VF.get_spaceV(0));
        const ActiveMesh<M> &Th(Vh.get_mesh());
        std::vector<int> idxV = Vh.idxAllElementFromBackMesh(kb, -1);

        const Cut_Part<Element> cutK(Th.get_cut_part(idxV[0], itq));
        measCut[0] = cutK.measure();
        measCut[1] = measCut[0];
        if (idxV.size() == 2) {
            const Cut_Part<Element> cutK(Th.get_cut_part(idxV[1], itq));
            measCut[1] = cutK.measure();
        }
    }

    const Rd linear_normal(-interface.normal(ifac));

    // GET THE QUADRATURE RULE
    const QFB &qfb(this->get_quadrature_formular_cutFace());

    for (int l = 0; l < VF.size(); ++l) {
        // if(!VF[l].on(domain)) continue;

        // FINITE ELEMENT SPACES && ELEMENTS
        const FESpace &Vhv(VF.get_spaceV(l));
        const FESpace &Vhu(VF.get_spaceU(l));
        bool same = (VF.isRHS() || (&Vhu == &Vhv));

        std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
        std::vector<int> idxU = (same) ? idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());

        int kv = VF[l].onWhatElementIsTestFunction(idxV);
        int ku = VF[l].onWhatElementIsTrialFunction(idxU);

        const FElement &FKu(Vhu[ku]);
        const FElement &FKv(Vhv[kv]);
        int domu = FKu.get_domain();
        int domv = FKv.get_domain();
        this->initIndex(FKu, FKv);

        // BF MEMORY MANAGEMENT -
        int lastop = getLastop(VF[l].du, VF[l].dv);
        RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
        RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
        What_d Fop = Fwhatd(lastop);

        // COMPUTE COEFFICIENT && NORMAL
        double coef = VF[l].computeCoefInterface(h, meas, measK, measCut, {domu, domv});
        // coef *= VF[l].computeCoefFromNormal(linear_normal);

        // // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {

            typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
            const Rd mip     = interface.mapToPhysicalFace(ifac, (RdHatBord)ip);
            const Rd face_ip = K.mapToReferenceElement(mip);
            double Cint      = meas * ip.getWeight() * cst_time;

            const Rd map_mip = mapping.map(kb, mip);

            mapping.computeInverseJacobian(kb, mip, invJ);
            Rd normal   = invJ * linear_normal;
            double DetJ = 1. / determinant(invJ);

            Cint *= coef * VF[l].c * DetJ * normal.norm();
            normal = normal / normal.norm();

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, face_ip, fv);
            if (!same)
                FKu.BF(Fop, face_ip, fu);

            mapping.transform(FKu, fu, invJ);
            if (!same)
                mapping.transform(FKv, fv, invJ);

            Cint *= VF[l].computeCoefFromNormal(normal);
            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domu, domv), map_mip,
                                                           tid, normal);

            if (In) {
                if (VF.isRHS())
                    this->addToRHS(VF[l], *In, FKv, fv, Cint);
                else
                    this->addToMatrix(VF[l], *In, FKu, FKv, fu, fv, Cint);
            } else {
                if (VF.isRHS()) {
                    // std::cout << Cint << std::endl;
                    // std::cout << fv << std::endl;
                    this->addToRHS(VF[l], FKv, fv, Cint);
                } else
                    this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
            }
        }
    }
}

// LAGRANGE MULTIPLIER
template <typename Mesh> void BaseFEM<Mesh>::addLagrangeMultiplier(const itemVFlist_t &VF, double val, const Mesh &Th) {
    assert(VF.isRHS());
    int ndf = this->rhs_.size();
    this->rhs_.resize(ndf + 1);
    this->rhs_(ndf) = val;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        BaseFEM<Mesh>::addLagrangeContribution(VF, k, nullptr, 0, 1);
        this->addLocalContributionLagrange(ndf);
    }
}
// ADD LAGRANGE contribution

template <typename M>
void BaseFEM<M>::addLagrangeContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq,
                                         double cst_time) {

    // Compute parameter coonected to the mesh.
    // on can take the one from the first test function
    const FESpace &Vh(*VF[0].fespaceV);
    const FElement &FK(Vh[k]);
    const Element &K(FK.T);
    double meas = K.measure();
    double h    = K.get_h();
    int domain  = FK.get_domain();
    int kb      = Vh.idxElementInBackMesh(k);

    // GET THE QUADRATURE RULE
    const QF &qf(this->get_quadrature_formular_K());
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

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
        double coef = VF[l].computeCoefElement(h, meas, meas, meas, domain);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
            typename QF::QuadraturePoint ip(qf[ipq]);
            const Rd mip = K.map(ip);
            double Cint  = meas * ip.getWeight() * cst_time;
            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, ip, fv);
            //   VF[l].applyFunNL(fu,fv);

            // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
            // Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, 0.);
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

template <typename M>
void BaseFEM<M>::addLagrangeMultiplier(const itemVFlist_t &VF, double val, const Interface<M> &gamma) {
    assert(VF.isRHS());

    int ndf = this->rhs_.size();
    this->rhs_.resize(ndf + 1);
    this->rhs_(ndf) = val;
    // this->nb_dof_ += 1; //? added but doesn't make a difference

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {

        addLagrangeContribution(VF, gamma, iface);

        this->addLocalContributionLagrange(ndf);
    }
}

// // First attempt
// template <typename M>
// void BaseFEM<M>::addLagrangeContribution(const itemVFlist_t &VF, const Interface<M> &interface, const int ifac) {

//     typedef typename FElement::RdHatBord RdHatBord;

//     // GET IDX ELEMENT CONTAINING FACE ON backMes
//     const int kb = interface.idxElementOfFace(ifac);
//     const Element &K(interface.get_element(kb));
//     double measK = K.measure();
//     double h     = K.get_h();
//     double meas  = interface.measure(ifac);
//     std::array<double, 2> measCut;
//     // KNM<double> invJ(Rd::d, Rd::d);

//     // { // compute each cut part
//     //     const FESpace &Vh(VF.get_spaceV(0));
//     //     const ActiveMesh<M> &Th(Vh.get_mesh());
//     //     std::vector<int> idxV = Vh.idxAllElementFromBackMesh(kb, -1);

//     //     const Cut_Part<Element> cutK(Th.get_cut_part(idxV[0], 0));
//     //     measCut[0] = cutK.measure();
//     //     measCut[1] = measCut[0];
//     //     if (idxV.size() == 2) {
//     //         const Cut_Part<Element> cutK(Th.get_cut_part(idxV[1], 0));
//     //         measCut[1] = cutK.measure();
//     //     }
//     // }

//     const Rd linear_normal(-interface.normal(ifac));

//     // GET THE QUADRATURE RULE
//     const QFB &qfb(this->get_quadrature_formular_cutFace());

//     for (int l = 0; l < VF.size(); ++l) {
//         // if(!VF[l].on(domain)) continue;

//         // FINITE ELEMENT SPACES && ELEMENTS
//         const FESpace &Vhv(VF.get_spaceV(l));

//         std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
//         int kv                = VF[l].onWhatElementIsTestFunction(idxV);

//         const FElement &FKv(Vhv[kv]);
//         int domv = FKv.get_domain();
//         this->initIndex(FKv, FKv);

//         // BF MEMORY MANAGEMENT -
//         //int lastop = getLastop(VF[l].du, VF[l].dv);
//         int lastop = getLastop(VF[l].dv, 0);

//         RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
//         What_d Fop = Fwhatd(lastop);

//         // COMPUTE COEFFICIENT && NORMAL
//         double coef = VF[l].computeCoefInterface(h, meas, measK, measCut, {domv, domv});
//         // coef *= VF[l].computeCoefFromNormal(linear_normal);

//         // // LOOP OVER QUADRATURE IN SPACE
//         for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {

//             typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
//             const Rd mip     = interface.mapToPhysicalFace(ifac, (RdHatBord)ip);
//             const Rd face_ip = K.mapToReferenceElement(mip);
//             double Cint      = meas * ip.getWeight();

//             //std::cout << Cint << "\n";

//             // const Rd map_mip = mapping.map(kb, mip);

//             // mapping.computeInverseJacobian(kb, mip, invJ);
//             // Rd normal   = invJ * linear_normal;

//             // double DetJ = 1. / determinant(invJ);

//             // Cint *= coef * VF[l].c * DetJ * normal.norm();
//             Cint *= coef * VF[l].c * linear_normal.norm();
//             // normal = normal / normal.norm();

//             // EVALUATE THE BASIS FUNCTIONS
//             FKv.BF(Fop, face_ip, fv);

//             // Cint *= VF[l].computeCoefFromNormal(normal);
//             Cint *= VF[l].computeCoefFromNormal(linear_normal);
//             // Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), domv, normal);
//             Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domv, domv), mip,
//             0.,
//                                                            linear_normal);

//             this->addToMatrix(VF[l], FKv, fv, Cint);
//         }
//     }
// }

template <typename M>
void BaseFEM<M>::addLagrangeContribution(const itemVFlist_t &VF, const Interface<M> &interface, const int ifac) {

    typedef typename FElement::RdHatBord RdHatBord;

    // GET IDX ELEMENT CONTAINING FACE ON backMes
    const int kb = interface.idxElementOfFace(ifac);
    const Element &K(interface.get_element(kb));
    double measK = K.measure();
    double h     = K.get_h();
    double meas  = interface.measure(ifac);
    std::array<double, 2> measCut;

    { // compute each cut part
        const FESpace &Vh(VF.get_spaceV(0));
        const ActiveMesh<M> &Th(Vh.get_mesh());
        std::vector<int> idxV = Vh.idxAllElementFromBackMesh(kb, -1);

        const Cut_Part<Element> cutK(Th.get_cut_part(idxV[0], 0));
        measCut[0] = cutK.measure();
        measCut[1] = measCut[0];
        if (idxV.size() == 2) {
            const Cut_Part<Element> cutK(Th.get_cut_part(idxV[1], 0));
            measCut[1] = cutK.measure();
        }
    }

    const Rd normal(-interface.normal(ifac));

    // GET THE QUADRATURE RULE
    const QFB &qfb(this->get_quadrature_formular_cutFace());

    for (int l = 0; l < VF.size(); ++l) {
        // if(!VF[l].on(domain)) continue;

        // FINITE ELEMENT SPACES && ELEMENTS
        const FESpace &Vhv(VF.get_spaceV(l));
        const FESpace &Vhu(VF.get_spaceU(l));
        bool same = (VF.isRHS() || (&Vhu == &Vhv));

        std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
        std::vector<int> idxU = (same) ? idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());

        int kv = VF[l].onWhatElementIsTestFunction(idxV);
        int ku = VF[l].onWhatElementIsTrialFunction(idxU);

        const FElement &FKu(Vhu[ku]);
        const FElement &FKv(Vhv[kv]);
        int domu = FKu.get_domain();
        int domv = FKv.get_domain();
        this->initIndex(FKu, FKv);

        // BF MEMORY MANAGEMENT -
        int lastop = getLastop(VF[l].du, VF[l].dv);
        RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
        RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
        What_d Fop = Fwhatd(lastop);

        // COMPUTE COEFFICIENT && NORMAL
        double coef = VF[l].computeCoefInterface(h, meas, measK, measCut, {domu, domv});
        coef *= VF[l].computeCoefFromNormal(normal);

        // // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {

            typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
            const Rd mip     = interface.mapToPhysicalFace(ifac, (RdHatBord)ip);
            const Rd face_ip = K.mapToReferenceElement(mip);
            double Cint      = meas * ip.getWeight();

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, face_ip, fv);
            if (!same)
                FKu.BF(Fop, face_ip, fu);
            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domu, domv), mip, 0.,
                                                           normal);
            Cint *= coef * VF[l].c;

            this->addToMatrix(VF[l], FKv, fv, Cint);
        }
    }
}

template <typename M>
void BaseFEM<M>::addLagrangeBorderContribution(const itemVFlist_t &VF, const Element &K, const BorderElement &BE,
                                               int ifac, const TimeSlab *In, int itq, double cst_time) {

    typedef typename FElement::RdHatBord RdHatBord;

    // Compute parameter connected to the mesh.
    double measK = K.measure();
    double meas  = K.mesureBord(ifac);
    double h     = K.get_h();
    Rd normal    = K.N(ifac);

    // U and V HAS TO BE ON THE SAME MESH
    const FESpace &Vh(VF.get_spaceV(0));
    int kb                = Vh.Th(K);
    std::vector<int> idxK = Vh.idxAllElementFromBackMesh(kb, -1);
    assert(idxK.size() == 1);
    // if(idxK.size() != 1) return;
    int k = VF[0].onWhatElementIsTestFunction(idxK);

    // FINTE ELEMENT SPACES && ELEMENTS
    const FESpace &Vhv(VF.get_spaceV(0));
    const FElement &FKv(Vhv[k]);
    int domain = FKv.get_domain();
    this->initIndex(FKv, FKv);

    // GET THE QUADRATURE RULE
    const QFB &qfb(this->get_quadrature_formular_dK());
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
    for (int l = 0; l < VF.size(); ++l) {
        if (!VF[l].on(domain))
            continue;

        // BF MEMORY MANAGEMENT -
        int lastop = getLastop(VF[l].du, VF[l].dv);
        RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop);
        What_d Fop = Fwhatd(lastop);

        // COMPUTE COEFFICIENT
        double coef = VF[l].computeCoefElement(h, meas, measK, measK, domain);
        coef *= VF[l].computeCoefFromNormal(normal);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
            typename QFB::QuadraturePoint ip(qfb[ipq]);
            const Rd mip     = BE.mapToPhysicalElement((RdHatBord)ip);
            const Rd ip_edge = K.mapToReferenceElement(mip);
            double Cint      = meas * ip.getWeight() * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, ip_edge, fv);

            // Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid,
            // normal);
            Cint *= coef * VF[l].c;

            if (In) {
                // if(VF.isRHS()) this->addToRHS(   VF[l], *In, FKv, fv, Cint);
                // else
                this->addToMatrix(VF[l], *In, FKv, fv, Cint);
            } else {
                // if(VF.isRHS()) this->addToRHS(VF[l], FKv, fv, Cint);
                // else
                this->addToMatrix(VF[l], FKv, fv, Cint);
            }
        }
    }
}

///
