// Cut integrals

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addElementContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq,
                                                    double cst_time) {

    // Get finite element space and element, active mesh, and mesh element
    const fespace_t &Vh(VF.get_spaceV(0));
    const ActiveMesh<M> &Th(Vh.get_mesh());
    const FElement &FK(Vh[k]);
    const Element &K(FK.T);

    // double meas = K.measure();
    // double h    = K.get_h();
    int domain = FK.get_domain();
    int kb     = Vh.idxElementInBackMesh(k);

#ifdef USE_OMP
    int iam = omp_get_thread_num();
#else
    int iam = 0;
#endif

    // Get coordinates of current quadrilateral
    const auto &V0(K.at(0)); // vertex 0
    const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

    algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
    algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

    // Get current time
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    phi.t = tid; // update time in level set function

    // Get quadrature rule for the intersection between the element K and the negative part of the level set function
    algoim::QuadratureRule<2> q =
        algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, quadrature_order);

    assert(q.nodes.size() != 0); // assert quadrature rule is not empty

    // Loop over the variational formulation items
    for (int l = 0; l < VF.size(); ++l) {
        if (!VF[l].on(domain))
            continue;

        // Finite element spaces and elements
        const fespace_t &Vhv(VF.get_spaceV(l));
        const fespace_t &Vhu(VF.get_spaceU(l));
        const FElement &FKv(Vhv[k]);
        const FElement &FKu(Vhu[k]);
        this->initIndex(FKu, FKv);

        // Basis functions memory management
        bool same  = (&Vhu == &Vhv);
        int lastop = getLastop(VF[l].du, VF[l].dv);

        long offset = iam * this->offset_bf_;
        RNMK_ fv(this->databf_ + offset, FKv.NbDoF(), FKv.N,
                 lastop); //  the value for basic function
        RNMK_ fu(this->databf_ + offset + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N,
                 lastop); //  the value for basic function
        What_d Fop = Fwhatd(lastop);

        // std::cout << "fu = " << fu << "\n";

        // Loop over quadrature in space
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            Rd cut_ip = K.mapToReferenceElement(mip); // map the quadrature points in the cut part to reference element
            const R weight = q.nodes.at(ipq).w;

            double Cint = weight * cst_time;

            // Evaluate the basis functions
            FKv.BF(Fop, cut_ip, fv);
            if (!same)
                FKu.BF(Fop, cut_ip, fu);

            // std::cout << "mip: " << mip
            //           << ",\t feval = " << VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid) << "\n";

            // Find and compute all the coefficients and parameters
            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);
            Cint *= VF[l].c;

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
    // getchar();
}

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addBilinearAlgoim(const itemVFlist_t &VF, const ActiveMesh<M> &Th) {
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
            // if (Th.isCut(k, 0)) {
            addElementContribution(VF, k, nullptr, 0, 1.);
            // } else {
            //     BaseFEM<M>::addElementContribution(VF, k, nullptr, 0, 1.);
            // }
            this->addLocalContribution();
        }
        bar.end();
    }
    // std::cout << " real time " << MPIcf::Wtime() - t0 << std::endl;
}

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addLinearAlgoim(const itemVFlist_t &VF, const ActiveMesh<M> &Th) {
    assert(VF.isRHS());
    progress bar(" Add Linear CutMesh", Th.last_element(), globalVariable::verbose);

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        bar += Th.next_element();

        // if (Th.isCut(k, 0)) {
        addElementContribution(VF, k, nullptr, 0, 1.);
        // } else {
        //     BaseFEM<M>::addElementContribution(VF, k, nullptr, 0, 1.);
        // }
        // if(Th.isCut(k, 0))  BaseCutFEM<M>::addElementContribution(VF,
        // k,nullptr, 0, 1.); else BaseFEM<M>::addElementContribution(VF,
        // k,nullptr, 0, 1.);
    }
    bar.end();
}

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addBilinearAlgoim(const itemVFlist_t &VF, const ActiveMesh<M> &Th, const TimeSlab &In) {

    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {

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
                // if (Th.isCut(k, itq))
                // If the element is cut, add contribution to local matrix using addElementContribution function
                // from BaseCutFEM class
                addElementContribution(VF, k, &In, itq, cst_time);
                // else
                //     // If the element is not cut, add contribution to local matrix using addElementContribution
                //     function
                //     // from BaseFEM class
                //     BaseFEM<M>::addElementContribution(VF, k, &In, itq, cst_time);

                // Add the contribution to the local matrix
                this->addLocalContribution();
            }
            bar.end();
        }
        // End parallel region
    }
}

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addLinearAlgoim(const itemVFlist_t &VF, const ActiveMesh<M> &Th, const TimeSlab &In) {
    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
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

            // if (Th.isCut(k, itq))
            addElementContribution(VF, k, &In, itq, cst_time);
            // else
            //     BaseFEM<M>::addElementContribution(VF, k, &In, itq, cst_time);
        }
    }
}

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addBilinearAlgoim(const itemVFlist_t &VF, const ActiveMesh<M> &Th, const int itq,
                                               const TimeSlab &In) {

    // Assert that the input is not a RHS
    assert(!VF.isRHS());

    // Get the quadrature time for iteration "itq"
    auto tq = this->get_quadrature_time(itq);

    // Calculate the time using the map function of the TimeSlab "In"
    double tid = In.map(tq);

    // Allocate memory for the time-dependent basis functions
    //! KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1); What does this do?
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
        // if (Th.isCut(k, itq))
        addElementContribution(VF, k, &In, itq, 1.);
        // Else, add its contribution using BaseFEM
        // else
        //     BaseFEM<M>::addElementContribution(VF, k, &In, itq, 1.);

        // Add the local contribution
        this->addLocalContribution();
    }

    // End the progress bar
    bar.end();
}

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addLinearAlgoim(const itemVFlist_t &VF, const ActiveMesh<M> &Th, int itq,
                                             const TimeSlab &In) {
    assert(VF.isRHS());
    auto tq    = this->get_quadrature_time(itq);
    double tid = In.map(tq);
    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time);
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (Th.isInactive(k, itq))
            continue;

        // if (Th.isCut(k, itq))
        addElementContribution(VF, k, &In, itq, 1.);
        // else
        //     BaseFEM<M>::addElementContribution(VF, k, &In, itq, 1.);

        this->addLocalContribution();
    }
}

// Interface integrals

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addInterfaceContribution(const itemVFlist_t &VF, const Interface<M> &interface, int ifac,
                                                      double tid, const TimeSlab *In, double cst_time, int itq) {
    typedef typename FElement::RdHatBord RdHatBord;

    phi.t = tid; // update time in level set function

    //  GET IDX ELEMENT CONTAINING FACE ON backMes
    const int kb = interface.idxElementOfFace(ifac);
    const Element &K(interface.get_element(kb));

    // const Rd normal(-interface.normal(ifac));

    const auto &V0(K.at(0)); // vertex 0
    const auto &V2(K.at(2)); // vertex 2   diagonally opposed

    algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
    algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

    algoim::QuadratureRule<2> q =
        algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

    assert(q.nodes.size() != 0); // assert quadrature rule not empty

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
        // double coef = VF[l].computeCoefInterface(h, meas, measK, measCut, {domu, domv});
        // coef *= VF[l].computeCoefFromNormal(normal);

        // double coef = VF[l].computeCoefFromNormal(normal);

        // LOOP OVER QUADRATURE IN SPACE

        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            // typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
            // const Rd mip     = interface.mapToPhysicalFace(ifac, (RdHatBord)ip);

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight   = q.nodes.at(ipq).w;
            const Rd face_ip = K.mapToReferenceElement(mip);
            double Cint      = weight * cst_time;

            const Rd normal(phi.normal(mip));

            assert(fabs(normal.norm() - 1) < 1e-14);
            double coef = VF[l].computeCoefFromNormal(normal);

            FKv.BF(Fop, face_ip, fv);

            if (!same)
                FKu.BF(Fop, face_ip, fu);

            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domu, domv), mip, tid,
                                                           normal);
            Cint *= coef * VF[l].c;

            // std::cout << "mip = " << mip << "\t" << "feval = "
            //           << VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domu, domv),
            //           mip,
            //                                                     tid, normal)
            //           << "\n";

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
    // getchar();
}

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addLagrangeContribution(const itemVFlist_t &VF, const Interface<mesh_t> &interface,
                                                     const int iface) {

    typedef typename FElement::RdHatBord RdHatBord;

    const fespace_t &Vh(VF.get_spaceV(0));

    const int kb = interface.idxElementOfFace(iface); // idx element in background mesh
    const Element &K(interface.get_element(kb));      //(interface.get_element(kb));
    // const FElement &FK(Vh[kb]);
    // int domain  = FK.get_domain();

    // const typename Interface::FaceIdx& face = interface[iface];  // the face
    // const R meas = interface.computeDx(face).norm();
    // const double h = meas;

    // const Rd linear_normal(-interface.normal(iface));
    // assert(fabs(linear_normal.norm() - 1) < 1e-14);

    int nend = this->rhs_.size() - 1;

    // GET THE QUADRATURE RULE
    // const QFB &qfb(this->get_quadrature_formular_cutFace());

    const auto &V0(K.at(0)); // vertex 0
    const auto &V2(K.at(2)); // vertex 2   diagonally opposed

    algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
    algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

    algoim::QuadratureRule<2> q =
        algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

    for (int l = 0; l < VF.size(); ++l) {

        // Finite element spaces and elements
        const fespace_t &Vhv(VF.get_spaceV(l));

        std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
        int kv                = VF[l].onWhatElementIsTestFunction(idxV);

        const FElement &FKv(Vhv[kv]);
        this->initIndex(FKv, FKv);

        // Basis function memory management
        int lastop = getLastop(VF[l].dv, 0);
        RNMK_ fv(this->databf_, FKv.NbDoF(), FKv.N, lastop); //  the value for basic fonction
        What_d Fop = Fwhatd(lastop);

        assert(q.nodes.size() != 0);
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const double weight = q.nodes.at(ipq).w;
            const Rd face_ip    = K.mapToReferenceElement(mip);
            double Cint         = weight;

            const Rd normal(phi.normal(mip));
            assert(fabs(normal.norm() - 1) < 1e-14);
            double coef = VF[l].computeCoefFromNormal(normal);
            // std::cout << VF[l].c << "\n";

            // QuadraturePoint ip(qfb[ipq]); // integration point
            // const Rd mip = interface.mapToFace(face,(RdHatBord)ip);
            // mapping.computeInverseJacobian(km, mip, invJ);
            // Rd normal = invJ*linear_normal;
            // double DetJ = 1./determinant(invJ);
            // const R Cint = meas * ip.getWeight()*DetJ*normal.norm();

            FKv.BF(Fop, face_ip, fv); // need point in local reference element

            // mapping.transform(FKv, fv, invJ);
            Cint *= coef * VF[l].c;

            this->addToMatrix(VF[l], FKv, fv, Cint);
            // for(int j = FKv.dfcbegin(VF[l].cv); j < FKv.dfcend(VF[l].cv); ++j) {
            //     (*this)(nend, FKv.loc2glb(j)) +=  Cint *  VF[l].c *fv(j,VF[l].cv,VF[l].dv);
            //     (*this)(FKv.loc2glb(j), nend) +=  Cint *  VF[l].c *fv(j,VF[l].cv,VF[l].dv);

            // }
        }

        // this->resetIndex();
    }
}

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addBilinearAlgoim(const itemVFlist_t &VF, const Interface<M> &gamma,
                                               std::list<int> label) {
    assert(!VF.isRHS());
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

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addLinearAlgoim(const itemVFlist_t &VF, const Interface<M> &gamma, std::list<int> label) {
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

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addBilinearAlgoim(const itemVFlist_t &VF, const AlgoimInterface<M, L> &gamma) {
    assert(!VF.isRHS());
    
    progress bar(" Add Bilinear Interface", gamma.last_element(), globalVariable::verbose);

    for (const auto &[kb, quadrature_rule] : gamma.get_cut_elements()) {

        assert(quadrature_rule.nodes.size() != 0);

        const Element &K(gamma.get_element(kb));

        for (int l = 0; l < VF.size(); ++l) {
            // FINITE ELEMENT SPACES && ELEMENTS
            const fespace_t &Vhv(VF.get_spaceV(l));
            const fespace_t &Vhu(VF.get_spaceU(l));
            bool same = (VF.isRHS() || (&Vhu == &Vhv));

            std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
            std::vector<int> idxU =
                (same) ? idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());

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

            for (int ipq = 0; ipq < quadrature_rule.nodes.size(); ++ipq) {

                Rd mip(quadrature_rule.nodes.at(ipq).x(0), quadrature_rule.nodes.at(ipq).x(1));
                const R weight   = quadrature_rule.nodes.at(ipq).w;
                const Rd face_ip = K.mapToReferenceElement(mip);
                double Cint      = weight;

                const Rd normal(phi.normal(mip));

                assert(fabs(normal.norm() - 1) < 1e-14);
                double coef = VF[l].computeCoefFromNormal(normal);

                FKv.BF(Fop, face_ip, fv);

                if (!same)
                    FKu.BF(Fop, face_ip, fu);

                Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domu, domv), mip,
                                                               0., normal);
                Cint *= coef * VF[l].c;

                this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
            }
        }
    
    }


    bar.end();
}

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addLinearAlgoim(const itemVFlist_t &VF, const AlgoimInterface<M, L> &gamma) {
    assert(VF.isRHS());
    
    progress bar(" Add Linear Interface", gamma.last_element(), globalVariable::verbose);

    for (const auto &[kb, quadrature_rule] : gamma.get_cut_elements()) {

        assert(quadrature_rule.nodes.size() != 0);

        const Element &K(gamma.get_element(kb));

        for (int l = 0; l < VF.size(); ++l) {
            // FINITE ELEMENT SPACES && ELEMENTS
            const fespace_t &Vhv(VF.get_spaceV(l));
            const fespace_t &Vhu(VF.get_spaceU(l));
            bool same = (VF.isRHS() || (&Vhu == &Vhv));

            std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
            std::vector<int> idxU =
                (same) ? idxV : Vhu.idxAllElementFromBackMesh(kb, VF[l].get_domain_trial_function());

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

            for (int ipq = 0; ipq < quadrature_rule.nodes.size(); ++ipq) {

                Rd mip(quadrature_rule.nodes.at(ipq).x(0), quadrature_rule.nodes.at(ipq).x(1));
                const R weight   = quadrature_rule.nodes.at(ipq).w;
                const Rd face_ip = K.mapToReferenceElement(mip);
                double Cint      = weight;

                const Rd normal(phi.normal(mip));

                assert(fabs(normal.norm() - 1) < 1e-14);
                double coef = VF[l].computeCoefFromNormal(normal);

                FKv.BF(Fop, face_ip, fv);

                if (!same)
                    FKu.BF(Fop, face_ip, fu);

                Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domu, domv), mip,
                                                               0., normal);
                Cint *= coef * VF[l].c;

                this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
            }
        }
    }


    bar.end();
}

// template <typename M, typename L>
// template <typename Fct>
// void AlgoimBaseCutFEM<M, L>::addInterfaceContribution(const Fct &f, const itemVFlist_t &VF, const Interface<M>
// &interface, int ifac,
//                                           double tid, const TimeSlab *In, double cst_time, int itq) {

// }