
template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addElementContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq,
                                                    double cst_time) {

    // Get finite element space and element, active mesh, and mesh element
    const fespace_t &Vh(VF.get_spaceV(0));
    const ActiveMesh<M> &Th(Vh.get_mesh());
    const FElement &FK(Vh[k]);
    const Element &K(FK.T);

    double meas = K.measure();
    double h    = K.get_h();
    int domain  = FK.get_domain();
    int kb      = Vh.idxElementInBackMesh(k);

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
}

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
    //getchar();
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

// template <typename M, typename L>
// template <typename Fct>
// void AlgoimBaseCutFEM<M, L>::addInterfaceContribution(const Fct &f, const itemVFlist_t &VF, const Interface<M>
// &interface, int ifac,
//                                           double tid, const TimeSlab *In, double cst_time, int itq) {

// }