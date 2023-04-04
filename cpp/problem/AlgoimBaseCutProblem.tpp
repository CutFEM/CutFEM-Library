
template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addElementContribution(const itemVFlist_t &VF, const int k, const TimeSlab *In, int itq,
                                                    double cst_time) {

    // GET CUT AND COMPUTE PARAMETERS
    const fespace_t &Vh(VF.get_spaceV(0));
    const ActiveMesh<M> &Th(Vh.get_mesh());
    // const Cut_Part<Element> cutK(Th.get_cut_part(k, itq));
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
    // GET THE QUADRATURE RULE
    // const QF &qf(this->get_quadrature_formular_cutK());

    const auto &V0(K.at(0)); // vertex 0
    const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

    algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
    algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

    // Get quadrature rule for the intersection between the element K and the negative part of the level set function
    algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, quadrature_order);

    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    phi.t = tid; // update time in level set function

    // LOOP OVER ELEMENTS IN THE CUT
    // for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

    // double meas_cut = cutK.measure(it);

    // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
    for (int l = 0; l < VF.size(); ++l) {
        if (!VF[l].on(domain))
            continue;

        // FINTE ELEMENT SPACES && ELEMENTS
        const fespace_t &Vhv(VF.get_spaceV(l));
        const fespace_t &Vhu(VF.get_spaceU(l));
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
        // double coef = VF[l].computeCoefElement(h, meas_cut, meas, meas_cut, domain);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            // typename QF::QuadraturePoint ip(qf[ipq]);
            // Rd mip      = cutK.mapToPhysicalElement(it, ip); // to the physical cut part
            // Rd cut_ip   = K.mapToReferenceElement(mip);      // back to the cut part in reference element

            const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            Rd cut_ip = K.mapToReferenceElement(mip); // map the quadrature points in the cut part to reference element
            const R weight = q.nodes.at(ipq).w;
            // double Cint = meas_cut * ip.getWeight() * cst_time;
            double Cint    = weight * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, cut_ip, fv);
            if (!same)
                FKu.BF(Fop, cut_ip, fu);
            //   VF[l].applyFunNL(fu,fv);

            // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);
            // Cint *= coef * VF[l].c;
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
    //}
}

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addInterfaceContribution(const itemVFlist_t &VF, const Interface<M> &interface, int ifac,
                                                      double tid, const TimeSlab *In, double cst_time, int itq) {
    typedef typename FElement::RdHatBord RdHatBord;

    phi.t = tid; // update time in level set function

    // GET IDX ELEMENT CONTAINING FACE ON backMes
    const int kb = interface.idxElementOfFace(ifac);
    const Element &K(interface.get_element(kb));

    const Rd normal(-interface.normal(ifac));

    const auto &V0(K.at(0)); // vertex 0
    const auto &V2(K.at(2)); // vertex 2   diagonally opposed

    algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
    algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

    algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

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

        double coef = VF[l].computeCoefFromNormal(normal);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            // typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
            // const Rd mip     = interface.mapToPhysicalFace(ifac, (RdHatBord)ip);

            const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight   = q.nodes.at(ipq).w;
            const Rd face_ip = K.mapToReferenceElement(mip);
            // double Cint      = meas * ip.getWeight() * cst_time;
            double Cint      = weight * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, face_ip, fv);
            if (!same)
                FKu.BF(Fop, face_ip, fu);
            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domu, domv), mip, tid,
                                                           normal);
            // Cint *= coef * VF[l].c;
            Cint *= VF[l].c;

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


template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addBilinear(const itemVFlist_t &VF, const Interface<M> &gamma) {
    assert(!VF.isRHS());
    progress bar(" Add Bilinear Interface", gamma.last_element(), globalVariable::verbose);

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        bar += gamma.next_element();
        
        addInterfaceContribution(VF, gamma, iface, 0., nullptr, 1., 0);
        this->addLocalContribution();
    }
    bar.end();
}

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addLinear(const itemVFlist_t &VF, const Interface<M> &gamma) {
    assert(VF.isRHS());
    progress bar(" Add Linear Interface", gamma.last_element(), globalVariable::verbose);

    for (int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {
        bar += gamma.next_element();

        const typename Interface<M>::Face &face = gamma[iface]; // the face

        addInterfaceContribution(VF, gamma, iface, 0., nullptr, 1., 0);

    }

    bar.end();
}

// template <typename M, typename L>
// template <typename Fct>
// void AlgoimBaseCutFEM<M, L>::addLinear(const Fct &f, const itemVFlist_t &VF, const Interface<M> &gamma) {
//     // TODO
// }

// template <typename M, typename L>
// template <typename Fct>
// void AlgoimBaseCutFEM<M, L>::addInterfaceContribution(const Fct &f, const itemVFlist_t &VF, const Interface<M> &interface, int ifac) {
//     // TODO
// }

template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addLagrangeMultiplier(const itemVFlist_t &VF, double val, const Interface<mesh_t> &gamma) {
    assert(VF.isRHS());
    
    int ndf = this->rhs_.size();
    this->rhs_.resize(ndf + 1);
    this->rhs_(ndf) = val;

    for(int iface = gamma.first_element(); iface < gamma.last_element(); iface += gamma.next_element()) {

        AlgoimBaseCutFEM<M, L>::addLagrangeContribution(VF, gamma, iface);
        this->addLocalContributionLagrange(ndf);

    }
}


template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addLagrangeContribution(const itemVFlist_t& VF, const Interface<mesh_t> &interface, const int iface) {

    //typedef typename QFB::QuadraturePoint QuadraturePoint;
    typedef typename FElement::RdHatBord RdHatBord;

    //KNM<double> invJ(Rd::d, Rd::d);
    const fespace_t &Vh(VF.get_spaceV(0));
    const int kb = interface.idxElementOfFace(iface);   // idx element in background mesh

    const Element &K(interface.get_element(kb));//(interface.get_element(kb));
    //const FElement &FK(Vh[kb]);
    //int domain  = FK.get_domain();

    //const typename Interface::FaceIdx& face = interface[iface];  // the face
    //const R meas = interface.computeDx(face).norm();
    //const double h = meas;
    const Rd linear_normal(-interface.normal(iface));
    assert(fabs(linear_normal.norm()-1)<1e-14);
    int nend = this->rhs_.size()-1;

    // GET THE QUADRATURE RULE
    // const QFB &qfb(this->get_quadrature_formular_cutFace());

    const auto &V0(K.at(0)); // vertex 0
    const auto &V2(K.at(2)); // vertex 2   diagonally opposed

    algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
    algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

    algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

    for(int l = 0; l < VF.size(); ++l) {
        
        // Finite element spaces and elements
        const fespace_t &Vhv(VF.get_spaceV(l));

        std::vector<int> idxV = Vhv.idxAllElementFromBackMesh(kb, VF[l].get_domain_test_function());
        int kv = VF[l].onWhatElementIsTestFunction(idxV);
        
        const FElement &FKv(Vhv[kv]);
        this->initIndex(FKv, FKv);

        // Basis function memory management
        int lastop = getLastop(VF[l].dv, 0);
        RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
        What_d Fop = Fwhatd(lastop);


        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {
            
            const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight   = q.nodes.at(ipq).w;
            // QuadraturePoint ip(qfb[ipq]); // integration point
            // const Rd mip = interface.mapToFace(face,(RdHatBord)ip);
            // mapping.computeInverseJacobian(km, mip, invJ);
            // Rd normal = invJ*linear_normal;
            // double DetJ = 1./determinant(invJ);
            // const R Cint = meas * ip.getWeight()*DetJ*normal.norm();
            const R Cint     = weight;

            FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element

            //mapping.transform(FKv, fv, invJ);

            this->addToMatrix(VF[l], FKv, fv, Cint);
            // for(int j = FKv.dfcbegin(VF[l].cv); j < FKv.dfcend(VF[l].cv); ++j) {
            //     (*this)(nend, FKv.loc2glb(j)) +=  Cint *  VF[l].c *fv(j,VF[l].cv,VF[l].dv);
            //     (*this)(FKv.loc2glb(j), nend) +=  Cint *  VF[l].c *fv(j,VF[l].cv,VF[l].dv);
                

            // }
        }
        //this->resetIndex();
    }
}
