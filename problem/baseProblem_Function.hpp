template <typename M>
void BaseFEM<M>::addToMatrix(const ItemVF<Rd::d> &VFi, const FElement &FKu, const FElement &FKv,
                             const RNMK_ &fu, const RNMK_ &fv, double Cint) {

    for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
        for (int j = FKu.dfcbegin(VFi.cu); j < FKu.dfcend(VFi.cu); ++j) {
            this->addToLocalContribution(FKv.loc2glb(i), FKu.loc2glb(j)) +=
                Cint * fv(i, VFi.cv, VFi.dv) * fu(j, VFi.cu, VFi.du);
        }
    }
}
template <typename M>
void BaseFEM<M>::addToMatrix(const ItemVF<Rd::d> &VFi, const TimeSlab &In, const FElement &FKu,
                             const FElement &FKv, const RNMK_ &fu, const RNMK_ &fv, double Cint) {
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    for (int it = In.dfcbegin(0); it < In.dfcend(0); ++it) {
        for (int jt = In.dfcbegin(0); jt < In.dfcend(0); ++jt) {
            const R tval = bf_time(it, 0, VFi.dtv) * bf_time(jt, 0, VFi.dtu);
            for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
                for (int j = FKu.dfcbegin(VFi.cu); j < FKu.dfcend(VFi.cu); ++j) {
                    this->addToLocalContribution(FKv.loc2glb(i, it), FKu.loc2glb(j, jt)) +=
                        Cint * tval * fv(i, VFi.cv, VFi.dv) * fu(j, VFi.cu, VFi.du);
                }
            }
        }
    }
}
template <typename M>
void BaseFEM<M>::addToMatrix(const ItemVF<Rd::d> &VFi, const FElement &FKv, const RNMK_ &fv,
                             double Cint) {

    for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
        this->addToLocalContribution(FKv.loc2glb(i), 0) += Cint * fv(i, VFi.cv, VFi.dv);
    }
}
template <typename M>
void BaseFEM<M>::addToMatrix(const ItemVF<Rd::d> &VFi, const TimeSlab &In, const FElement &FKv,
                             const RNMK_ &fv, double Cint) {
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    for (int it = In.dfcbegin(0); it < In.dfcend(0); ++it) {
        for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {
            this->addToLocalContribution(FKv.loc2glb(i, it), 0) +=
                Cint * fv(i, VFi.cv, VFi.dv) * bf_time(it, 0, VFi.dtv);
        }
    }
}
template <typename M>
void BaseFEM<M>::addToMatrix_Opt(const ItemVF<Rd::d> &VFi, const FElement &FK, const RNMK_ &fv,
                                 double Cint) {

    for (int i = FK.dfcbegin(VFi.cv); i < FK.dfcend(VFi.cv); ++i) {
        for (int j = FK.dfcbegin(VFi.cu); j < FK.dfcend(VFi.cu); ++j) {
            this->addToLocalContribution_Opt(i, j) +=
                Cint * fv(i, VFi.cv, VFi.dv) * fv(j, VFi.cu, VFi.du);
        }
    }
}

template <typename M>
void BaseFEM<M>::addToRHS(const ItemVF<Rd::d> &VFi, const FElement &FKv, const RNMK_ &fv,
                          double Cint) {
    for (int i = FKv.dfcbegin(VFi.cv); i < FKv.dfcend(VFi.cv); ++i) {

        (*this)(FKv.loc2glb(i)) += Cint * fv(i, VFi.cv, VFi.dv);
    }
}

template <typename M>
void BaseFEM<M>::addToRHS(const ItemVF<Rd::d> &VFi, const TimeSlab &In, const FElement &FKv,
                          const RNMK_ &fv, double Cint) {
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
template <typename Mesh>
void BaseFEM<Mesh>::addBilinear(const ListItemVF<Rd::d> &VF, const Mesh &Th) {
    assert(!VF.isRHS());
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        BaseFEM<Mesh>::addElementContribution(VF, k, nullptr, 0, 1.);

        this->addLocalContribution();
    }
}
template <typename Mesh>
void BaseFEM<Mesh>::addLinear(const ListItemVF<Rd::d> &VF, const Mesh &Th) {
    assert(VF.isRHS());
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        BaseFEM<Mesh>::addElementContribution(VF, k, nullptr, 0, 1.);
    }
}
template <typename M>
void BaseFEM<M>::addElementContribution(const ListItemVF<Rd::d> &VF, const int k,
                                        const TimeSlab *In, int itq, double cst_time) {

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
                    this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
            }
        }
    }
}

template <typename M>
void BaseFEM<M>::addElementContribution_Opt(const ListItemVF<Rd::d> &VF, const int k,
                                            const TimeSlab *In, int itq, double cst_time) {

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

    // GET THE QUADRATURE RULE
    const QF &qf(this->get_quadrature_formular_K());
    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

    // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
    this->initIndex(FK, FK);
    // BF MEMORY MANAGEMENT -
    int lastop = VF.get_lastOp();
    RNMK_ fv(this->databf_, FK.NbDoF(), FK.N,
             lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);
    this->loc_mat.init(FK.NbDoF(), FK.NbDoF());
    this->loc_mat = 0.;
    for (int l = 0; l < VF.size(); ++l) {

        // if(!VF[l].on(domain)) continue;
        // BF MEMORY MANAGEMENT -
        // int lastop = getLastop(VF[l].du, VF[l].dv);
        // RNMK_ fv(this->databf_,FK.NbDoF(),FK.N,lastop); //  the value for
        // basic fonction What_d Fop = Fwhatd(lastop);

        // COMPUTE COEFFICIENT
        double coef = VF[l].computeCoefElement(h, meas, meas, meas, domain);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
            typename QF::QuadraturePoint ip(qf[ipq]);
            const Rd mip = K.map(ip);
            double Cint  = meas * ip.getWeight() * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
            FK.BF(Fop, ip, fv);

            // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);

            Cint *= coef * VF[l].c;

            if (In) {
                if (VF.isRHS())
                    this->addToRHS(VF[l], *In, FK, fv, Cint);
                else
                    this->addToMatrix(VF[l], *In, FK, FK, fv, fv, Cint);
            } else {
                if (VF.isRHS())
                    this->addToRHS(VF[l], FK, fv, Cint);
                else
                    this->addToMatrix_Opt(VF[l], FK, fv, Cint);
            }
        }
    }
}

// INTEGRATION ON INNER FACE
template <typename Mesh>
void BaseFEM<Mesh>::addBilinear(const ListItemVF<Rd::d> &VF, const Mesh &Th, const CFacet &b) {
    assert(!VF.isRHS());
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
        this->addLocalContribution();
    }
}
template <typename Mesh>
void BaseFEM<Mesh>::addLinear(const ListItemVF<Rd::d> &VF, const Mesh &Th, const CFacet &b) {
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
void BaseFEM<M>::addFaceContribution(const ListItemVF<Rd::d> &VF, const std::pair<int, int> &e1,
                                     const std::pair<int, int> &e2, const TimeSlab *In, int itq,
                                     double cst_time) {

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
        RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N,
                 lastop);
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

            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(
                std::make_pair(kbu, kbv), std::make_pair(domain, domain), mip, tid, normal);
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

// INTEGRATION ON BOUNDARY
template <typename Mesh>
void BaseFEM<Mesh>::addBilinear(const ListItemVF<Rd::d> &VF, const Mesh &Th, const CBorder &b,
                                list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);

    for (int idx_be = Th.first_boundary_element(); idx_be < Th.last_boundary_element();
         idx_be += Th.next_boundary_element()) {

        int ifac;
        const int kb = Th.BoundaryElement(idx_be, ifac);
        const Element &K(Th[kb]);
        const BorderElement &BE(Th.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {
            BaseFEM<Mesh>::addBorderContribution(VF, K, BE, ifac, nullptr, 0, 1.);
            this->addLocalContribution();
        }
    }
}
template <typename Mesh>
void BaseFEM<Mesh>::addLinear(const ListItemVF<Rd::d> &VF, const Mesh &Th, const CBorder &b,
                              list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);

    for (int idx_be = Th.first_boundary_element(); idx_be < Th.last_boundary_element();
         idx_be += Th.next_boundary_element()) {

        int ifac;
        const int kb = Th.BoundaryElement(idx_be, ifac);

        const Element &K(Th[kb]);
        const BorderElement &BE(Th.be(idx_be));
        if (util::contain(label, BE.lab) || all_label) {

            // CHECK IF IT IS A CUT EDGE
            BaseFEM<Mesh>::addBorderContribution(VF, K, BE, ifac, nullptr, 0, 1.);
        }
    }
}
template <typename M>
void BaseFEM<M>::addBorderContribution(const ListItemVF<Rd::d> &VF, const Element &K,
                                       const BorderElement &BE, int iiifac, const TimeSlab *In,
                                       int itq, double cst_time) {

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
        RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N,
                 lastop);
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
void BaseFEM<Mesh>::setDirichlet(const FunFEM<Mesh> &, const Mesh &, list<int> label) {

    // bool all_label = (label.size() == 0);
    // std::map<int, double> df2set;
    //
    // for(int idx_be=Th.first_boundary_element();
    // idx_be<Th.last_boundary_element(); idx_be+= Th.next_boundary_element()) {
    //
    //   int ifac;
    //   const int kb = Th.BoundaryElement(idx_be, ifac);
    //   const Element & K(Th[kb]);
    //   const BorderElement & BE(Th.be(idx_be));
    //   if(util::contain(label, BE.lab) || all_label) {
    //
    //     // df2rm.insert({dfu, val});
    //
    //   }
}

// INTEGRATION ON INTERFACE
// On Faces
template <typename M>
void BaseFEM<M>::addBilinear(const ListItemVF<Rd::d> &VF, const Interface<M> &gamma,
                             list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);

    for (int iface = gamma.first_element(); iface < gamma.last_element();
         iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, gamma, iface, 0., nullptr, 1.);
            this->addLocalContribution();
        }
    }
}

template <typename M>
void BaseFEM<M>::addBilinear(const ListItemVF<Rd::d> &VF, const Interface<M> &gamma,
                             const TimeSlab &In, int itq, list<int> label) {
    assert(!VF.isRHS());
    bool all_label = (label.size() == 0);
    auto tq        = this->get_quadrature_time(itq);
    double tid     = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions

    for (int iface = gamma.first_element(); iface < gamma.last_element();
         iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, gamma, iface, tid, &In, 1.);
            this->addLocalContribution();
        }
    }
}

template <typename M>
void BaseFEM<M>::addBilinear(const ListItemVF<Rd::d> &VF, const TimeInterface<M> &gamma,
                             const TimeSlab &In, list<int> label) {

    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addBilinear(VF, gamma, In, itq, label);
    }
}

template <typename M>
void BaseFEM<M>::addBilinear(const ListItemVF<Rd::d> &VF, const TimeInterface<M> &gamma,
                             const TimeSlab &In, int itq, list<int> label) {
    assert(!VF.isRHS());
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

            addInterfaceContribution(VF, *gamma[itq], iface, tid, &In, cst_time);
            this->addLocalContribution();
        }
    }
}

template <typename M>
void BaseFEM<M>::addLinear(const ListItemVF<Rd::d> &VF, const Interface<M> &gamma,
                           list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);

    for (int iface = gamma.first_element(); iface < gamma.last_element();
         iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, gamma, iface, 0., nullptr, 1.);
        }
    }
}

template <typename M>
void BaseFEM<M>::addLinear(const ListItemVF<Rd::d> &VF, const Interface<M> &gamma,
                           const TimeSlab &In, int itq, list<int> label) {
    assert(VF.isRHS());
    bool all_label = (label.size() == 0);
    auto tq        = this->get_quadrature_time(itq);
    double tid     = In.map(tq);

    KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
    RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
    In.BF(tq.x, bf_time); // compute time basic funtions

    for (int iface = gamma.first_element(); iface < gamma.last_element();
         iface += gamma.next_element()) {
        const typename Interface<M>::Face &face = gamma[iface]; // the face
        if (util::contain(label, face.lab) || all_label) {

            addInterfaceContribution(VF, gamma, iface, tid, &In, 1.);
        }
        this->addLocalContribution();
    }
}

template <typename M>
void BaseFEM<M>::addLinear(const ListItemVF<Rd::d> &VF, const TimeInterface<M> &gamma,
                           const TimeSlab &In, list<int> label) {
    assert(VF.isRHS());

    for (int itq = 0; itq < this->get_nb_quad_point_time(); ++itq) {
        addLinear(VF, gamma, In, itq, label);
    }
}

template <typename M>
void BaseFEM<M>::addLinear(const ListItemVF<Rd::d> &VF, const TimeInterface<M> &gamma,
                           const TimeSlab &In, int itq, list<int> label) {
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

            addInterfaceContribution(VF, *gamma[itq], iface, tid, &In, cst_time);
        }
    }
}

template <typename M>
void BaseFEM<M>::addInterfaceContribution(const ListItemVF<Rd::d> &VF,
                                          const Interface<M> &interface, int ifac, double tid,
                                          const TimeSlab *In, double cst_time) {
                                            
    typedef typename FElement::RdHatBord RdHatBord;

    // GET IDX ELEMENT CONTAINING FACE ON backMes
    const int kb = interface.idxElementOfFace(ifac);
    const Element &K(interface.get_element(kb));
    double measK = K.measure();
    double h     = K.get_h();
    double meas  = interface.measure(ifac);

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
        RNMK_ fu(this->databf_ + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N,
                 lastop);
        What_d Fop = Fwhatd(lastop);

        // COMPUTE COEFFICIENT && NORMAL
        double coef = VF[l].computeCoefInterface(h, meas, measK, measK, domu, domv);
        coef *= VF[l].computeCoefFromNormal(normal);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {

            typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
            const Rd mip     = interface.mapToPhysicalFace(ifac, (RdHatBord)ip);
            const Rd face_ip = K.mapToReferenceElement(mip);
            double Cint      = meas * ip.getWeight() * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, face_ip, fv);
            if (!same)
                FKu.BF(Fop, face_ip, fu);

            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(
                std::make_pair(kb, kb), std::make_pair(domu, domv), mip, tid, normal);
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

// LAGRANGE MULTIPLIER
template <typename Mesh>
void BaseFEM<Mesh>::addLagrangeMultiplier(const ListItemVF<Rd::d> &VF, double val, const Mesh &Th) {
    assert(VF.isRHS());
    int ndf = this->rhs_.size();
    this->rhs_.resize(ndf + 1);
    this->rhs_(ndf) = val;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        BaseFEM<Mesh>::addLagrangeContribution(VF, k);
        this->addLocalContributionLagrange(ndf);
    }
}
// ADD LAGRANGE contribution

template <typename M>
void BaseFEM<M>::addLagrangeContribution(const ListItemVF<Rd::d> &VF, const int k) {

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
            double Cint  = meas * ip.getWeight();
            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, ip, fv);
            //   VF[l].applyFunNL(fu,fv);

            // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
            // Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip,
            // 0.);
            Cint *= coef * VF[l].c;

            // Cint = coef;
            this->addToMatrix(VF[l], FKv, fv, Cint);
        }
    }
}

template <typename M>
void BaseFEM<M>::addLagrangeBorderContribution(const ListItemVF<Rd::d> &VF, const Element &K,
                                               const BorderElement &BE, int ifac,
                                               const TimeSlab *In, int itq, double cst_time) {

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

            // Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip,
            // tid, normal);
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
