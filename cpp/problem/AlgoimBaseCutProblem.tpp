
template <typename M, typename L>
void AlgoimBaseCutFEM<M, L>::addElementContribution(const ListItemVF<Rd::d> &VF, const int k, const TimeSlab *In, int itq,
                                           double cst_time) {

    // GET CUT AND COMPUTE PARAMETERS
    const FESpace &Vh(VF.get_spaceV(0));
    const ActiveMesh<M> &Th(Vh.get_mesh());
    //const Cut_Part<Element> cutK(Th.get_cut_part(k, itq));
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
    //const QF &qf(this->get_quadrature_formular_cutK());

	const auto &V0(K.at(0));	// vertex 0
	const auto &V2(K.at(2));	// vertex 2 (diagonally opposed)

	algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
	algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

	// Get quadrature rule for the intersection between the element K and the negative part of the level set function
	algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, 5);

    auto tq    = this->get_quadrature_time(itq);
    double tid = (In) ? (double)In->map(tq) : 0.;

	phi.t = tid;	// update time in level set function

    // LOOP OVER ELEMENTS IN THE CUT
    //for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

        //double meas_cut = cutK.measure(it);
        
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
		//double coef = VF[l].computeCoefElement(h, meas_cut, meas, meas_cut, domain);


		// LOOP OVER QUADRATURE IN SPACE
		for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

			// typename QF::QuadraturePoint ip(qf[ipq]);
			// Rd mip      = cutK.mapToPhysicalElement(it, ip); // to the physical cut part
			// Rd cut_ip   = K.mapToReferenceElement(mip);      // back to the cut part in reference element

			const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
			Rd cut_ip   = K.mapToReferenceElement(mip);      // map the quadrature points in the cut part to reference element
			const R weight = q.nodes.at(ipq).w;
			//double Cint = meas_cut * ip.getWeight() * cst_time;
			double Cint = weight * cst_time;

			// EVALUATE THE BASIS FUNCTIONS
			FKv.BF(Fop, cut_ip, fv);
			if (!same)
				FKu.BF(Fop, cut_ip, fu);
			//   VF[l].applyFunNL(fu,fv);

			// FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
			Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);
			//Cint *= coef * VF[l].c;
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
void AlgoimBaseCutFEM<M, L>::addInterfaceContribution(const ListItemVF<Rd::d> &VF, const Interface<M> &interface, int ifac,
                                          double tid, const TimeSlab *In, double cst_time, int itq) {
    typedef typename FElement::RdHatBord RdHatBord;

	phi.t = tid;	// update time in level set function

    // GET IDX ELEMENT CONTAINING FACE ON backMes
    const int kb = interface.idxElementOfFace(ifac);
    const Element &K(interface.get_element(kb));
    // double measK = K.measure();
    // double h     = K.get_h();
    // double meas  = interface.measure(ifac);
    //std::array<double, 2> measCut;

    // { // compute each cut part
    //     const FESpace &Vh(VF.get_spaceV(0));
    //     const ActiveMesh<M> &Th(Vh.get_mesh());
    //     std::vector<int> idxV = Vh.idxAllElementFromBackMesh(kb, -1);

    //     const Cut_Part<Element> cutK(Th.get_cut_part(idxV[0], itq));
    //     measCut[0] = cutK.measure();
    //     measCut[1] = measCut[0];
    //     if (idxV.size() == 2) {
    //         const Cut_Part<Element> cutK(Th.get_cut_part(idxV[1], itq));
    //         measCut[1] = cutK.measure();
    //     }
    // }

    const Rd normal(-interface.normal(ifac));

    // GET THE QUADRATURE RULE
    //const QFB &qfb(this->get_quadrature_formular_cutFace());

	const auto &V0(K.at(0)); // vertex 0
	const auto &V2(K.at(2)); // vertex 2   diagonally opposed

	algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
	algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

	algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, 1);

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
        //double coef = VF[l].computeCoefInterface(h, meas, measK, measCut, {domu, domv});
        //coef *= VF[l].computeCoefFromNormal(normal);
		
		double coef = VF[l].computeCoefFromNormal(normal);

        // LOOP OVER QUADRATURE IN SPACE
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {
			
            // typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
            // const Rd mip     = interface.mapToPhysicalFace(ifac, (RdHatBord)ip);

			const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
			const R weight = q.nodes.at(ipq).w;
            const Rd face_ip = K.mapToReferenceElement(mip);
            //double Cint      = meas * ip.getWeight() * cst_time;
			double Cint      = weight * cst_time;

            // EVALUATE THE BASIS FUNCTIONS
            FKv.BF(Fop, face_ip, fv);
            if (!same)
                FKu.BF(Fop, face_ip, fu);
            Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kb, kb), std::make_pair(domu, domv), mip, tid,
                                                           normal);
            //Cint *= coef * VF[l].c;
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
