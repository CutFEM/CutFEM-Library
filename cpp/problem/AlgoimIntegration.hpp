
#ifndef ALGOIM_INTEGRATION_HPP
#define ALGOIM_INTEGRATION_HPP

template <typename Mesh, typename L>
double L2_norm_surface_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(const R2, int i),
                    const Interface<Mesh> &interface, L& phi) {
    
    using fespace_t     = GFESpace<Mesh>;
    using FElement      = typename fespace_t::FElement;
    using Element       = typename Mesh::Element;
    using Rd            = typename FElement::Rd;

    // typedef typename FElement::QFB QFB;
    // typedef typename QFB::QuadraturePoint QuadraturePoint;

    // const QFB &qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
    What_d Fop = Fwhatd(op_id);

    double val = 0.;
    int quadrature_order = 5;

    for (int iface = interface.first_element(); iface < interface.last_element(); iface +=
    interface.next_element()) {

        const int kb = interface.idxElementOfFace(iface); // idx on backMesh
        const Element &K(interface.get_element(kb));
        //const R meas = interface.measure(iface);

        const auto &V0(K.at(0)); // vertex 0
        const auto &V2(K.at(2)); // vertex 2   diagonally opposed

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        algoim::QuadratureRule<2> q = algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight   = q.nodes.at(ipq).w;
            const R Cint = weight;

            double a = fh->evalOnBackMesh(kb, 0, mip) - fex(mip, fh->cu);
            val += Cint * a * a;
        }
    }
    double val_receive = 0;
#ifdef USE_MPI
    MPIcf::AllReduce(val, val_receive, MPI_SUM);
#else
    val_receive = val;
#endif

    return val_receive;
}


template <typename Mesh, typename L>
double L2_norm_surface(const FunFEM<Mesh> &fh, R(fex)(const R2, int i), const Interface<Mesh> &interface, L& phi, int c0,
                  int num_comp) {

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_surface_2(ui, fex, interface, phi);
    }
    return sqrt(val);
}


template <typename L, typename fct_t>
double integral_saye(fct_t &fh, const Interface<MeshQuad2> &interface, int cu, L phi) {

	using mesh_t        = MeshQuad2;
    using fespace_t     = GFESpace<mesh_t>;
    using itemVFlist_t  = ListItemVF<mesh_t>;
    using FElement      = typename fespace_t::FElement;
    using Rd            = typename FElement::Rd;
    using QF            = typename FElement::QF;
    using QFB           = typename FElement::QFB;
    using Element       = typename mesh_t::Element;
    using BorderElement = typename mesh_t::BorderElement;

	double val = 0.;
    int quadrature_order = 5;

    for (int iface = interface.first_element(); iface < interface.last_element();
            iface += interface.next_element()) {

    const int kb = interface.idxElementOfFace(iface); // idx on backMesh

    const auto &T(interface.get_element(kb));
    const auto &V0(T.at(0)); // vertex 0
    const auto &V2(T.at(2)); // vertex 2   diagonally opposed

    algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
    algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

    algoim::QuadratureRule<2> q =
        algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

    for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

        const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
        const R weight = q.nodes.at(ipq).w;

        const R Cint = weight;
        if constexpr (std::is_same_v<fct_t, FunFEM<MeshQuad2>>) {
            val += Cint * fh.evalOnBackMesh(kb, 0, mip, cu, 0);
        } else {
            val += Cint * fh->evalOnBackMesh(kb, 0, mip, nullptr);
        }
    }
    }
                
	double val_receive = 0;
	#ifdef USE_MPI
	MPIcf::AllReduce(val, val_receive, MPI_SUM);
	#else
	val_receive = val;
	#endif

	return val_receive;
}

#endif