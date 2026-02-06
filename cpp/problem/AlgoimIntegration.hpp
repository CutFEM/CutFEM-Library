
#ifndef ALGOIM_INTEGRATION_HPP
#define ALGOIM_INTEGRATION_HPP

#include "../common/AlgoimInterface.hpp"

#include "../common/AlgoimInterface.hpp"

const int quadrature_order_integration = 4;
const int bernstein_deg                = 2;

template <typename L>
double L2_norm_surface_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(const R2, int i),
                         const Interface<MeshQuad2> &interface, L &phi, const int order_space = quadrature_order_integration) {

    using fespace_t = GFESpace<MeshQuad2>;
    using FElement  = typename fespace_t::FElement;
    using Element   = typename MeshQuad2::Element;
    using Rd        = typename FElement::Rd;

    // typedef typename FElement::QFB QFB;
    // typedef typename QFB::QuadraturePoint QuadraturePoint;

    // const QFB &qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
    What_d Fop = Fwhatd(op_id);

    double val = 0.;

    for (int iface = interface.first_element(); iface < interface.last_element(); iface += interface.next_element()) {

        const int kb = interface.idxElementOfFace(iface); // idx on backMesh
        const Element &K(interface.get_element(kb));
        // const R meas = interface.measure(iface);

        const auto &V0(K.at(0)); // vertex 0
        const auto &V2(K.at(2)); // vertex 2   diagonally opposed

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        algoim::QuadratureRule<2> q =
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, order_space);

        // assert((q.nodes.size() == 1 * quadrature_order_integration) ||
        //    (q.nodes.size() == 2 * quadrature_order_integration));
        assert(q.nodes.size() != 0);
        assert((order_space <= q.nodes.size()) || (q.nodes.size() <= 2 * order_space));
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight = q.nodes.at(ipq).w;
            const R Cint   = weight;

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

template <typename L>
double L2_norm_surface_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(const R2, int i, double t),
                         const Interface<MeshQuad2> &interface, double tt, L &phi, const int order_space = quadrature_order_integration) {

    using fespace_t = GFESpace<MeshQuad2>;
    using FElement  = typename fespace_t::FElement;
    using Element   = typename MeshQuad2::Element;
    using Rd        = typename FElement::Rd;

    // typedef typename FElement::QFB QFB;
    // typedef typename QFB::QuadraturePoint QuadraturePoint;

    phi.t = tt;

    // const QFB &qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
    What_d Fop = Fwhatd(op_id);

    double val = 0.;

    for (int iface = interface.first_element(); iface < interface.last_element(); iface += interface.next_element()) {

        const int kb = interface.idxElementOfFace(iface); // idx on backMesh
        const Element &K(interface.get_element(kb));
        // const R meas = interface.measure(iface);

        const auto &V0(K.at(0)); // vertex 0
        const auto &V2(K.at(2)); // vertex 2   diagonally opposed

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        algoim::QuadratureRule<2> q =
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, order_space);

        // assert((q.nodes.size() == 1 * quadrature_order_integration) ||
        //    (q.nodes.size() == 2 * quadrature_order_integration));
        assert((order_space <= q.nodes.size()) || (q.nodes.size() <= 2 * order_space));
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight = q.nodes.at(ipq).w;
            const R Cint   = weight;

            double a = fh->evalOnBackMesh(kb, 0, mip, tt) - fex(mip, fh->cu, tt);
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

template <typename L>
double L2_norm_surface_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, double t),
                         const Interface<MeshQuad2> &interface, double tt, L &phi, const int order_space = quadrature_order_integration) {

    using fespace_t = GFESpace<MeshQuad2>;
    using FElement  = typename fespace_t::FElement;
    using Element   = typename MeshQuad2::Element;
    using Rd        = typename FElement::Rd;

    // typedef typename FElement::QFB QFB;
    // typedef typename QFB::QuadraturePoint QuadraturePoint;

    phi.t = tt;

    // const QFB &qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
    What_d Fop = Fwhatd(op_id);

    double val = 0.;

    for (int iface = interface.first_element(); iface < interface.last_element(); iface += interface.next_element()) {

        const int kb = interface.idxElementOfFace(iface); // idx on backMesh
        const Element &K(interface.get_element(kb));
        // const R meas = interface.measure(iface);

        const auto &V0(K.at(0)); // vertex 0
        const auto &V2(K.at(2)); // vertex 2   diagonally opposed

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        algoim::QuadratureRule<2> q =
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, order_space);

        assert(q.nodes.size() != 0);
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight = q.nodes.at(ipq).w;
            const R Cint   = weight;

            double a = fh->evalOnBackMesh(kb, 0, mip, tt) - fex(mip, fh->cu, tt);
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

template <meshTriag Mesh, typename L, typename FEX>
double L2_norm_surface_2(const std::shared_ptr<ExpressionVirtual> &fh, const FEX &fex, const Interface<Mesh> &interface, L &phi, const int order_space = quadrature_order_integration) {

    using fespace_t = GFESpace<Mesh>;
    using FElement  = typename fespace_t::FElement;
    using Element   = typename Mesh::Element;
    using Rd        = typename FElement::Rd;

    // typedef typename FElement::QFB QFB;
    // typedef typename QFB::QuadraturePoint QuadraturePoint;

    double val = 0.;

    for (int iface = interface.first_element(); iface < interface.last_element(); iface += interface.next_element()) {

        const int kb = interface.idxElementOfFace(iface); // idx on backMesh
        const Element &K(interface.get_element(kb));
        // const R meas = interface.measure(iface);

        using Vec2 = algoim::uvector<double,2>;
        std::array<Vec2,3> verts = {
            Vec2(K.at(0)[0], K.at(0)[1]),
            Vec2(K.at(1)[0], K.at(1)[1]),
            Vec2(K.at(2)[0], K.at(2)[1])
        };

        const auto q = cuttri_multipoly::compute_cut_triangle_quadrature(
            verts,
            phi,
            bernstein_deg,   // bernstein degree
            quadrature_order_integration    // algoim 1D quadrature order
        );

        // assert(q.x_surf.size() != 0);
        for (int ipq = 0; ipq < q.x_surf.size(); ++ipq) {

            Rd mip(q.x_surf.at(ipq)[0], q.x_surf.at(ipq)[1]);
            
            const R Cint   = q.w_surf.at(ipq);

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


template <meshTriag Mesh, typename L, typename FEX>
double L2_norm_surface_2(const std::shared_ptr<ExpressionVirtual> &fh, const FEX &fex, const Interface<Mesh> &interface, double tt, L &phi, const int order_space = quadrature_order_integration) {
    
    using fespace_t = GFESpace<Mesh>;
    using FElement  = typename fespace_t::FElement;
    using Element   = typename Mesh::Element;
    using Rd        = typename FElement::Rd;

    // typedef typename FElement::QFB QFB;
    // typedef typename QFB::QuadraturePoint QuadraturePoint;

    phi.t = tt;

    double val = 0.;

    for (int iface = interface.first_element(); iface < interface.last_element(); iface += interface.next_element()) {

        const int kb = interface.idxElementOfFace(iface); // idx on backMesh
        const Element &K(interface.get_element(kb));
        // const R meas = interface.measure(iface);

        using Vec2 = algoim::uvector<double,2>;
        std::array<Vec2,3> verts = {
            Vec2(K.at(0)[0], K.at(0)[1]),
            Vec2(K.at(1)[0], K.at(1)[1]),
            Vec2(K.at(2)[0], K.at(2)[1])
        };

        const auto q = cuttri_multipoly::compute_cut_triangle_quadrature(
            verts,
            phi,
            bernstein_deg,   // bernstein degree
            quadrature_order_integration    // algoim 1D quadrature order
        );

        // assert(q.x_surf.size() != 0);
        for (int ipq = 0; ipq < q.x_surf.size(); ++ipq) {

            Rd mip(q.x_surf.at(ipq)[0], q.x_surf.at(ipq)[1]);
            
            const R Cint   = q.w_surf.at(ipq);

            double a = fh->evalOnBackMesh(kb, 0, mip, tt) - fex(mip, fh->cu, tt);
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


template <typeMesh M, typename L, typename FEX>
double L2_norm_surface(const FunFEM<M> &fh, const FEX &fex, const Interface<M> &interface, L &phi, int c0, int num_comp, const int order_space = quadrature_order_integration) {
    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_surface_2(ui, fex, interface, phi, order_space);
    }
    return sqrt(val);
}

//! REMOVE BELOW IF THE ONE ABOVE WORKS
template <typeMesh M, typename L>
double L2_norm_surface(const FunFEM<M> &fh, R(fex)(const R2, int i), const Interface<M> &interface, L &phi,
                       int c0, int num_comp, const int order_space = quadrature_order_integration) {

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_surface_2(ui, fex, interface, phi, order_space);
    }
    return sqrt(val);
}

template <typeMesh M, typename L, typename FEX>
double L2_norm_surface(const FunFEM<M> &fh, const FEX &fex, const Interface<M> &interface, double tt, L &phi, int c0, int num_comp, const int order_space = quadrature_order_integration) {
    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_surface_2(ui, fex, interface, tt, phi, order_space);
    }
    return sqrt(val);
}

//! REMOVE BELOW IF THE ONE ABOVE WORKS
template <typeMesh M, typename L>
double L2_norm_surface(const FunFEM<M> &fh, R(fex)(const R2, int i, double t), const Interface<M> &interface,
                       double tt, L &phi, int c0, int num_comp, const int order_space = quadrature_order_integration) {

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_surface_2(ui, fex, interface, tt, phi, order_space);
    }
    return sqrt(val);
}

template <typeMesh M, typename L>
double L2_norm_surface(const FunFEM<M> &fh, R(fex)(double *, int i, double t), const Interface<M> &interface,
                       double tt, L &phi, int c0, int num_comp, const int order_space = quadrature_order_integration) {

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_surface_2(ui, fex, interface, tt, phi, order_space);
    }
    return sqrt(val);
}

// L2(In x Gamma)
template <typename L, typename fct_t>
double L2L2_norm_surf(const FunFEM<MeshQuad2> &fh, const fct_t &f, const TimeInterface<MeshQuad2> &gamma, const TimeSlab &In,
                 const QuadratureFormular1d &qTime, L &phi, const int order_space = quadrature_order_integration) {
    using mesh_t    = MeshQuad2;
    using fespace_t = GFESpace<mesh_t>;
    using FElement  = typename fespace_t::FElement;
    using Rd        = typename FElement::Rd;
    using Element   = typename mesh_t::Element;

    double val = 0.;

    const int domain = 0; // do only for main domain

    // Loop in time
    for (int itq = 0; itq < qTime.n; ++itq) {

        // Get quadrature points in time
        GQuadraturePoint<R1> tq((qTime)[itq]);
        const double t = In.mapToPhysicalElement(tq);
        phi.t          = t;

        double weight_time = In.T.measure() * tq.a;

        // Loop in space
        double weight_space = 0.;
        for (int iface = gamma[itq]->first_element(); iface < gamma[itq]->last_element(); iface += gamma[itq]->next_element()) {

            const int kb = (*gamma[itq]).idxElementOfFace(iface);
            const Element &K((*gamma[itq]).get_element(kb));

            // Get coordinates of current quadrilateral
            const auto &V0(K.at(0)); // vertex 0
            const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

            algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
            algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

            // Get quadrature rule
            algoim::QuadratureRule<2> q =
                algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, order_space);

            // Loop over quadrature in space
            assert(q.nodes.size() != 0);

            double weight_K = 0.;
            for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

                Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));

                const R weight = q.nodes.at(ipq).w;

                assert(weight > 0);

                double err = fh.evalOnBackMesh(kb, 0, mip, t, 0, 0, 0) - f(mip, 0, t);

                weight_K += weight * err * err;
            }
            weight_space += weight_K;
        }
        val += weight_space * weight_time;
    }

    return val; // return \int_{I_n} ||u(t)-uh(t)||_{L^2(Omega(t))}^2 dt
}


// L2(L2(Omega(t)), 0, T)
template <typename L, typename fct_t>
double L2L2_norm(const FunFEM<MeshQuad2> &fh, const fct_t &f, const ActiveMesh<MeshQuad2> &Th, const TimeSlab &In,
                 const QuadratureFormular1d &qTime, L &phi, const int order_space = quadrature_order_integration) {
    using mesh_t    = MeshQuad2;
    using fespace_t = GFESpace<mesh_t>;
    using FElement  = typename fespace_t::FElement;
    using Rd        = typename FElement::Rd;
    using Element   = typename mesh_t::Element;

    double val = 0.;

    const int domain = 0; // do only for main domain

    // Loop in time
    for (int itq = 0; itq < qTime.n; ++itq) {

        // Get quadrature points in time
        GQuadraturePoint<R1> tq((qTime)[itq]);
        const double t = In.mapToPhysicalElement(tq);
        phi.t          = t;

        double weight_time = In.T.measure() * tq.a;

        // Loop in space
        double weight_space = 0.;
        for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

            if (domain != Th.get_domain_element(k))
                continue;

            if (Th.isInactive(k, itq))
                continue;

            const Element &K(Th[k]);
            int kb = Th.idxElementInBackMesh(k);

            // Get coordinates of current quadrilateral
            const auto &V0(K.at(0)); // vertex 0
            const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

            algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
            algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

            // Get quadrature rule
            algoim::QuadratureRule<2> q =
                algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, order_space);

            // Loop over quadrature in space
            assert(q.nodes.size() != 0);

            double weight_K = 0.;
            for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

                Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));

                const R weight = q.nodes.at(ipq).w;

                double err = fh.evalOnBackMesh(kb, domain, mip, t, 0, 0, 0) - f(mip, 0, domain, t);

                weight_K += weight * err * err;
            }
            weight_space += weight_K;
        }
        val += weight_space * weight_time;
    }

    return val; // return \int_{I_n} ||u(t)-uh(t)||_{L^2(Omega(t))}^2 dt
}


// L2(H1(Omega(t)), 0, T)
template <typename L>
double L2H1_norm(const FunFEM<MeshQuad2> &fh, const FunFEM<MeshQuad2> &f, const ActiveMesh<MeshQuad2> &Th, const TimeSlab &In,
                 const QuadratureFormular1d &qTime, L &phi, const int order_space = quadrature_order_integration) {
    using mesh_t    = MeshQuad2;
    using fespace_t = GFESpace<mesh_t>;
    using FElement  = typename fespace_t::FElement;
    using Rd        = typename FElement::Rd;
    using Element   = typename mesh_t::Element;

    //! REQUIRES F TO BE SCALAR

    double val = 0.;

    const int domain = 0; // do only for main domain

    const auto &dfhx = dx(fh.expr());
    const auto &dfhy = dy(fh.expr());
    const auto &dfx = dx(f.expr());
    const auto &dfy = dy(f.expr());

    // Loop in time
    for (int itq = 0; itq < qTime.n; ++itq) {

        // Get quadrature points in time
        GQuadraturePoint<R1> tq((qTime)[itq]);
        const double t = In.mapToPhysicalElement(tq);
        phi.t          = t;

        double weight_time = In.T.measure() * tq.a;

        // Loop in space
        double weight_space = 0.;
        for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

            if (domain != Th.get_domain_element(k))
                continue;

            if (Th.isInactive(k, itq))
                continue;

            const Element &K(Th[k]);
            int kb = Th.idxElementInBackMesh(k);

            // Get coordinates of current quadrilateral
            const auto &V0(K.at(0)); // vertex 0
            const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

            algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
            algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

            // Get quadrature rule
            algoim::QuadratureRule<2> q =
                algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, order_space);

            // Loop over quadrature in space
            assert(q.nodes.size() != 0);

            double weight_K = 0.;
            for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

                Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));

                const R weight = q.nodes.at(ipq).w;

        
                //double derrx = dfhx->eval(kb, mip, nullptr) - dfx->eval(kb, mip, nullptr);
                //double derry = dfhy->eval(kb, mip, nullptr) - dfy->eval(kb, mip, nullptr);
                double derrx = dfhx->evalOnBackMesh(kb, domain, mip, t, nullptr) - dfx->evalOnBackMesh(kb, domain, mip, t, nullptr);
                double derry = dfhy->evalOnBackMesh(kb, domain, mip, t, nullptr) - dfy->evalOnBackMesh(kb, domain, mip, t, nullptr);

                weight_K += weight * (derrx * derrx + derry * derry);

            }

            weight_space += weight_K;
            
        }
        val += weight_space * weight_time;

        
    }

    return val; // return \int_{I_n} ||u(t)-uh(t)||_{L^2(Omega(t))}^2 dt
}

// Time-dependent bulk L2 norm
template <typename Mesh, typename L>
double L2_norm_cut(const FunFEM<Mesh> &fh, R(fex)(double *, int i, int dom, double tt), const TimeSlab &In,
                   const QuadratureFormular1d &qTime, const int itq, L &phi, int c0, int num_comp, const int order_space = quadrature_order_integration) {

    const GFESpace<Mesh> &Vh(*fh.Vh);
    const ActiveMesh<Mesh> &Th(Vh.get_mesh());
    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_cut_2(ui, fex, Th, In, qTime, itq, phi, order_space);
    }
    return sqrt(val);
}

template <typeMesh M, typename L>
double L2_norm_cut_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom, double tt),
                     const ActiveMesh<M> &Th, const TimeSlab &In, const QuadratureFormular1d &qTime, const int itq,
                     L &phi, const int order_space = quadrature_order_integration) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for (int i = 0; i < nb_dom; ++i) {
        val += L2_norm_cut_2(fh, fex, i, Th, In, qTime, itq, phi, order_space);
    }
    return val;
}

template <typename L>
double L2_norm_cut_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom, double tt),
                     int domain, const ActiveMesh<MeshQuad2> &Th, const TimeSlab &In, const QuadratureFormular1d &qTime,
                     const int itq, L &phi, const int order_space = quadrature_order_integration) {
    typedef GFESpace<MeshQuad2> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename ActiveMesh<MeshQuad2>::Element Element;
    typedef typename FElement::Rd Rd;

    double val = 0.;

    GQuadraturePoint<R1> tq((qTime)[itq]);
    const double t = In.mapToPhysicalElement(tq);
    phi.t          = t;
    // std::cout << "t = " << t << "\n";
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (domain != Th.get_domain_element(k))
            continue;

        // Integrate only on Omega(t_itq)
        if (Th.isInactive(k, itq))
            continue;

        const Element &K(Th[k]);
        int kb = Th.idxElementInBackMesh(k);

        int kk = k;

        // Get coordinates of current quadrilateral
        const auto &V0(K.at(0)); // vertex 0
        const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        algoim::QuadratureRule<2> q =
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, order_space);

        // std::cout << "K = [" << V0[0] << ", " << V2[0] << "] x [" << V0[1] << ", " << V2[1] << "]\n";
        // std::cout << phi(xymin) << "\n";

        // Loop over quadrature in space
        assert(q.nodes.size() != 0);
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight = q.nodes.at(ipq).w;

            double a = fh->eval(kk, mip) - fex(mip, fh->cu, domain, t);

            val += weight * a * a;
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

// Stationary bulk L2 norm
template <typename Mesh, typename L>
double L2_norm_cut(const FunFEM<Mesh> &fh, R(fex)(double *, int i, int dom), L &phi, int c0, int num_comp, const int order_space = quadrature_order_integration) {

    const GFESpace<Mesh> &Vh(*fh.Vh);
    const ActiveMesh<Mesh> &Th(Vh.get_mesh());
    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_cut_2(ui, fex, Th, phi, order_space);
    }
    return sqrt(val);
}

template <typename M, typename L>
double L2_norm_cut_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom),
                     const ActiveMesh<M> &Th, L &phi, const int order_space = quadrature_order_integration) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for (int i = 0; i < nb_dom; ++i) {
        val += L2_norm_cut_2(fh, fex, i, Th, phi, order_space);
    }
    return val;
}

template <typename L>
double L2_norm_cut_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom), int domain,
                     const ActiveMesh<MeshQuad2> &Th, L &phi, const int order_space = quadrature_order_integration) {
    typedef GFESpace<MeshQuad2> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename ActiveMesh<MeshQuad2>::Element Element;
    typedef typename FElement::Rd Rd;

    double val = 0.;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (domain != Th.get_domain_element(k))
            continue;

        // const Cut_Part<Element> cutK(Th.get_cut_part(k, 0));
        const Element &K(Th[k]);
        int kb = Th.idxElementInBackMesh(k);

        int kk = k;
        // if(macro){
        //   if(!macro->isRootFat(k)) {
        //     kk = macro->getIndexRootElement(k);
        //   }
        // }

        // Get coordinates of current quadrilateral
        const auto &V0(K.at(0)); // vertex 0
        const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        algoim::QuadratureRule<2> q =
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, order_space);

        // Loop over quadrature in space
        assert(q.nodes.size() != 0);
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight = q.nodes.at(ipq).w;
            assert(weight > 0);

            double a = fh->eval(kk, mip) - fex(mip, fh->cu, domain);

            val += weight * a * a;
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

template <meshTriag M, typename L, typename FEX>
double L2_norm_cut_2(const std::shared_ptr<ExpressionVirtual> &fh, const FEX &fex, int domain, const ActiveMesh<M> &Th, L &phi, const int order_space = quadrature_order_integration) {
    typedef GFESpace<M> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename ActiveMesh<M>::Element Element;
    typedef typename FElement::Rd Rd;
    typedef typename FElement::QF QF;
    typedef typename QF::QuadraturePoint QuadraturePoint;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(quadrature_order_integration));

    double val = 0.;
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (domain != Th.get_domain_element(k))
            continue;

        if (Th.isInactive(k, 0))
            continue;

        const Element &K(Th[k]);
        int kb = Th.idxElementInBackMesh(k);
        
        int kk = k;

        if (Th.isCut(k, 0)) {
            using Vec2 = algoim::uvector<double,2>;
            std::array<Vec2,3> verts = {
                Vec2(K.at(0)[0], K.at(0)[1]),
                Vec2(K.at(1)[0], K.at(1)[1]),
                Vec2(K.at(2)[0], K.at(2)[1])
            };

            const auto q = cuttri_multipoly::compute_cut_triangle_quadrature(
                verts,
                phi,
                bernstein_deg,   // bernstein degree
                quadrature_order_integration              // algoim 1D quadrature order
            );

            // Loop over quadrature in space
            // assert(q.x_vol.size() != 0);
            for (int ipq = 0; ipq < q.x_vol.size(); ++ipq) {

                Rd mip(q.x_vol.at(ipq)[0], q.x_vol.at(ipq)[1]);
                const R weight = q.w_vol.at(ipq);
                assert(weight > 0);

                double a = fh->eval(kk, mip) - fex(mip, fh->cu, domain);

                val += weight * a * a;
            } 
        } else {
            double meas = K.measure();
            
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                QuadraturePoint ip(qf[ipq]);
                Rd mip = K.mapToPhysicalElement(ip);
                double weight  = meas * ip.getWeight();
                assert(weight > 0);

                double a = fh->eval(kk, mip) - fex(mip, fh->cu, domain);
                
                val += weight * a * a;
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


template <meshTriag mesh_t, typename L>
double max_norm_cut(const std::shared_ptr<ExpressionVirtual> &fh, const ActiveMesh<mesh_t> &Th, int domain, L& phi, const int order_space = quadrature_order_integration) {

    using fespace_t = GFESpace<mesh_t>;
    using fe_t      = typename fespace_t::FElement;
    using e_t       = typename mesh_t::Element;
    using QF        = typename fe_t::QF;
    using QFB       = typename fe_t::QFB;
    using v_t       = typename fe_t::Rd;
    using qp_t      = typename QF::QuadraturePoint;

    const QF &qf(*QF_Simplex<typename fe_t::RdHat>(order_space));
    const QFB &qfb(*QF_Simplex<typename fe_t::RdHatBord>(order_space));

    What_d Fop = Fwhatd(op_id);
    double val = 0.;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (domain != Th.get_domain_element(k))
            continue;
        const e_t &K(Th[k]);
        const Cut_Part<e_t> cutK(Th.get_cut_part(k, 0));
        int kb = Th.idxElementInBackMesh(k);

        for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                auto ip(qf[ipq]); // integration point
                v_t mip = cutK.mapToPhysicalElement(it, ip);
                val     = std::max(val, fabs(fh->eval(k, mip)));
            }

            for (int ifac = 0; ifac < e_t::nea; ++ifac) {
                for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
                    auto ip(qfb[ipq]); // integration point
                    auto ipf = K.mapToReferenceElement(ip, ifac);
                    v_t mip  = cutK.mapToPhysicalElement(it, ipf);
                    val      = std::max(val, fabs(fh->eval(k, mip)));
                }
            }
        }
    }
    double val_receive = 0;
#ifdef USE_MPI
    MPIcf::AllReduce(val, val_receive, MPI_MAX);
#else
    val_receive = val;
#endif

    return val_receive;
}


template <typename L, typename fct_t>
double integral_algoim(fct_t &fh, const Interface<MeshQuad2> &interface, int cu, L &phi, const int order_space = quadrature_order_integration) {

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

    for (int iface = interface.first_element(); iface < interface.last_element(); iface += interface.next_element()) {

        const int kb = interface.idxElementOfFace(iface); // idx on backMesh

        const auto &T(interface.get_element(kb));
        const auto &V0(T.at(0)); // vertex 0
        const auto &V2(T.at(2)); // vertex 2   diagonally opposed

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        algoim::QuadratureRule<2> q =
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, order_space);

        assert(q.nodes.size() != 0);
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
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

template <typename L, typename fct_t>
double integral_algoim(fct_t &fh, const Interface<MeshQuad2> &interface, int cu, L &phi, const double t, const int order_space = quadrature_order_integration) {

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

    phi.t = t;

    for (int iface = interface.first_element(); iface < interface.last_element(); iface += interface.next_element()) {

        const int kb = interface.idxElementOfFace(iface); // idx on backMesh

        const auto &T(interface.get_element(kb));
        const auto &V0(T.at(0)); // vertex 0
        const auto &V2(T.at(2)); // vertex 2   diagonally opposed

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        algoim::QuadratureRule<2> q =
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, order_space);

        assert(q.nodes.size() != 0);
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
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

template <typename L, typename fct_t>
double integral_algoim(fct_t &fh, const Interface<MeshQuad2> &interface, int cu, L &phi, const TimeSlab &In,
                       const QuadratureFormular1d &qTime, const int itq, const int order_space = quadrature_order_integration) {

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

    GQuadraturePoint<R1> tq((qTime)[itq]);
    const double t = In.mapToPhysicalElement(tq);
    phi.t          = t;

    for (int iface = interface.first_element(); iface < interface.last_element(); iface += interface.next_element()) {

        const int kb = interface.idxElementOfFace(iface); // idx on backMesh

        const auto &T(interface.get_element(kb));
        const auto &V0(T.at(0)); // vertex 0
        const auto &V2(T.at(2)); // vertex 2   diagonally opposed

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        algoim::QuadratureRule<2> q =
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, order_space);

        assert(q.nodes.size() != 0);
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight = q.nodes.at(ipq).w;

            const R Cint = weight;
            if constexpr (std::is_same_v<fct_t, FunFEM<MeshQuad2>>) {
                val += Cint * fh.evalOnBackMesh(kb, 0, mip, t, cu, 0, 0);
            } else if constexpr (std::is_same_v<fct_t, std::shared_ptr<ExpressionVirtual>>) {
                val += Cint * fh->evalOnBackMesh(kb, 0, mip, t, nullptr);
            }

            else {
                val += Cint * fh(mip, 0, t);
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

/**
 * @brief Integrate stationary function over stationary domain
 *
 * @tparam L
 * @tparam fct_t
 * @param fh
 * @param Th
 * @param phi
 * @return double
 */
template <typename L, typename fct_t>
double integral_algoim(fct_t &fh, const ActiveMesh<MeshQuad2> &Th, L &phi, int c0, const int order_space = quadrature_order_integration) {

    assert(Th.get_nb_domain() == 1);

    using mesh_t        = MeshQuad2;
    using fespace_t     = GFESpace<mesh_t>;
    using itemVFlist_t  = ListItemVF<mesh_t>;
    using FElement      = typename fespace_t::FElement;
    using Rd            = typename FElement::Rd;
    using QF            = typename FElement::QF;
    using QFB           = typename FElement::QFB;
    using Element       = typename mesh_t::Element;
    using BorderElement = typename mesh_t::BorderElement;

    double val       = 0.;
    const int domain = 0;
    const int cu     = 1; // number of components

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (domain != Th.get_domain_element(k))
            continue;

        const Element &K(Th[k]);
        int kb = Th.idxElementInBackMesh(k);

        // Get coordinates of current quadrilateral
        const auto &V0(K.at(0)); // vertex 0
        const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        // Get quadrature rule
        algoim::QuadratureRule<2> q =
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, order_space);

        // if (Th.isCut(k, 0)) {
        //     std::cout << "This element is not cut: here comes the quadrature: "
        //               << "\n";
        //     for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {
        //         Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
        //         std::cout << "mip = " << mip[0] << ", " << mip[1] << "\n";
        //     }
        // }
        // Loop over quadrature in space
        assert(q.nodes.size() != 0);
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight = q.nodes.at(ipq).w;

            const R Cint = weight;
            if constexpr (std::is_same_v<fct_t, FunFEM<MeshQuad2>>) {
                val += Cint * fh.evalOnBackMesh(kb, domain, mip, 0, 0);
            } else {
                val += Cint * fh->evalOnBackMesh(kb, domain, mip);
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

/**
 * @brief Integrate function over time-dependent cut domain over In
 */
template <typename L, typename fct_t>
double integral_algoim(fct_t &fh, const int cu, const ActiveMesh<MeshQuad2> &Th, L &phi, const TimeSlab &In,
                       const QuadratureFormular1d &qTime, const int order_space = quadrature_order_integration) {

    const int number_of_domains = Th.get_nb_domain();
    double val                  = 0.;

    for (int i = 0; i < number_of_domains; ++i) {
        val += integral_algoim(fh, cu, Th, i, phi, In, qTime, order_space);
    }

    return val;
}

template <typename L, typename fct_t>
double integral_algoim(fct_t &fh, const int cu, const ActiveMesh<MeshQuad2> &Th, const int domain, L &phi,
                       const TimeSlab &In, const QuadratureFormular1d &qTime, const int order_space = quadrature_order_integration) {

    using mesh_t    = MeshQuad2;
    using fespace_t = GFESpace<mesh_t>;
    using FElement  = typename fespace_t::FElement;
    using Rd        = typename FElement::Rd;
    using Element   = typename mesh_t::Element;

    double val = 0.;

    for (int itq = 0; itq < qTime.n; ++itq) {

        GQuadraturePoint<R1> tq((qTime)[itq]);
        const double t = In.mapToPhysicalElement(tq);
        phi.t          = t;

        for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

            if (domain != Th.get_domain_element(k))
                continue;
            // if (Th.isInactive(k, itq))
            //     continue;

            const Element &K(Th[k]);
            int kb = Th.idxElementInBackMesh(k);

            // Get coordinates of current quadrilateral
            const auto &V0(K.at(0)); // vertex 0
            const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

            algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
            algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

            // Get quadrature rule
            algoim::QuadratureRule<2> q =
                algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, order_space);

            // Loop over quadrature in space
            //assert(q.nodes.size() != 0);
            if (q.nodes.size()==0)
                continue;

            for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

                Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
                const R weight = q.nodes.at(ipq).w;

                assert(weight > 0);

                const R Cint = weight * In.T.measure() * tq.a;
                if constexpr (std::is_same_v<fct_t, FunFEM<MeshQuad2>>) {
                    val += Cint * fh.evalOnBackMesh(kb, domain, mip, t, cu, 0, 0);

                } else if constexpr (std::is_same_v<fct_t, std::shared_ptr<ExpressionVirtual>>) {
                    val += Cint * fh->evalOnBackMesh(kb, domain, mip, t, 0);
                } else {
                    val += Cint * fh(mip, domain, t);
                }
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

/**
 * @brief Integrate function over time-dependent cut domains in specific quadrature point
 *
 * @tparam L Algoim level set function struct
 * @tparam fct_t Integrand type
 * @param fh Integrand
 * @param Th Active Mesh
 * @param phi Algoim level set function
 * @param In Time slab
 * @param qTime Time quadrature rule
 * @param itq Time quadrature point
 * @return double Integral of the integrand over the geometry in time quadrature point itq of the time slab In
 *
 * @note This integral does NOT scale with dT. It loops over all defined subdomains in the active mesh.
 */
template <typename L, typename fct_t>
double integral_algoim(const fct_t &fh, const ActiveMesh<MeshQuad2> &Th, L &phi, const TimeSlab &In,
                       const QuadratureFormular1d &qTime, const int itq, const int order_space = quadrature_order_integration) {

    const int number_of_domains = Th.get_nb_domain();
    double val                  = 0.;
    for (int i = 0; i < number_of_domains; ++i) {
        val += integral_algoim(fh, Th, i, phi, In, qTime, itq, order_space);
    }

    return val;
}
/**
 * @brief Integrate function over time-dependent cut domain in specific quadrature point.
 * @note For specific subdomain, or if the active mesh only contains one subdomain.
 * @param domain Domain index
 *
 */
template <typename L, typename fct_t>
double integral_algoim(const fct_t &fh, const ActiveMesh<MeshQuad2> &Th, const int domain, L &phi, const TimeSlab &In,
                       const QuadratureFormular1d &qTime, const int itq, const int order_space = quadrature_order_integration) {

    using mesh_t    = MeshQuad2;
    using fespace_t = GFESpace<mesh_t>;
    using FElement  = typename fespace_t::FElement;
    using Rd        = typename FElement::Rd;
    using Element   = typename mesh_t::Element;

    double val = 0.;

    GQuadraturePoint<R1> tq((qTime)[itq]);
    const double t = In.mapToPhysicalElement(tq);
    phi.t          = t;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (domain != Th.get_domain_element(k))
            continue;

        if (Th.isInactive(k, itq))
            continue;

        const Element &K(Th[k]);
        int kb = Th.idxElementInBackMesh(k);

        // Get coordinates of current quadrilateral
        const auto &V0(K.at(0)); // vertex 0
        const auto &V2(K.at(2)); // vertex 2 (diagonally opposed)

        algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
        algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

        // Get quadrature rule
        algoim::QuadratureRule<2> q =
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, order_space);

        // Loop over quadrature in space
        assert(q.nodes.size() != 0);
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight = q.nodes.at(ipq).w;

            assert(weight > 0);

            if constexpr (std::is_same_v<fct_t, FunFEM<MeshQuad2>>) {
                val += weight * fh.evalOnBackMesh(kb, domain, mip, 0, 0);
            } else if constexpr ((std::is_same_v<fct_t, std::shared_ptr<ExpressionVirtual>>) || (std::is_same_v<fct_t, std::shared_ptr<ExpressionSum>>)) {
                val += weight * fh->evalOnBackMesh(kb, domain, mip, t, 0);
            } else {
                val += weight * fh(mip, domain, t);
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

/**
 * @brief Integrate FunFEM function over interface and over In
 *
 * @tparam M Mesh
 * @tparam L Algoim level set struct
 * @param fh FunFEM object
 * @param In Time Slab
 * @param gamma TimeInterface
 * @param phi L object
 * @param cu Component number
 * @return double int_In int_Gamma(t) fh ds dt
 * @note This integral scales with dT.
 */
// template <typename M, typename L>
// double integral_algoim(FunFEM<M> &fh, const TimeSlab &In, const TimeInterface<M> &gamma, L &phi, int cu) {

//     using mesh_t    = MeshQuad2;
//     using fespace_t = GFESpace<mesh_t>;
//     using FElement  = typename fespace_t::FElement;
//     using Rd        = typename FElement::Rd;
//     using Element   = typename mesh_t::Element;

//     double val = 0.;

//     // Loop over time quadrature points
//     for (int it = 0; it < gamma.size(); ++it) {
//         const Interface<mesh_t> &interface(*gamma(it));
//         const QuadratureFormular1d *qTime(gamma.get_quadrature_time());
//         GQuadraturePoint<R1> tq((*qTime)[it]);
//         const double t = In.mapToPhysicalElement(tq);
//         phi.t          = t;

//         for (int iface = interface.first_element(); iface < interface.last_element();
//              iface += interface.next_element()) {

//             const int kb = interface.idxElementOfFace(iface); // idx on backMesh
//             // const R meas = interface.measure(iface);

//             const auto &T(interface.get_element(kb));
//             const auto &V0(T.at(0)); // vertex 0
//             const auto &V2(T.at(2)); // vertex 2   diagonally opposed

//             algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
//             algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

//             algoim::QuadratureRule<2> q =
//                 algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order_integration);

//             assert(q.nodes.size() != 0);
//             for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

//                 Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
//                 const R weight             = q.nodes.at(ipq).w;
//                 // const R Cint = meas * ip.getWeight() * In.T.mesure() * tq.a;
//                 const R Cint               = weight * In.T.measure() * tq.a;
//                 const int domain_interface = 0;
//                 val += Cint * fh.evalOnBackMesh(kb, domain_interface, mip, t, cu, 0, 0);
//             }
//         }
//     }
//     double val_receive = 0;
// #ifdef USE_MPI
//     MPIcf::AllReduce(val, val_receive, MPI_SUM);
// #else
//     val_receive = val;
// #endif

//     return val_receive;
// }

template <typename L, typename fct_t>
double integral_algoim(fct_t &fh, const TimeSlab &In, const TimeInterface<MeshQuad2> &gamma, L &phi, int cu, const int order_space = quadrature_order_integration) {

    using mesh_t    = MeshQuad2;
    using fespace_t = GFESpace<mesh_t>;
    using FElement  = typename fespace_t::FElement;
    using Rd        = typename FElement::Rd;
    using Element   = typename mesh_t::Element;

    double val = 0.;

    // Loop over time quadrature points
    for (int it = 0; it < gamma.size(); ++it) {
        const Interface<mesh_t> &interface(*gamma(it));
        const QuadratureFormular1d *qTime(gamma.get_quadrature_time());
        GQuadraturePoint<R1> tq((*qTime)[it]);
        const double t = In.mapToPhysicalElement(tq);
        phi.t          = t;

        for (int iface = interface.first_element(); iface < interface.last_element();
             iface += interface.next_element()) {

            const int kb = interface.idxElementOfFace(iface); // idx on backMesh
            // const R meas = interface.measure(iface);

            const auto &T(interface.get_element(kb));
            const auto &V0(T.at(0)); // vertex 0
            const auto &V2(T.at(2)); // vertex 2   diagonally opposed

            algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
            algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y

            algoim::QuadratureRule<2> q =
                algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, order_space);

            assert(q.nodes.size() != 0);
            for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

                Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
                const R weight             = q.nodes.at(ipq).w;
                // const R Cint = meas * ip.getWeight() * In.T.mesure() * tq.a;
                const R Cint               = weight * In.T.measure() * tq.a;
                const int domain_interface = 0;
                // val += Cint * fh.evalOnBackMesh(kb, domain_interface, mip, t, cu, 0, 0);

                if constexpr (std::is_same_v<fct_t, FunFEM<MeshQuad2>>) {
                    val += Cint * fh.evalOnBackMesh(kb, domain_interface, mip, t, cu, 0, 0);
                } else if constexpr (std::is_same_v<fct_t, std::shared_ptr<ExpressionVirtual>>) {
                    val += Cint * fh->evalOnBackMesh(kb, domain_interface, mip, t, cu);
                } else {
                    val += Cint * fh(mip, 0, t);
                }
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







// integrals over faceas
template <typename L, typename fct_t>
double integral_algoim(fct_t &fh, const ActiveMesh<MeshQuad2> &Th, L &phi, const CFacet &b) {
    
    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        
        int kb  = Th.idxElementInBackMesh(k, 0);    // element index in background mesh

        for (int ifac = 0; ifac < MeshQuad2::Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            
            int kn  = Th.ElementAdj(k, jfac);
            int kbn = Th.idxElementInBackMesh(kn, 0);   // element index of neighbor in background mesh


            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kbn == -1 || kbn < kb) 
                continue;
            
            std::cout << "\n";
            std::cout << "k = " << k << ", kn = " << kn << "\n";
            std::cout << "kb = " << kb << ", kbn = " << kbn << "\n";
            std::cout << "ifac = " << ifac << ", jfac = " << jfac << "\n";
            

            // // std::pair<int, int> e1 = std::make_pair(kb, ifac);
            // // std::pair<int, int> e2 = std::make_pair(kbn, jfac);
            
            // // int kbi = e1.first, ifac = e1.second;
            // // int kbj = e2.first, jfac = e2.second;
            
            // // int ki = Th.idxElementFromBackMesh(kbi);
            // // int kj = Th.idxElementFromBackMesh(kbj);

            const MeshQuad2::Element &K(Th[k]);
            const MeshQuad2::Element &Kn(Th[kn]);



            // Get coordinates of current quadrilateral
            const auto &V0(K[0]); // vertex 0
            const auto &V2(K[2]); // vertex 2 (diagonally opposed)

            algoim::uvector<double, 2> xymin{V0[0], V0[1]}; // min x and y
            algoim::uvector<double, 2> xymax{V2[0], V2[1]}; // max x and y


            // print coordinates of the quadrilateral
            std::cout << "V0 = " << V0[0] << ", " << V0[1] << "\n";
            std::cout << "V2 = " << V2[0] << ", " << V2[1] << "\n\n";

            //getchar();
            int dim = -1, side = -1;
            
            if (ifac % 2 == 0) {
                dim = 1;
                side = 1;    // vertical direction, top face
            }
            else {
                dim = 0;
                side = 1;    // horizontal direction, right face 
            }

            assert(dim != -1 || side != -1);

            std::cout << "dim = " << dim << ", side = " << side << "\n";
            algoim::QuadratureRule<2> q =
                    algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), dim, side, 5);

            for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {
                R2 mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
                std::cout << "mip = " << mip[0] << ", " << mip[1] << "\n";
            }

            getchar();

            // const Cut_Part<Element> cutKi(Th.get_cut_part(ki, itq));
            // const Cut_Part<Element> cutKj(Th.get_cut_part(kj, itq));
            // typename Element::Face face;
            // const Cut_Part<typename Element::Face> cutFace(Th.get_cut_face(face, ki, ifac, itq));

            // double measK = cutKi.measure() + cutKj.measure();
            // double measF = Ki.mesureBord(ifac);
            // double h     = 0.5 * (Ki.get_h() + Kj.get_h());
            // int domain   = FKi.get_domain();
            // Rd normal    = Ki.N(ifac);

            // int kbi = Vh.idxElementInBackMesh(ki);
            // int kbj = Vh.idxElementInBackMesh(kj);
            

            // // LOOP OVER ELEMENTS IN THE CUT
            // for (auto it = cutFace.element_begin(); it != cutFace.element_end(); ++it) {

            //     for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
            //         typename QFB::QuadraturePoint ip(qfb[ipq]);
            //         Rd mip = cutFace.mapToPhysicalElement(it, (RdHatBord)ip);
            //         double Cint  = meas * ip.getWeight() * cst_time;

            //         // EVALUATE THE BASIS FUNCTIONS
            //         FKv.BF(Fop, FKv.T.mapToReferenceElement(mip), fv);
            //         if (!same)
            //             FKu.BF(Fop, FKu.T.mapToReferenceElement(mip), fu);
            //         //   VF[l].applyFunNL(fu,fv);

            //         Cint *= VF[l].evaluateFunctionOnBackgroundMesh(std::make_pair(kbu, kbv), std::make_pair(domain, domain),
            //                                                     mip, tid, normal);
            //         Cint *= coef * VF[l].c;
            //         if (In) {
            //             if (VF.isRHS())
            //                 this->addToRHS(VF[l], *In, FKv, fv, Cint);
            //             else
            //                 this->addToMatrix(VF[l], *In, FKu, FKv, fu, fv, Cint);
            //         } else {
            //             if (VF.isRHS())
            //                 this->addToRHS(VF[l], FKv, fv, Cint);
            //             else
            //                 this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint);
            //         }
            //     }
            
            // }
                
        }
        
    }

    return 0;
}


#endif