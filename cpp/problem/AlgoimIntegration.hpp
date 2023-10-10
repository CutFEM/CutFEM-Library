
#ifndef ALGOIM_INTEGRATION_HPP
#define ALGOIM_INTEGRATION_HPP

#include "../common/AlgoimInterface.hpp"

#include "../common/AlgoimInterface.hpp"

const int quadrature_order = 5;

template <typename Mesh, typename L>
double L2_norm_surface_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(const R2, int i),
                         const Interface<Mesh> &interface, L &phi) {

    using fespace_t = GFESpace<Mesh>;
    using FElement  = typename fespace_t::FElement;
    using Element   = typename Mesh::Element;
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
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

        // assert((q.nodes.size() == 1 * quadrature_order) ||
        //    (q.nodes.size() == 2 * quadrature_order));
        assert(q.nodes.size() != 0);
        assert((quadrature_order <= q.nodes.size()) || (q.nodes.size() <= 2 * quadrature_order));
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
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

template <typename Mesh, typename L>
double L2_norm_surface_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(const R2, int i, double t),
                         const Interface<Mesh> &interface, double tt, L &phi) {

    using fespace_t = GFESpace<Mesh>;
    using FElement  = typename fespace_t::FElement;
    using Element   = typename Mesh::Element;
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
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

        // assert((q.nodes.size() == 1 * quadrature_order) ||
        //    (q.nodes.size() == 2 * quadrature_order));
        assert((quadrature_order <= q.nodes.size()) || (q.nodes.size() <= 2 * quadrature_order));
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            const Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
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

template <typename Mesh, typename L>
double L2_norm_surface_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, double t),
                         const Interface<Mesh> &interface, double tt, L &phi) {

    using fespace_t = GFESpace<Mesh>;
    using FElement  = typename fespace_t::FElement;
    using Element   = typename Mesh::Element;
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
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

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

template <typename Mesh, typename L>
double L2_norm_surface(const FunFEM<Mesh> &fh, R(fex)(const R2, int i), const Interface<Mesh> &interface, L &phi,
                       int c0, int num_comp) {

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_surface_2(ui, fex, interface, phi);
    }
    return sqrt(val);
}

template <typename Mesh, typename L>
double L2_norm_surface(const FunFEM<Mesh> &fh, R(fex)(const R2, int i, double t), const Interface<Mesh> &interface,
                       double tt, L &phi, int c0, int num_comp) {

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_surface_2(ui, fex, interface, tt, phi);
    }
    return sqrt(val);
}

template <typename Mesh, typename L>
double L2_norm_surface(const FunFEM<Mesh> &fh, R(fex)(double *, int i, double t), const Interface<Mesh> &interface,
                       double tt, L &phi, int c0, int num_comp) {

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_surface_2(ui, fex, interface, tt, phi);
    }
    return sqrt(val);
}

// L2(In x Omega)
template <typename L, typename fct_t>
double L2_norm_T(const FunFEM<MeshQuad2> &fh, const fct_t &f, const ActiveMesh<MeshQuad2> &Th, const TimeSlab &In,
                 const QuadratureFormular1d &qTime, L &phi) {
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
                algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, quadrature_order);

            // Loop over quadrature in space
            assert(q.nodes.size() != 0);

            double weight_K = 0.;
            for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

                Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
                
                const R weight = q.nodes.at(ipq).w;

                double err = fh.evalOnBackMesh(kb, domain, mip, t, 0, 0, 0) - f(mip, 0, domain, t);

                weight_K += weight*err*err;

            }
            weight_space += weight_K;
        }
        val += weight_space*weight_time;

    }

    return sqrt(val);
}

// Time-dependent bulk L2 norm
template <typename Mesh, typename L>
double L2_norm_cut(const FunFEM<Mesh> &fh, R(fex)(double *, int i, int dom, double tt), const TimeSlab &In,
                   const QuadratureFormular1d &qTime, const int itq, L &phi, int c0, int num_comp) {

    const GFESpace<Mesh> &Vh(*fh.Vh);
    const ActiveMesh<Mesh> &Th(Vh.get_mesh());
    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_cut_2(ui, fex, Th, In, qTime, itq, phi);
    }
    return sqrt(val);
}

template <typename M, typename L>
double L2_norm_cut_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom, double tt),
                     const ActiveMesh<M> &Th, const TimeSlab &In, const QuadratureFormular1d &qTime, const int itq,
                     L &phi) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for (int i = 0; i < nb_dom; ++i) {
        val += L2_norm_cut_2(fh, fex, i, Th, In, qTime, itq, phi);
    }
    return val;
}

template <typename Mesh, typename L>
double L2_norm_cut_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom, double tt),
                     int domain, const ActiveMesh<Mesh> &Th, const TimeSlab &In, const QuadratureFormular1d &qTime,
                     const int itq, L &phi) {
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename ActiveMesh<Mesh>::Element Element;
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
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, quadrature_order);

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
double L2_norm_cut(const FunFEM<Mesh> &fh, R(fex)(double *, int i, int dom), L &phi, int c0, int num_comp) {

    const GFESpace<Mesh> &Vh(*fh.Vh);
    const ActiveMesh<Mesh> &Th(Vh.get_mesh());
    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2_norm_cut_2(ui, fex, Th, phi);
    }
    return sqrt(val);
}

template <typename M, typename L>
double L2_norm_cut_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom),
                     const ActiveMesh<M> &Th, L &phi) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for (int i = 0; i < nb_dom; ++i) {
        val += L2_norm_cut_2(fh, fex, i, Th, phi);
    }
    return val;
}

template <typename Mesh, typename L>
double L2_norm_cut_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom), int domain,
                     const ActiveMesh<Mesh> &Th, L &phi) {
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename ActiveMesh<Mesh>::Element Element;
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
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, quadrature_order);

        // Loop over quadrature in space
        assert(q.nodes.size() != 0);
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight = q.nodes.at(ipq).w;

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

template <typename L, typename fct_t>
double integral_algoim(fct_t &fh, const Interface<MeshQuad2> &interface, int cu, L &phi) {

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
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

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
double integral_algoim(fct_t &fh, const Interface<MeshQuad2> &interface, int cu, L &phi, const double t) {

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
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

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
                       const QuadratureFormular1d &qTime, const int itq) {

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
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

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
double integral_algoim(fct_t &fh, const ActiveMesh<MeshQuad2> &Th, L &phi, int c0) {

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
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, quadrature_order);

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
                       const QuadratureFormular1d &qTime) {

    const int number_of_domains = Th.get_nb_domain();
    double val                  = 0.;

    for (int i = 0; i < number_of_domains; ++i) {
        val += integral_algoim(fh, cu, Th, i, phi, In, qTime);
    }

    return val;
}

template <typename L, typename fct_t>
double integral_algoim(fct_t &fh, const int cu, const ActiveMesh<MeshQuad2> &Th, const int domain, L &phi,
                       const TimeSlab &In, const QuadratureFormular1d &qTime) {

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
                algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, quadrature_order);

            // Loop over quadrature in space
            assert(q.nodes.size() != 0);
            for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

                Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
                const R weight = q.nodes.at(ipq).w;

                const R Cint = weight * In.T.measure() * tq.a;
                if constexpr (std::is_same_v<fct_t, FunFEM<MeshQuad2>>) {
                    val += Cint * fh.evalOnBackMesh(kb, domain, mip, t, cu, 0, 0);

                } else if constexpr (std::is_same_v<fct_t, std::shared_ptr<ExpressionVirtual>>) {
                    val += Cint * fh->evalOnBackMesh(kb, domain, mip, t);
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
double integral_algoim(fct_t &fh, const ActiveMesh<MeshQuad2> &Th, L &phi, const TimeSlab &In,
                       const QuadratureFormular1d &qTime, const int itq) {

    const int number_of_domains = Th.get_nb_domain();
    double val                  = 0.;
    for (int i = 0; i < number_of_domains; ++i) {
        val += integral_algoim(fh, Th, i, phi, In, qTime, itq);
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
double integral_algoim(fct_t &fh, const ActiveMesh<MeshQuad2> &Th, const int domain, L &phi, const TimeSlab &In,
                       const QuadratureFormular1d &qTime, const int itq) {

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
            algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), -1, -1, quadrature_order);

        // Loop over quadrature in space
        assert(q.nodes.size() != 0);
        for (int ipq = 0; ipq < q.nodes.size(); ++ipq) {

            Rd mip(q.nodes.at(ipq).x(0), q.nodes.at(ipq).x(1));
            const R weight = q.nodes.at(ipq).w;

            if constexpr (std::is_same_v<fct_t, FunFEM<MeshQuad2>>) {
                val += weight * fh.evalOnBackMesh(kb, domain, mip, 0, 0);
            } else if constexpr (std::is_same_v<fct_t, std::shared_ptr<ExpressionVirtual>>) {
                val += weight * fh->evalOnBackMesh(kb, domain, mip);
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
//                 algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

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
double integral_algoim(fct_t &fh, const TimeSlab &In, const TimeInterface<MeshQuad2> &gamma, L &phi, int cu) {

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
                algoim::quadGen<2>(phi, algoim::HyperRectangle<double, 2>(xymin, xymax), 2, -1, quadrature_order);

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

#endif