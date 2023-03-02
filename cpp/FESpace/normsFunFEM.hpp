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

template <typename M>
double L2normCut(const FunFEM<M> &fh, R(fex)(double *, int i, int dom, double tt), double t, int c0, int num_comp,
                 const MacroElement<M> *macro = nullptr) {

    const GFESpace<M> &Vh(*fh.Vh);
    const ActiveMesh<M> &Th(Vh.get_mesh());
    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2normCut_2(ui, fex, Th, t, macro);
    }
    return sqrt(val);
}

template <typeMesh mesh_t, FunctionDomain fct>
double L2normCut(const FunFEM<mesh_t> &fh, fct fex, int c0, int num_comp, const MacroElement<mesh_t> *macro = nullptr) {

    const GFESpace<mesh_t> &Vh(*fh.Vh);
    const ActiveMesh<mesh_t> &Th(Vh.get_mesh());
    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2normCut_2(ui, fex, Th, macro);
    }
    return sqrt(val);
}

template <typeMesh mesh_t, FunctionDomain fct>
double L2normCut(const std::shared_ptr<ExpressionVirtual> &fh, fct fex, const ActiveMesh<mesh_t> &Th,
                 const MacroElement<mesh_t> *macro = nullptr) {
    double val = L2normCut_2(fh, fex, Th, macro);
    return sqrt(val);
}

template <typename M>
double L2normCut(const std::shared_ptr<ExpressionVirtual> &fh, const ActiveMesh<M> &Th,
                 const MacroElement<M> *macro = nullptr) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for (int i = 0; i < nb_dom; ++i) {
        val += L2normCut_2(fh, i, Th, macro);
    }
    return sqrt(val);
}

template <typename M>
double L2normCut(const std::shared_ptr<ExpressionVirtual> &fh, const ActiveMesh<M> &Th, int dom,
                 const MacroElement<M> *macro = nullptr) {
    double val = L2normCut_2(fh, dom, Th, macro);
    return sqrt(val);
}

template <typeMesh mesh_t, FunctionDomain fct_t>
double L2normCut_2(const std::shared_ptr<ExpressionVirtual> &fh, fct_t fex, const ActiveMesh<mesh_t> &Th,
                   const MacroElement<mesh_t> *macro) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for (int i = 0; i < nb_dom; ++i) {
        val += L2normCut_2(fh, fex, i, Th, macro);
    }
    return val;
}

template <typename M>
double L2normCut_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom, double tt),
                   const ActiveMesh<M> &Th, double t, const MacroElement<M> *macro) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for (int i = 0; i < nb_dom; ++i) {
        val += L2normCut_2(fh, fex, i, Th, t, macro);
    }
    return val;
}

template <typename Mesh>
double L2normCut_2(const std::shared_ptr<ExpressionVirtual> &fh, int domain, const ActiveMesh<Mesh> &Th,
                   const MacroElement<Mesh> *macro) {
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename ActiveMesh<Mesh>::Element Element;
    typedef typename FElement::QF QF;
    typedef typename FElement::Rd Rd;
    typedef typename QF::QuadraturePoint QuadraturePoint;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(5));

    double val = 0.;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (domain != Th.get_domain_element(k))
            continue;

        const Cut_Part<Element> cutK(Th.get_cut_part(k, 0));
        int kb = Th.idxElementInBackMesh(k);

        int kk = k;
        // if(macro){
        //   if(!macro->isRootFat(k)) {
        //     kk = macro->getIndexRootElement(k);
        //   }
        // }

        for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

            const R meas = cutK.measure(it);

            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

                QuadraturePoint ip(qf[ipq]); // integration point
                Rd mip       = cutK.mapToPhysicalElement(it, ip);
                const R Cint = meas * ip.getWeight();
                // std::cout << " before "  << std::endl;

                double a = fh->eval(kk, mip);
                // std::cout << " after "  << std::endl;

                val += Cint * a * a;
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

template <typeMesh mesh_t, FunctionDomain fct_t>
double L2normCut_2(const std::shared_ptr<ExpressionVirtual> &fh, fct_t fex, int domain, const ActiveMesh<mesh_t> &Th,
                   const MacroElement<mesh_t> *macro) {
    using fespace_t = GFESpace<mesh_t>;
    using fe_t      = typename fespace_t::FElement;
    using e_t       = typename mesh_t::Element;
    using QF        = typename fe_t::QF;
    using v_t       = typename fe_t::Rd;
    using qp_t      = typename QF::QuadraturePoint;

    const QF &qf(*QF_Simplex<typename fe_t::RdHat>(5));

    double val = 0.;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (domain != Th.get_domain_element(k))
            continue;

        const Cut_Part<e_t> cutK(Th.get_cut_part(k, 0));
        int kb = Th.idxElementInBackMesh(k);

        int kk = k;

        for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

            double meas = cutK.measure(it);
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

                qp_t ip(qf[ipq]);
                v_t mip           = cutK.mapToPhysicalElement(it, ip);
                const double Cint = meas * ip.getWeight();
                double a          = fh->eval(kk, mip) - fex(mip, fh->cu, domain);

                val += Cint * a * a;
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

template <typename Mesh>
double L2normCut_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom, double tt),
                   int domain, const ActiveMesh<Mesh> &Th, double t, const MacroElement<Mesh> *macro) {
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename FElement::Rd Rd;
    typedef typename QF::QuadraturePoint QuadraturePoint;
    typedef typename ActiveMesh<Mesh>::Element Element;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(5));

    double val = 0.;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (domain != Th.get_domain_element(k))
            continue;

        const Cut_Part<Element> cutK(Th.get_cut_part(k, 0));
        int kb = Th.idxElementInBackMesh(k);

        int kk = k;
        // if(macro){
        //   if(!macro->isRootFat(k)) {
        //     kk = macro->getIndexRootElement(k);
        //   }
        // }

        for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
            const R meas = cutK.measure(it);

            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

                QuadraturePoint ip(qf[ipq]); // integration point
                Rd mip       = cutK.mapToPhysicalElement(it, ip);
                const R Cint = meas * ip.getWeight();

                double a = fh->eval(kk, mip) - fex(mip, fh->cu, domain, t);

                val += Cint * a * a;
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

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

template <typename Mesh>
double L2normSurf_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i),
                    const Interface<Mesh> &interface) {
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QFB QFB;
    typedef typename FElement::Rd Rd;
    typedef typename QFB::QuadraturePoint QuadraturePoint;

    const QFB &qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
    What_d Fop = Fwhatd(op_id);

    double val = 0.;

    for (int iface = interface.first_element(); iface < interface.last_element(); iface += interface.next_element()) {

        const int kb = interface.idxElementOfFace(iface); // idx on backMesh
        const R meas = interface.measure(iface);

        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {

            QuadraturePoint ip(qfb[ipq]); // integration point
            const Rd mip = interface.mapToPhysicalFace(iface, (typename FElement::RdHatBord)ip);
            const R Cint = meas * ip.getWeight();

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
template <typename Mesh>
double L2normSurf_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, double t),
                    const Interface<Mesh> &interface, double tt) {
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QFB QFB;
    typedef typename FElement::Rd Rd;
    // typedef GenericInterface<Mesh> Interface;
    typedef typename QFB::QuadraturePoint QuadraturePoint;

    const QFB &qfb(*QF_Simplex<typename FElement::RdHatBord>(5));
    What_d Fop = Fwhatd(op_id);

    double val = 0.;

    for (int iface = interface.first_element(); iface < interface.last_element(); iface += interface.next_element()) {

        const int kb = interface.idxElementOfFace(iface); // idx on backMesh
        const R meas = interface.measure(iface);

        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {

            QuadraturePoint ip(qfb[ipq]); // integration point
            Rd mip       = interface.mapToPhysicalFace(iface,
                                                       (typename FElement::RdHatBord)ip); // mip - map quadrature point
            const R Cint = meas * ip.getWeight();

            // double a = fh->eval(k, mip, tt) - fex(mip, fh->cu, tt);
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
template <typename Mesh>
double L2normSurf(const FunFEM<Mesh> &fh, R(fex)(double *, int i, double t), const Interface<Mesh> &interface,
                  double tt, int c0, int num_comp) {

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2normSurf_2(ui, fex, interface, tt);
    }
    return sqrt(val);
}
template <typename Mesh>
double L2normSurf(const FunFEM<Mesh> &fh, R(fex)(double *, int i), const Interface<Mesh> &interface, int c0,
                  int num_comp) {

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2normSurf_2(ui, fex, interface);
    }
    return sqrt(val);
}
template <typename Mesh>
double L2normSurf(const FunFEM<Mesh> &fh, R(fex)(double *, int i, double t), const Interface<Mesh> *interface,
                  double tt, int c0, int num_comp) {

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2normSurf_2(ui, fex, *interface, tt);
    }
    return sqrt(val);
}
template <typename Mesh>
double L2normSurf(const FunFEM<Mesh> &fh, R(fex)(double *, int i), const Interface<Mesh> *interface, int c0,
                  int num_comp) {

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2normSurf_2(ui, fex, *interface);
    }
    return sqrt(val);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

template <typename M>
double L2norm_2(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i), const M &Th) {
    typedef M Mesh;
    typedef GFESpace<Mesh> FESpace;
    typedef typename Mesh::Element Element;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename FElement::Rd Rd;
    typedef typename QF::QuadraturePoint QuadraturePoint;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(5));
    What_d Fop = Fwhatd(op_id);
    double val = 0.;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        const Element &K(Th[k]);

        const R meas = K.mesure();
        for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

            QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip       = K.mapToPhysicalElement(ip);
            const R Cint = meas * ip.getWeight();
            double a     = fh->eval(k, mip) - fex(mip, fh->cu);

            std::cout << " ----------- " << std::endl;
            std::cout << mip << std::endl;
            std::cout << fh->eval(k, mip) << "\t" << fex(mip, fh->cu) << std::endl;
            val += Cint * a * a;
        }
        getchar();
    }
    double val_receive = 0;
#ifdef USE_MPI
    MPIcf::AllReduce(val, val_receive, MPI_SUM);
#else
    val_receive = val;
#endif

    return val_receive;
}

template <typename M> double L2norm_2(const std::shared_ptr<ExpressionVirtual> &fh, const M &Th) {
    typedef Mesh2 Mesh;
    typedef GFESpace<Mesh> FESpace;
    typedef typename Mesh::Element Element;
    typedef typename FESpace::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename FElement::Rd Rd;
    typedef typename QF::QuadraturePoint QuadraturePoint;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(5));
    What_d Fop = Fwhatd(op_id);
    double val = 0.;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        const Element &K(Th[k]);

        const R meas = K.mesure();
        for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

            QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip       = K.mapToPhysicalElement(ip);
            const R Cint = meas * ip.getWeight();
            double a     = fh->eval(k, mip);
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

template <typename M> double L2norm(const FunFEM<M> &fh, R(fex)(double *, int i), int c0, int num_comp) {

    const GFESpace<M> &Vh(*fh.Vh);
    const M &Th(Vh.Th);

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        auto ui = fh.expr(i);
        val += L2norm_2(ui, fex, Th);
    }
    return sqrt(val);
}

template <typename M>
double L2norm(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i), const M &Th) {

    double val = L2norm_2(fh, fex, Th);

    return sqrt(val);
}

template <typename M> double L2norm(const FunFEM<M> &fh, int c0, int num_comp) {

    const GFESpace<M> &Vh(*fh.Vh);
    const M &Th(Vh.Th);

    double val = 0;
    for (int i = c0; i < num_comp + c0; ++i) {
        ExpressionFunFEM<M> ui(fh, i, op_id);
        val += L2norm_2(ui, Th);
    }
    return sqrt(val);
}

template <typename M> double L2norm(const ExpressionVirtual &fh, const M &Th) {

    double val = L2norm_2(fh, Th);

    return sqrt(val);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

//
// //
// template<typename M>
// double maxNorm(const ExpressionVirtual& fh,R (fex)(const typename
// GFESpace<M>::FElement::Rd, int i),const GFESpace<M>& Vh) {
//
//   typedef M Mesh;
//   typedef GFESpace<Mesh> FESpace;
//   typedef typename FESpace::FElement FElement;
//   typedef typename FElement::QF QF;
//   typedef typename FElement::Rd Rd;
//   typedef typename Mesh::Partition Partition;
//   typedef typename QF::QuadraturePoint QuadraturePoint;
//   typedef typename TypeCutData<Rd::d>::CutData CutData;
//
//   const QF& qf(*QF_Simplex<typename FElement::RdHat>(0));
//   What_d Fop = Fwhatd(op_id);
//   double val = 0.;
//
//   for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
//
//     const FElement& FK(Vh[k]);
//     const R meas = FK.getMeasure();
//     for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//
//       QuadraturePoint ip(qf[ipq]); // integration point
//       Rd mip = FK.map(ip);
//       const R Cint = meas * ip.getWeight();
//
//       val = std::max(val, fabs(fh->eval(k, mip)-fex(mip, fh->cu)));
//     }
//   }
//
//   double val_receive = 0;
//   MPIcf::AllReduce(val, val_receive, MPI_MAX);
// #ifdef USE_MPI
// MPIcf::AllReduce(val, val_receive, MPI_MAX);
// #else
// val_receive = val;
// #endif //

//   return val_receive;
// }
template <typename M> double maxNormCut(const std::shared_ptr<ExpressionVirtual> &fh, const ActiveMesh<M> &Th) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for (int i = 0; i < nb_dom; ++i) {
        val += std::max(val, maxNormCut(fh, Th, i));
    }
    return val;
}
template <typename Mesh>
double maxNormCut(const std::shared_ptr<ExpressionVirtual> &fh, const ActiveMesh<Mesh> &Th, int domain) {

    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename ActiveMesh<Mesh>::Element Element;
    typedef typename FElement::QF QF;
    typedef typename FElement::Rd Rd;
    typedef typename QF::QuadraturePoint QuadraturePoint;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(3));
    What_d Fop = Fwhatd(op_id);

    double val = 0.;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {

        if (domain != Th.get_domain_element(k))
            continue;

        const Cut_Part<Element> cutK(Th.get_cut_part(k, 0));
        int kb = Th.idxElementInBackMesh(k);

        for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

                QuadraturePoint ip(qf[ipq]); // integration point
                Rd mip = cutK.mapToPhysicalElement(it, ip);

                val = std::max(val, fabs(fh->eval(k, mip)));
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

template <typeMesh mesh_t, FunctionDomain fct_t>
double maxNormCut(const std::shared_ptr<ExpressionVirtual> &fh, fct_t fex, const ActiveMesh<mesh_t> &Th) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for (int i = 0; i < nb_dom; ++i) {
        val += std::max(val, maxNormCut(fh, fex, Th, i));
    }
    return val;
}

template <typeMesh mesh_t, FunctionDomain fct_t>
double maxNormCut(const std::shared_ptr<ExpressionVirtual> &fh, fct_t fex, const ActiveMesh<mesh_t> &Th, int domain) {

    using fespace_t = GFESpace<mesh_t>;
    using fe_t      = typename fespace_t::FElement;
    using e_t       = typename mesh_t::Element;
    using QF        = typename fe_t::QF;
    using QFB       = typename fe_t::QFB;
    using v_t       = typename fe_t::Rd;
    using qp_t      = typename QF::QuadraturePoint;

    const QF &qf(*QF_Simplex<typename fe_t::RdHat>(5));
    const QFB &qfb(*QF_Simplex<typename fe_t::RdHatBord>(5));

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
                val     = std::max(val, fabs(fh->eval(k, mip) - fex(mip, fh->cu, domain)));
            }

            for (int ifac = 0; ifac < e_t::nea; ++ifac) {
                for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
                    auto ip(qfb[ipq]); // integration point
                    auto ipf = K.mapToReferenceElement(ip, ifac);
                    v_t mip  = cutK.mapToPhysicalElement(it, ipf);
                    val      = std::max(val, fabs(fh->eval(k, mip) - fex(mip, fh->cu, domain)));
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

template <typename M>
double maxNormCut(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom),
                  const ActiveMesh<M> &Th, const std::vector<R2> &sample_node) {
    int nb_dom = Th.get_nb_domain();
    double val = 0.;
    for (int i = 0; i < nb_dom; ++i) {
        val += std::max(val, maxNormCut(fh, fex, Th, i, sample_node));
    }
    return val;
}

template <typename Mesh>
double maxNormCut(const std::shared_ptr<ExpressionVirtual> &fh, R(fex)(double *, int i, int dom),
                  const ActiveMesh<Mesh> &cutTh, int domain, const std::vector<R2> &sample_node) {

    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename ActiveMesh<Mesh>::Element Element;
    typedef typename FElement::Rd Rd;

    What_d Fop = Fwhatd(op_id);
    double val = 0.;
    int k0     = 0;

    const auto &Th(cutTh.Th);

    progress bar(" Max noorm in Omega_i", sample_node.size(), globalVariable::verbose);
    int ii = 0;
    ;
    for (auto mip : sample_node) {
        bar++;
        int kb = geometry::find_triangle_contenant_p(Th, mip, k0);
        k0     = kb;
        int k  = cutTh.idxElementFromBackMesh(kb, domain);
        if (k == -1)
            continue;
        const auto cutK(cutTh.get_cut_part(k, 0));
        for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {

            if (!geometry::p_dans_triangle(mip, cutK.get_vertex(it, 0), cutK.get_vertex(it, 1), cutK.get_vertex(it, 2)))
                continue;
            val = std::max(val, fabs(fh->eval(k, mip) - fex(mip, fh->cu, domain)));
        }
    }
    bar.end();
    double val_receive = 0;
#ifdef USE_MPI
    MPIcf::AllReduce(val, val_receive, MPI_MAX);
#else
    val_receive = val;
#endif
    return val_receive;
}

template <typename Mesh> double maxNorm(const std::shared_ptr<ExpressionVirtual> &fh, const Mesh &Th) {

    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename ActiveMesh<Mesh>::Element Element;
    typedef typename FElement::QF QF;
    typedef typename FElement::Rd Rd;
    typedef typename QF::QuadraturePoint QuadraturePoint;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(0));
    What_d Fop = Fwhatd(op_id);

    double val = 0.;

    for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
        const Element &K(Th[k]);

        for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {

            QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip = K.mapToPhysicalElement(ip);

            val = std::max(val, fabs(fh->eval(k, mip)));
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

/*---------------------------- cut local L2 norm ----------------------------
Integrates only over elements whose edges are not part of the stabilised edges,
thereby neglecting elements messed up by the mcdonald stab
*/
// template<typename M>
// double L2normCutLoc_2(const ExpressionVirtual& fh,R (fex)(const typename
// GFESpace<M>::FElement::Rd, int i, int dom),int domain, const GFESpace<M>& Vh,
// const MacroElement<M>& bigMac ) {
//   typedef Mesh2 Mesh;
//   typedef GFESpace<Mesh> FESpace;
//   typedef typename FESpace::FElement FElement;
//   typedef typename FElement::QF QF;
//   typedef typename FElement::Rd Rd;
//   typedef typename Mesh::Partition Partition;
//   typedef typename Mesh::Element Element; // [Needed for neighbor element]
//   typedef typename QF::QuadraturePoint QuadraturePoint;
//   typedef typename TypeCutData<Rd::d>::CutData CutData;
//
//
//   const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
//   What_d Fop = Fwhatd(op_id);
//   double val = 0.;
//
//   vector<int> elements_to_integrate(Vh.NbElement(),1);
//   // std::vector<std::pair<int,int>>::iterator it;
//   // 0. Mark all elements as to integrate
//   // 1. Loop over all elements
//   // 2. Remove from elements_to_integrate when isSmall or neighbor
//   // 3. Standard L2 error, but only over elements_to_integrate
//
//
// assert(0);
//   // for (auto it=bigMac.edges_to_stabilize.begin();
//   it<bigMac.edges_to_stabilize.end(); it++) {
//   //
//   //   int k = it->first;
//   //   int ifac = it->second;
//   //   //
//   //   int ifacn = ifac;
//   //   // for(int ifac=0;ifac<3;++ifac){
//   //     // int ifacn = ifac;
//   //
//   //     int k_back  = Vh.Th(Vh[k].T);
//   //     int kn_back = Vh.Th.ElementAdj(k_back,ifacn);
//   //     if(kn_back == -1) kn_back = k_back;
//   //     int kn = Vh.idxElementFromBackMesh(kn_back, domain);
//   //     if(kn == -1) kn = k;
//   //     elements_to_integrate[k] = 0;
//   // }
//
//
//
//   for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
//     const FElement& FK(Vh[k]);
//     if(domain != FK.whichDomain()) continue;
//
//     if(elements_to_integrate[k]==0) continue;
//
//     const int kb = Vh.Th(FK.T);
//     CutData cutData(Vh.getInterface(0).getCutData(kb));     // get the cut
//     data const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
//     ElementSignEnum the_part = cutK.what_part(domain);
//     double locV = 0;
//
//     for(typename Partition::const_element_iterator it =
//     cutK.element_begin(the_part); it != cutK.element_end(the_part); ++it){
//
//       const R meas = cutK.mesure(it);
//       for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//
//         QuadraturePoint ip(qf[ipq]); // integration point
//         Rd mip = cutK.toK(it, ip);
//         const R Cint = meas * ip.getWeight();
//
//         double a = fh->eval(k, mip) - fex(mip, fh->cu, domain);
//
//         val += Cint * a * a;
//         locV += Cint * a * a;
//       }
//     }
//     // std::cout << k << "\t" << locV << "\t"<< val << std::endl;
//   }
//   double val_receive = 0;
// #ifdef USE_MPI
// MPIcf::AllReduce(val, val_receive, MPI_SUM);
// #else
// val_receive = val;
// #endif
//
//   return val_receive;
// }
//
// template<typename M>
// double L2normCutLoc_2(const ExpressionVirtual& fh,R (fex)(const typename
// GFESpace<M>::FElement::Rd, int i, int dom), const GFESpace<M>& Vh, const
// MacroElement<M>& bigMac) {
//   return L2normCutLoc_2(fh,fex,0,Vh,bigMac) +
//   L2normCutLoc_2(fh,fex,1,Vh,bigMac);
// }
//
// // Several components, ie vector valued L2norm
// template<typename M>
// double L2normCutLoc( const FunFEM<M>& fh,R (fex)(const typename
// GFESpace<M>::FElement::Rd, int i, int dom), int c0, int num_comp , const
// MacroElement<M>& bigMac) {
//
//   const GFESpace<M>& Vh(*fh.Vh);
//
//   double val = 0;
//   for(int i=c0;i<num_comp+c0;++i) {
//     ExpressionFunFEM<M> ui(fh, i, op_id);
//     val += L2normCutLoc_2(ui,fex,Vh,bigMac);
//   }
//   return sqrt(val);
// }
//
// // Only one component
// template<typename M>
// double L2normCutLoc( const ExpressionVirtual& fh, R (fex)(const typename
// GFESpace<M>::FElement::Rd, int i, int dom), const GFESpace<M>& Vh, const
// MacroElement<M>& bigMac) {
//   double val = L2normCutLoc_2(fh,fex,Vh,bigMac);
//   return sqrt(val);
// }
