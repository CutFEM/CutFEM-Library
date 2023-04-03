/**
 * @file cpp/FESpace/limiter.hpp
 * @brief
 *
 * @copyright

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

#ifndef CUTFEM_LIMITER_HPP
#define CUTFEM_LIMITER_HPP

#include "FESpace.hpp"
#include "macroElement.hpp"
#include "../num/print_container.hpp"
#include <tuple>

namespace limiter {

inline KNM<double> inv(const KNM<double> &M) {

    KNM<double> R(M.N(), M.M());
    if (M.N() == 2) {
        double Delta = M(0, 0) * M(1, 1) - M(1, 0) * M(0, 1);
        R(0, 0)      = M(1, 1) / Delta;
        R(1, 1)      = M(0, 0) / Delta;
        R(0, 1)      = -M(0, 1) / Delta;
        R(1, 0)      = -M(1, 0) / Delta;
        return R;
    }

    if (M.N() == 3) {
        double Delta = M(0, 0) * (M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2)) -
                       M(0, 1) * (M(1, 0) * M(2, 2) - M(2, 0) * M(1, 2)) +
                       M(0, 2) * (M(1, 0) * M(2, 1) - M(2, 0) * M(1, 1));
        R(0, 0) = (M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2)) / Delta;
        R(1, 0) = -(M(1, 0) * M(2, 2) - M(2, 0) * M(1, 2)) / Delta;
        R(2, 0) = (M(1, 0) * M(2, 1) - M(2, 0) * M(1, 1)) / Delta;
        R(0, 1) = -(M(0, 1) * M(2, 2) - M(2, 1) * M(0, 2)) / Delta;
        R(1, 1) = (M(0, 0) * M(2, 2) - M(2, 0) * M(0, 2)) / Delta;
        R(2, 1) = -(M(0, 0) * M(2, 1) - M(2, 0) * M(0, 1)) / Delta;
        R(0, 2) = (M(0, 1) * M(1, 2) - M(1, 1) * M(0, 2)) / Delta;
        R(1, 2) = -(M(0, 0) * M(1, 2) - M(1, 0) * M(0, 2)) / Delta;
        R(2, 2) = (M(0, 0) * M(1, 1) - M(1, 0) * M(0, 1)) / Delta;
        return R;
    }

    if (M.N() > 3) {
        std::cout << "matrix is too big, cannot be inverted" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return R;
}

// check mean value to satisfy maximuim principle
template <typename Mesh> void check_maximum_principle(std::map<int, double> &u_mean, double min_u, double max_u) {

    for (auto &p : u_mean) {
        double val = p.second;
        if (min_u > val || max_u < val) {
            std::cout << "element \t" << p.first << "\t" << min_u << "\t" << val << "\t" << max_u << std::endl;
        }
    }
}

namespace CutFEM {

// -------------------------------------------------------------------------
// COMPUTE THE MIN AND THE MAX VALUE OF A FE FUNCTION ON A CUT MESH
// DEPENDING OF THE POLYNOMIAL ORDER.
// BY DEFAULT IT USES LINEAR ELEMENTS

template <typename Mesh> std::tuple<double, double> findMinAndMaxValue_P0(const FunFEM<Mesh> &uh) {
    typedef typename Mesh::Rd Rd;
    typedef typename Mesh::Element Element;
    typedef typename GFESpace<Mesh>::FElement FElement;
    const GFESpace<Mesh> &Wh(*uh.Vh);
    const ActiveMesh<Mesh> &Kh(Wh.get_mesh());

    double min_val, max_val;
    double min_val_loc = 1e300;
    double max_val_loc = -1e300;
    int k_min, k_max;

    for (int k = Wh.first_element(); k < Wh.last_element(); k += Wh.next_element()) {
        const FElement &FK(Wh[k]);

        Rd mip     = FK.T.barycenter();
        double val = uh.eval(k, mip);
        if (val > max_val_loc) {
            max_val_loc = val;
            k_max       = k;
        }
        if (val < min_val_loc) {
            min_val_loc = val;
            k_min       = k;
        }
    }
#ifdef USE_MPI
    MPIcf::AllReduce(min_val_loc, min_val, MPI_MIN);
    MPIcf::AllReduce(max_val_loc, max_val, MPI_MAX);
#else
    min_val = min_val_loc;
    max_val = max_val_loc;
#endif
    return {min_val, max_val};
}

template <typename Mesh> std::tuple<double, double> findMinAndMaxValue_P1(const FunFEM<Mesh> &uh) {
    typedef typename Mesh::Rd Rd;
    typedef typename Mesh::Element Element;
    typedef typename GFESpace<Mesh>::FElement FElement;
    int nv_loc = Mesh::Rd::d + 1;
    const GFESpace<Mesh> &Wh(*uh.Vh);
    const ActiveMesh<Mesh> &Kh(Wh.get_mesh());

    double min_val, max_val;
    double min_val_loc = 1e300;
    double max_val_loc = -1e300;
    // Apply limiter in all element
    int k_min, k_max;
    // for(int k=Wh.first_element(); k<Wh.last_element();++k) {
    for (int k = Wh.first_element(); k < Wh.last_element(); k += Wh.next_element()) {
        const FElement &FK(Wh[k]);
        const Cut_Part<Element> cutK(Kh.get_cut_part(k, 0));

        for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
            // get the local std::min and std::max on K
            for (int ipq = 0; ipq < nv_loc; ++ipq) {
                Rd mip     = cutK.get_vertex(it, ipq);
                // double val = uh.eval(k, mip);
                double val = uh.eval(k, mip);
                if (val > max_val_loc) {
                    max_val_loc = val;
                    k_max       = k;
                }
                if (val < min_val_loc) {
                    min_val_loc = val;
                    k_min       = k;
                }
                // max_val  = std::max(max_val, val);
                // min_val  = std::min(min_val, val);
            }
        }
    }
#ifdef USE_MPI
    MPIcf::AllReduce(min_val_loc, min_val, MPI_MIN);
    MPIcf::AllReduce(max_val_loc, max_val, MPI_MAX);
#else
    min_val = min_val_loc;
    max_val = max_val_loc;
#endif

    return {min_val, max_val};
}

template <typename Mesh> std::tuple<double, double> findMinAndMaxValue(const FunFEM<Mesh> &uh) {
    const int fctOrder = uh.getPolynomialOrder();
    if (fctOrder == 0) {
        return findMinAndMaxValue_P0(uh);
    } else if (fctOrder == 1) {
        return findMinAndMaxValue_P1(uh);
    } else {
        std::cout << " Min and max compute using linear elements. Not exact for "
                     "this polynomial order"
                  << std::endl;
        return findMinAndMaxValue_P1(uh);
    }
}

// -----------------------------------------------------------------------------
template <typename Mesh>
std::map<int, double> computeMeanValue(const FunFEM<Mesh> &uh, const MacroElement<Mesh> &macro) {

    using QF      = typename GFElement<Mesh>::QF;
    using CutMesh = ActiveMesh<Mesh>;
    const QF &qf(*QF_Simplex<typename GFElement<Mesh>::RdHat>(3));

    std::map<int, double> u_mean;
    const auto &Vh(*uh.Vh);
    const CutMesh &Th(Vh.get_mesh());

    // 1) Compute on all elements
    // for (int k = Vh.first_element(); k < Vh.last_element(); k +=
    // Vh.next_element()) {

    //     const auto cutK = Th.get_cut_part(k, 0);
    //     double meas_K   = 0.;
    //     double mean_val = 0.;
    //     for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it)
    //     {
    //         double meas_cut = cutK.measure(it);
    //         meas_K += meas_cut;

    //         for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
    //             auto ip  = qf[ipq]; // integration point
    //             auto mip = cutK.mapToPhysicalElement(it, ip);
    //             mean_val += meas_cut * ip.getWeight() * uh.eval(k, mip);
    //         }
    //     }
    //     u_mean[k] = mean_val / meas_K;
    // }

    // 2) Compute on the macro elements
    for (auto &p : macro.macro_element) {

        const auto &MK(p.second);
        double mean_val = 0.;

        for (auto k : MK.idx_element) {
            const auto FK   = Vh[k];
            const auto cutK = macro.Th_.get_cut_part(k, 0);

            for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
                double meas = cutK.measure(it);

                for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                    auto ip  = qf[ipq];
                    auto mip = cutK.mapToPhysicalElement(it, ip);
                    mean_val += meas * ip.getWeight() * uh.eval(k, mip);
                }
            }
        }
        mean_val                    = mean_val / MK.area_total_;
        u_mean[MK.idx_root_element] = mean_val;
    }
    return u_mean;
}

// Only for discontinuous Lagrange
template <typename Mesh> std::vector<double> extendToMacro_P0(const FunFEM<Mesh> &uh, const MacroElement<Mesh> &macro) {
    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename Mesh::Rd Rd;
    std::vector<double> u_new(uh.array().begin(), uh.array().end());
    const GFESpace<Mesh> &Wh(*uh.Vh);
    for (auto &p : macro.macro_element) {

        // compute the value of the dof of the elements in the
        // macro element.
        // u_{h,M} = {sum_{K in M} |K \cap \Omega|*u_{h}^K }{sum_{K in M} |K
        // \cap \Omega|}
        const MElement &MK(p.second);
        int n_element = MK.size();

        // loop over the elements that has to be changed
        // assert(0);
        // need to fix this with u_mean
        for (auto s : MK.idx_element) {
            const FElement &FK(Wh[s]);
            assert(FK.NbDoF() == FK.NbNode());
            double area_total = MK.area_total_;

            // loop over the dof (in that case just node evaluation)
            for (int df = 0; df < FK.NbNode(); ++df) {
                Rd P       = FK.Pt(df);
                int df_glb = FK.loc2glb(df);
                double val = 0;

                for (auto k : MK.idx_element) {
                    double s = macro.get_area(k);
                    val += s * uh.eval(k, P);
                }
                val           = val / area_total;
                u_new[df_glb] = val;
            }
        }
    }
    return u_new;
}

// Only for discontinuous Lagrange
template <typename Mesh>
std::vector<double> extendToMacro_P1(const FunFEM<Mesh> &uh, const std::map<int, double> &u_mean,
                                     const MacroElement<Mesh> &macro) {
    typedef typename Mesh::Element Element;
    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename Mesh::Rd Rd;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(3));
    const GFESpace<Mesh> &Wh(*uh.Vh);

    std::vector<double> u_new(uh.array().begin(), uh.array().end());

    FunFEM<Mesh> fun_u_new(Wh, u_new);

    for (auto &p : macro.macro_element) {

        // compute the value of the dof of the elements in the
        // macro element.
        // u_{h,M} = {sum_{K in M} |K \cap \Omega|*u_{h}^K }{sum_{K in M} |K
        // \cap \Omega|}
        const MElement &MK(p.second);
        int n_element         = MK.size();
        double area_total     = MK.area_total_;
        int idx_root          = MK.idx_root_element;
        double mean_val_macro = u_mean.find(idx_root)->second;

        // loop over the elements that has to be changed
        for (auto s : MK.idx_element) {
            const FElement &FK(Wh[s]);

            // loop over the dof (in that case just node evaluation)
            for (int df = 0; df < FK.NbDoF(); ++df) {
                Rd P       = FK.Pt(df);
                int df_glb = FK.loc2glb(df);
                double val = 0;
                for (auto k : MK.idx_element) {
                    double ss = macro.get_area(k);
                    val += ss * uh.eval(k, P);
                }
                val           = val / area_total;
                u_new[df_glb] = val;
            }
        }
        double C0 = 0;
        for (auto k : MK.idx_element) {
            const FElement &FK(Wh[k]);
            const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k, 0));
            for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
                double meas = cutK.measure(it);

                for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                    typename QF::QuadraturePoint ip(qf[ipq]); // integration point
                    Rd mip = cutK.mapToPhysicalElement(it, ip);
                    C0 -= meas * ip.getWeight() * fun_u_new.eval(k, mip);

                    // double Cint    = meas * ip.getWeight();
                    // double val_old = uh.eval(k, mip);
                    // C0 += Cint * (val_old - fun_u_new.eval(k, mip));
                }
            }
        }
        C0 /= area_total;
        C0 += mean_val_macro;

        //                   Rd mip         = cutK.mapToPhysicalElement(it, ip);
        //            double Cint    = meas * ip.getWeight();
        //            double val_old = uh.eval(k, mip);
        //            mean_val += Cint * val_old;
        //            C0 += Cint * (val_old - fun_u_new.eval(k, mip));

        //   C0               = C0 / area_total;
        //   mean_val         = mean_val / area_total;

        for (auto k : MK.idx_element) {
            const FElement &FK(Wh[k]);
            for (int df = 0; df < FK.NbDoF(); ++df) {
                int df_glb = FK.loc2glb(df);
                u_new[df_glb] += C0;
            }
        }
    }
    return u_new;
}

template <typename Mesh>
std::vector<double> extendToMacro(const FunFEM<Mesh> &uh, const std::map<int, double> &u_mean,
                                  const MacroElement<Mesh> &macro) {
    const BasisFctType basisFctType = uh.getBasisFctType();
    if (basisFctType == BasisFctType::P0) {
        return extendToMacro_P0(uh, macro);
    } else if (basisFctType == BasisFctType::P1dc) {
        return extendToMacro_P1(uh, u_mean, macro);
    } else {
        std::cout << " Extension to macro element not implemented for those element" << std::endl;
        assert(0);
        return std::vector<double>();
    }
}

template <typename Mesh>
void extendToMacro_P0(const FunFEM<Mesh> &uh, Rn &u_new, std::map<int, double> &u_mean,
                      const MacroElement<Mesh> &macro) {
    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename Mesh::Rd Rd;
    u_new = uh.v;
    const GFESpace<Mesh> &Wh(*uh.Vh);
    for (auto &p : macro.macro_element) {

        // compute the value of the dof of the elements in the
        // macro element.
        // u_{h,M} = {sum_{K in M} |K \cap \Omega|*u_{h}^K }{sum_{K in M} |K
        // \cap \Omega|}
        const MElement &MK(p.second);
        int n_element = MK.size();

        // loop over the elements that has to be changed
        // assert(0);
        // need to fix this with u_mean
        for (auto s : MK.idx_element) {
            const FElement &FK(Wh[s]);
            assert(FK.NbDoF() == FK.NbNode());
            double area_total = MK.area_total_;

            // loop over the dof (in that case just node evaluation)
            for (int df = 0; df < FK.NbNode(); ++df) {
                Rd P       = FK.Pt(df);
                int df_glb = FK.loc2glb(df);
                double val = 0;

                for (auto k : MK.idx_element) {
                    double s = macro.get_area(k);
                    val += s * uh.eval(k, P);
                }
                val           = val / area_total;
                u_new[df_glb] = val;
            }
        }
    }
    // return u_new;
}

// Only for discontinuous Lagrange
template <typename Mesh>
void extendToMacro_P1(const FunFEM<Mesh> &uh, Rn &u_new, std::map<int, double> &u_mean,
                      const MacroElement<Mesh> &macro) {
    typedef typename Mesh::Element Element;
    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename Mesh::Rd Rd;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(3));
    const GFESpace<Mesh> &Wh(*uh.Vh);

    // std::vector<double> u_new(uh.array().begin(), uh.array().end());

    FunFEM<Mesh> fun_u_new(Wh, u_new);

    for (auto &p : macro.macro_element) {

        // compute the value of the dof of the elements in the
        // macro element.
        // u_{h,M} = {sum_{K in M} |K \cap \Omega|*u_{h}^K }{sum_{K in M} |K
        // \cap \Omega|}
        const MElement &MK(p.second);
        int n_element         = MK.size();
        double area_total     = MK.area_total_;
        int idx_root          = MK.idx_root_element;
        double mean_val_macro = u_mean.find(idx_root)->second;

        // loop over the elements that has to be changed
        for (auto s : MK.idx_element) {
            const FElement &FK(Wh[s]);

            // loop over the dof (in that case just node evaluation)
            for (int df = 0; df < FK.NbDoF(); ++df) {
                Rd P       = FK.Pt(df);
                int df_glb = FK.loc2glb(df);
                double val = 0;
                for (auto k : MK.idx_element) {
                    double ss = macro.get_area(k);
                    val += ss * uh.eval(k, P);
                }
                val           = val / area_total;
                u_new[df_glb] = val;
            }
        }
        double C0       = 0;
        double mean_val = 0.;
        for (auto k : MK.idx_element) {
            const FElement &FK(Wh[k]);
            const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k, 0));
            for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
                double meas = cutK.measure(it);

                for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                    typename QF::QuadraturePoint ip(qf[ipq]); // integration point
                    Rd mip = cutK.mapToPhysicalElement(it, ip);
                    // C0 -= meas * ip.getWeight() * fun_u_new.eval(k, mip);

                    double Cint    = meas * ip.getWeight();
                    double val_old = uh.eval(k, mip);
                    mean_val += Cint * val_old;
                    C0 += Cint * (val_old - fun_u_new.eval(k, mip));
                }
            }
        }
        C0 /= area_total;
        C0 += mean_val_macro;

        //                   Rd mip         = cutK.mapToPhysicalElement(it, ip);
        //            double Cint    = meas * ip.getWeight();
        //            double val_old = uh.eval(k, mip);
        //            mean_val += Cint * val_old;
        //            C0 += Cint * (val_old - fun_u_new.eval(k, mip));

        //   C0               = C0 / area_total;
        mean_val         = mean_val / area_total;
        u_mean[idx_root] = mean_val;
        for (auto k : MK.idx_element) {
            const FElement &FK(Wh[k]);
            for (int df = 0; df < FK.NbDoF(); ++df) {
                int df_glb = FK.loc2glb(df);
                u_new[df_glb] += C0;
            }
        }
    }
}

template <typename Mesh>
void extendToMacro(const FunFEM<Mesh> &uh, Rn &u_new, std::map<int, double> &u_mean, const MacroElement<Mesh> &macro) {
    const BasisFctType basisFctType = uh.getBasisFctType();
    if (basisFctType == BasisFctType::P0) {
        extendToMacro_P0(uh, u_new, u_mean, macro);
    } else if (basisFctType == BasisFctType::P1dc) {
        extendToMacro_P1(uh, u_new, u_mean, macro);
    } else {
        std::cout << " Extension to macro element not implemented for those element" << std::endl;
    }
}

template <typename Mesh> void extendToMacro(const FunFEM<Mesh> &uh, Rn &u_new, const MacroElement<Mesh> &macro) {
    const BasisFctType basisFctType = uh.getBasisFctType();
    std::map<int, double> u_mean;
    if (basisFctType == BasisFctType::P0) {
        extendToMacro_P0(uh, u_new, u_mean, macro);
    } else if (basisFctType == BasisFctType::P1dc) {
        extendToMacro_P1(uh, u_new, u_mean, macro);
    } else {
        std::cout << " Extension to macro element not implemented for those element" << std::endl;
    }
}

template <typename Mesh>
void boundPreservingLimiter_P1(const FunFEM<Mesh> &uh, Rn &u_new, double min_val, double max_val,
                               //   const std::map<int, double> map_mean_value,
                               const MacroElement<Mesh> &macro) {
    typedef typename Mesh::Element Element;
    typedef typename Mesh::Rd Rd;
    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename FElement::QF QF;
    const QF &qf(*QF_Simplex<typename FElement::RdHat>(2));

    int nv_loc = Mesh::Rd::d + 1;
    const GFESpace<Mesh> &Wh(*uh.Vh);
    const ActiveMesh<Mesh> &Kh(macro.Th_);

    // const auto &uM = fun_uM.array();
    // std::vector<double> u_new(uM.begin(), uM.end());
    u_new = uh.v;
    std::map<int, double> map_mean_value;
    Rn uM(u_new);
    // std::vector<double> uMv(u_new.begin(), u_new.end());
    //=
    // Rn_ uM(uMv);
    extendToMacro(uh, uM, map_mean_value, macro);

    FunFEM<Mesh> fun_uM(Wh, uM);

    u_new = uM;

    // Apply limiter in all element
    for (int k = Wh.first_element(); k < Wh.last_element(); ++k) {

        const FElement &FK(Wh[k]);
        const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k, 0));
        double u_bar_K = 0.;
        double max_K = -1e300, min_K = 1e300;
        double areaCut = 0.;
        for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
            double meas = cutK.measure(it);
            areaCut += meas;
            // get the mean value on K
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                typename QF::QuadraturePoint ip(qf[ipq]); // integration point
                Rd mip      = cutK.mapToPhysicalElement(it, ip);
                double Cint = meas * ip.getWeight();
                // u_bar_K += Cint*uh.eval(k, mip);
                u_bar_K += Cint * fun_uM.eval(k, mip);
            }

            // get the local std::min and max on K
            for (int ipq = 0; ipq < nv_loc; ++ipq) {
                Rd mip     = cutK.get_vertex(it, ipq);
                // double val = uh.eval(k, mip);
                double val = fun_uM.eval(k, mip);

                max_K = std::max(max_K, val);
                min_K = std::min(min_K, val);
            }
        }
        u_bar_K = u_bar_K / areaCut;

        // 3) compute theta
        double v1 = fabs((min_val - u_bar_K) / (min_K - u_bar_K));
        if ((min_K - u_bar_K) == 0.)
            v1 = 1.; //
        double v2 = fabs((max_val - u_bar_K) / (max_K - u_bar_K));
        if ((max_K - u_bar_K) == 0.)
            v2 = 1.;
        double theta = std::min(std::min(v1, v2), 1.);
        // if( fabs(theta) < Epsilon) theta = 0.;

        // 4) replace the dof
        for (int df = 0; df < FK.NbDoF(); ++df) {
            int df_glb    = FK.loc2glb(df);
            u_new[df_glb] = theta * (uM(df_glb) - u_bar_K) + u_bar_K;
        }
    }

    // loop over macro element
    for (auto &p : macro.macro_element) {
        const MElement &MK(p.second);
        int idx_root = MK.idx_root_element;

        // 1) the mean value
        double u_bar_M = map_mean_value.find(idx_root)->second;
        // 2) get local min and std::max
        double max_M = -1e300, min_M = 1e300;
        for (auto k : MK.idx_element) {
            const FElement &FK(Wh[k]);
            const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k, 0));
            for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
                for (int ipq = 0; ipq < nv_loc; ++ipq) {
                    Rd mip     = cutK.get_vertex(it, ipq);
                    double val = fun_uM.eval(k, mip);
                    max_M      = std::max(max_M, val);
                    min_M      = std::min(min_M, val);
                }
            }
        }

        // 3) compute theta
        double v1    = fabs((min_val - u_bar_M) / (min_M - u_bar_M - globalVariable::Epsilon));
        double v2    = fabs((max_val - u_bar_M) / (max_M - u_bar_M + globalVariable::Epsilon));
        double theta = std::min(std::min(v1, v2), 1.);
        // std::cout << theta << std::endl;
        if (fabs(theta) < globalVariable::Epsilon)
            theta = 0.;
        // 4) replace the dof
        for (auto k : MK.idx_element) {
            const FElement &FK(Wh[k]);
            for (int df = 0; df < FK.NbDoF(); ++df) {
                int df_glb    = FK.loc2glb(df);
                u_new[df_glb] = theta * (uM[df_glb] - u_bar_M) + u_bar_M;
            }
        }
    }

    // return u_new;
}

// -----------------------------------------------------------------------------
template <typename Mesh>
std::vector<double> boundPreservingLimiter_P1(const FunFEM<Mesh> &fun_uM, double min_val, double max_val,
                                              const std::map<int, double> map_mean_value,
                                              const MacroElement<Mesh> &macro) {
    typedef typename Mesh::Element Element;
    typedef typename Mesh::Rd Rd;
    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename FElement::QF QF;
    const QF &qf(*QF_Simplex<typename FElement::RdHat>(2));

    int nv_loc = Mesh::Rd::d + 1;
    const GFESpace<Mesh> &Wh(*fun_uM.Vh);
    const ActiveMesh<Mesh> &Kh(macro.Th_);

    const auto &uM = fun_uM.array();
    std::vector<double> u_new(uM.begin(), uM.end());

    // Apply limiter in all element
    for (int k = Wh.first_element(); k < Wh.last_element(); ++k) {

        const FElement &FK(Wh[k]);
        const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k, 0));
        double u_bar_K = 0.;
        double max_K = -1e300, min_K = 1e300;
        double areaCut = 0.;
        for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
            double meas = cutK.measure(it);
            areaCut += meas;
            // get the mean value on K
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                typename QF::QuadraturePoint ip(qf[ipq]); // integration point
                Rd mip      = cutK.mapToPhysicalElement(it, ip);
                double Cint = meas * ip.getWeight();
                // u_bar_K += Cint*uh.eval(k, mip);
                u_bar_K += Cint * fun_uM.eval(k, mip);
            }

            // get the local std::min and max on K
            for (int ipq = 0; ipq < nv_loc; ++ipq) {
                Rd mip     = cutK.get_vertex(it, ipq);
                // double val = uh.eval(k, mip);
                double val = fun_uM.eval(k, mip);

                max_K = std::max(max_K, val);
                min_K = std::min(min_K, val);
            }
        }
        u_bar_K = u_bar_K / areaCut;

        // 3) compute theta
        double v1 = fabs((min_val - u_bar_K) / (min_K - u_bar_K));
        if ((min_K - u_bar_K) == 0.)
            v1 = 1.; //
        double v2 = fabs((max_val - u_bar_K) / (max_K - u_bar_K));
        if ((max_K - u_bar_K) == 0.)
            v2 = 1.;
        double theta = std::min(std::min(v1, v2), 1.);
        // if( fabs(theta) < Epsilon) theta = 0.;

        // 4) replace the dof
        for (int df = 0; df < FK.NbDoF(); ++df) {
            int df_glb    = FK.loc2glb(df);
            u_new[df_glb] = theta * (uM(df_glb) - u_bar_K) + u_bar_K;
        }
    }

    // loop over macro element
    for (auto &p : macro.macro_element) {
        const MElement &MK(p.second);
        int idx_root = MK.idx_root_element;

        // 1) the mean value
        double u_bar_M = map_mean_value.find(idx_root)->second;
        // 2) get local min and std::max
        double max_M = -1e300, min_M = 1e300;
        for (auto k : MK.idx_element) {
            const FElement &FK(Wh[k]);
            const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k, 0));
            for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
                for (int ipq = 0; ipq < nv_loc; ++ipq) {
                    Rd mip     = cutK.get_vertex(it, ipq);
                    double val = fun_uM.eval(k, mip);
                    max_M      = std::max(max_M, val);
                    min_M      = std::min(min_M, val);
                }
            }
        }

        // 3) compute theta
        double v1    = fabs((min_val - u_bar_M) / (min_M - u_bar_M - globalVariable::Epsilon));
        double v2    = fabs((max_val - u_bar_M) / (max_M - u_bar_M + globalVariable::Epsilon));
        double theta = std::min(std::min(v1, v2), 1.);
        // std::cout << theta << std::endl;
        if (fabs(theta) < globalVariable::Epsilon)
            theta = 0.;
        // 4) replace the dof
        for (auto k : MK.idx_element) {
            const FElement &FK(Wh[k]);
            for (int df = 0; df < FK.NbDoF(); ++df) {
                int df_glb    = FK.loc2glb(df);
                u_new[df_glb] = theta * (uM[df_glb] - u_bar_M) + u_bar_M;
            }
        }
    }

    return u_new;
}

template <typename Mesh>
std::vector<double> applyBoundPreservingLimiter(const FunFEM<Mesh> &uh, double min_val, double max_val,
                                                const std::map<int, double> u_mean, const MacroElement<Mesh> &macro) {
    if (uh.getBasisFctType() == BasisFctType::P1dc) {
        return boundPreservingLimiter_P1(uh, min_val, max_val, u_mean, macro);
    } else {
        std::cout << " No boundary preserving limiter implemented for this elements" << std::endl;
        exit(EXIT_FAILURE);
        return std::vector<double>();
    }
}

template <typename Mesh>
std::vector<double> applyBoundPreservingLimiter(const FunFEM<Mesh> &uh, double min_val, double max_val,
                                                const MacroElement<Mesh> &macro) {
    const std::map<int, double> u_mean = limiter::CutFEM::computeMeanValue(uh, macro);
    std::vector<double> uh_M           = limiter::CutFEM::extendToMacro(uh, u_mean, macro);
    FunFEM<Mesh> fun_uM(*uh.Vh, uh_M);

    if (uh.getBasisFctType() == BasisFctType::P1dc) {
        return boundPreservingLimiter_P1(fun_uM, min_val, max_val, u_mean, macro);
    } else {
        std::cout << " No boundary preserving limiter implemented for this elements" << std::endl;
        exit(EXIT_FAILURE);
        return std::vector<double>();
    }
}

template <typename Mesh>
void applyBoundPreservingLimiter(const FunFEM<Mesh> &uh, Rn &u_new, double min_val, double max_val,
                                 const MacroElement<Mesh> &macro) {
    if (uh.getBasisFctType() == BasisFctType::P1dc) {
        boundPreservingLimiter_P1(uh, u_new, min_val, max_val, macro);
    } else {
        std::cout << " No boundary preserving limiter implemented for this elements" << std::endl;
    }
}

template <typename Mesh>
void check_mean_value(const FunFEM<Mesh> &uh, const MacroElement<Mesh> &macro, double minV, double maxV, Rn &u_mean) {
    typedef typename Mesh::Element Element;
    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename Mesh::Rd Rd;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(5));
    const GFESpace<Mesh> &Wh(*uh.Vh);

    if (u_mean.size() != Wh.get_nb_element()) {
        u_mean.resize(Wh.get_nb_element());
    }
    u_mean = 0.;

    // compute  on each element
    for (int k = Wh.first_element(); k < Wh.last_element(); ++k) {

        const FElement &FK(Wh[k]);
        const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k, 0));
        double u_bar_K = 0.;
        double areaCut = 0.;
        for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
            double meas = cutK.measure(it);
            areaCut += meas;
            // get the mean value on K
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                typename QF::QuadraturePoint ip(qf[ipq]); // integration point
                Rd mip      = cutK.mapToPhysicalElement(it, ip);
                double Cint = meas * ip.getWeight();
                u_bar_K += Cint * uh.eval(k, mip);
            }
        }
        u_bar_K   = u_bar_K / areaCut;
        u_mean[k] = u_bar_K;
    }

    // compute mean value of the macro elements
    for (auto &p : macro.macro_element) {
        const MElement &MK(p.second);
        int idx_root      = MK.idx_root_element;
        double area_total = MK.area_total_;
        double areaCut    = 0.;
        double u_bar_K    = 0.;
        for (auto k : MK.idx_element) {
            const FElement &FK(Wh[k]);
            const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k, 0));
            for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
                double meas = cutK.measure(it);
                areaCut += meas;
                for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                    typename QF::QuadraturePoint ip(qf[ipq]); // integration point
                    Rd mip      = cutK.mapToPhysicalElement(it, ip);
                    double Cint = meas * ip.getWeight();
                    u_bar_K += Cint * uh.eval(k, mip);
                }
            }
        }
        assert(fabs(area_total - areaCut) < 1e-14);
        u_bar_K = u_bar_K / areaCut;
        for (auto k : MK.idx_element) {
            u_mean[k] = u_bar_K;
        }
    }
    for (int i = 0; i < u_mean.size(); ++i) {
        double u_bar_K = u_mean(i);
        if (u_bar_K < minV - globalVariable::Epsilon || maxV + globalVariable::Epsilon < u_bar_K) {
            std::cout << " Mean value does not satisfies maximum principle " << std::endl;
            std::cout << "[m, M] = [ " << minV << "\t" << u_bar_K << "\t" << maxV << " ]" << std::endl;
        }
    }
}

template <typename Mesh> std::vector<double> computeTroubleCellIndicator(const FunFEM<Mesh> &uh) {

    using Element = typename Mesh::Element;
    using CutMesh = ActiveMesh<Mesh>;

    const auto &qf(*QF_Simplex<typename Mesh::Rd>(5));
    const GFESpace<Mesh> &Vh(*uh.Vh);
    const CutMesh &Th(Vh.get_mesh());
    int nt = Vh.get_nb_element();
    std::vector<double> indicator(nt);
    std::vector<double> indicator_loc(nt);

    for (int k = Vh.first_element(); k < Vh.last_element(); k += Vh.next_element()) {
        const auto cutK = Th.get_cut_part(k, 0);
        const auto FK(Vh[k]);
        double s     = 0.;
        double maxPj = 0.;

        // compute mean of p_j j=0,1,2,3 on K_0 and mean of p_j j=1,2,3 on K_j
        for (int ifac = 0; ifac < Mesh::Element::nea; ++ifac) {

            int ifacn = ifac;
            int kn    = Vh.getNeighborElement(k, ifacn);
            if (kn == -1)
                continue;
            auto FKj(Vh[kn]);

            double Pj      = 0.;
            double sj      = 0.;
            double meas_K0 = 0.;
            // compute mean of p_j j=0,1,2,3 on K_0
            for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
                double meas_cut = cutK.measure(it);
                meas_K0 += meas_cut;
                for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                    auto ip(qf[ipq]);
                    auto mip  = cutK.mapToPhysicalElement(it, ip);
                    double p0 = meas_cut * ip.getWeight() * uh.eval(k, mip, 0, op_id);
                    double pj = meas_cut * ip.getWeight() * uh.eval(kn, mip, 0, op_id);
                    sj += p0 - pj;
                    Pj += p0;
                }
            }
            Pj /= meas_K0;
            sj /= meas_K0;
            s += std::fabs(sj);
            maxPj            = std::max(maxPj, std::fabs(Pj));
            Pj               = 0.;
            double meas_Kj   = 0.;
            // compute mean of p_j j=1,2,3 on K_j
            const auto cutKj = Th.get_cut_part(kn, 0);
            for (auto it = cutKj.element_begin(); it != cutKj.element_end(); ++it) {
                double meas_cut = cutKj.measure(it);
                meas_Kj += meas_cut;
                for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                    auto ip(qf[ipq]);
                    auto mip = cutKj.mapToPhysicalElement(it, ip);
                    Pj += meas_cut * ip.getWeight() * uh.eval(kn, mip, 0, op_id);
                }
            }
            Pj /= meas_Kj;
            maxPj = std::max(maxPj, std::fabs(Pj));
        }
        // std::cout << s << "\t" << maxPj << "\t" << s / (maxPj +
        // globalVariable::Epsilon) << std::endl; getchar();
        indicator_loc[k] = s / (maxPj + globalVariable::Epsilon);
    }
#ifdef USE_MPI
    MPIcf::AllReduce(indicator_loc.data(), indicator.data(), nt, MPI_MAX);
    return indicator;
#else
    return indicator_loc;
#endif
}

template <typename Mesh>
std::tuple<double, double> minAndMaxAverageNeighbor(const GFElement<Mesh> &FK, const FunFEM<Mesh> &uh) {
    using QF      = typename GFElement<Mesh>::QF;
    using CutMesh = ActiveMesh<Mesh>;

    const QF &qf(*QF_Simplex<typename GFElement<Mesh>::RdHat>(2));
    const CutMesh &Th(FK.Vh.get_mesh());

    double m = 1e300, M = -1e300;
    for (int ifac = 0; ifac < Mesh::Element::nea; ++ifac) { // loop over the edges /faces
        int ifacn = ifac;
        int kn    = FK.Vh.getNeighborElement(FK.index(), ifacn);
        if (kn == -1)
            continue; // border edge
        const GFElement<Mesh> &FKn(FK.Vh[kn]);
        const auto cutKn = Th.get_cut_part(kn, 0);

        double val     = 0.;
        double meas_Kn = 0.;
        for (auto it = cutKn.element_begin(); it != cutKn.element_end(); ++it) {
            double meas_cut = cutKn.measure(it);
            meas_Kn += meas_cut;
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                typename QF::QuadraturePoint ip(qf[ipq]); // integration point
                auto mip = cutKn.mapToPhysicalElement(it, ip);
                val += meas_cut * ip.getWeight() * uh.eval(kn, mip, 0, op_id);
            }
        }
        val /= meas_Kn;
        m = std::min(m, val);
        M = std::max(M, val);
    }
    return {m, M};
}

template <typename Mesh> double averageOnKj(const GFElement<Mesh> &FK, const FunFEM<Mesh> &uh) {
    using QF      = typename GFElement<Mesh>::QF;
    using CutMesh = ActiveMesh<Mesh>;

    const QF &qf(*QF_Simplex<typename GFElement<Mesh>::RdHat>(2));
    const CutMesh &Th(FK.Vh.get_mesh());

    int k           = FK.index();
    const auto cutK = Th.get_cut_part(k, 0);

    double Uj     = 0.;
    double meas_K = 0.;
    for (auto it = cutK.element_begin(); it != cutK.element_end(); ++it) {
        double meas_cut = cutK.measure(it);
        meas_K += meas_cut;
        for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
            typename QF::QuadraturePoint ip(qf[ipq]); // integration point
            auto mip = cutK.mapToPhysicalElement(it, ip);
            Uj += meas_cut * ip.getWeight() * uh.eval(FK.index(), mip, 0, op_id);
        }
    }
    return Uj / meas_K;
}

template <typename Mesh>
double computeAlpha(const GFElement<Mesh> &FK, const FunFEM<Mesh> &uh, double Uj, double m, double M) {
    using QFB     = typename GFElement<Mesh>::QFB;
    using CutMesh = ActiveMesh<Mesh>;

    const QFB &qfb(*QF_Simplex<typename GFElement<Mesh>::RdHatBord>(4));
    const CutMesh &Th(FK.Vh.get_mesh());

    int k           = FK.index();
    const auto cutK = Th.get_cut_part(k, 0);

    double min_alpha_i = 1e300;

    for (int ifac = 0; ifac < Mesh::Element::nea; ++ifac) {
        int ifacn = ifac;
        int kn    = FK.Vh.getNeighborElement(k, ifacn);
        if (kn == -1)
            continue;
        double alpha_i = 1.;

        if (Th.isCutFace(k, ifac, 0)) {
            typename Mesh::Element::Face face;
            const Cut_Part<typename Mesh::Element::Face> cutFace(Th.get_cut_face(face, k, ifac, 0));
            for (auto it = cutFace.element_begin(); it != cutFace.element_end(); ++it) {
                for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
                    typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point

                    const R2 mip = cutFace.mapToPhysicalElement(it, (R1)ip);
                    // compute U_j(x_i^*)
                    double ux    = uh.eval(FK.index(), mip, 0, op_id);
                    if (ux - M > 0) {
                        alpha_i = (M - Uj) / (ux - Uj);
                    } else if (ux - m < 0) {
                        alpha_i = (m - Uj) / (ux - Uj);
                    } else {
                        alpha_i = 1.;
                    }
                    alpha_i     = std::max(alpha_i, 0.);
                    min_alpha_i = std::min(min_alpha_i, alpha_i);
                    assert(min_alpha_i < 1e300);
                }
            }
        } else {
            for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
                typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
                const R2 mip = FK.T(FK.T.toKref((R1)ip, ifac));
                // compute U_j(x_i^*)
                double ux    = uh.eval(FK.index(), mip, 0, op_id);
                if (ux - M > 0) {
                    alpha_i = (M - Uj) / (ux - Uj);
                } else if (ux - m < 0) {
                    alpha_i = (m - Uj) / (ux - Uj);
                } else {
                    alpha_i = 1.;
                }
                alpha_i     = std::max(alpha_i, 0.);
                min_alpha_i = std::min(min_alpha_i, alpha_i);
                assert(min_alpha_i < 1e300);
            }
        }
    }

    return min_alpha_i;
}

template <typename Mesh>
std::vector<double> applySlopeLimiter(const FunFEM<Mesh> &uh, const std::vector<double> &cell_indicator,
                                      const double tol) {

    std::vector<double> u_new(uh.array().begin(), uh.array().end());

    KNM<double> B(3, 3), invB(3, 3);
    KN<double> U_K(3);
    KN<double> dof_Taylor(3);
    const auto &bf_Taylor(DataFE<Mesh2>::P1dcTaylor);
    std::vector<R2> xref_i = {R2(0., 0.), R2(1., 0.), R2(0., 1.)};
    KNMK<double> phi(3, 1, 1);

    const auto &Vh(*uh.Vh);
    if (cell_indicator.size() != Vh.get_nb_element()) {
        std::cout << " Cell indicator array not of the required size " << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int k = Vh.first_element(); k < Vh.last_element(); k += Vh.next_element()) {

        if (cell_indicator[k] < tol)
            continue;

        const auto FK(Vh[k]);

        // compute m and M , min max average of neighbor
        auto [m, M]  = minAndMaxAverageNeighbor(FK, uh);
        // compute solution average on K_j
        double Uj    = averageOnKj(FK, uh);
        // std::cout << "min max \t" << m << "\t" << M << std::endl;
        // std::cout << " average Uj \t" << Uj << std::endl;
        double alpha = computeAlpha(FK, uh, Uj, m, M);

        const auto &T(FK.T);
        // 1) Compute B and invB and U_K
        for (int i = 0; i < T.nv; ++i) {
            bf_Taylor.FB(Fop_D0, T, xref_i[i], phi);
            B(i, '.') = phi('.', 0, op_id);
            U_K(i)    = uh.eval(k, (R2)T[i], 0, op_id);
        }
        invB = inv(B);

        // 2) Get the dof of the Taylor FE
        dof_Taylor = invB * U_K;

        // 3) Modify U_K
        for (int i = 0; i < T.nv; ++i) {
            u_new[FK.loc2glb(i)] =
                dof_Taylor(0) * B(i, 0) + alpha * (dof_Taylor(1) * B(i, 1) + dof_Taylor(2) * B(i, 2));
        }
    }
    return u_new;
}

} // namespace CutFEM

namespace FEM {

// -------------------------------------------------------------------------
// COMPUTE THE MIN AND THE MAX VALUE OF A FE FUNCTION ON A  MESH
// DEPENDING OF THE POLYNOMIAL ORDER.
// BY DEFAULT IT USES LINEAR ELEMENTS

template <typename Mesh> std::tuple<double, double> findMinAndMaxValue_P0(const FunFEM<Mesh> &uh) {
    typedef typename Mesh::Rd Rd;
    typedef typename Mesh::Element Element;
    typedef typename GFESpace<Mesh>::FElement FElement;
    const GFESpace<Mesh> &Wh(*uh.Vh);
    const Mesh &Kh(Wh.Th);
    double min_val, max_val;
    double min_val_loc = 1e300;
    double max_val_loc = -1e300;
    int k_min, k_max;

    for (int k = Wh.first_element(); k < Wh.last_element(); k += Wh.next_element()) {
        const FElement &FK(Wh[k]);

        Rd mip     = FK.T.barycenter();
        double val = uh.eval(k, mip);
        if (val > max_val_loc) {
            max_val_loc = val;
            k_max       = k;
        }
        if (val < min_val_loc) {
            min_val_loc = val;
            k_min       = k;
        }
    }
#ifdef USE_MPI
    MPIcf::AllReduce(min_val_loc, min_val, MPI_MIN);
    MPIcf::AllReduce(max_val_loc, max_val, MPI_MAX);
#else
    min_val = min_val_loc;
    max_val = max_val_loc;
#endif
    return {min_val, max_val};
}

template <typename Mesh> std::tuple<double, double> findMinAndMaxValue_P1(const FunFEM<Mesh> &uh) {
    typedef typename Mesh::Rd Rd;
    typedef typename Mesh::Element Element;
    typedef typename GFESpace<Mesh>::FElement FElement;
    int nv_loc = Mesh::Rd::d + 1;
    const GFESpace<Mesh> &Wh(*uh.Vh);
    const Mesh &Kh(Wh.Th);

    double min_val, max_val;
    double min_val_loc = 1e300;
    double max_val_loc = -1e300;
    // Apply limiter in all element
    int k_min, k_max;
    for (int k = Wh.first_element(); k < Wh.last_element(); k += Wh.next_element()) {
        const FElement &FK(Wh[k]);

        // get the local min and std::max on K
        for (int ipq = 0; ipq < nv_loc; ++ipq) {
            Rd P       = Rd::KHat[ipq];
            Rd mip     = FK.mapToPhysicalElement(P);
            double val = uh.eval(k, mip);
            if (val > max_val_loc) {
                max_val_loc = val;
                k_max       = k;
            }
            if (val < min_val_loc) {
                min_val_loc = val;
                k_min       = k;
            }
        }
    }
#ifdef USE_MPI
    MPIcf::AllReduce(min_val_loc, min_val, MPI_MIN);
    MPIcf::AllReduce(max_val_loc, max_val, MPI_MAX);
#else
    min_val = min_val_loc;
    max_val = max_val_loc;
#endif
    return {min_val, max_val};
}

template <typename Mesh> std::tuple<double, double> findMinAndMaxValue(const FunFEM<Mesh> &uh) {
    const int fctOrder = uh.getPolynomialOrder();
    if (fctOrder == 0) {
        return findMinAndMaxValue_P0(uh);
    } else if (fctOrder == 1) {
        return findMinAndMaxValue_P1(uh);
    } else {
        std::cout << " Min and std::max compute using linear elements. Not exact for "
                     "this polynomial order"
                  << std::endl;
        return findMinAndMaxValue_P1(uh);
    }
}

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

template <typename Mesh> std::vector<double> computeMeanValue(const FunFEM<Mesh> &uh, double minV, double maxV) {
    typedef typename Mesh::Element Element;
    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename Mesh::Rd Rd;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(5));
    const GFESpace<Mesh> &Wh(*uh.Vh);
    const Mesh &Kh(Wh.Th);

    std::vector<double> u_mean(Wh.get_nb_element(), 0.);

    for (int k = Wh.first_element(); k < Wh.last_element(); ++k) {

        const FElement &FK(Wh[k]);
        const Element &K(FK.T);

        double u_bar_K = 0.;
        double meas    = K.measure();
        // get the mean value on K
        for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
            typename QF::QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip      = K.mapToPhysicalElement(ip);
            double Cint = ip.getWeight();
            u_bar_K += Cint * uh.eval(k, mip);
        }
        u_mean[k] = u_bar_K;
        if (u_bar_K < minV || maxV < u_bar_K) {
            std::cout << " Mean value does not satisfies maximum principle " << std::endl;
            std::cout << "[m, M] = [ " << minV << "\t" << u_bar_K << "\t" << maxV << " ]" << std::endl;
        }
    }
    return u_mean;
}

template <typename Mesh> void checkMeanValue(const FunFEM<Mesh> &uh, double minV, double maxV) {
    typedef typename Mesh::Element Element;
    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename FElement::QF QF;
    typedef typename Mesh::Rd Rd;

    const QF &qf(*QF_Simplex<typename FElement::RdHat>(5));
    const GFESpace<Mesh> &Wh(*uh.Vh);
    const Mesh &Kh(Wh.Th);

    for (int k = Wh.first_element(); k < Wh.last_element(); ++k) {

        const FElement &FK(Wh[k]);
        const Element &K(FK.T);

        double u_bar_K = 0.;
        double meas    = K.measure();
        // get the mean value on K
        for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
            typename QF::QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip      = K.mapToPhysicalElement(ip);
            double Cint = ip.getWeight();
            u_bar_K += Cint * uh.eval(k, mip);
        }
        if (u_bar_K < minV || maxV < u_bar_K) {
            std::cout << " Mean value does not satisfies maximum principle " << std::endl;
            std::cout << "[m, M] = [ " << minV << "\t" << u_bar_K << "\t" << maxV << " ]" << std::endl;
            getchar();
        }
    }
}

template <typename Mesh>
std::vector<double> boundPreservingLimiter_P1(const FunFEM<Mesh> &uh, double min_val, double max_val) {
    typedef typename Mesh::Element Element;
    typedef typename Mesh::Rd Rd;
    typedef typename GFESpace<Mesh>::FElement FElement;
    typedef typename FElement::QF QF;
    const QF &qf(*QF_Simplex<typename FElement::RdHat>(2));
    std::vector<double> u_new(uh.array().begin(), uh.array().end());

    int nv_loc = Mesh::Rd::d + 1;
    const GFESpace<Mesh> &Wh(*uh.Vh);
    const Mesh &Kh(Wh.Th);

    // Apply limiter in all element
    for (int k = Wh.first_element(); k < Wh.last_element(); ++k) {

        const FElement &FK(Wh[k]);
        const Element &K(FK.T);

        double u_bar_K = 0.;
        double max_K = -1e300, min_K = 1e300;
        double meas = K.measure();
        // get the mean value on K
        for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
            typename QF::QuadraturePoint ip(qf[ipq]); // integration point
            Rd mip      = K.mapToPhysicalElement(ip);
            double Cint = ip.getWeight();
            u_bar_K += Cint * uh.eval(k, mip);
        }
        // if(u_bar_K < min_val ) u_bar_K = min_val;

        // get the local min and std::max on K
        for (int ipq = 0; ipq < nv_loc; ++ipq) {
            Rd mip     = K[ipq];
            double val = uh.eval(k, mip);

            max_K = std::max(max_K, val);
            min_K = std::min(min_K, val);
        }

        // 3) compute theta
        double v1 = fabs((min_val - u_bar_K) / (min_K - u_bar_K));
        if ((min_K - u_bar_K) == 0.)
            v1 = 1.; //
        double v2 = fabs((max_val - u_bar_K) / (max_K - u_bar_K));
        if ((max_K - u_bar_K) == 0.)
            v2 = 1.;
        double theta = std::min(std::min(v1, v2), 1.);

        // 4) replace the dof
        for (int df = 0; df < FK.NbDoF(); ++df) {
            int df_glb    = FK.loc2glb(df);
            u_new[df_glb] = theta * (uh(df_glb) - u_bar_K) + u_bar_K;

            // if(u_new(df_glb) - min_val < 0 && u_new(df_glb) - min_val >
            // -Epsilon ) u_new(df_glb) = min_val; if(u_new(df_glb) - max_val >
            // 0 && u_new(df_glb) - max_val < Epsilon ) u_new(df_glb) = max_val;

            // if(u_new(df_glb) < min_val) {
            //   std::cout << " ----------------------- " << std::endl;
            //   std::cout << " Element " << k << "\n"
            //             << " theta   " << theta << "\n"
            //             << " std::min K   " << min_K << "\n"
            //             << " minK-u  " << min_K - u_bar_K << "\n"
            //             << " minv-u  " << min_val - u_bar_K << "\n"
            //             << " u_bar_K " << u_bar_K << "\n"
            //             << " uh      " << uh(df_glb) << "\n"
            //             << " u_new   " << u_new(df_glb) << std::endl;
            //             getchar();
            //
            // }
        }
    }
    return u_new;
}

template <typename Mesh> std::tuple<double, double> minmaxP1(const FunFEM<Mesh> &uh) {
    typedef typename Mesh::Rd Rd;
    typedef typename Mesh::Element Element;
    typedef typename GFESpace<Mesh>::FElement FElement;
    int nv_loc = Mesh::Rd::d + 1;
    const GFESpace<Mesh> &Wh(*uh.Vh);
    const Mesh &Kh(Wh.Th);

    double min_val = 1e300;
    double max_val = -1e300;

    // Apply limiter in all element
    for (int k = Wh.first_element(); k < Wh.last_element(); ++k) {
        const FElement &FK(Wh[k]);

        // get the local std::min and std::max on K
        for (int ipq = 0; ipq < nv_loc; ++ipq) {
            Rd mip     = FK.T[ipq];
            // double val = uh.eval(k, mip);
            double val = uh.eval(k, mip);
            max_val    = std::max(max_val, val);
            min_val    = std::min(min_val, val);
        }
    }
    return {min_val, max_val};
}

template <typename Mesh>
std::vector<double> applyBoundPreservingLimiter(const FunFEM<Mesh> &uh, double min_val, double max_val) {
    if (uh.getBasisFctType() == BasisFctType::P1dc) {
        return boundPreservingLimiter_P1(uh, min_val, max_val);
    } else {
        std::cout << " No boundary preserving limiter implemented for this elements" << std::endl;
        exit(EXIT_FAILURE);
        return std::vector<double>();
    }
}

template <typename Mesh> std::vector<double> computeTroubleCellIndicator(const FunFEM<Mesh> &uh) {

    const auto &qf(*QF_Simplex<typename Mesh::Rd>(5));
    const auto &Vh(*uh.Vh);
    int nt = Vh.get_nb_element();
    std::vector<double> indicator(nt);
    std::vector<double> indicator_loc(nt);

    for (int k = Vh.first_element(); k < Vh.last_element(); k += Vh.next_element()) {

        auto FK(Vh[k]);
        // double meas_K0 = FK.get_measure();
        double s     = 0.;
        double maxPj = 0.;

        // compute mean of p_j j=0,1,2,3 on K_0 and mean of p_j j=1,2,3 on K_j
        for (int ifac = 0; ifac < Mesh::Element::nea; ++ifac) {

            int ifacn = ifac;
            int kn    = Vh.getNeighborElement(k, ifacn);
            if (kn == -1)
                continue;
            auto FKj(Vh[kn]);

            // double meas_Kj = FKj.get_measure();
            double Pj = 0.;
            double sj = 0.;
            // compute mean of p_j j=0,1,2,3 on K_0
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                auto ip(qf[ipq]);
                auto mip = FK.T((typename Mesh::Rd)ip);

                double p0 = ip.getWeight() * uh.eval(k, mip, 0, op_id);
                double pj = ip.getWeight() * uh.eval(kn, mip, 0, op_id);
                sj += p0 - pj;
                Pj += p0;
            }
            s += std::fabs(sj);
            maxPj = std::max(maxPj, std::fabs(Pj));
            Pj    = 0.;

            // compute mean of p_j j=1,2,3 on K_j
            for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
                auto ip(qf[ipq]);
                auto mip = FKj.T((typename Mesh::Rd)ip);
                Pj += ip.getWeight() * uh.eval(kn, mip, 0, op_id);
            }
            maxPj = std::max(maxPj, std::fabs(Pj));
        }
        indicator_loc[k] = s / (maxPj + globalVariable::Epsilon);
    }
#ifdef USE_MPI
    MPIcf::AllReduce(indicator_loc.data(), indicator.data(), nt, MPI_MAX);
    return indicator;
#else
    return indicator_loc;
#endif
}

template <typename Mesh>
std::tuple<double, double> minAndMaxAverageNeighbor(const GFElement<Mesh> &FK, const FunFEM<Mesh> &uh) {
    typedef typename GFElement<Mesh>::QF QF;
    const QF &qf(*QF_Simplex<typename GFElement<Mesh>::RdHat>(2));

    double m = 1e300, M = -1e300;
    for (int ifac = 0; ifac < Mesh::Element::nea; ++ifac) { // loop over the edges /faces
        int ifacn = ifac;
        int kn    = FK.Vh.getNeighborElement(FK.index(), ifacn);
        if (kn == -1)
            continue; // border edge
        const GFElement<Mesh> &FKn(FK.Vh[kn]);
        double val = 0.;
        for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
            typename QF::QuadraturePoint ip(qf[ipq]); // integration point
            const R2 mip = FKn.T((R2)ip);
            val += ip.getWeight() * uh.eval(kn, mip, 0, op_id);
        }
        m = std::min(m, val);
        M = std::max(M, val);
    }
    return {m, M};
}

template <typename Mesh> double averageOnKj(const GFElement<Mesh> &FK, const FunFEM<Mesh> &uh) {
    typedef typename GFElement<Mesh>::QF QF;
    const QF &qf(*QF_Simplex<typename GFElement<Mesh>::RdHat>(2));

    double Uj = 0.;
    for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
        typename QF::QuadraturePoint ip(qf[ipq]); // integration point
        const R2 mip = FK.T((R2)ip);
        Uj += ip.getWeight() * uh.eval(FK.index(), mip, 0, op_id);
    }
    return Uj;
}

template <typename Mesh>
double computeAlpha(const GFElement<Mesh> &FK, const FunFEM<Mesh> &uh, double Uj, double m, double M) {
    typedef typename GFElement<Mesh>::QFB QFB;
    const QFB &qfb(*QF_Simplex<typename GFElement<Mesh>::RdHatBord>(4));

    double min_alpha_i = 1e300;

    for (int ifac = 0; ifac < Mesh::Element::nea; ++ifac) {
        int ifacn = ifac;
        int kn    = FK.Vh.getNeighborElement(FK.index(), ifacn);
        if (kn == -1)
            continue;
        double alpha_i = 1.;
        for (int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq) {
            typename QFB::QuadraturePoint ip(qfb[ipq]); // integration point
            const R2 mip = FK.T(FK.T.toKref((R1)ip, ifac));
            // compute U_j(x_i^*)
            double ux    = uh.eval(FK.index(), mip, 0, op_id);
            if (ux - M > 0) {
                alpha_i = (M - Uj) / (ux - Uj);
            } else if (ux - m < 0) {
                alpha_i = (m - Uj) / (ux - Uj);
            } else {
                alpha_i = 1.;
            }
            alpha_i     = std::max(alpha_i, 0.);
            min_alpha_i = std::min(min_alpha_i, alpha_i);
            assert(min_alpha_i < 1e300);
        }
    }

    return min_alpha_i;
}

template <typename Mesh>
std::vector<double> applySlopeLimiter(const FunFEM<Mesh> &uh, const std::vector<double> &cell_indicator,
                                      const double tol) {

    std::vector<double> u_new(uh.array().begin(), uh.array().end());
    // std::vector<double> u_new(uh.size());

    KNM<double> B(3, 3), invB(3, 3);
    KN<double> U_K(3);
    KN<double> dof_Taylor(3);
    const auto &bf_Taylor(DataFE<Mesh2>::P1dcTaylor);
    std::vector<R2> xref_i = {R2(0., 0.), R2(1., 0.), R2(0., 1.)};
    KNMK<double> phi(3, 1, 1);

    const auto &Vh(*uh.Vh);
    if (cell_indicator.size() != Vh.get_nb_element()) {
        std::cout << " Cell indicator array not of the required size " << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int k = Vh.first_element(); k < Vh.last_element(); k += Vh.next_element()) {

        if (cell_indicator[k] < tol)
            continue;

        const auto FK(Vh[k]);

        // compute m and M , min max average of neighbor
        auto [m, M] = minAndMaxAverageNeighbor(FK, uh);
        // compute solution average on K_j
        double Uj   = averageOnKj(FK, uh);
        std::cout << "min max \t" << m << "\t" << M << std::endl;
        std::cout << " average Uj \t" << Uj << std::endl;
        double alpha = computeAlpha(FK, uh, Uj, m, M);

        const auto &T(FK.T);
        // 1) Compute B and invB and U_K
        for (int i = 0; i < T.nv; ++i) {
            bf_Taylor.FB(Fop_D0, T, xref_i[i], phi);
            B(i, '.') = phi('.', 0, op_id);
            U_K(i)    = uh.eval(k, (R2)T[i], 0, op_id);
        }
        invB = inv(B);

        // 2) Get the dof of the Taylor FE
        dof_Taylor = invB * U_K;

        // 3) Modify U_K
        for (int i = 0; i < T.nv; ++i) {
            u_new[FK.loc2glb(i)] =
                dof_Taylor(0) * B(i, 0) + alpha * (dof_Taylor(1) * B(i, 1) + dof_Taylor(2) * B(i, 2));
        }
    }
    return u_new;
}

}; // namespace FEM
}; // namespace limiter

// class Limiter {
//
//   typedef typename Mesh2::Element Element;
//   typedef typename FESpace2::FElement FElement;
//   typedef typename FElement::QF QF;
//   typedef typename FElement::QFB QFB;
//
//   const QF&  qf ;
//   const QFB& qfb;
//
// public:
//   Rn indicator;
//
//   Limiter() : qf(*QF_Simplex<R2>(3)), qfb(*QF_Simplex<R1>(3)){}
//
//   void KXRCF_indicator(const Fun2_h& uh, const Fun2_h& flux);
//
//   void limiting(const Fun2_h& uh, Rn& lh);
//   void min_and_max_average_neighbor(const FElement& FK, const Fun2_h&uh,
//   double&m, double&M); double average_on_Kj(const FElement& FK, const
//   Fun2_h&uh); double compute_alpha(const FElement& FK, const Fun2_h&uh,
//   double Uj, double m, double M);
// };
//
//
// class Box_Side {
// public:
//   // a side is the part that satisfies ax+by+c < 0
//   double a, b, c;
//   Box_Side(double aa, double bb, double cc) : a(aa), b(bb), c(cc) {}
//
//   double f(const R2 A) const {
//     return (a*A.x+b*A.y+c);
//   }
//   bool is_negative(const R2 P) const {
//     return (f(P) <=0);
//   }
//   bool is_positive(const R2 P) const {
//     return (f(P) > 0);
//   }
//   bool cut_edge(const R2 A, const R2 B) const {
//     return (f(A)*f(B)< 0);
//   }
//   R2 get_cut_node(const R2 A, const R2 B) const {
//     double t = -f(A)/(f(B)-f(A));
//     return (1-t) * A + t * B;
//   }
//
// };
//
// class Filter_Box {
//   typedef typename Mesh2::Element Element;
// public:
//   static const int nv = 3;                         // triangle
//   typedef SortArray<Ubyte, nv> TriaIdx;           // the vertices of a
//   triangle of the cut:
//
//   std::vector<Box_Side> sides_;
//   std::vector<R2>      vertices_;
//   std::vector<TriaIdx> triangles_;               // idx of the vertices
//
//   R2 minXY, maxXY;
//
//   Filter_Box(const R2 mminXY, const R2 mmaxXY) : minXY(mminXY),
//   maxXY(mmaxXY){
//     sides_.push_back(Box_Side( 0., -1.,  minXY.y));
//     sides_.push_back(Box_Side( 1.,  0., -maxXY.x));
//     sides_.push_back(Box_Side( 0.,  1., -maxXY.y));
//     sides_.push_back(Box_Side(-1.,  0.,  minXY.x));
//   }
//
//   bool contain(const R2 P) const{
//     for(auto it=sides_.begin(); it != sides_.end();++it) {
//       if(it->is_positive(P)) return false;
//     }
//     return true;
//   }
//
//   bool intersect_with(const typename Mesh2::Element& K) const{
//     // if we find one node inside the box it is enough to know
//     for(int i=0; i<3;++i) {
//       if(this->contain((R2) K[i])) return true;
//     }
//     return false;
//   }
//
//   void cut_triangle(const typename Mesh2::Element& K) ;
//
//   void print(std::string filename = "newTriangle.dat") const ;
//
// };
//
// class Filter {
// public:
//
//
// };
// template<typename Mesh>
// void check_mean_value(const FunFEM<Mesh>& uh, const MacroElement<Mesh>&
// macro, double minV, double maxV, Rn& u_mean) {
//   typedef typename Mesh::Element Element;
//   typedef typename GFESpace<Mesh>::FElement FElement;
//   typedef typename FElement::QF QF;
//   typedef typename Mesh::Rd Rd;
//
//   const QF& qf(*QF_Simplex<typename FElement::RdHat>(5));
//   const GFESpace<Mesh>& Wh(*uh.Vh);
//
//   if(u_mean.size() != Wh.get_nb_element()){
//     u_mean.resize(Wh.get_nb_element());
//   }
//   u_mean = 0.;
//
//   // compute mean value on each element
//   for(int k=Wh.first_element(); k<Wh.last_element();++k) {
//
//     const FElement& FK(Wh[k]);
//     const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k,0));
//     double u_bar_K = 0.;
//     double areaCut = 0.;
//     for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
//       double meas = cutK.measure(it);
//       areaCut += meas;
//       // get the mean value on K
//       for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//         typename QF::QuadraturePoint ip(qf[ipq]); // integration point
//         Rd mip = cutK.mapToPhysicalElement(it, ip);
//         double Cint = meas * ip.getWeight();
//         u_bar_K += Cint*uh.eval(k, mip);
//       }
//     }
//     u_bar_K = u_bar_K / areaCut;
//     u_mean[k] = u_bar_K;
//   }
//
//   // compute mean value of the macro elements
//   for (auto& p : macro.macro_element) {
//     const MElement& MK(p.second);
//     int idx_root = MK.idx_root_element;
//     double area_total = MK.area_total_;
//     double areaCut = 0.;
//     double u_bar_K = 0.;
//     for(auto k : MK.idx_element ) {
//       const FElement& FK(Wh[k]);
//       const Cut_Part<Element> cutK(macro.Th_.get_cut_part(k,0));
//       for(auto it = cutK.element_begin();it != cutK.element_end(); ++it){
//         double meas = cutK.measure(it);
//         areaCut += meas;
//         for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {
//           typename QF::QuadraturePoint ip(qf[ipq]); // integration point
//           Rd mip = cutK.mapToPhysicalElement(it, ip);
//           double Cint = meas * ip.getWeight();
//           u_bar_K += Cint*uh.eval(k, mip);
//         }
//       }
//     }
//     assert(fabs(area_total - areaCut) < 1e-14);
//     u_bar_K = u_bar_K / areaCut;
//   }
// }

#endif
