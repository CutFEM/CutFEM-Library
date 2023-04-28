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
/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef _GTYPE_OF_FE_HPP
#define _GTYPE_OF_FE_HPP

#include <iostream>
#include <cassert>

#include "../num/util.hpp"
#include "../common/RNM.hpp"
#include "../common/Mesh1dn.hpp"
#include "operationTypeFE.hpp"
// #include "transformation.hpp"

enum class BasisFctType {
    P0,
    P1,
    P2,
    P3,
    P4,
    P1dc,
    P2dc,
    P3dc,
    P0poly,
    P1poly,
    P2poly,
    P2BR,
    RT0,
    RT1,
    RT2,
    BDM1,
    BDM2,
    UNDEFINED
};

/*
 *   Basis class for the FE
 *   Contains divers information
 */
class dataTypeOfFE {

  private:
    const int *data; // array to save all the data needed for the FE
    const int *dataalloc;

  public:
    const int ndfonVertex; // number of dof on points
    const int ndfonEdge;   // number of dof on egdes
    const int ndfonFace;   // number of dof on faces
    const int ndfonVolume; // number of dof in volume
    const int nbDoF;       // dof per element
    const int nbNode;      // node per element = dof so far
    int N;                 // dim espace arrivée N = 1
    int nbOfFE;
    bool constDfPerNode = true;

    int const *const nbNodeOnWhat;
    int const *const DFOnWhat;
    int const *const DFOfNode;
    int const *const NodeOfDF;

    const int *ndfOn() const { return &ndfonVertex; }

    dataTypeOfFE(const std::vector<int> &nitemdim, const KN<dataTypeOfFE const *> &);
    dataTypeOfFE(const std::vector<int> &nitemdim, std::vector<dataTypeOfFE const *> &&);

    dataTypeOfFE(const std::vector<int> &nitem, const int *Data, int nbdf, int NN);

    virtual ~dataTypeOfFE() {
        if (dataalloc)
            delete[] dataalloc;
    }
};

template <class RdHat> class InterpolationMatrix;
template <class Mesh> class GbaseFElement;
template <class Mesh> class GFElement;

struct IPJ {
    int i, p, j; // i is DoF, p is Point, j is componante
    IPJ(int ii = 0, int pp = 0, int jj = 0) : i(ii), p(pp), j(jj) {}
};

/*
 *   Main class for the FE structure
 *
 */
template <class Mesh> class GTypeOfFE : public dataTypeOfFE {
  public:
    typedef typename Mesh::Element Element;
    typedef Element E;
    typedef typename Element::RdHat RdHat;
    typedef typename Element::Rd Rd;
    typedef GFElement<Mesh> FElement;
    typedef double (*bf_type)(double *, int);

    BasisFctType basisFctType = BasisFctType::UNDEFINED;
    int polynomialOrder;

    int NbPtforInterpolation;   // Nb of interpolation points per elemen
    int NbcoefforInterpolation; // Nb of interpolation points per element

    KN<IPJ> ipj_Pi_h;   //?
    KN<RdHat> Pt_Pi_h;  //?
    double *coef_Pi_h;  //?
    // KN<int> begin_coef_Pi_h, end_coef_Pi_h;

    KN<GTypeOfFE<Mesh> *> Sub_ToFE;
    KN<int> begin_dfcomp, end_dfcomp;

    virtual void get_Coef_Pi_h(const GbaseFElement<Mesh> &K, KN_<double> &v) const {
        assert(coef_Pi_h);
        // std::cout << "get coef in GTypeOfFE" << std::endl;
        v = KN_<double>(coef_Pi_h, ipj_Pi_h.N());
    }

    virtual bf_type referenceBasisFunction(int i) const {
        return [](double *x, int c0) { return 0.; };
    };
    /*
     * Simple constructor for simple finite element
     * (not for union of FE)
     */
    GTypeOfFE(const int nbdf, const int NN, const int *data, int kPi, int npPi, double *coef_Pi_h_a = 0)
        : dataTypeOfFE(Element::itemTopology(), data, nbdf, NN),

          NbPtforInterpolation(npPi), NbcoefforInterpolation(kPi), ipj_Pi_h(kPi), Pt_Pi_h(npPi), coef_Pi_h(coef_Pi_h_a),
          // begin_coef_Pi_h, end_coef_Pi_h;
          Sub_ToFE(nbOfFE), begin_dfcomp(N, 0), end_dfcomp(N, nbdf)
    // begin_dfcomp(data+4*nbdf+4+N),
    // end_dfcomp(data+4*nbdf+4+2*N)
    {
        Sub_ToFE = this;
    }

    /*
     *  Constructor for GTypeOfFE_Sum
     *
     */
    GTypeOfFE(const KN<GTypeOfFE<Mesh> const *> &t)
        : dataTypeOfFE(Element::itemTopology(), t), NbPtforInterpolation(this->nbNode),
          NbcoefforInterpolation(this->nbNode), Sub_ToFE(nbOfFE), begin_dfcomp(N, 0), end_dfcomp(N, this->nbDoF)

    {

        Sub_ToFE = this;
    }

    GTypeOfFE(std::vector<GTypeOfFE<Mesh> const *> &t)
        : dataTypeOfFE(Element::itemTopology(), std::vector<dataTypeOfFE const *>(t.begin(), t.end())),
          NbPtforInterpolation(this->nbNode), NbcoefforInterpolation(this->nbNode), Sub_ToFE(nbOfFE),
          begin_dfcomp(N, 0), end_dfcomp(N, this->nbDoF)

    {

        Sub_ToFE = this;
    }

    //   /*
    //  *  Constructor for GTypeOfFE_Time
    //  *
    //  */
    // GTypeOfFE(const KN<GTypeOfFE<Mesh> const *> &t, const GTypeOfFE<Mesh1>*
    // tt)
    //   :
    //   dataTypeOfFE(Element::nitemdim,t, tt),
    //   NbPtforInterpolation(this->nbNode),
    //   NbcoefforInterpolation(this->nbNode),
    //   PtInterpolation(0),
    //   coefInterpolation(0,0),
    //   Sub_ToFE(nbOfFE),
    //   begin_dfcomp(N,0),
    //   end_dfcomp(N,this->nbDoF)

    // {
    //   Sub_ToFE = this;
    // }

    // virtual void init(InterpolationMatrix<RdHat> & M) const;
    virtual void FB(const What_d whatd, const Element &K, const Rd &P, KNMK_<R> &val) const = 0;
    virtual void FB(const What_d whatd, const Element &K, const Rd &P, KNMK_<R> &val, const KNM_<R> &J) const {
        assert(0);
    };

    ~GTypeOfFE() {}

  private:
    GTypeOfFE(const GTypeOfFE &);
    void operator=(const GTypeOfFE &);
};

// template<class Mesh>
// void GTypeOfFE<Mesh>::init(InterpolationMatrix<RdHat> & M) const
//   {
//     assert(M.np==NbPtforInterpolation);
//     assert(M.ncoef==NbcoefforInterpolation);
//
//     M.P=PtInterpolation;
//     M.coef=coefInterpolation;
//   }

/*
 *   Structure that will contain the different FE
 */
template <class mesh> struct DataFE {
    static GTypeOfFE<mesh> &P0;
    static GTypeOfFE<mesh> &P1;
    static GTypeOfFE<mesh> &P2;
    static GTypeOfFE<mesh> &P3;
    static GTypeOfFE<mesh> &P4;

    static GTypeOfFE<mesh> &P0Poly;
    static GTypeOfFE<mesh> &P1Poly;
    static GTypeOfFE<mesh> &P2Poly;
    static GTypeOfFE<mesh> &P3Poly;

    static GTypeOfFE<mesh> &RT0;
    static GTypeOfFE<mesh> &RT0m;
    static GTypeOfFE<mesh> &RT1;
    static GTypeOfFE<mesh> &RT2;
    static GTypeOfFE<mesh> &BDM1;
    static GTypeOfFE<mesh> &BDM2;
    static GTypeOfFE<mesh> &P2BR;

    static GTypeOfFE<mesh> &P1dc;
    static GTypeOfFE<mesh> &P2dc;
    static GTypeOfFE<mesh> &P3dc;

    static GTypeOfFE<mesh> &P0sc;
    static GTypeOfFE<mesh> &P1dcsc;

    static GTypeOfFE<mesh> &P1dcTaylor;
};

template <int D> void jacobianLinearTransformation(KNM<R> &J, const Triangle2 &T);

template <int D> void inverseJacobian(const KNM<R> &J, const double inv_detJ, KNM<R> &invJ);

template <int D> double inverseDeterminant(const KNM<R> &J);

template <int D>
void piolatTransformation(const KN_<R> &Ji, double inv_detJ, const std::vector<std::vector<double>> &bf, RN_ &bfMat);

template <int D> double piolatTransformation(const KN_<R> &Ji, double inv_detJ, std::vector<double> &phi);

template <int D>
void piolatTransformationGradient(int df, const KNM_<R> &J, const KNM_<R> &invJ, double inv_det_J,
                                  const KNM_<R> &Dphi_ref, KNMK_<R> &bfMat);

template <> inline void jacobianLinearTransformation<2>(KNM<R> &J, const Triangle2 &T) {
    R2 xv0  = T.at(0);
    R2 xv1  = T.at(1);
    R2 xv2  = T.at(2);
    J(0, 0) = xv1[0] - xv0[0];
    J(0, 1) = xv2[0] - xv0[0];
    J(1, 0) = xv1[1] - xv0[1];
    J(1, 1) = xv2[1] - xv0[1];
}

template <> inline void inverseJacobian<2>(const KNM<R> &J, const double inv_detJ, KNM<R> &invJ) {
    invJ(0, 0) = inv_detJ * J(1, 1);
    invJ(0, 1) = -inv_detJ * J(0, 1);
    invJ(1, 0) = -inv_detJ * J(1, 0);
    invJ(1, 1) = inv_detJ * J(0, 0);
}

template <> inline double inverseDeterminant<2>(const KNM<R> &J) {
    return 1. / (J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0));
}

// template <>
// inline void
// piolatTransformation<2>(const KN_<R> &Ji, double inv_detJ,
//                         const std::vector<std::vector<double>> &bf_ref,
//                         RN_ &bfMat) {
//     int df = 0;
//     for (const auto &phi : bf_ref) {
//         bfMat[df++] =
//             inv_detJ * (phi[0] * std::abs(Ji[0]) + phi[1] * std::abs(Ji[1]));
//     }
//     bfMat[9]  = inv_detJ * (bf_ref[9][0] * (Ji[0]) + bf_ref[9][1] * (Ji[1]));
//     bfMat[10] = inv_detJ * (bf_ref[10][0] * (Ji[0]) + bf_ref[10][1] *
//     (Ji[1])); bfMat[11] = inv_detJ * (bf_ref[11][0] * (Ji[0]) + bf_ref[11][1]
//     * (Ji[1]));
// }
template <>
inline void piolatTransformation<2>(const KN_<R> &Ji, double inv_detJ, const std::vector<std::vector<double>> &bf_ref,
                                    RN_ &bfMat) {
    int df = 0;
    for (const auto &phi : bf_ref) {
        bfMat[df++] = inv_detJ * (phi[0] * (Ji[0]) + phi[1] * (Ji[1]));
    }
}

template <> inline double piolatTransformation<2>(const KN_<R> &Ji, double inv_detJ, std::vector<double> &phi) {
    return inv_detJ * (phi[0] * Ji[0] + phi[1] * Ji[1]);
}

template <>
inline void piolatTransformationGradient<2>(int df, const KNM_<R> &J, const KNM_<R> &invJ, double inv_det_J,
                                            const KNM_<R> &Dphi_ref, KNMK_<R> &bfMat) {

    std::vector<std::vector<double>> J_Dphi = {
        {J(0, 0) * Dphi_ref(0, 0) + J(0, 1) * Dphi_ref(1, 0), J(0, 0) * Dphi_ref(0, 1) + J(0, 1) * Dphi_ref(1, 1)},
        {J(1, 0) * Dphi_ref(0, 0) + J(1, 1) * Dphi_ref(1, 0), J(1, 0) * Dphi_ref(0, 1) + J(1, 1) * Dphi_ref(1, 1)}};

    bfMat(df, 0, op_dx) = inv_det_J * (J_Dphi[0][0] * invJ(0, 0) + J_Dphi[0][1] * invJ(1, 0));
    bfMat(df, 0, op_dy) = inv_det_J * (J_Dphi[0][0] * invJ(0, 1) + J_Dphi[0][1] * invJ(1, 1));
    bfMat(df, 1, op_dx) = inv_det_J * (J_Dphi[1][0] * invJ(0, 0) + J_Dphi[1][1] * invJ(1, 0));
    bfMat(df, 1, op_dy) = inv_det_J * (J_Dphi[1][0] * invJ(0, 1) + J_Dphi[1][1] * invJ(1, 1));
}

#endif
