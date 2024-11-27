#ifndef __NUMBER_SPACE_HPP__
#define __NUMBER_SPACE_HPP__

#include "GTypeOfFE.hpp"

class TypeOfFE_NumberSpace2d : public GTypeOfFE<Mesh2> {

    typedef Mesh2 Mesh;
    typedef typename Mesh::Element E;
    static const int nbNodeOnItem[4];

  public:
    static const int k   = 0;
    static const int ndf = 1;
    static int Data[];
    static double alpha_Pi_h[];

    TypeOfFE_NumberSpace2d() : GTypeOfFE<Mesh2>(ndf, 1, Data, 1, 1, alpha_Pi_h) {

        GTypeOfFE<Mesh2>::basisFctType   = BasisFctType::P0;
        GTypeOfFE<Mesh>::polynomialOrder = 0;

        static const R2 Pt[1] = {R2(1. / 3, 1. / 3)};

        for (int i = 0; i < ndf; ++i) {
            Pt_Pi_h[i]  = Pt[i];
            ipj_Pi_h[i] = IPJ(i, i, 0);
        }

        // for (int i = 0; i < ndf; ++i) {
        //     Pt_Pi_h[i]  = Pt[i];
        //     ipj_Pi_h[i] = IPJ(i, i, 0);
        // }
    }

    // void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
    //   for (int i = 0; i < 3; ++i) v[i] = 1;
    // }

    void FB(const What_d, const Element &, const Rd &, RNMK_ &) const;
};

// const int TypeOfFE_NumberSpace2d::nbNodeOnItem[4] = {0, 0, 0, 0};
// int TypeOfFE_NumberSpace2d::Data[]                = {
//     0,          // the support number  of the node of the df
//     0,          // the number of the df on  the node
//     0,          // the node of the df
//     0,          // which are de df on sub FE
//     0, 0, 0, 0, // nb node on what
//     0,          // for each compontant $j=0,N-1$ it give the sub FE associated
//     0,          // begin_dfcomp
//     0           // end_dfcomp
// };
// double TypeOfFE_NumberSpace2d::alpha_Pi_h[] = {1.};

// void TypeOfFE_NumberSpace2d::FB(const What_d whatd, const Element &K, const R2 &P, RNMK_ &val) const {
//     assert(val.N() >= 1);
//     assert(val.M() == 1);

//     val = 0;
//     RN_ f0(val('.', 0, op_id));

//     if (whatd & Fop_D0) {
//         f0[0] = 1.;
//     }
// }


// static TypeOfFE_NumberSpace2d Number_Space;
// GTypeOfFE<Mesh2> &PNumberSpace(Number_Space);
// template <> GTypeOfFE<Mesh2> &DataFE<Mesh2>::NumberSpace = Number_Space;

#endif // __NUMBER_SPACE_HPP__