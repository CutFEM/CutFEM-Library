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
#ifndef INTERPOLATION_HPP_
#define INTERPOLATION_HPP_

#include "../common/RNM.hpp"
#include "../common/global.hpp"
#include "FESpace.hpp"

// #include "../parallel/cfmpi.hpp"

/*
Interpolate f : Rd->R    on space Vh
- output : fh contains the values
*/
template <Space F, FunctionScalar fct_t> void interpolate(const F &Mh, KN_<double> fh, fct_t f) {
    // std::cout << " need to double check this interpolate function and add "
    //           << std::endl;
    // getchar();
    // << std::endl; assert(0);
    typedef typename F::Rd Rd;
    typedef typename F::Element::RdHat RdHat;
    // Rn fhSend(Mh.nbDoF); fhSend = 1e+50;
    fh.init(Mh.nbDoF);
    const int d   = 1;
    const int nve = Mh.MaxNbNodePerElement;
    KNM<R> Vpf(nve, 1);              // value of f at the interpolation points
    KN<R> ggf(Mh.MaxNbDFPerElement); // stock the values of the dof of the interpolate

    progress bar(" Interpolating", Mh.NbElement(), globalVariable::verbose);
    // for (int t=Mh.first_element();t<Mh.last_element();
    //      t+= Mh.next_element()) {      // loop over element
    for (int t = 0; t < Mh.NbElement(); t += 1) {
        bar += 1;
        typename F::FElement K(Mh[t]);
        const int nbdf = K.NbDoF(); // nof local

        for (int p = 0; p < K.tfe->NbPtforInterpolation; p++) { // all interpolation points
            Rd P(K.Pt(p));                                      // the coordinate of P in K hat
            Vpf(p, 0) = f(P);
        }

        K.Pi_h(Vpf, ggf);
        for (int df = 0; df < nbdf; df++) {
            // fhSend[K(df)] =  ggf[df] ;
            fh[K(df)] = ggf[df];
        }
    }

    bar.end();
    // MPIcf::AllReduce(fhSend, fh, MPI_MIN);
}

/*
Interpolate f : Rd->R    on space Vh
- output : fh contains the values
*/

template <Space F, FunctionLevelSet fct_t> void interpolate(const F &Mh, KN_<double> fh, fct_t f) {
    // std::cout << " need to double check this interpolate function and add MPI"
    // << std::endl; assert(0);
    typedef typename F::Rd Rd;
    typedef typename F::Element::RdHat RdHat;
    // Rn fhSend(Mh.nbDoF); fhSend = 1e+50;
    assert(fh.size() == Mh.nbDoF);
    // fh.init(Mh.nbDoF);
    const int d   = Mh.N;
    const int nve = Mh.TFE(0)->NbPtforInterpolation;
    KNM<R> Vpf(nve, d);              // value of f at the interpolation points
    KN<R> ggf(Mh.MaxNbDFPerElement); // stock the values of the dof of the interpolate

    progress bar(" Interpolating", Mh.NbElement(), globalVariable::verbose);

    // for (int t=Mh.first_element();t<Mh.last_element();
    //      t+= Mh.next_element()) {      // loop over element
    for (int t = 0; t < Mh.NbElement(); t += 1) {
        bar += 1;
        typename F::FElement K(Mh[t]);
        const int nbdf = K.NbDoF(); // nof local

        for (int p = 0; p < K.tfe->NbPtforInterpolation; p++) { // all interpolation points
            Rd P(K.Pt(p));                                      // the coordinate of P in K hat

            for (int i = 0; i < d; ++i) {
                Vpf(p, i) = f(P, i);
            }
        }

        K.Pi_h(Vpf, ggf);
        for (int df = 0; df < nbdf; df++) {
            // fhSend[K(df)] =  ggf[df] ;
            fh[K(df)] = ggf[df];
        }
    }
    // MPIcf::AllReduce(fhSend, fh, MPI_MIN);
    bar.end();
}

/*
Interpolate f : Rd->R    on space Vh
- output : fh contains the values
*/
template <Space F, FunctionDomain fct> void interpolate(const F &Mh, KN_<double> fh, fct f) {
    typedef typename F::Rd Rd;
    typedef typename F::Element::RdHat RdHat;
    assert(fh.size() == Mh.nbDoF);

    const int d   = Mh.N;
    const int nve = Mh.TFE(0)->NbPtforInterpolation;
    KNM<R> Vpf(nve, d);              // value of f at the interpolation points
    KN<R> ggf(Mh.MaxNbDFPerElement); // stock the values of the dof of the interpolate
    progress bar(" Interpolating", Mh.NbElement(), globalVariable::verbose);

    for (int t = 0; t < Mh.NbElement(); t += 1) {
        bar += 1;
        typename F::FElement K(Mh[t]);
        const int nbdf   = K.NbDoF(); // nof local
        const int domain = K.get_domain();

        for (int p = 0; p < K.tfe->NbPtforInterpolation; p++) { // all interpolation points
            Rd P = K.Pt(p);                                     // the coordinate of P in K hat
            for (int i = 0; i < d; ++i) {
                Vpf(p, i) = f((double *)P, i, domain);
            }
        }

        K.Pi_h(Vpf, ggf);

        for (int df = 0; df < nbdf; df++) {
            fh[K(df)] = ggf[df];
        }
    }

    bar.end();
}

/*
Interpolate f : Rd->R    on space Vh
- output : fh contains the values
*/
template <Space F, FunctionDomainTime fct> void interpolate(const F &Mh, KN_<double> fh, fct f, R tid) {
    // std::cout << " need to double check this interpolate function and add MPI"
    // << std::endl; assert(0);
    typedef typename F::Rd Rd;
    typedef typename F::Element::RdHat RdHat;
    // Rn fhSend(Mh.nbDoF); fhSend = 1e+50;
    assert(fh.size() == Mh.nbDoF);
    // fh.init(Mh.nbDoF);
    const int d   = Mh.N;
    const int nve = Mh.TFE(0)->NbPtforInterpolation;
    KNM<R> Vpf(nve, d);              // value of f at the interpolation points
    KN<R> ggf(Mh.MaxNbDFPerElement); // stock the values of the dof of the interpolate
    progress bar(" Interpolating", Mh.NbElement(), globalVariable::verbose);

    // for (int t=Mh.first_element();t<Mh.last_element();
    //      t+= Mh.next_element()) {      // loop over element
    for (int t = 0; t < Mh.NbElement(); t += 1) {
        bar += 1;
        typename F::FElement K(Mh[t]);
        const int nbdf   = K.NbDoF(); // nof local
        const int domain = K.get_domain();

        for (int p = 0; p < K.tfe->NbPtforInterpolation; p++) { // all interpolation points
            Rd P(K.Pt(p));                                      // the coordinate of P in K hat
            for (int i = 0; i < d; ++i) {
                Vpf(p, i) = f(P, i, domain, tid);
            }
        }

        K.Pi_h(Vpf, ggf);

        for (int df = 0; df < nbdf; df++) {
            fh[K(df)] = ggf[df];
        }
    }
    bar.end();
    // MPIcf::AllReduce(fhSend, fh, MPI_MIN);
}

/*
Interpolate f : Rd->R    on space Vh
- output : fh contains the values
*/
template <Space F> void interpolate(const F &Mh, KN_<double> fh, R (*f)(double *, int, R), R tid) {
    // std::cout << " need to double check this interpolate function and add MPI"
    // << std::endl; assert(0);
    typedef typename F::Rd Rd;
    typedef typename F::Element::RdHat RdHat;
    // Rn fhSend(Mh.nbDoF); fhSend = 1e+50;
    assert(fh.size() == Mh.nbDoF);
    // fh.init(Mh.nbDoF);
    const int d   = Mh.N;
    const int nve = Mh.TFE(0)->NbPtforInterpolation;
    KNM<R> Vpf(nve, d);              // value of f at the interpolation points
    KN<R> ggf(Mh.MaxNbDFPerElement); // stock the values of the dof of the interpolate
    progress bar(" Interpolating", Mh.NbElement(), globalVariable::verbose);

    // for (int t=Mh.first_element();t<Mh.last_element();
    //      t+= Mh.next_element()) {      // loop over element
    for (int t = 0; t < Mh.NbElement(); t += 1) {
        bar += 1;
        typename F::FElement K(Mh[t]);
        const int nbdf = K.NbDoF(); // nof local

        for (int p = 0; p < K.tfe->NbPtforInterpolation; p++) { // all interpolation points
            Rd P(K.Pt(p));                                      // the coordinate of P in K hat
            for (int i = 0; i < d; ++i) {
                Vpf(p, i) = f(P, i, tid);
            }
        }

        K.Pi_h(Vpf, ggf);

        for (int df = 0; df < nbdf; df++) {
            // fhSend[K(df)] =  ggf[df] ;
            fh[K(df)] = ggf[df];
        }
    }
    bar.end();
    // MPIcf::AllReduce(fhSend, fh, MPI_MIN);
}

/*
Interpolate f : Rd->R    on space time Vh
- output : fh contains the values
*/
template <Space F> void interpolate(const F &Mh, const TimeSlab &In, KN_<double> fh, R (*f)(double *, int, R)) {
    // std::cout << " need to double check this interpolate function and add MPI"
    // << std::endl; assert(0);
    typedef typename F::Rd Rd;
    typedef typename F::Element::RdHat RdHat;
    // Rn fhSend(Mh.nbDoF); fhSend = 1e+50;
    assert(fh.size() == Mh.nbDoF * In.NbDoF());
    // fh.init(Mh.nbDoF);
    const int d   = Mh.N;
    const int nve = Mh.TFE(0)->NbPtforInterpolation;
    const int nvt = In.tfe->NbPtforInterpolation;
    KNM<R> Vpt(nvt, 1);
    KNM<R> ggft(Mh.MaxNbDFPerElement,
                In.NbDoF()); // stock the values of the dof of the interpolate

    KNMK<R> Vpft(nve, d, nvt); // value of f at the interpolation points
    progress bar(" Interpolating", Mh.NbElement(), globalVariable::verbose);

    // for (int t=Mh.first_element();t<Mh.last_element();
    //      t+= Mh.next_element()) {      // loop over element
    for (int t = 0; t < Mh.NbElement(); t += 1) {
        bar += 1;
        typename F::FElement K(Mh[t]);
        const int nbdf = K.NbDoF(); // nof local

        // compute all the value , space and time
        for (int it = 0; it < In.tfe->NbPtforInterpolation; ++it) {
            const R1 &tq(In.Pt(it));
            for (int p = 0; p < K.tfe->NbPtforInterpolation; p++) { // all interpolation points
                Rd P(K.Pt(p));                                      // the coordinate of P in K hat
                for (int i = 0; i < d; ++i) {
                    Vpft(p, i, it) = f(P, i, tq);
                }
            }
        }

        // perfom the spacxe interpolation for each time dof
        for (int it = 0; it < In.tfe->NbPtforInterpolation; ++it) {
            KN_<R> ggg(ggft('.', it));
            K.Pi_h(Vpft('.', '.', it), ggg);
        }

        // perfom time interpolation and save/replacve the value in the matrix by
        // the good one
        for (int df = 0; df < nbdf; df++) {
            for (int it = 0; it < In.tfe->NbPtforInterpolation; ++it) {
                Vpt(it, 0) = ggft(df, it);
            }
            KN_<R> ggg(ggft(df, '.'));
            In.Pi_h(Vpt, ggg);
        }

        for (int it = 0; it < In.tfe->NbPtforInterpolation; ++it) {
            for (int df = 0; df < nbdf; df++) {
                fh[K.loc2glb(df, it)] = ggft(df, it); //[df] ;
            }
        }
    }
    bar.end();
    // MPIcf::AllReduce(fhSend, fh, MPI_MIN);
}

template <typename Mesh>
void interpolateOnBackGroundMesh(FunFEM<Mesh> &uh, const FunFEM<Mesh> &fh, const FunFEM<Mesh> &ls) {

    using cutmesh_t = ActiveMesh<Mesh>;
    using Rd        = typename Mesh::Rd;

    const auto &Vh_cut = *fh.Vh;
    const auto &cutTh  = Vh_cut.get_mesh();
    const auto &Vh     = *uh.Vh;

    uh.v = 0.;
    for (int k = 0; k < Vh.NbElement(); ++k) {
        const auto &FK(Vh[k]);

        auto idx_K = cutTh.idxAllElementFromBackMesh(k, -1);
        if (idx_K.size() > 1) {
            for (int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
                Rd x     = FK.Pt(i);
                R lsval  = ls.eval(k, x);
                int kcut = (lsval > 0) ? idx_K[0] : idx_K[1];
                for (int ci = 0; ci < Rd::d; ++ci) {
                    uh(FK(i + FK.dfcbegin(ci))) = fh.eval(kcut, x, ci);
                }
            }
        } else {
            int kcut = idx_K[0];
            for (int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
                Rd x = FK.Pt(i);
                for (int ci = 0; ci < Rd::d; ++ci) {
                    uh(FK(i + FK.dfcbegin(ci))) = fh.eval(kcut, x, ci);
                }
            }
        }
    }
}

template <typename Mesh>
void interpolateOnBackGroundMesh(FunFEM<Mesh> &uh, const FunFEM<Mesh> &fh, const FunFEM<Mesh> &ls, double tt) {

    using cutmesh_t = ActiveMesh<Mesh>;
    using Rd        = typename Mesh::Rd;

    const auto &Vh_cut = *fh.Vh;
    const auto &cutTh  = Vh_cut.get_mesh();
    const auto &Vh     = *uh.Vh;

    uh.v = 0.;
    for (int k = 0; k < Vh.NbElement(); ++k) {
        const auto &FK(Vh[k]);

        auto idx_K = cutTh.idxAllElementFromBackMesh(k, -1);
        if (idx_K.size() > 1) {
            for (int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
                Rd x     = FK.Pt(i);
                R lsval  = ls.eval(k, x);
                int kcut = (lsval > 0) ? idx_K[0] : idx_K[1];
                for (int ci = 0; ci < Rd::d; ++ci) {
                    uh(FK(i + FK.dfcbegin(ci))) = fh.eval(kcut, x, tt, ci, op_id);
                }
            }
        } else {
            int kcut = idx_K[0];
            for (int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
                Rd x = FK.Pt(i);
                for (int ci = 0; ci < Rd::d; ++ci) {
                    uh(FK(i + FK.dfcbegin(ci))) = fh.eval(kcut, x, tt, ci, op_id);
                }
            }
        }
    }
}

#endif
