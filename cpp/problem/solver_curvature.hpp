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
#ifndef CUTFEM_CPP_PROBLEM_SOLVER_CURVATURE_HPP
#define CUTFEM_CPP_PROBLEM_SOLVER_CURVATURE_HPP

#include "baseProblem.hpp"

namespace solver {

namespace cutfem {

namespace curvature {
template <typename mesh_t> std::vector<double> solve(const CutFESpace<mesh_t> &Vh, const Interface<mesh_t> &interface) {

    const int D = mesh_t::Rd::d;
    CutFEM<mesh_t> problem(Vh);
    const auto &Kh = Vh.get_mesh();

    TestFunction<mesh_t> H(Vh, D), v(Vh, D);
    Normal n;
    Rnm Id(D, D);
    Id = 0.;

    for (int i = 0; i < D; ++i)
        Id(i, i) = 1.;

    problem.addBilinear((H, v)
                        //+ (grad(H) * n, grad(v) * n) * 1e-2
                        ,
                        interface);
    problem.addLinear(-contractProduct(Id, gradS(v)), interface);
    problem.addBilinear((jump(grad(H) * n), jump(grad(v) * n)) * 1e-1, Kh, INTEGRAL_INNER_FACET);

    // if(deg == 2) {
    //   TestFunction<Rd::d> grad2un = grad(gradun)*n;
    //   this->addEdgeIntegral(innerProduct(1e-2*h*h*jump(grad2un),
    //   jump(grad2un)));
    // }
    problem.solve();
    return problem.rhs_;
}

} // namespace curvature

} // namespace cutfem

} // namespace solver

// template <typename M> class Curvature {
//   public:
//     typedef M Mesh;
//     typedef CutFESpace<Mesh> CutSpace;
//     typedef typename CutSpace::FElement FElement;
//     typedef typename FElement::Rd Rd;
//     typedef Mapping<Mesh> IsoMapping;
//     typedef FunFEM<Mesh> Fun_h;
//     typedef TestFunction<Mesh> FunTest;

//     static const int D = Rd::d;
//     const CutSpace &Vh;
//     const ActiveMesh<Mesh> &Kh;
//     const Interface<Mesh> &interface;

//     Curvature(const CutSpace &vh, const Interface<Mesh> *inter) : Vh(vh), Kh(vh.cutTh), interface(*inter) {}
//     Curvature(const CutSpace &vh, const Interface<Mesh> &inter) : Vh(vh), Kh(vh.cutTh), interface(inter) {}

//   public:
//     std::vector<double> solve() {

//         CutFEM<Mesh2> problem(Vh);
//         FunTest H(Vh, D), v(Vh, D);
//         Normal n;
//         Rnm Id(D, D);
//         Id = 0.;

//         for (int i = 0; i < D; ++i)
//             Id(i, i) = 1.;

//         problem.addBilinear((H, v)
//                             //+ (grad(H) * n, grad(v) * n) * 1e-2
//                             ,
//                             interface);
//         problem.addLinear(-contractProduct(Id, gradS(v)), interface);
//         problem.addBilinear((jump(grad(H) * n), jump(grad(v) * n)) * 1e-1, Kh, INTEGRAL_INNER_FACET);

//         // if(deg == 2) {
//         //   TestFunction<Rd::d> grad2un = grad(gradun)*n;
//         //   this->addEdgeIntegral(innerProduct(1e-2*h*h*jump(grad2un),
//         //   jump(grad2un)));
//         // }
//         problem.solve();
//         return problem.rhs_;
//     }
//     std::vector<double> solve(const ExpressionVirtual &w) {

//         CutFEM<Mesh2> problem(Vh);
//         FunTest H(Vh, D), v(Vh, D);
//         Normal n;
//         Rnm Id(D, D);
//         Id = 0.;
//         for (int i = 0; i < D; ++i)
//             Id(i, i) = 1.;
//         //
//         // //a(u,v)_Gamma
//         // double t0 = CPUtime();
//         problem.addBilinear((H, v) //+ (gradun,gradun)*1e-2
//                             ,
//                             interface
//                             // , {}
//                             // , mapping
//         );
//         // l(v)_Omega
//         problem.addLinear(-contractProduct(Id, w * gradS(v)), interface
//                           // , {}
//                           // , mapping
//         );
//         //
//         // ListItemVF<Rd::d> Sh = (jump(gradun),jump(gradun))*1e-2;
//         // this->addEdgeIntegral(Sh);
//         problem.addBilinear((jump(grad(H) * n), jump(grad(v) * n)) * 1e-2, Kh, INTEGRAL_INNER_FACET);

//         // if(deg == 2) {
//         //   TestFunction<Rd::d> grad2un = grad(gradun)*n;
//         //   this->addEdgeIntegral(innerProduct(1e-2*h*h*jump(grad2un),
//         //   jump(grad2un)));
//         // }
//         problem.solve();
//         return problem.rhs_;
//     }
//     std::vector<double> solve(const IsoMapping &mapping) {

//         CutFEM<Mesh2> problem(Vh);
//         FunTest H(Vh, D), v(Vh, D);
//         Normal n;
//         Rnm Id(D, D);
//         Id = 0.;
//         for (int i = 0; i < D; ++i)
//             Id(i, i) = 1.;
//         //
//         // //a(u,v)_Gamma
//         // double t0 = CPUtime();
//         problem.addBilinear((H, v) //+ (gradun,gradun)*1e-2
//                             ,
//                             interface, mapping);
//         // l(v)_Omega
//         problem.addLinear(-contractProduct(Id, gradS(v)), interface, mapping

//         );
//         //
//         // ListItemVF<Rd::d> Sh = (jump(gradun),jump(gradun))*1e-2;
//         // this->addEdgeIntegral(Sh);
//         std::cout << " ADD THE H SCALING IN STABILIZATION " << std::endl;
//         if (Vh.polynomialOrder == 1) {
//             problem.addBilinear((jump(grad(H) * n), jump(grad(v) * n)) * 1e-1, Kh, INTEGRAL_INNER_FACET);
//         }
//         if (Vh.polynomialOrder == 2) {
//             problem.addBilinear((jump(grad(H) * n), jump(grad(v) * n)) * 1e-1 +
//                                     (jump(grad(grad(H) * n) * n), jump(grad(grad(v) * n) * n)) * 1e-2,
//                                 Kh, INTEGRAL_INNER_FACET);
//         }

//         // if(deg == 2) {
//         //   TestFunction<Rd::d> grad2un = grad(gradun)*n;
//         //   this->addEdgeIntegral(innerProduct(1e-2*h*h*jump(grad2un),
//         //   jump(grad2un)));
//         // }
//         problem.solve();
//         return problem.rhs_;
//     }
//     std::vector<double> solve(const ExpressionVirtual &w, const IsoMapping &mapping) {

//         CutFEM<Mesh2> problem(Vh);
//         FunTest H(Vh, D), v(Vh, D);
//         Normal n;
//         Rnm Id(D, D);
//         Id = 0.;
//         for (int i = 0; i < D; ++i)
//             Id(i, i) = 1.;
//         //
//         // //a(u,v)_Gamma
//         // double t0 = CPUtime();
//         problem.addBilinear((H, v) //+ (gradun,gradun)*1e-2
//                             ,
//                             interface, mapping);
//         // l(v)_Omega
//         problem.addLinear(-contractProduct(Id, w * gradS(v)), interface, mapping

//         );
//         std::cout << " ADD THE H SCALING IN STABILIZATION " << std::endl;
//         if (Vh.polynomialOrder == 1) {
//             problem.addBilinear((jump(grad(H) * n), jump(grad(v) * n)) * 1e-1, Kh, INTEGRAL_INNER_FACET);
//         }
//         if (Vh.polynomialOrder == 2) {
//             problem.addBilinear((jump(grad(H) * n), jump(grad(v) * n)) * 1e-1 +
//                                     (jump(grad(grad(H) * n) * n), jump(grad(grad(v) * n) * n)) * 1e-2,
//                                 Kh, INTEGRAL_INNER_FACET);
//         }

//         problem.solve();
//         return problem.rhs_;
//     }
// };

#endif
