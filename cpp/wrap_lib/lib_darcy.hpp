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
#ifndef LIB_DARCY_HPP
#define LIB_DARCY_HPP
#include "lib_cutfem.hpp"

class pyDarcy : public pyProblem {

    using mesh_t       = Mesh2;
    using cutmesh_t    = ActiveMeshT2;
    using fespace_t    = FESpace2;
    using cutfespace_t = CutFESpaceT2;
    using Rd           = R2;
    using interface_t  = InterfaceLevelSet<mesh_t>;
    using macro_t      = MacroElement<mesh_t>;
    using fct_t        = FunFEM<mesh_t>;

    enum { classic = 0, mixed = 1 };

  public:
    void delete_obj() override {}
    void init_space(double (*f)(double *), std::string FE_type) override {

        fespace_t Lh(*Kh_p, DataFE<mesh_t>::P1);

        fct_t levelSet(Lh, f); // fun_levelSet);
        inter_p = std::make_unique<InterfaceLevelSet<mesh_t>>(*Kh_p, levelSet, 0);
        const auto &interface(*inter_p);

        if (FE_type == "RT0") {
            Vh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::RT0);
            Qh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P0);
        } else if (FE_type == "BDM1") {
            Vh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::BDM1);
            Qh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P0);
        } else if (FE_type == "RT1") {
            Vh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::RT1);
            Qh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P1dc);
        }

        Khi_p = std::make_unique<ActiveMesh<mesh_t>>(*Kh_p, *inter_p);
        if (globalVariable::verbose > 0) {
            Khi_p->info();
        }
        Wh_p = std::make_unique<cutfespace_t>(*Khi_p, *Vh_p);
        if (globalVariable::verbose > 0) {
            std::cout << " Velocity space : " << std::endl;
            Wh_p->info();
        }
        Ph_p = std::make_unique<cutfespace_t>(*Khi_p, *Qh_p);
        if (globalVariable::verbose > 0) {
            std::cout << " Pressure space : " << std::endl;
            Ph_p->info();
        }
        problem.initSpace(*Wh_p);
        problem.add(*Ph_p);
    }

    void add_bulk_integral(double (*f)(double *, int, int)) override {

        TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);
        fct_t fq(*Ph_p, f);
        problem.addBilinear(innerProduct(u, v) - innerProduct(p, div(v)) + innerProduct(div(u), q), *Khi_p);

        problem.addLinear(innerProduct(fq.expr(), q), *Khi_p);
    }
    void add_interface_integral(double (*f)(double *, int, int)) override {

        double xi   = 3. / 4;
        double xi0  = (xi - 0.5) / 2.;
        double mu_G = 2. / 3 * 0.250001;
        fct_t phat(*Ph_p, f);

        Normal n;
        TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);

        problem.addBilinear(innerProduct(mu_G * average(u * n), average(v * n)) +
                                innerProduct(xi0 * mu_G * jump(u * n), jump(v * n)),
                            *inter_p);
        problem.addLinear(-innerProduct(phat.expr(), jump(v * n)), *inter_p);
    }
    void add_natural_BC(double (*f)(double *, int, int)) {

        fespace_t Lh(*Kh_p, DataFE<mesh_t>::P2);
        fct_t p0(Lh, f);

        Normal n;
        TestFunction<2> v(*Wh_p, 2);
        problem.addLinear(-innerProduct(p0.expr(), v * n) // Only on Gamma_N (pressure)
                          ,
                          *Khi_p, INTEGRAL_BOUNDARY);
    }

    // void add_lagrange_multiplier_classic() override {}
    // void add_lagrange_multiplier_mixed() override {}
    void post_processing_pressure(double (*f)(double *, int, int)) override {}

    void add_macro_stabilization(double dlt_i, int stab_method) {
        double Cu  = stab_param.Cu;
        double Cp  = stab_param.Cp;
        double h_i = 1. / (mesh_param.nx[0] - 1);
        MacroElement<mesh_t> macro(*Khi_p, dlt_i);

        Normal n;
        TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);
        TestFunction<2> grad2un = grad(grad(u) * n) * n;

        problem.addFaceStabilization(innerProduct(Cu * h_i * jump(u),
                                                  jump(v)) // [Method 1: Remove jump in vel]
                                         + innerProduct(Cu * pow(h_i, 3) * jump(grad(u) * n), jump(grad(v) * n)) +
                                         innerProduct(Cu * pow(h_i, 5) * jump(grad2un), jump(grad2un)) -
                                         innerProduct(Cp * h_i * jump(p), jump(div(v))) +
                                         innerProduct(Cp * h_i * jump(div(u)), jump(q)) -
                                         innerProduct(Cp * pow(h_i, 3) * jump(grad(p)), jump(grad(div(v)))) +
                                         innerProduct(Cp * pow(h_i, 3) * jump(grad(div(v))), jump(grad(q))),
                                     *Khi_p, macro);
    }
    void add_full_stabilization(int stab_method) {
        double Cu  = stab_param.Cu;
        double Cp  = stab_param.Cp;
        double h_i = 1. / (mesh_param.nx[0] - 1);

        Normal n;
        TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);
        TestFunction<2> grad2un = grad(grad(u) * n) * n;

        problem.addFaceStabilization(innerProduct(Cu * h_i * jump(u),
                                                  jump(v)) // [Method 1: Remove jump in vel]
                                         + innerProduct(Cu * pow(h_i, 3) * jump(grad(u) * n), jump(grad(v) * n)) +
                                         innerProduct(Cu * pow(h_i, 5) * jump(grad2un), jump(grad2un)) -
                                         innerProduct(Cp * h_i * jump(p), jump(div(v))) +
                                         innerProduct(Cp * h_i * jump(div(u)), jump(q)) -
                                         innerProduct(Cp * pow(h_i, 3) * jump(grad(p)), jump(grad(div(v)))) +
                                         innerProduct(Cp * pow(h_i, 3) * jump(grad(div(v))), jump(grad(q))),
                                     *Khi_p);
    }
};

typedef pyDarcy Darcy2;

#endif
