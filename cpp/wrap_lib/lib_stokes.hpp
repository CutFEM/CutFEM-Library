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
#ifndef LIB_STOKES_HPP
#define LIB_STOKES_HPP

#include "lib_cutfem.hpp"

class pyFictitiousStokesRT : public pyProblem {

    using mesh_t       = Mesh2;
    using cutmesh_t    = ActiveMeshT2;
    using fespace_t    = FESpace2;
    using cutfespace_t = CutFESpaceT2;
    using Rd           = R2;
    using interface_t  = InterfaceLevelSet<mesh_t>;
    using macro_t      = MacroElement<mesh_t>;
    using fct_t        = FunFEM<mesh_t>;

    double mu{1.};
    double sigma{1.e-2};

    std::unique_ptr<Lagrange2> FEvelocity_p;
    std::unique_ptr<fespace_t> Velh_p, Ph_ex_p;

    enum { classic = 0, mixed = 1 };

  public:
    void delete_obj() {}
    void init_mu(double mmu) { mu = mmu; }
    void init_sigma(double s) { sigma = s; }
    void init_space(double (*f)(double *), std::string FE_type) override {

        Lh_p    = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P1);
        Ph_ex_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P4);

        FEvelocity_p = std::make_unique<Lagrange2>(4);
        Velh_p       = std::make_unique<fespace_t>(*Kh_p, *FEvelocity_p);

        fct_t levelSet(*Lh_p, f); // fun_levelSet);
        inter_p = std::make_unique<InterfaceLevelSet<mesh_t>>(*Kh_p, levelSet, 0);
        const auto &interface(*inter_p);
        if (FE_type == "BDM1") {
            Vh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::BDM1);
            Qh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P0);
        } else if (FE_type == "RT1") {
            Vh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::RT1);
            Qh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P1dc);
        } else {
            std::cout << " Cannot use this element pair " << std::endl;
            exit(EXIT_FAILURE);
        }

        Khi_p = std::make_unique<ActiveMesh<mesh_t>>(*Kh_p);
        Khi_p->truncate(*inter_p, 1);

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
        double hi = 1. / (mesh_param.nx[0] - 1);
        fct_t fh(*Velh_p, f);

        TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);

        problem.addBilinear(contractProduct(mu * grad(u), grad(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                            *Khi_p);
        problem.addLinear(innerProduct(fh.exprList(), v), *Khi_p);
        problem.addBilinear(-innerProduct(mu * average(grad(u * t) * n, 0.5, 0.5), jump(v * t)) +
                                innerProduct(jump(u * t), mu * average(grad(v * t) * n, 0.5, 0.5)) +
                                innerProduct(1. / hi * sigma * jump(u * t), jump(v * t)),
                            *Khi_p, INTEGRAL_INNER_EDGE_2D);
    }
    void add_interface_integral(double (*f)(double *, int, int)) override {
        double hi           = 1. / (mesh_param.nx[0] - 1);
        double penaltyParam = penalty_param.gamma;

        fct_t gh(*Velh_p, f);
        TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);

        problem.addBilinear(-innerProduct(mu * grad(u) * n, v)                // natural
                                + innerProduct(u, mu * grad(v) * n)           // symmetry
                                + innerProduct(1. / hi * penaltyParam * u, v) // stability
                                + innerProduct(p, v * n)                      // natural
                            ,
                            *inter_p);
        problem.addLinear(+innerProduct(gh.exprList(), mu * grad(v) * n) +
                              innerProduct(gh.exprList(), 1. / hi * penaltyParam * v),
                          *inter_p);
    }

    void add_lagrange_multiplier() override {

        TestFunction<2> p(*Ph_p, 1), v(*Wh_p, 2);

        CutFEM<mesh_t> lagr(*Wh_p);
        lagr.add(*Ph_p);
        lagr.addLinear(innerProduct(1., p), *Khi_p);
        std::vector<double> lag_row(lagr.rhs_.begin(), lagr.rhs_.end());
        lagr.rhs_ = 0.;
        lagr.addLinear(innerProduct(1, v * n), *inter_p);
        problem.addLagrangeVecToRowAndCol(lag_row, lagr.rhsSpan(), 0);
    }

    void add_macro_stabilization(double dlt_i, int stab_method) override {
        double Cu = stab_param.Cu;
        double Cp = stab_param.Cp;
        double hi = 1. / (mesh_param.nx[0] - 1);
        MacroElement<mesh_t> macro(*Khi_p, dlt_i);

        Normal n;
        TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);
        TestFunction<2> grad2un = grad(grad(u) * n) * n;

        problem.addFaceStabilization(innerProduct(Cu * pow(hi, -1) * jump(u), jump(v)) +
                                         innerProduct(Cu * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
                                         innerProduct(Cu * pow(hi, 3) * jump(grad2un), jump(grad2un)) -
                                         innerProduct(Cp * pow(hi, 1) * jump(p), jump(div(v))) +
                                         innerProduct(Cp * pow(hi, 1) * jump(div(u)), jump(q)) -
                                         innerProduct(Cp * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                                         innerProduct(Cp * pow(hi, 3) * jump(grad(div(v))), jump(grad(q))),
                                     *Khi_p, macro);
    }
    void add_full_stabilization(int stab_method) override {
        double Cu = stab_param.Cu;
        double Cp = stab_param.Cp;
        double hi = 1. / (mesh_param.nx[0] - 1);

        Normal n;
        TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);
        TestFunction<2> grad2un = grad(grad(u) * n) * n;

        problem.addFaceStabilization(innerProduct(Cu * pow(hi, -1) * jump(u), jump(v)) +
                                         innerProduct(Cu * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
                                         innerProduct(Cu * pow(hi, 3) * jump(grad2un), jump(grad2un)) -
                                         innerProduct(Cp * pow(hi, 1) * jump(p), jump(div(v))) +
                                         innerProduct(Cp * pow(hi, 1) * jump(div(u)), jump(q)) -
                                         innerProduct(Cp * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                                         innerProduct(Cp * pow(hi, 3) * jump(grad(div(v))), jump(grad(q))),
                                     *Khi_p);
    }

    void post_processing_pressure(double (*f)(double *, int, int)) override {

        int nb_flux_dof           = Wh_p->get_nb_dof();
        std::span<double> data_ph = std::span<double>(problem.rhs_.data() + nb_flux_dof, Ph_p->get_nb_dof());
        fct_t ph(*Ph_p, data_ph);
        fct_t exactp(*Ph_ex_p, f);

        double meanP    = integral(*Khi_p, exactp.expr(), 0);
        double meanPfem = integral(*Khi_p, ph.expr(), 0);

        TestFunction<2> q(*Ph_p, 1);

        CutFEM<mesh_t> post(*Ph_p);
        post.addLinear(innerProduct(1, q), *Khi_p);
        double area = post.rhs_.sum();
        ph.v -= meanPfem / area;
        ph.v += meanP / area;
    }
};

class pyFictitiousStokesVorticity : public pyProblem {

    using mesh_t       = Mesh2;
    using cutmesh_t    = ActiveMeshT2;
    using fespace_t    = FESpace2;
    using cutfespace_t = CutFESpaceT2;
    using Rd           = R2;
    using interface_t  = InterfaceLevelSet<mesh_t>;
    using macro_t      = MacroElement<mesh_t>;
    using fct_t        = FunFEM<mesh_t>;

    double mu{1.};
    double sigma{1.e-2};

    std::unique_ptr<Lagrange2> FEvelocity_p;
    std::unique_ptr<fespace_t> Velh_p, Ph_ex_p;
    std::unique_ptr<fespace_t> Bh_p, Uh_p;

    enum { primal = 0, mixed = 1 };

  public:
    void delete_obj() { std::cout << "Delete Vorticity Stokes object " << std::endl; }
    void init_mu(double mmu) { mu = mmu; }
    void init_sigma(double s) { sigma = s; }

    void init_space(double (*f)(double *), std::string FE_type) override {

        Lh_p    = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P1);
        Ph_ex_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P4);

        FEvelocity_p = std::make_unique<Lagrange2>(4);
        Velh_p       = std::make_unique<fespace_t>(*Kh_p, *FEvelocity_p);

        fct_t levelSet(*Lh_p, f); // fun_levelSet);
        inter_p = std::make_unique<InterfaceLevelSet<mesh_t>>(*Kh_p, levelSet, 0);
        const auto &interface(*inter_p);
        if (FE_type == "RT0") {
            Bh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P1);
            Vh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::RT0);
            Qh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P0);

        } else if (FE_type == "RT1") {
            Bh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P2);
            Vh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::RT1);
            Qh_p = std::make_unique<fespace_t>(*Kh_p, DataFE<mesh_t>::P1dc);
        } else {
            std::cout << " Cannot use this element pair " << std::endl;
            exit(EXIT_FAILURE);
        }

        Khi_p = std::make_unique<ActiveMesh<mesh_t>>(*Kh_p);
        Khi_p->truncate(*inter_p, 1);

        if (globalVariable::verbose > 0) {
            Khi_p->info();
        }
        Uh_p = std::make_unique<cutfespace_t>(*Khi_p, *Bh_p);
        if (globalVariable::verbose > 0) {
            std::cout << " Vorticity space : " << std::endl;
            Uh_p->info();
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
        problem.initSpace(*Uh_p);
        problem.add(*Wh_p);
        problem.add(*Ph_p);
    }

    void add_bulk_integral(double (*f)(double *, int, int)) override {
        const auto &Khi(*Khi_p.get());
        double hi = 1. / (mesh_param.nx[0] - 1);

        fct_t fh(*Velh_p, f);

        TestFunction<2> w(*Uh_p, 1, 0), tau(*Uh_p, 1, 0), p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);

        problem.addBilinear( // w = curl u
            innerProduct(1. / mu * w, tau) - innerProduct(u, rotgrad(tau)), Khi);
        problem.addBilinear( // mu Delta u + grad p
            innerProduct(rotgrad(w), v) - innerProduct(p, div(v)), Khi);
        problem.addLinear(+innerProduct(fh.exprList(2), v), Khi);
        problem.addBilinear(+innerProduct(div(u), q), Khi);
    }

    void add_interface_integral(double (*f)(double *, int, int)) override {
        double hi           = 1. / (mesh_param.nx[0] - 1);
        double penaltyParam = penalty_param.gamma;

        fct_t u0(*Velh_p, f);
        TestFunction<2> w(*Uh_p, 1, 0), tau(*Uh_p, 1, 0), p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);

        problem.addBilinear(+innerProduct(p, v * n) + innerProduct(1. / hi * penaltyParam * u * n, v * n), *inter_p);
        problem.addLinear(+innerProduct(u0 * t, tau) + innerProduct(u0 * n, 1. / hi * penaltyParam * v * n), *inter_p);
    }

    void add_lagrange_multiplier() override {

        TestFunction<2> p(*Ph_p, 1), v(*Wh_p, 2);

        CutFEM<mesh_t> lagr(*Uh_p);
        lagr.add(*Wh_p);
        lagr.add(*Ph_p);
        lagr.addLinear(innerProduct(1., p), *Khi_p);
        std::vector<double> lag_row(lagr.rhs_.begin(), lagr.rhs_.end());
        lagr.rhs_ = 0.;
        lagr.addLinear(innerProduct(1, v * n), *inter_p);
        problem.addLagrangeVecToRowAndCol(lag_row, lagr.rhsSpan(), 0);
    }

    void add_macro_stabilization(double dlt_i, int stab_method) override {
        double Cu, Cpu, Cp, Cw, Cwu;
        if (stab_method == primal) {
            Cu  = stab_param.Cu;
            Cpu = stab_param.Cpu;
            Cwu = 0.;
            Cp  = 0.;
            Cw  = stab_param.Cw;
        } else if (stab_method == mixed) {
            Cu  = 0.;
            Cwu = stab_param.Cwu;
            Cpu = stab_param.Cpu;
            Cp  = 0.;
            Cw  = 0.;
        } else {
            std::cout << " Stabilization method not found" << std::endl;
            return;
        }

        double hi = 1. / (mesh_param.nx[0] - 1);
        MacroElement<mesh_t> macro(*Khi_p, dlt_i);
        Normal n;
        TestFunction<2> w(*Uh_p, 1, 0), tau(*Uh_p, 1, 0), p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);
        TestFunction<2> grad2un = grad(grad(u) * n) * n;
        TestFunction<2> grad2wn = grad(grad(w) * n) * n;
        problem.addFaceStabilization(innerProduct(Cw * pow(hi, 3) * jump(grad(w) * n), jump(grad(tau) * n)) +
                                         innerProduct(Cw * pow(hi, 5) * jump(grad2wn), jump(grad2wn)) +

                                         +innerProduct(Cwu * pow(hi, 1) * jump(rotgrad(w)), jump(v)) -
                                         innerProduct(Cwu * pow(hi, 1) * jump(u), jump(rotgrad(tau))) +
                                         innerProduct(Cwu * pow(hi, 3) * jump(grad(rotgrad(w))), jump(grad(v))) -
                                         innerProduct(Cwu * pow(hi, 3) * jump(grad(u)), jump(grad(rotgrad(tau))))

                                         + innerProduct(Cu * pow(hi, 1) * jump(u), jump(v)) +
                                         innerProduct(Cu * pow(hi, 3) * jump(grad(u) * n), jump(grad(v) * n)) +
                                         innerProduct(Cu * pow(hi, 5) * jump(grad2un), jump(grad2un))

                                         - innerProduct(Cpu * pow(hi, 1) * jump(p), jump(div(v))) +
                                         innerProduct(Cpu * pow(hi, 1) * jump(div(u)), jump(q)) -
                                         innerProduct(Cpu * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                                         innerProduct(Cpu * pow(hi, 3) * jump(grad(div(u))), jump(grad(q))),
                                     *Khi_p, macro);
    }

    void add_full_stabilization(int stab_method) override {
        double Cu, Cpu, Cp, Cw, Cwu;
        if (stab_method == primal) {
            Cu  = stab_param.Cu;
            Cpu = stab_param.Cpu;
            Cwu = 0.;
            Cp  = 0.;
            Cw  = stab_param.Cw;
        } else if (stab_method == mixed) {
            Cu  = 0.;
            Cwu = stab_param.Cwu;
            Cpu = stab_param.Cpu;
            Cp  = 0.;
            Cw  = 0.;
        } else {
            std::cout << " Stabilization method not found" << std::endl;
            return;
        }

        double hi = 1. / (mesh_param.nx[0] - 1);
        Normal n;
        TestFunction<2> w(*Uh_p, 1, 0), tau(*Uh_p, 1, 0), p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);
        TestFunction<2> grad2un = grad(grad(u) * n) * n;
        TestFunction<2> grad2wn = grad(grad(w) * n) * n;
        problem.addFaceStabilization(innerProduct(Cw * pow(hi, 3) * jump(grad(w) * n), jump(grad(tau) * n)) +
                                         innerProduct(Cw * pow(hi, 5) * jump(grad2wn), jump(grad2wn)) +

                                         +innerProduct(Cwu * pow(hi, 1) * jump(rotgrad(w)), jump(v)) -
                                         innerProduct(Cwu * pow(hi, 1) * jump(u), jump(rotgrad(tau))) +
                                         innerProduct(Cwu * pow(hi, 3) * jump(grad(rotgrad(w))), jump(grad(v))) -
                                         innerProduct(Cwu * pow(hi, 3) * jump(grad(u)), jump(grad(rotgrad(tau))))

                                         + innerProduct(Cu * pow(hi, 1) * jump(u), jump(v)) +
                                         innerProduct(Cu * pow(hi, 3) * jump(grad(u) * n), jump(grad(v) * n)) +
                                         innerProduct(Cu * pow(hi, 5) * jump(grad2un), jump(grad2un))

                                         - innerProduct(Cpu * pow(hi, 1) * jump(p), jump(div(v))) +
                                         innerProduct(Cpu * pow(hi, 1) * jump(div(u)), jump(q)) -
                                         innerProduct(Cpu * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                                         innerProduct(Cpu * pow(hi, 3) * jump(grad(div(u))), jump(grad(q))),
                                     *Khi_p);
    }

    void post_processing_pressure(double (*f)(double *, int, int)) override {

        int nb_vort_dof = Uh_p->get_nb_dof();
        int nb_flux_dof = Vh_p->get_nb_dof();
        std::span<double> data_ph =
            std::span<double>(problem.rhs_.data() + nb_flux_dof + nb_vort_dof, Ph_p->get_nb_dof());
        fct_t ph(*Ph_p, data_ph);
        fct_t p0(*Ph_ex_p, f);
        double meanP    = integral(*Khi_p, p0.expr(), 0);
        double meanPfem = integral(*Khi_p, ph.expr(), 0);
        TestFunction<2> q(*Ph_p, 1, 0);
        CutFEM<mesh_t> post(*Ph_p.get());

        post.addLinear(innerProduct(1, q), *Khi_p);
        double area = post.rhs_.sum();
        ph.v -= meanPfem / area;
        ph.v += meanP / area;
    }

    void write_vtk_file(std::string filename) override {
        int ndof_w = Uh_p->get_nb_dof();
        int ndof_u = Wh_p->get_nb_dof();
        int ndof_p = Ph_p->get_nb_dof();

        Rn_ data_wh = sub_array(problem.rhs_, 0, ndof_w);
        Rn_ data_uh = sub_array(problem.rhs_, ndof_w, ndof_u);
        Rn_ data_ph = sub_array(problem.rhs_, ndof_w + ndof_u, ndof_p);

        fct_t wh(*Uh_p.get(), data_wh);
        fct_t uh(*Wh_p.get(), data_uh);
        fct_t ph(*Ph_p.get(), data_ph);
        auto femSol_0dx = dx(uh.expr(0));
        auto femSol_1dy = dy(uh.expr(1));

        Paraview<mesh_t> writer(*Khi_p, filename);

        writer.add(uh, "velocity", 0, 2);
        writer.add(ph, "pressure", 0, 1);
        writer.add(wh, "vorticity", 0, 1);
        writer.add(femSol_0dx + femSol_1dy, "divergence");
    }

    double L2error_div(double (*f)(double *, int, int)) override {

        Rn_ data_uh = sub_array(problem.rhs_, Uh_p->get_nb_dof(), Wh_p->get_nb_dof());
        fct_t uh(*Wh_p, data_uh);
        auto femSol_0dx = dx(uh.expr(0));
        auto femSol_1dy = dy(uh.expr(1));

        return L2normCut(femSol_0dx + femSol_1dy, f, *Khi_p);
    }

    double L2error_u(double (*f)(double *, int, int)) override {
        Rn_ data_uh = sub_array(problem.rhs_, Uh_p->get_nb_dof(), Wh_p->get_nb_dof());
        fct_t uh(*Wh_p, data_uh);
        return L2normCut(uh, f, 0, 2);
    }

    double L2error_p(double (*f)(double *, int, int)) override {
        Rn_ data_ph = sub_array(problem.rhs_, Uh_p->get_nb_dof() + Wh_p->get_nb_dof(), Ph_p->get_nb_dof());
        fct_t ph(*Ph_p.get(), data_ph);
        return L2normCut(ph, f, 0, 1);
    }
};

typedef pyFictitiousStokesRT FictitiousStokesRT_2;
typedef pyFictitiousStokesVorticity FictitiousStokesVorticity_2;
#endif