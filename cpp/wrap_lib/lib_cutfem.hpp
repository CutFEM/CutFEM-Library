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
#ifndef LIB_CUTFEM_HPP
#define LIB_CUTFEM_HPP

#include "../tool.hpp"

class pySolverTool {
    using mesh_t = Mesh2;

  public:
    CutFEM<mesh_t> problem;

    void solve_umfpack() { problem.solve("umfpack"); }
    int get_nz() { return problem.mat_[0].size(); }
    int get_size_A() { return problem.rhs_.size(); }
    void get_CSR_data(int32_t *r, int32_t *c, double *v, double *b) {
        int n = problem.rhs_.size();
        buil_CSR_array(n, problem.mat_[0], r, c, v);
        for (int i = 0; i < n; ++i)
            b[i] = problem.rhs_[i];
    }
    void get_COO_data(int32_t *r, int32_t *c, double *v, double *b) {
        int i = 0;
        int n = problem.rhs_.size();
        for (const auto &aij : problem.mat_[0]) {
            r[i] = static_cast<int32_t>(aij.first.first);
            c[i] = static_cast<int32_t>(aij.first.second);
            v[i] = aij.second;
            i++;
        }
        for (int i = 0; i < n; ++i)
            b[i] = problem.rhs_[i];
    }
    void get_back_sol(double *x) {
        for (int i = 0; i < problem.rhs_.size(); ++i)
            problem.rhs_[i] = x[i];
    }
};

struct pyMixProblemTool {

    using mesh_t       = Mesh2;
    using cutmesh_t    = ActiveMeshT2;
    using fespace_t    = FESpace2;
    using cutfespace_t = CutFESpaceT2;
    using Rd           = R2;
    using interface_t  = InterfaceLevelSet<mesh_t>;
    using macro_t      = MacroElement<mesh_t>;
    using fct_t        = FunFEM<mesh_t>;

    std::unique_ptr<mesh_t> Kh_p;
    std::unique_ptr<cutmesh_t> Khi_p;
    std::unique_ptr<fespace_t> Lh_p, Vh_p, Qh_p;
    std::unique_ptr<cutfespace_t> Wh_p, Ph_p;
    std::unique_ptr<interface_t> inter_p;
    std::unique_ptr<macro_t> macro_p;

    Normal n;
    Tangent t;

    struct {
        std::vector<int> nx{11, 11};
        std::vector<double> orx{0., 0.};
        std::vector<double> lx{1., 1.};
    } mesh_param;

    struct {
        int type{1};
        double Cu{0.};
        double Cp{0.};
        double Cpu{0.};
        double Cw{0.};
        double Cwu{0.};
    } stab_param;

    struct {
        double gamma = 4e3;
    } penalty_param;

    void build_mesh(std::span<int> nnx, std::span<double> oorx, std::span<double> llx) {
        for (int i = 0; i < 2; ++i) {
            mesh_param.nx[i]  = nnx[i];
            mesh_param.orx[i] = oorx[i];
            mesh_param.lx[i]  = llx[i];
        }
        Kh_p = std::make_unique<mesh_t>(mesh_param.nx[0], mesh_param.nx[1], mesh_param.orx[0], mesh_param.orx[1],
                                        mesh_param.lx[0], mesh_param.lx[1]);
    }

    void set_stabilization_penalty(double ccu, double ccp) {
        stab_param.Cu = ccu;
        stab_param.Cp = ccp;
    }

    void set_stab_Cu(double l) { stab_param.Cu = l; }
    void set_stab_Cp(double l) { stab_param.Cp = l; }
    void set_stab_Cpu(double l) { stab_param.Cpu = l; }
    void set_stab_Cw(double l) { stab_param.Cw = l; }
    void set_stab_Cwu(double l) { stab_param.Cwu = l; }
    void set_penalty_param(double l) { penalty_param.gamma = l; }
};

class pyProblem : public pyMixProblemTool, public pySolverTool {
    using mesh_t       = Mesh2;
    using cutmesh_t    = ActiveMeshT2;
    using fespace_t    = FESpace2;
    using cutfespace_t = CutFESpaceT2;
    using Rd           = R2;
    using interface_t  = InterfaceLevelSet<mesh_t>;
    using macro_t      = MacroElement<mesh_t>;
    using fct_t        = FunFEM<mesh_t>;

  public:
    virtual double L2error_div(double (*f)(double *, int, int)) {
        Rn_ data_uh = sub_array(problem.rhs_, 0, Wh_p->get_nb_dof());
        fct_t uh(*Wh_p, data_uh);
        auto femSol_0dx = dx(uh.expr(0));
        auto femSol_1dy = dy(uh.expr(1));

        return L2normCut(femSol_0dx + femSol_1dy, f, *Khi_p);
    }
    virtual double L2error_u(double (*f)(double *, int, int)) {
        Rn_ data_uh = sub_array(problem.rhs_, 0, Wh_p->get_nb_dof());
        fct_t uh(*Wh_p, data_uh);
        return L2normCut(uh, f, 0, 2);
    }
    virtual double L2error_p(double (*f)(double *, int, int)) {
        Rn_ data_ph = sub_array(problem.rhs_, Wh_p->get_nb_dof(), Ph_p->get_nb_dof());
        fct_t ph(*Ph_p.get(), data_ph);
        return L2normCut(ph, f, 0, 1);
    }

    virtual void write_vtk_file(std::string filename) {
        int ndof_u  = Wh_p->get_nb_dof();
        int ndof_p  = Ph_p->get_nb_dof();
        Rn_ data_uh = sub_array(problem.rhs_, 0, ndof_u);
        Rn_ data_ph = sub_array(problem.rhs_, ndof_u, ndof_p);

        fct_t uh(*Wh_p.get(), data_uh);
        fct_t ph(*Ph_p.get(), data_ph);
        auto femSol_0dx = dx(uh.expr(0));
        auto femSol_1dy = dy(uh.expr(1));

        Paraview<mesh_t> writer(*Khi_p, filename);

        writer.add(uh, "velocity", 0, 2);
        writer.add(ph, "pressure", 0, 1);
        writer.add(femSol_0dx + femSol_1dy, "divergence");
    }

    virtual void add_lagrange_multiplier_mixed() {

        TestFunction<2> p(*Ph_p, 1), v(*Wh_p, 2);

        CutFEM<mesh_t> lagr(*Wh_p);
        lagr.add(*Ph_p);
        lagr.addLinear(innerProduct(1., p), *Khi_p);
        std::vector<double> lag_row(lagr.rhs_.begin(), lagr.rhs_.end());
        lagr.rhs_ = 0.;
        lagr.addLinear(innerProduct(1, v * n), *inter_p);
        problem.addLagrangeVecToRowAndCol(lag_row, lagr.rhsSpan(), 0);
    }
    virtual void add_lagrange_multiplier_classic(double val) {
        TestFunction<2> p(*Ph_p, 1);
        problem.addLagrangeMultiplier(innerProduct(1., p), val, *Khi_p);
    }

    virtual void init_space(double (*f)(double *), std::string FE_type)    = 0;
    virtual void add_bulk_integral(double (*f)(double *, int, int))        = 0;
    virtual void add_interface_integral(double (*f)(double *, int, int))   = 0;
    virtual void add_macro_stabilization(double dlt_i, int stab_method)    = 0;
    virtual void add_full_stabilization(int stab_method)                   = 0;
    // virtual void add_lagrange_multiplier_classic()                         = 0;
    // virtual void add_lagrange_multiplier_mixed()                           = 0;
    virtual void post_processing_pressure(double (*f)(double *, int, int)) = 0;
    virtual void delete_obj()                                              = 0;
    virtual ~pyProblem()                                                   = default;
};

#endif