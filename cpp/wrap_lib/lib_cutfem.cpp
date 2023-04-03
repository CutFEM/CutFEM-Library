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
#include "lib_cutfem.hpp"

extern "C" {

void build_mesh(pyProblem *darcy, int nx, int ny, R orx, R ory, R lx, R ly) {
    std::vector<int> nnx{nx, ny};
    std::vector<double> oorx{orx, ory};
    std::vector<double> llx{lx, ly};
    darcy->build_mesh(nnx, oorx, llx);
}

void delete_object(pyProblem *darcy) { darcy->delete_obj(); }

void init_space(pyProblem *darcy, double (*f)(double *), char *st_type) { darcy->init_space(f, std::string(st_type)); }

void set_stabilization_penalty(pyProblem *darcy, double cu, double cp) { darcy->set_stabilization_penalty(cu, cp); }

void set_stab_Cu(pyProblem *problem, double l) { problem->set_stab_Cu(l); }
void set_stab_Cp(pyProblem *problem, double l) { problem->set_stab_Cp(l); }
void set_stab_Cpu(pyProblem *problem, double l) { problem->set_stab_Cpu(l); }
void set_stab_Cwu(pyProblem *problem, double l) { problem->set_stab_Cwu(l); }
void set_stab_Cw(pyProblem *problem, double l) { problem->set_stab_Cw(l); }
void set_penalty_param(pyProblem *problem, double l) { problem->set_penalty_param(l); }

void add_bulk_integral(pyProblem *darcy, double (*f)(double *, int, int)) { darcy->add_bulk_integral(f); }

void add_interface_integral(pyProblem *darcy, double (*f)(double *, int, int)) { darcy->add_interface_integral(f); }

void add_macro_stabilization(pyProblem *darcy, double di, int i) { darcy->add_macro_stabilization(di, i); }

void add_full_stabilization(pyProblem *darcy, int i) { darcy->add_full_stabilization(i); }

void add_lagrange_multiplier_mixed(pyProblem *stokes) { stokes->add_lagrange_multiplier_mixed(); }
void add_lagrange_multiplier_classic(pyProblem *stokes, double val) { stokes->add_lagrange_multiplier_classic(val); }

void post_processing_pressure(pyProblem *stokes, double (*f)(double *, int, int)) {
    stokes->post_processing_pressure(f);
}

int get_nz(pyProblem *darcy) { return darcy->get_nz(); }

int get_size_A(pyProblem *darcy) { return darcy->get_size_A(); }

void get_COO_data(pyProblem *darcy, int32_t *r, int32_t *c, double *v, double *b) { darcy->get_COO_data(r, c, v, b); }

void get_CSR_data(pyProblem *darcy, int32_t *r, int32_t *c, double *v, double *b) { darcy->get_CSR_data(r, c, v, b); }

void give_back_sol(pyProblem *darcy, double *x) { darcy->get_back_sol(x); }

void write_vtk_file(pyProblem *darcy, char *st_type) { darcy->write_vtk_file(std::string(st_type)); }

double L2error_div(pyProblem *darcy, double (*f)(double *, int, int)) { return darcy->L2error_div(f); }

double L2error_u(pyProblem *darcy, double (*f)(double *, int, int)) { return darcy->L2error_u(f); }

double L2error_p(pyProblem *darcy, double (*f)(double *, int, int)) { return darcy->L2error_p(f); }

void solve_umfpack(pyProblem *darcy) { return darcy->solve_umfpack(); }

void set_verbose(int s) { globalVariable::verbose = s; }
}