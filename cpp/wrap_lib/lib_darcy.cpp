#include "lib_darcy.hpp"


extern "C"
{

  Darcy2* Darcy2_new(Darcy2* darcy){ return new Darcy2();}
  void Darcy2_build_mesh(Darcy2* darcy, int nx, int ny, R orx, R ory, R lx, R ly) {
    std::vector<int> nnx{nx,ny};
    std::vector<double> oorx{orx,ory};
    std::vector<double> llx{lx,ly};
    darcy->build_mesh(nnx.data(),oorx.data(),llx.data());
  }
  void Darcy2_delete(Darcy2* darcy) {
    darcy->delete_obj();
  }
  void Darcy2_init_space(Darcy2* darcy, char* st_type ) {
    darcy->init_space(std::string(st_type));
  }
  void Darcy2_add_bulk_integral(Darcy2* darcy, double (*f)(double *, int, int) ) {
    darcy->add_bulk_integral(f);
  }
  void Darcy2_add_interface_integral(Darcy2* darcy, double (*f)(double *, int, int) ) {
    darcy->add_interface_integral(f);
  }
  void Darcy2_add_natural_BC(Darcy2* darcy, double (*f)(double *, int, int) ) {
    darcy->add_natural_BC(f);
  }
  void Darcy2_set_stabilization_penalty(Darcy2* darcy, double cu, double cp){
    darcy->set_stabilization_penalty(cu,cp);
  }
  void Darcy2_add_macro_stabilization(Darcy2* darcy, double di){
    darcy->add_macro_stabilization(di);
  }
  void Darcy2_add_full_stabilization(Darcy2* darcy){
    darcy->add_full_stabilization();
  }
  int Darcy2_get_nz(Darcy2* darcy) {
    return darcy->get_nz();
  }
  int Darcy2_get_size_A(Darcy2* darcy) {
    return darcy->get_size_A();
  }
  void Darcy2_get_COO_data(Darcy2* darcy, int32_t* r, int32_t* c, double* v, double* b) {
    darcy->get_COO_data(r,c,v,b);
  }
  void Darcy2_get_CSR_data(Darcy2* darcy, int32_t* r, int32_t* c, double* v, double* b) {
    darcy->get_CSR_data(r,c,v,b);
  }
  void Darcy2_give_back_sol(Darcy2* darcy, double* x) {
    darcy->get_back_sol(x);
  }

  void Darcy2_write_vtk_file(Darcy2* darcy, char* st_type ) {
    darcy->write_vtk_file(std::string(st_type));
  }
  double Darcy2_L2error_div(Darcy2* darcy, double (*f)(double *, int, int) ) {
    return darcy->L2error_div(f);
  }
  double Darcy2_L2error_u(Darcy2* darcy, double (*f)(double *, int, int) ) {
    return darcy->L2error_u(f);
  }
  double Darcy2_L2error_p(Darcy2* darcy, double (*f)(double *, int, int) ) {
    return darcy->L2error_p(f);
  }
  void Darcy2_solve_umfpack(Darcy2* darcy) { return darcy->solve_umfpack();}
  void set_verbose(int s) { globalVariable::verbose = s;}
}
