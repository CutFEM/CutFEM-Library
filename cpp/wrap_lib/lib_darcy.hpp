#ifndef LIB_DARCY_HPP
#define LIB_DARCY_HPP
#include <memory>

#include "../problem/baseProblem.hpp"
#include "../FESpace/paraview.hpp"
#include "../num/matlab.hpp"

struct KernelDarcy2 {
   typedef Mesh2 Mesh;
   typedef ActiveMeshT2 CutMesh;
   typedef FESpace2 Space;
   typedef CutFESpaceT2 CutSpace;
   typedef R2 Rd;

   static const int D = 2;
};

template <typename Kernel> class Darcy {

   typedef typename Kernel::Mesh Mesh;
   typedef typename Kernel::Space Space;
   typedef typename Kernel::CutSpace CutSpace;
   typedef typename Kernel::Rd Rd;
   typedef InterfaceLevelSet<Mesh> Interface;

   const int D = Kernel::D;

   std::shared_ptr<Mesh> Kh_p;
   std::shared_ptr<ActiveMesh<Mesh>> Khi_p;
   std::shared_ptr<Space> Vh_p, Qh_p;
   std::shared_ptr<CutSpace> Wh_p, Ph_p;
   std::shared_ptr<Interface> inter_p;
   std::shared_ptr<MacroElement<Mesh>> macro_p;

   CutFEM<Mesh> darcy;
   Kernel kernel;

   struct {
      int type           = 1;
      double Cu          = 1.;
      double Cp          = 1.;
      double delta_macro = 0.25;
   } stab_param;

   struct {
      std::vector<int> nx{11, 11, 11};
      std::vector<double> orx{0., 0., 0.};
      std::vector<double> lx{1., 1., 1.};
   } mesh_param;

 public:
   void build_mesh(int *nnx, double *oorx, double *llx) {
      for (int i = 0; i < D; ++i) {
         mesh_param.nx[i]  = nnx[i];
         mesh_param.orx[i] = oorx[i];
         mesh_param.lx[i]  = llx[i];
      }
      this->init_mesh();
   }
   void init_mesh();
   void delete_obj() {}
   void set_stabilization_type(std::string st_type) {
      if (st_type == "macro") {
         stab_param.type = 1;
         return;
      } else if (st_type == "full") {
         stab_param.type = 2;
         return;
      } else {
         stab_param.type = 0;
         return;
      }
   }

   void init_space(double (*f)(double *), std::string FE_type) {

      Space Lh(*Kh_p, DataFE<Mesh>::P1);

      FunFEM<Mesh> levelSet(Lh, f); // fun_levelSet);
      inter_p = std::make_shared<Interface>(*Kh_p, levelSet);
      const auto &interface(*inter_p);

      if (FE_type == "RT0") {
         Vh_p = std::make_shared<Space>(*Kh_p, DataFE<Mesh>::RT0);
         Qh_p = std::make_shared<Space>(*Kh_p, DataFE<Mesh>::P0);
      } else if (FE_type == "BDM1") {
         Vh_p = std::make_shared<Space>(*Kh_p, DataFE<Mesh>::BDM1);
         Qh_p = std::make_shared<Space>(*Kh_p, DataFE<Mesh>::P0);
      } else if (FE_type == "RT1") {
         Vh_p = std::make_shared<Space>(*Kh_p, DataFE<Mesh>::RT1);
         Qh_p = std::make_shared<Space>(*Kh_p, DataFE<Mesh>::P1dc);
      }

      Khi_p = std::make_shared<ActiveMesh<Mesh>>(*Kh_p, *inter_p);
      if (globalVariable::verbose > 0) {
         Khi_p->info();
      }
      Wh_p = std::make_shared<CutSpace>(*Khi_p, *Vh_p);
      if (globalVariable::verbose > 0)
         Wh_p->info();
      Ph_p = std::make_shared<CutSpace>(*Khi_p, *Qh_p);
      if (globalVariable::verbose > 0)
         Ph_p->info();

      darcy.initSpace(*Wh_p);
      darcy.add(*Ph_p);
   }

   void add_bulk_integral(double (*f)(double *, int, int)) {

      TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);
      FunFEM<Mesh> fq(*Ph_p, f);
      darcy.addBilinear(innerProduct(u, v) - innerProduct(p, div(v)) +
                            innerProduct(div(u), q),
                        *Khi_p);

      darcy.addLinear(innerProduct(fq.expression(), q), *Khi_p);
   }
   void add_interface_integral(double (*f)(double *, int, int)) {

      double xi   = 3. / 4;
      double xi0  = (xi - 0.5) / 2.;
      double mu_G = 2. / 3 * 0.250001;
      FunFEM<Mesh> phat(*Ph_p, f);

      Normal n;
      TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);

      darcy.addBilinear(innerProduct(mu_G * average(u * n), average(v * n)) +
                            innerProduct(xi0 * mu_G * jump(u * n), jump(v * n)),
                        *inter_p);
      darcy.addLinear(-innerProduct(phat.expression(), jump(v * n)), *inter_p);
   }
   void add_natural_BC(double (*f)(double *, int, int)) {

      Space Lh(*Kh_p, DataFE<Mesh>::P2);
      FunFEM<Mesh> p0(Lh, f);

      Normal n;
      TestFunction<2> v(*Wh_p, 2);
      darcy.addLinear(
          -innerProduct(p0.expression(), v * n) // Only on Gamma_N (pressure)
          ,
          *Khi_p, INTEGRAL_BOUNDARY);
   }

   void set_stabilization_penalty(double ccu, double ccp) {
      stab_param.Cu = ccu;
      stab_param.Cp = ccp;
   }
   void add_macro_stabilization(double dlt_i) {
      stab_param.delta_macro = dlt_i;
      double Cu              = stab_param.Cu;
      double Cp              = stab_param.Cp;
      double h_i             = 1. / (mesh_param.nx[0] - 1);
      MacroElement<Mesh> macro(*Khi_p, stab_param.delta_macro);

      Normal n;
      TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);
      TestFunction<2> grad2un = grad(grad(u) * n) * n;

      darcy.addFaceStabilization(
          innerProduct(Cu * h_i * jump(u),
                       jump(v)) // [Method 1: Remove jump in vel]
              + innerProduct(Cu * pow(h_i, 3) * jump(grad(u) * n),
                             jump(grad(v) * n)) +
              innerProduct(Cu * pow(h_i, 5) * jump(grad2un), jump(grad2un)) -
              innerProduct(Cp * h_i * jump(p), jump(div(v))) +
              innerProduct(Cp * h_i * jump(div(u)), jump(q)) -
              innerProduct(Cp * pow(h_i, 3) * jump(grad(p)),
                           jump(grad(div(v)))) +
              innerProduct(Cp * pow(h_i, 3) * jump(grad(div(v))),
                           jump(grad(q))),
          *Khi_p, macro);
   }
   void add_full_stabilization() {
      double Cu  = stab_param.Cu;
      double Cp  = stab_param.Cp;
      double h_i = 1. / (mesh_param.nx[0] - 1);

      Normal n;
      TestFunction<2> p(*Ph_p, 1), q(*Ph_p, 1), u(*Wh_p, 2), v(*Wh_p, 2);
      TestFunction<2> grad2un = grad(grad(u) * n) * n;

      darcy.addFaceStabilization(
          innerProduct(Cu * h_i * jump(u),
                       jump(v)) // [Method 1: Remove jump in vel]
              + innerProduct(Cu * pow(h_i, 3) * jump(grad(u) * n),
                             jump(grad(v) * n)) +
              innerProduct(Cu * pow(h_i, 5) * jump(grad2un), jump(grad2un)) -
              innerProduct(Cp * h_i * jump(p), jump(div(v))) +
              innerProduct(Cp * h_i * jump(div(u)), jump(q)) -
              innerProduct(Cp * pow(h_i, 3) * jump(grad(p)),
                           jump(grad(div(v)))) +
              innerProduct(Cp * pow(h_i, 3) * jump(grad(div(v))),
                           jump(grad(q))),
          *Khi_p);
   }

   void solve_umfpack() { darcy.solve("umfpack"); }
   int get_nz() { return darcy.mat_[0].size(); }
   int get_size_A() { return darcy.rhs_.size(); }
   void get_CSR_data(int32_t *r, int32_t *c, double *v, double *b) {
      int n = darcy.rhs_.size();
      buil_CSR_array(n, darcy.mat_[0], r, c, v);
      for (int i = 0; i < n; ++i)
         b[i] = darcy.rhs_[i];
   }
   void get_COO_data(int32_t *r, int32_t *c, double *v, double *b) {
      int i = 0;
      int n = darcy.rhs_.size();
      for (const auto &aij : darcy.mat_[0]) {
         r[i] = static_cast<int32_t>(aij.first.first);
         c[i] = static_cast<int32_t>(aij.first.second);
         v[i] = aij.second;
         i++;
      }
      for (int i = 0; i < n; ++i)
         b[i] = darcy.rhs_[i];
   }
   void get_back_sol(double *x) {
      for (int i = 0; i < darcy.rhs_.size(); ++i)
         darcy.rhs_[i] = x[i];
   }

   double L2error_div(double (*f)(double *, int, int)) {
      Rn_ data_uh = sub_array(darcy.rhs_, 0, Wh_p->get_nb_dof());
      FunFEM<Mesh> uh(*Wh_p, data_uh);
      ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
      ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

      return L2normCut(femSol_0dx + femSol_1dy, f, *Khi_p);
   }
   double L2error_u(double (*f)(double *, int, int)) {
      Rn_ data_uh = sub_array(darcy.rhs_, 0, Wh_p->get_nb_dof());
      FunFEM<Mesh> uh(*Wh_p, data_uh);
      return L2normCut(uh, f, 0, 2);
   }
   double L2error_p(double (*f)(double *, int, int)) {
      Rn_ data_ph =
          sub_array(darcy.rhs_, Wh_p->get_nb_dof(), Ph_p->get_nb_dof());
      FunFEM<Mesh> ph(*Ph_p, data_ph);
      return L2normCut(ph, f, 0, 1);
   }
   void write_vtk_file(std::string filename) {
      int ndof_u  = Wh_p->get_nb_dof();
      int ndof_p  = Ph_p->get_nb_dof();
      Rn_ data_uh = sub_array(darcy.rhs_, 0, ndof_u);
      Rn_ data_ph = sub_array(darcy.rhs_, ndof_u, ndof_p);

      FunFEM<Mesh> uh(*Wh_p, data_uh);
      FunFEM<Mesh> ph(*Ph_p, data_ph);
      ExpressionFunFEM<Mesh> femSol_0dx(uh, 0, op_dx);
      ExpressionFunFEM<Mesh> femSol_1dy(uh, 1, op_dy);

      Paraview<Mesh> writer(*Khi_p, filename);

      writer.add(uh, "velocity", 0, 2);
      writer.add(ph, "pressure", 0, 1);
      writer.add(femSol_0dx + femSol_1dy, "divergence");
   }
};

template <> void Darcy<KernelDarcy2>::init_mesh() {
   Kh_p = std::make_shared<Mesh>(mesh_param.nx[0], mesh_param.nx[1],
                                 mesh_param.orx[0], mesh_param.orx[1],
                                 mesh_param.lx[0], mesh_param.lx[1]);
}

typedef Darcy<KernelDarcy2> Darcy2;

#endif
