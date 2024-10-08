#include "cpp/cutfem.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using paraview_t = Paraview<mesh_t>;

int main(int argc, char **argv) {

    globalVariable::verbose = 0;
    MPIcf cfMPI(argc, argv);

    int nx = 2;

    for (int i = 0; i < 5; ++i) {
        mesh_t Kh(nx, nx, 0., 0., 1., 1.);
        Kh.info();

        mesh_t Kh_ref = refine_barycentric(Kh);

        nx = 2 * (nx - 1) + 1;

        paraview_t writer(Kh, "meshing" + std::to_string(i) + ".vtk");
        paraview_t writer2(Kh_ref, "meshing_ref" + std::to_string(i) + ".vtk");
    }
}
