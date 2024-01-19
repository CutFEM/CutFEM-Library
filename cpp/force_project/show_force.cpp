#include "../force_project/util.hpp"
#include "../force_project/csvfile.hpp"
#include "../force_project/yaml_reader.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;
using gamma_t    = InterfaceLevelSet<mesh_t>;
using paraview_t = Paraview<mesh_t>;

using namespace force_project;

int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);
    CutFEMLogger::initialize("../cpp/force_project/data/force_test_log.txt");

    YamlReader input_data("../cpp/force_project/data/output.yaml");

    YamlReaderNode Xc_yaml = input_data["Xc"];
    YamlReaderNode Yc_yaml = input_data["Yc"];

    std::vector<double> Xc(Xc_yaml.num_children());
    std::vector<double> Yc(Yc_yaml.num_children());

    for (int i = 0; auto &ei : Xc)
        Xc_yaml[i++] >> ei;

    for (int i = 0; auto &ei : Yc)
        Yc_yaml[i++] >> ei;

    MeshDimensions box;
    int nx = 50, ny = 50;
    int number_of_time_step = 120;

    double r_0 = 0.25;
    // double xc  = 0.1;
    // double yc  = 0.5;
    double dt  = 0.01;

    // Define Mesh
    mesh_t Kh(nx, ny, box.x_min, box.y_min, box.lx, box.ly);
    space_t Lh(Kh, DataFE<Mesh2>::P1);

    int n = Xc.size();

    for (int i = 0, ifig = 0; i < n; ++i) {

        double xc = Xc[i];
        double yc = Yc[i];

        auto ls = [xc, yc, r_0](R2 P) { return sqrt((P.x - xc) * (P.x - xc) + (P.y - yc) * (P.y - yc)) - r_0; };
        fct_t levelSet(Lh, ls);

        if (i % 2 == 0) {
            paraview_t writer(Kh, "AI_geometry" + std::to_string(ifig++) + ".vtk");
            writer.add(levelSet, "levelSet0", 0, 1);
        }
    }
}
