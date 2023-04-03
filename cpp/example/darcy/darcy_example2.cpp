#include "../tool.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t::D>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;

double sq_SW   = -1. + 1e-10;    // 1.1 for all convg tables
double sq_LGTH = 2 * 1. + 2e-10; // 1.1 for all convg tables

double c1   = 100;
double c2   = 42;
double c3   = 200;
double a3   = 50;
double c    = 10.;
double c0   = 0; // c= 10., 7. ; c0 = 1e2,7e2
double mu_G = 2. * c / 1e2;
double fun_levelSet(R2 P, ) {
    double x = P[0], y = P[1];
    return (c1 * pow(x, 6) + 20 * x * x - c2 * y * y + c3 * pow(y, 4) - c);
}
double fun_dxg(R2 P) {
    double x = P[0], y = P[1];
    return (6 * c1 * pow(x, 5) + 20 * 2 * x);
}
double fun_dyg(R2 P) {
    double x = P[0], y = P[1];
    return (-2 * c2 * y + 4 * c3 * y * y * y);
}
double fun_normgradg(R2 P) {
    double x = P[0], y = P[1];
    double g = fun_levelSet(P, 0);
    return sqrt(pow(fun_dxg(P), 2) + pow(fun_dyg(P), 2) + c0 * g * g);
}
double fun_normal(R2 P, const int i) {
    double x = P[0], y = P[1];
    double dxg       = fun_dxg(P);
    double dyg       = fun_dyg(P);
    double normgradg = fun_normgradg(P);
    return -((i == 0) * dxg / normgradg + (i == 1) * dyg / normgradg);
}

double fun_exact_p(R2 P, int compInd, int dom) {
    return (dom == 0) * (fun_levelSet(P, 0) + c) / 1e2; // 1,0
}
double fun_force(R2 P, int compInd, int dom) {
    double x = P[0], y = P[1];
    return (dom == 0) * ((compInd == 0) * (fun_dxg(P) / 1e2 + fun_normal(P, 0)) +
                         (compInd == 1) * (fun_dyg(P) / 1e2 + fun_normal(P, 1)));
}
double fun_neumann(R2 P, int compInd, int dom) {
    return (dom == 0) * ((compInd == 0) * (fun_normal(P, 0)) + (compInd == 1) * (fun_normal(P, 1)));
}
double fun_exact_u(R2 P, int compInd, int dom) {
    return (dom == 0) * ((compInd == 0) * (fun_normal(P, 0)) + (compInd == 1) * (fun_normal(P, 1)));
}
double fun_div(R2 P, int compInd, int dom) {
    double x = P[0], y = P[1];
    double val =
        (-12 * c3 * y * y + 2 * c2) / sqrt(pow(6 * c1 * pow(x, 5) + 40 * x, 2) +
                                           c0 * pow(c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) - c2 * y * y - c, 2) +
                                           pow(-4 * c3 * y * y * y + 2 * c2 * y, 2)) -
        (30 * c1 * pow(x, 4) + 40) / sqrt(pow(6 * c1 * pow(x, 5) + 40 * x, 2) +
                                          c0 * pow(c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) - c2 * y * y - c, 2) +
                                          pow(-4 * c3 * y * y * y + 2 * c2 * y, 2)) -
        ((-4 * c3 * y * y * y + 2 * c2 * y) * (2 * (-4 * c3 * y * y * y + 2 * c2 * y) * (-12 * c3 * y * y + 2 * c2) -
                                               2 * c0 * (-4 * c3 * y * y * y + 2 * c2 * y) *
                                                   (c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) - c2 * y * y - c))) /
            (2 * pow(sqrt(pow(6 * c1 * pow(x, 5) + 40 * x, 2) +
                          c0 * pow(c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) - c2 * y * y - c, 2) +
                          pow(-4 * c3 * y * y * y + 2 * c2 * y, 2)),
                     3)) +
        ((2 * (6 * c1 * pow(x, 5) + 40 * x) * (30 * c1 * pow(x, 4) + 40) +
          2 * c0 * (6 * c1 * pow(x, 5) + 40 * x) * (c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) - c2 * y * y - c)) *
         (6 * c1 * pow(x, 5) + 40 * x)) /
            (2 * pow(sqrt(pow(6 * c1 * pow(x, 5) + 40 * x, 2) +
                          c0 * pow(c1 * pow(x, 6) + 20 * x * x + c3 * pow(y, 4) - c2 * y * y - c, 2) +
                          pow(-4 * c3 * y * y * y + 2 * c2 * y, 2)),
                     3));

    return (dom == 0) * val;
}
double fun_interfacePr(R2 P) { return (-1. / 8 * mu_G + c / (2. * 1e2)); }

int main(int argc, char **argv) {

    globalVariable::verbose = 0;
#ifdef USE_MPI
    MPIcf cfMPI(argc, argv);
#endif
