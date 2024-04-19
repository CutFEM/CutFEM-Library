#include "../cutfem.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;

double shift        = 0.5;
double interfaceRad = sqrt(0.25);

double fun_levelSet(R2 P) { return sqrt((P.x - shift) * (P.x - shift) + (P.y - shift) * (P.y - shift)) - interfaceRad; }

double fun_rhs(R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 40 * x * x * x - 40 * x * y * y - 32 * y + 16;
    else
        return -40 * x * x * y + 32 * x + 40 * y * y * y - 16;
}

double fun_exact_u(R2 P, int i, int dom) {
    double x = P.x;
    double y = P.y;
    if (i == 0)
        return 2 * (x * x - x + 0.25 + y * y - y) * (2 * y - 1);
    else
        return -2 * (x * x - x + 0.25 + y * y - y) * (2 * x - 1);
}

double fun_exact_p(R2 P, int i, int dom) {
    double x = P.x;
    double y = P.y;
    return 10 * (x * x - y * y) * (x * x - y * y);
}

int main(int argc, char **argv) {

    MPIcf cfMPI(argc, argv);

    int nx              = 11;
    double penaltyParam = 4e3;
    double sigma        = 1e-2;
    double uPenParam    = 1e0;
    double pPenParam    = 1e0;
    double mu           = 1.;

    Normal n;
    Tangent t;

    std::vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    globalVariable::verbose = 2;

    for (int i = 0; i < 4; ++i) {

        double hi = 1. / (nx - 1);

        mesh_t Kh(nx, nx, -0.0000001, -0.00000001, 1.0000002, 1.000002);
        space_t Lh(Kh, DataFE<Mesh2>::P1);
        fct_t levelSet(Lh, fun_levelSet);

        InterfaceLevelSet<mesh_t> interface(Kh, levelSet);

        cutmesh_t Khi(Kh);
        Khi.truncate(interface, 1);

        MacroElement<mesh_t> macro(Khi, 1.);

        Lagrange2 FEvelocity(4);
        space_t VELh(Kh, FEvelocity);
        space_t SCAh(Kh, DataFE<mesh_t>::P2);

        space_t Wh(Kh, DataFE<mesh_t>::RT1);
        cutspace_t Vh(Khi, Wh);

        space_t Qh(Kh, DataFE<mesh_t>::P1dc);
        cutspace_t Ph(Khi, Qh);

        fct_t fh(VELh, fun_rhs);
        fct_t gh(VELh, fun_exact_u);
        fct_t exactp(SCAh, fun_exact_p);

        CutFEM<mesh_t> stokes(Vh);
        stokes.add(Ph);

        funtest_t u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);

        stokes.addBilinear(contractProduct(mu * grad(u), grad(v)) - innerProduct(p, div(v)) + innerProduct(div(u), q),
                           Khi);
        stokes.addLinear(innerProduct(fh.exprList(), v), Khi);
        // [NECESSARY FOR RT / BDM ELEMENTS]
        stokes.addBilinear(-innerProduct(mu * average(grad(u * t) * n, 0.5, 0.5), jump(v * t)) +
                               innerProduct(jump(u * t), mu * average(grad(v * t) * n, 0.5, 0.5)) +
                               innerProduct(1. / hi * sigma * jump(u * t), jump(v * t)),
                           Khi, INTEGRAL_INNER_EDGE_2D);

        stokes.addBilinear(-innerProduct(mu * grad(u) * n, v)                // natural
                               + innerProduct(u, mu * grad(v) * n)           // symmetry
                               + innerProduct(1. / hi * penaltyParam * u, v) // stability
                               + innerProduct(p, v * n)                      // natural
                           ,
                           interface);
        stokes.addLinear(+innerProduct(gh.exprList(), mu * grad(v) * n) +
                             innerProduct(gh.exprList(), 1. / hi * penaltyParam * v),
                         interface);

        // funtest_t grad2un = grad(grad(u) * n) * n;
        // stokes.addFaceStabilization(                                  // [h^(2k+1) h^(2k+1)]
        //     +innerProduct(uPenParam * pow(hi, -1) * jump(u), jump(v)) // [Method 1: Remove jump in vel]
        //         + innerProduct(uPenParam * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n))
        //         // +innerProduct(uPenParam*pow(hi,3)*jump(grad2un), jump(grad2un))
        //         - innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
        //         innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
        //         innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
        //         innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(v))), jump(grad(q))),
        //     Khi, macro);

        stokes.addPatchStabilization(+innerProduct(uPenParam * pow(hi, -2) * jump(u), jump(v)) -
                                         innerProduct(pPenParam * pow(hi, 0) * jump(p), jump(div(v))) +
                                         innerProduct(pPenParam * pow(hi, 0) * jump(div(u)), jump(q)),
                                     Khi);

        // [Sets uniqueness of the pressure]
        CutFEM<Mesh2> lagr(Vh);
        lagr.add(Ph);
        lagr.addLinear(innerProduct(1., p), Khi);
        std::vector<double> lag_row(lagr.rhs_.begin(), lagr.rhs_.end());
        std::fill(lagr.rhs_.begin(), lagr.rhs_.end(), 0.);

        lagr.addLinear(innerProduct(1, v * n), interface);

        stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhsSpan(), 0);

        matlab::Export(stokes.mat_[0], "mat" + std::to_string(i) + "Cut.dat");

        stokes.solve("umfpack");

        int nb_flux_dof           = Vh.get_nb_dof();
        std::span<double> data_uh = std::span<double>(stokes.rhs_.data(), nb_flux_dof);
        std::span<double> data_ph = std::span<double>(stokes.rhs_.data() + nb_flux_dof, Ph.get_nb_dof());
        fct_t uh(Vh, data_uh);
        fct_t ph(Ph, data_ph);

        // [Post process pressure]
        double meanP    = integral(Khi, exactp.expr(), 0);
        double meanPfem = integral(Khi, ph.expr(), 0);
        CutFEM<Mesh2> post(Ph);
        post.addLinear(innerProduct(1, q), Khi);
        // double area = post.rhs_.sum();
        double area = std::accumulate(post.rhs_.begin(), post.rhs_.end(), 0);

        ph.v -= meanPfem / area;
        ph.v += meanP / area;

        // Plotting
        // {
        //     Paraview<mesh_t> writer(Khi, "stokes_" + std::to_string(i) + ".vtk");
        //     writer.add(uh, "velocity", 0, 2);
        //     writer.add(ph, "pressure", 0, 1);
        //     writer.add(dx(uh.expr(0)) + dy(uh.expr(1)), "divergence");
        // }
        auto uh_0dx = dx(uh.expr(0));
        auto uh_1dy = dy(uh.expr(1));

        double errU      = L2normCut(uh, fun_exact_u, 0, 2);
        double errP      = L2normCut(ph, fun_exact_p, 0, 1);
        double errDiv    = L2normCut(uh_0dx + uh_1dy, Khi);
        double maxErrDiv = maxNormCut(uh_0dx + uh_1dy, Khi);

        ul2.push_back(errU);
        pl2.push_back(errP);
        divl2.push_back(errDiv);
        divmax.push_back(maxErrDiv);
        h.push_back(hi);
        if (i == 0) {
            convu.push_back(0);
            convp.push_back(0);
        } else {
            convu.push_back(log(ul2[i] / ul2[i - 1]) / log(h[i] / h[i - 1]));
            convp.push_back(log(pl2[i] / pl2[i - 1]) / log(h[i] / h[i - 1]));
        }

        nx = 2 * nx - 1;
    }

    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ') << "err_p"
              << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ') << "err u"
              << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15) << std::setfill(' ')
                  << pl2[i] << std::setw(15) << std::setfill(' ') << convp[i] << std::setw(15) << std::setfill(' ')
                  << ul2[i] << std::setw(15) << std::setfill(' ') << convu[i] << std::setw(15) << std::setfill(' ')
                  << divl2[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << divmax[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
