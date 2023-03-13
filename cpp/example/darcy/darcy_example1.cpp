#include "../tool.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t::D>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;

double sq_SW        = 0. + 1e-10;
double sq_LGTH      = 1. + 2e-10;
double interfaceRad = 0.250001;
double mu_G         = 2. / 3 * interfaceRad;
R2 shift(0.5, 0.5);

double fun_radius2(R2 P) { return norm2_2(P - shift); }
double fun_levelSet(R2 P, const int i) { return norm2(P - shift) - interfaceRad; }
double fun_dirichlet(R2 P, int compInd) { return 0; }

double fun_neumann(R2 P, int compInd, int dom) {
    double r2      = fun_radius2(P);
    double radius2 = interfaceRad * interfaceRad;
    return r2 / (2 * radius2) + 3. / 2.;
}
double fun_force(R2 P, int compInd) { return 0; }

double fun_div(R2 P, int compInd, int dom) {
    double radius2 = interfaceRad * interfaceRad;
    if (dom == 0)
        return -2. / radius2;
    else
        return -4. / radius2;
}
double fun_exact_u(R2 P, int compInd, int dom) {
    DA<double, 2> X(P[0], 0), Y(P[1], 1);
    DA<double, 2> r2  = (X - shift.x) * (X - shift.x) + (Y - shift.y) * (Y - shift.y);
    double radius2    = interfaceRad * interfaceRad;
    double cst        = (dom == 0) * 3. / 2;
    double mul        = (dom == 0) * 2 + (dom == 1) * 1;
    DA<double, 2> val = r2 / (mul * radius2) + cst;
    return -val.d[compInd];
}

double fun_exact_p(R2 P, int compInd, int dom) {
    DA<double, 2> X(P[0], 0), Y(P[1], 1);
    DA<double, 2> r2  = (X - shift.x) * (X - shift.x) + (Y - shift.y) * (Y - shift.y);
    double radius2    = interfaceRad * interfaceRad;
    double cst        = (dom == 0) * 3. / 2;
    double mul        = (dom == 0) * 2 + (dom == 1) * 1;
    DA<double, 2> val = r2 / (mul * radius2) + cst;
    return val.val;
}
double fun_interfacePr(R2 P, int compInd) { return 19. / 12; }

int main(int argc, char **argv) {

    globalVariable::verbose = 0;
#ifdef USE_MPI
    MPIcf cfMPI(argc, argv);
#endif

    int nx = 11;

    std::vector<double> uPrint, pPrint, divPrint, divPrintLoc, maxDivPrint, h, convuPr, convpPr, convdivPr,
        convdivPrLoc, convmaxdivPr;
    std::vector<double> ratioCut1, ratioCut2;
    int iters = 3;

    for (int i = 0; i < iters; ++i) {

        mesh_t Kh(nx, nx, sq_SW, sq_SW, sq_LGTH, sq_LGTH);

        space_t Lh(Kh, DataFE<mesh_t>::P1);
        fct_t levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<mesh_t> interface(Kh, levelSet);

        space_t Vh(Kh, DataFE<mesh_t>::RT0);
        space_t Qh(Kh, DataFE<mesh_t>::P0);

        cutmesh_t Kh_i(Kh, interface);
        cutspace_t Wh(Kh_i, Vh);
        cutspace_t Ph(Kh_i, Qh);

        CutFEM<mesh_t> darcy(Wh);
        darcy.add(Ph);

        const double h_i  = 2. / (nx - 1);
        const double invh = 1. / h_i;
        MacroElement<mesh_t> macro(Kh_i, 0.25);

        double xi  = 3. / 4;
        double xi0 = (xi - 0.5) / 2.;

        Lagrange2 FEvelocity(4);
        space_t V2h(Kh, FEvelocity);
        cutspace_t W2h(Kh_i, V2h);
        space_t Q2h(Kh, DataFE<mesh_t>::P2dc);
        cutspace_t P2h(Kh_i, Q2h);

        fct_t fv(W2h, fun_force);
        fct_t fq(P2h, fun_div);
        fct_t u0(W2h, fun_exact_u);
        fct_t p0(P2h, fun_exact_p);
        fct_t phat(Ph, fun_interfacePr);
        fct_t pex(Ph, fun_exact_p);

        Normal n;
        Tangent t;
        funtest_t p(Ph, 1), q(Ph, 1), u(Wh, 2), v(Wh, 2);

        double uPenParam = 1e-1;
        double pPenParam = 1e-1;
        double jumpParam = 1e0;

        darcy.addBilinear(innerProduct(u, v) - innerProduct(p, div(v)) + innerProduct(div(u), q), Kh_i);

        darcy.addBilinear(innerProduct(mu_G * average(u * n), average(v * n)) +
                              innerProduct(xi0 * mu_G * jump(u * n), jump(v * n)),
                          interface);

        darcy.addLinear(-innerProduct(phat.expr(), jump(v * n)), interface);

        darcy.addLinear(innerProduct(fv.exprList(2), v) + innerProduct(fq.expr(), q), Kh_i);

        darcy.addLinear(-innerProduct(p0.expr(), v * n), Kh_i, INTEGRAL_BOUNDARY);

        funtest_t grad2un = grad(grad(u) * n) * n;
        darcy.addFaceStabilization(innerProduct(uPenParam * h_i * jump(u), jump(v)) +
                                       innerProduct(uPenParam * pow(h_i, 3) * jump(grad(u) * n), jump(grad(v) * n)) +
                                       innerProduct(uPenParam * pow(h_i, 5) * jump(grad2un), jump(grad2un)) -
                                       innerProduct(pPenParam * h_i * jump(p), jump(div(v))) +
                                       innerProduct(pPenParam * h_i * jump(div(u)), jump(q)) -
                                       innerProduct(pPenParam * pow(h_i, 3) * jump(grad(p)), jump(grad(div(v)))) +
                                       innerProduct(pPenParam * pow(h_i, 3) * jump(grad(div(v))), jump(grad(q))),
                                   Kh_i, macro);

        darcy.solve("umfpack");

        // EXTRACT SOLUTION
        int idx0_s = Wh.get_nb_dof();
        std::span<double> data_uh{std::span(darcy.rhs_.data(), Wh.get_nb_dof())};
        std::span<double> data_ph{std::span(darcy.rhs_.data() + Wh.get_nb_dof(), Ph.get_nb_dof())};
        fct_t uh(Wh, data_uh);
        fct_t ph(Ph, data_ph);
        auto femSol_0dx = dx(uh.expr(0));
        auto femSol_1dy = dy(uh.expr(1));

        // L2 norm vel
        double errU      = L2normCut(uh, fun_exact_u, 0, 2);
        double errP      = L2normCut(ph, fun_exact_p, 0, 1);
        double errDiv    = L2normCut(femSol_0dx + femSol_1dy, fun_div, Kh_i);
        double maxErrDiv = maxNormCut(femSol_0dx + femSol_1dy, fun_div, Kh_i);
        // [PLOTTING]
        // if (MPIcf::IamMaster())
        // {
        //    fct_t solh(Wh, fun_exact_u);
        //    // solh.v -= uh.v;
        //    // solh.v.map(fabs);
        //    fct_t divSolh(Ph, fun_div);
        //    ExpressionFunFEM<mesh_t> femDiv(divSolh, 0, op_id);

        //    // Paraview<mesh_t> writer(Kh_i, "darcyRT2_"+to_string(i)+".vtk");
        //    Paraview<mesh_t> writer(Kh_i, "darcyscotti" + to_string(i) + ".vtk");

        //    writer.add(uh, "velocity", 0, 2);
        //    writer.add(ph, "pressure", 0, 1);
        //    writer.add(femSol_0dx + femSol_1dy, "divergence");
        //    writer.add(solh, "velocityExact", 0, 2);
        //    writer.add(fabs((femSol_0dx + femSol_1dy) - femDiv),
        //               "divergenceError");
        //    writer.add(fq, "divu", 0, 1);
        //    writer.add(fqq, "divuex", 0, 1);
        // }

        pPrint.push_back(errP);
        uPrint.push_back(errU);
        divPrint.push_back(errDiv);
        maxDivPrint.push_back(maxErrDiv);
        h.push_back(h_i);

        if (i == 0) {
            convpPr.push_back(0);
            convuPr.push_back(0);
            convdivPr.push_back(0);
            convdivPrLoc.push_back(0);
            convmaxdivPr.push_back(0);
        } else {
            convpPr.push_back(log(pPrint[i] / pPrint[i - 1]) / log(h[i] / h[i - 1]));
            convuPr.push_back(log(uPrint[i] / uPrint[i - 1]) / log(h[i] / h[i - 1]));
            convdivPr.push_back(log(divPrint[i] / divPrint[i - 1]) / log(h[i] / h[i - 1]));
            convmaxdivPr.push_back(log(maxDivPrint[i] / maxDivPrint[i - 1]) / log(h[i] / h[i - 1]));
        }
        nx = 2 * nx - 1;
    }

    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ') << "err_p"
              << std::setw(15) << std::setfill(' ') << "conv p" << std::setw(15) << std::setfill(' ') << "err u"
              << std::setw(15) << std::setfill(' ') << "conv u" << std::setw(15) << std::setfill(' ') << "err divu"
              << std::setw(15) << std::setfill(' ')
              << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ') << "err maxdivu" << std::setw(15) << std::setfill(' ')
              << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < uPrint.size(); ++i) {
        std::cout << std::left << std::setprecision(5) << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15)
                  << std::setfill(' ') << pPrint[i] << std::setw(15) << std::setfill(' ') << convpPr[i] << std::setw(15)
                  << std::setfill(' ') << uPrint[i] << std::setw(15) << std::setfill(' ') << convuPr[i] << std::setw(15)
                  << std::setfill(' ') << divPrint[i] << std::setw(15) << std::setfill(' ')
                  << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << maxDivPrint[i]
                  // << std::setw(15) << std::setfill(' ') << ratioCut1[i]
                  // << std::setw(15) << std::setfill(' ') << ratioCut2[i]
                  << std::setw(15) << std::setfill(' ') << convmaxdivPr[i] << std::endl;
    }
}
