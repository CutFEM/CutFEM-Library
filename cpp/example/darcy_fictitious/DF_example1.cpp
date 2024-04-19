#include "cpp/cutfem.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;

// double sq_SW        = 0. + 1e-10;
// double sq_LGTH      = 1. + 2e-10;
double interfaceRad = 0.250001;
double mu_G         = 2. / 3 * interfaceRad;
R2 shift(0.5, 0.5);

double fun_levelSet(R2 P, const int i) { return -(norm2(P - shift) - interfaceRad); }
double fun_div(R2 P, int compInd, int dom) { return 2.; }
double fun_exact_u(R2 P, int compInd, int dom) {
    if (compInd == 0) {
        return P.x - 0.5;
    } else {
        return P.y - 0.5;
    }
}
double fun_exact_p(R2 P, int compInd, int dom) { return -0.5 * (P.x * P.x - P.x + P.y * P.y - P.y); }

int main(int argc, char **argv) {

    globalVariable::verbose = 0;

    MPIcf cfMPI(argc, argv);

    int nx = 11;

    std::vector<double> uPrint, pPrint, divPrint, divPrintLoc, maxDivPrint, h, convuPr, convpPr, convdivPr,
        convdivPrLoc, convmaxdivPr;
    std::vector<double> ratioCut1, ratioCut2;
    int iters = 5;

    for (int i = 0; i < iters; ++i) {

        mesh_t Kh(nx, nx, 0., 0., 1., 1.);

        space_t Lh(Kh, DataFE<mesh_t>::P1);
        fct_t levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<mesh_t> interface(Kh, levelSet);

        space_t V2h(Kh, DataFE<mesh_t>::RT2);
        space_t Q2h(Kh, DataFE<mesh_t>::P2);
        space_t Vh(Kh, DataFE<mesh_t>::RT0);
        space_t Qh(Kh, DataFE<mesh_t>::P0);

        cutmesh_t Kh_i(Kh);
        Kh_i.truncate(interface, -1);

        cutspace_t W2h(Kh_i, V2h);
        cutspace_t P2h(Kh_i, Q2h);
        cutspace_t Wh(Kh_i, Vh);
        cutspace_t Ph(Kh_i, Qh);

        // SURFACE MESH
        cutmesh_t Kh_itf(Kh);
        Kh_itf.createSurfaceMesh(interface);
        cutspace_t Ph_itf(Kh_itf, Qh);

        CutFEM<mesh_t> darcy(Wh);
        darcy.add(Ph);
        darcy.add(Ph_itf);

        const double h_i  = 2. / (nx - 1);
        const double invh = 1. / h_i;

        MacroElement<mesh_t> macro(Kh_i, 1.);

        fct_t fq(Ph, fun_div);
        fct_t u0(W2h, fun_exact_u);
        fct_t p0(P2h, fun_exact_p);
        auto exactp = p0.expr(0);

        Normal n;
        Tangent t;
        funtest_t p(Ph, 1), q(Ph, 1), u(Wh, 2), v(Wh, 2);
        funtest_t p_itf(Ph_itf, 1), q_itf(Ph_itf, 1);

        double uPenParam   = 1e1;
        double pPenParam   = 1e1;
        double itfPenParam = 1e-1;

        darcy.addBilinear(innerProduct(u, v) - innerProduct(p, div(v)) + innerProduct(div(u), q), Kh_i);

        darcy.addLinear(innerProduct(fq.expr(), q), Kh_i);

        darcy.addBilinear(innerProduct(p_itf, v * n) + innerProduct(u * n, q_itf), interface);

        darcy.addLinear(innerProduct(u0.exprList(2), q_itf * n), interface);

        // [GHOST PENALTY]
        funtest_t grad2un = grad(grad(u) * n) * n;
        darcy.addFaceStabilization(
            // [h^(2k+1) h^(2k+1)]
            //  innerProduct(uPenParam*pow(h_i,1)*jump(u), jump(v))
            // +innerProduct(uPenParam*pow(h_i,3)*jump(grad(u)*n), jump(grad(v)*n))
            // +innerProduct(uPenParam*pow(h_i,5)*jump(grad2un), jump(grad2un))
            // +innerProduct(pPenParam*pow(h_i,1)*jump(p), jump(q))
            // +innerProduct(pPenParam*pow(h_i,3)*jump(grad(p)), jump(grad(q)))

            +innerProduct(uPenParam * pow(h_i, 1) * jump(u), jump(v)) +
                innerProduct(uPenParam * pow(h_i, 3) * jump(grad(u) * n), jump(grad(v) * n))
                // +innerProduct(uPenParam*pow(h_i,5)*jump(grad2un), jump(grad2un))
                - innerProduct(pPenParam * pow(h_i, 1) * jump(p), jump(div(v))) +
                innerProduct(pPenParam * pow(h_i, 1) * jump(div(u)), jump(q)) -
                innerProduct(pPenParam * pow(h_i, 3) * jump(grad(p)), jump(grad(div(v)))) +
                innerProduct(pPenParam * pow(h_i, 3) * jump(grad(div(v))), jump(grad(q))),
            Kh_i, macro);

        darcy.addFaceStabilizationMixed(
            // -innerProduct(itfPenParam * pow(h_i, 1) * jump(p_itf), jump(q_itf))
            //-
            //     innerProduct(itfPenParam * pow(h_i, 3) * jump(grad(p_itf)), jump(grad(q_itf))) // (1)
            -innerProduct(itfPenParam * pow(h_i, 0) * jump(p_itf), jump(div(v))) -
                innerProduct(itfPenParam * pow(h_i, 0) * jump(div(u)), jump(q_itf)),
            Kh_itf
            // , macro_itf
        );

        darcy.addLagrangeMultiplier(innerProduct(1, q), 0, Kh_i);

        // darcy.addLinear(-innerProduct(phat.expr(), jump(v * n)), interface);

        // darcy.addLinear(-innerProduct(p0.expr(), v * n), Kh_i, INTEGRAL_BOUNDARY);

        // funtest_t grad2un = grad(grad(u) * n) * n;
        // darcy.addFaceStabilization(innerProduct(uPenParam * h_i * jump(u), jump(v)) +
        //                                innerProduct(uPenParam * pow(h_i, 3) * jump(grad(u) * n), jump(grad(v) * n)) +
        //                                innerProduct(uPenParam * pow(h_i, 5) * jump(grad2un), jump(grad2un)) -
        //                                innerProduct(pPenParam * h_i * jump(p), jump(div(v))) +
        //                                innerProduct(pPenParam * h_i * jump(div(u)), jump(q)) -
        //                                innerProduct(pPenParam * pow(h_i, 3) * jump(grad(p)), jump(grad(div(v)))) +
        //                                innerProduct(pPenParam * pow(h_i, 3) * jump(grad(div(v))), jump(grad(q))),
        //                            Kh_i, macro);

        matlab::Export(darcy.mat_[0], "mat" + std::to_string(i) + "Cut.dat");

        darcy.solve("umfpack");

        double area = integral(Kh_i, 1.);

        // EXTRACT SOLUTION
        int idx0_s = Wh.get_nb_dof();
        std::span<double> data_uh{std::span(darcy.rhs_.data(), Wh.get_nb_dof())};
        std::span<double> data_ph{std::span(darcy.rhs_.data() + Wh.get_nb_dof(), Ph.get_nb_dof())};
        fct_t uh(Wh, data_uh);
        fct_t ph(Ph, data_ph);

        // [Post process pressure]
        auto fem_p      = ph.expr(0);
        double meanP    = integral(Kh_i, exactp, 0);
        double meanPfem = integral(Kh_i, fem_p, 0);
        ph.v += (meanP - meanPfem) / area;

        // L2 norm vel
        auto femSol_0dx  = dx(uh.expr(0));
        auto femSol_1dy  = dy(uh.expr(1));
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
              << "conv maxdivu" << "\n"
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
