#include "../tool.hpp"

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

    int nx              = 11;
    double penaltyParam = 4e3;
    double wPenParam1   = 1e1;
    double wPenParam2   = 0.;
    double uPenParam    = 1e1;
    double pPenParam    = 1e1;
    double sigma        = 1e-2;
    double mu           = 1.;

    Normal n;
    Tangent t;

    std::vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    globalVariable::verbose = 2;

    for (int i = 0; i < 4; ++i) {

        double hi = 1. / (nx - 1);

        mesh_t Kh(nx, nx, -0.127, -0.127, 1.333, 1.333);
        space_t Lh(Kh, DataFE<Mesh2>::P1);
        fct_t levelSet(Lh, fun_levelSet);

        InterfaceLevelSet<mesh_t> interface(Kh, levelSet);

        cutmesh_t Khi(Kh);
        Khi.truncate(interface, 1);

        MacroElement<mesh_t> macro(Khi, 0.8);

        Lagrange2 FEvelocity(4);
        space_t VELh(Kh, FEvelocity);
        space_t SCAh(Kh, DataFE<mesh_t>::P4);

        space_t Wh(Kh, DataFE<mesh_t>::RT1);
        cutspace_t Vh(Khi, Wh);

        space_t Qh(Kh, DataFE<mesh_t>::P1dc);
        cutspace_t Ph(Khi, Qh);

        space_t Uh_(Kh, DataFE<mesh_t>::P2);
        cutspace_t Uh(Khi, Uh_);

        fct_t fh(VELh, fun_rhs);
        fct_t u0(VELh, fun_exact_u);
        fct_t p0(SCAh, fun_exact_p);

        CutFEM<mesh_t> stokes(Uh);
        stokes.add(Vh);
        stokes.add(Ph);

        funtest_t w(Uh, 1, 0), tau(Uh, 1, 0), u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);

        // [Bulk]
        stokes.addBilinear( // w = curl u
            innerProduct(1. / mu * w, tau) - innerProduct(u, rotgrad(tau)), Khi);
        stokes.addBilinear( // mu Delta u + grad p
            innerProduct(rotgrad(w), v) - innerProduct(p, div(v)), Khi);
        stokes.addLinear(+innerProduct(fh.exprList(2), v), Khi);
        stokes.addBilinear(+innerProduct(div(u), q), Khi);

        //[interface]
        stokes.addBilinear(+innerProduct(p, v * n) + innerProduct(1. / hi * penaltyParam * u * n, v * n) // stability
                           ,
                           interface);
        stokes.addLinear(+innerProduct(u0 * t, tau) + innerProduct(u0 * n, 1. / hi * penaltyParam * v * n), interface);

        funtest_t grad2un = grad(grad(u) * n) * n;
        funtest_t grad2wn = grad(grad(w) * n) * n;
        stokes.addFaceStabilization(
            /* "Primal" stab (lw,0,la) */
            +innerProduct(wPenParam1 * pow(hi, 3) * jump(grad(w) * n), jump(grad(tau) * n)) +
                innerProduct(wPenParam1 * pow(hi, 5) * jump(grad2wn), jump(grad2wn)) +

                +innerProduct(wPenParam2 * pow(hi, 1) * jump(rotgrad(w)), jump(v)) -
                innerProduct(wPenParam2 * pow(hi, 1) * jump(u), jump(rotgrad(tau))) +
                innerProduct(wPenParam2 * pow(hi, 3) * jump(grad(rotgrad(w))), jump(grad(v))) -
                innerProduct(wPenParam2 * pow(hi, 3) * jump(grad(u)), jump(grad(rotgrad(tau))))

                + innerProduct(uPenParam * pow(hi, 1) * jump(u), jump(v)) +
                innerProduct(uPenParam * pow(hi, 3) * jump(grad(u) * n), jump(grad(v) * n)) +
                innerProduct(uPenParam * pow(hi, 5) * jump(grad2un), jump(grad2un))

                - innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
                innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
                innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(u))), jump(grad(q))),
            Khi, macro);

        // [Sets uniqueness of the pressure]

        CutFEM<mesh_t> lagr(Uh);
        lagr.add(Vh);
        lagr.add(Ph);
        lagr.addLinear(innerProduct(1, p), Khi);
        std::vector<double> lag_row(lagr.rhs_.begin(), lagr.rhs_.end());
        lagr.rhs_ = 0.;
        lagr.addLinear(innerProduct(1, v * n), interface);
        stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhsSpan(), 0);

        stokes.solve("umfpack");

        // EXTRACT SOLUTION
        int nb_vort_dof           = Uh.get_nb_dof();
        int nb_flux_dof           = Vh.get_nb_dof();
        std::span<double> data_wh = std::span<double>(stokes.rhs_.data(), nb_vort_dof);
        std::span<double> data_uh = std::span<double>(stokes.rhs_.data() + nb_vort_dof, nb_flux_dof);
        std::span<double> data_ph = std::span<double>(stokes.rhs_.data() + nb_flux_dof + nb_vort_dof, Ph.get_nb_dof());

        fct_t wh(Uh, data_wh);
        fct_t uh(Vh, data_uh);
        fct_t ph(Ph, data_ph);

        // [Post process pressure]
        double meanP    = integral(Khi, p0.expr(), 0);
        double meanPfem = integral(Khi, ph.expr(), 0);
        CutFEM<mesh_t> post(Ph);
        post.addLinear(innerProduct(1, q), Khi);
        double area = post.rhs_.sum();
        ph.v -= meanPfem / area;
        ph.v += meanP / area;

        {
            Paraview<mesh_t> writer(Khi, "stokesVorticity_" + std::to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(dx(uh.expr(0)) + dy(uh.expr(1)), "divergence");
        }
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
