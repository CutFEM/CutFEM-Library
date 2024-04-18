#include "../cutfem.hpp"

using mesh_t     = Mesh3;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;



double shift        = 0.5;
double interfaceRad = 0.25001;               
double mu_G         = 2. / 3 * interfaceRad; 

R fun_radius2(const R3 P) {
    return (P[0] - shift) * (P[0] - shift) + (P[1] - shift) * (P[1] - shift) + (P.z - shift) * (P.z - shift);
}
R fun_levelSet(const R3 P, const int i) {
    return sqrt((P[0] - shift) * (P[0] - shift) + (P[1] - shift) * (P[1] - shift) + (P.z - shift) * (P.z - shift)) -
           interfaceRad;
}

R fun_dirichlet(const R3 P, int compInd) { return 0; }
R fun_neumann(const R3 P, int compInd, int dom) {
    R r2      = fun_radius2(P);
    R radius2 = interfaceRad * interfaceRad;
    return r2 / (2 * radius2) + 3. / 2.;
}
R fun_force(const R3 P, int compInd) { return 0; }

R fun_div(const R3 P, int compInd, int dom) { // is also exact divergence
    R radius2 = interfaceRad * interfaceRad;
    if (dom == 0) // r2>radius2
        return -3. / radius2;
    else
        return -6. / radius2;
}
R fun_exact_u(const R3 P, int compInd, int dom) {
    DA<R, 3> X(P[0], 0), Y(P[1], 1), Z(P.z, 2);
    DA<R, 3> r2  = (X - shift) * (X - shift) + (Y - shift) * (Y - shift) + (Z - shift) * (Z - shift);
    R radius2      = interfaceRad * interfaceRad;
    R cst          = (dom == 0) * 3. / 2;
    R mul          = (dom == 0) * 2 + (dom == 1) * 1;
    DA<R, 3> val = r2 / (mul * radius2) + cst;
    return -val.d[compInd];

    // return P.norme();
}
R fun_exact_p(const R3 P, int compInd, int dom) {
    DA<R, 3> X(P[0], 0), Y(P[1], 1), Z(P.z, 2);
    DA<R, 3> r2  = (X - shift) * (X - shift) + (Y - shift) * (Y - shift) + (Z - shift) * (Z - shift);
    R radius2      = interfaceRad * interfaceRad;
    R cst          = (dom == 0) * 3. / 2;
    R mul          = (dom == 0) * 2 + (dom == 1) * 1;
    DA<R, 3> val = r2 / (mul * radius2) + cst;
    return val.val;
}
R fun_interfacePr(const R3 P, int compInd) { return 19. / 12; }

int main(int argc, char **argv) {


/// MUMPS allow to solve for lower mesh sizes but the divergence error inrcreses 
/// while umfpack gives good divergence error but fails for lower mesh sizes
    #ifndef USE_MPI
    assert(0 && " need mpi for mumps");
    #else
    MPIcf cfMPI(argc, argv);
    #endif
    int nx = 11;

    std::vector<double> uPrint, pPrint, divPrint, divPrintLoc, maxDivPrint, h, convuPr, convpPr, convdivPr,
        convdivPrLoc, convmaxdivPr;
    std::vector<double> ratioCut1, ratioCut2;
    int iters = 3;

    for (int i = 0; i < iters; ++i) {
        mesh_t Kh(nx, nx, nx, 0., 0., 0., 1., 1., 1.);

        space_t Lh(Kh, DataFE<mesh_t>::P1);
        fct_t levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<mesh_t> interface(Kh, levelSet);

        Lagrange3 FEvelocity(2);
        space_t Vh(Kh, DataFE<mesh_t>::RT0);
        space_t Qh(Kh, DataFE<mesh_t>::P0);

        ActiveMesh<mesh_t> Kh_i(Kh, interface);
        cutspace_t Wh(Kh_i, Vh);
        cutspace_t Ph(Kh_i, Qh);

        CutFEM<mesh_t> darcy(Wh);
        darcy.add(Ph);
        const R h_i  = 1. / (nx - 1);
        const R invh = 1. / h_i;
        R xi         = 3. / 4;
        R xi0        = (xi - 0.5) / 2.;

        // We define fh on the cutSpace
        fct_t fq(Ph, fun_div);
        fct_t p0(Lh, fun_neumann);
        fct_t phat(Ph, fun_interfacePr);


        Normal n;
        Tangent t;
        funtest_t p(Ph, 1), q(Ph, 1), u(Wh, 3), v(Wh, 3);

        double uPenParam = 1e0; 
        double pPenParam = 1e0; 
        double jumpParam = 1e0;

        darcy.addBilinear(innerProduct(u, v) - innerProduct(p, div(v)) + innerProduct(div(u), q), Kh_i);
        darcy.addBilinear(innerProduct(mu_G * average(u * n), average(v * n)) +
                              innerProduct(xi0 * mu_G * jump(u * n), jump(v * n)) 
                          ,interface);

        darcy.addLinear(-innerProduct(phat.expr(), jump(v * n)), interface);
        darcy.addLinear(innerProduct(fq.expr(), q), Kh_i);

        darcy.addLinear(
            -innerProduct(p0.expr(), v * n)
            ,
            Kh_i, INTEGRAL_BOUNDARY);

        funtest_t grad2un = grad(grad(u) * n) * n;
        darcy.addFaceStabilization( 
            innerProduct(uPenParam * pow(h_i, 1) * jump(u), jump(v)) 
                + innerProduct(uPenParam * pow(h_i, 3) * jump(grad(u) * n), jump(grad(v) * n))
                - innerProduct(pPenParam * h_i * jump(p), jump(div(v))) +
                innerProduct(pPenParam * h_i * jump(div(u)), jump(q))
            ,
            Kh_i
        );

        darcy.solve("mumps");

        // EXTRACT SOLUTION
        int idx0_s  = Wh.get_nb_dof();
        std::span<double> data_uh{std::span(darcy.rhs_.data(), Wh.get_nb_dof())};
        std::span<double> data_ph{std::span(darcy.rhs_.data() + Wh.get_nb_dof(), Ph.get_nb_dof())};
        fct_t uh(Wh, data_uh);
        fct_t ph(Ph, data_ph);

        auto femSol_0dx = dx(uh.expr(0));
        auto femSol_1dy = dy(uh.expr(1));
        auto femSol_1dz = dz(uh.expr(2));

        // L2 norm vel
        R errU      = L2normCut(uh, fun_exact_u, 0, 3);
        R errP      = L2normCut(ph, fun_exact_p, 0, 1);
        R errDiv    = L2normCut(femSol_0dx + femSol_1dy + femSol_1dz, fun_div, Kh_i);
        R maxErrDiv = maxNormCut(femSol_0dx + femSol_1dy + femSol_1dz, fun_div, Kh_i);

        // [PLOTTING]
        #ifdef USE_MPI
        if(MPIcf::IamMaster())
        #endif
        {
          fct_t solh(Wh, fun_exact_u);
          fct_t solph(Ph, fun_exact_p);
          fct_t divSolh(Ph, fun_div);
          auto femDiv = divSolh.expr();
          Paraview<mesh_t> writer(Kh_i, "darcy_example1_3D_"+std::to_string(i)+".vtk");
        
          writer.add(uh, "velocity" , 0, 3);
          writer.add(solh, "velocityExact" , 0, 3);
          writer.add(ph, "pressure" , 0, 1);
          writer.add(solph, "pressureExact" , 0, 1);
          writer.add(femSol_0dx+femSol_1dy, "divergence");
          writer.add(fabs((femSol_0dx+femSol_1dy+femSol_1dz)-femDiv),
          "divergenceError");
        }

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
              << std::left << std::setw(10) << std::setfill(' ') << "h" << std::setw(15) << std::setfill(' ')
              << "err_p"
              // << std::setw(15) << std::setfill(' ') << "conv p"
              << std::setw(15) << std::setfill(' ')
              << "err u"
              // << std::setw(15) << std::setfill(' ') << "conv u"
              << std::setw(15) << std::setfill(' ')
              << "err divu"
              // << std::setw(15) << std::setfill(' ') << "conv divu"
              // << std::setw(15) << std::setfill(' ') << "err_new divu"
              // << std::setw(15) << std::setfill(' ') << "convLoc divu"
              << std::setw(15) << std::setfill(' ')
              << "err maxdivu"
              // << std::setw(15) << std::setfill(' ') << "conv maxdivu"
              << "\n"
              << std::endl;
    for (int i = 0; i < uPrint.size(); ++i) {
        std::cout << std::left << std::setprecision(5) << std::setw(10) << std::setfill(' ') << h[i] << std::setw(15)
                  << std::setfill(' ')
                  << pPrint[i]
                  // << std::setw(15) << std::setfill(' ') << convpPr[i]
                  << std::setw(15) << std::setfill(' ')
                  << uPrint[i]
                  // << std::setw(15) << std::setfill(' ') << convuPr[i]
                  << std::setw(15) << std::setfill(' ')
                  << divPrint[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPr[i]
                  // << std::setw(15) << std::setfill(' ') << divPrintLoc[i]
                  // << std::setw(15) << std::setfill(' ') << convdivPrLoc[i]
                  << std::setw(15) << std::setfill(' ')
                  << maxDivPrint[i]
                  // << std::setw(15) << std::setfill(' ') << ratioCut1[i]
                  // << std::setw(15) << std::setfill(' ') << ratioCut2[i]
                  // << std::setw(15) << std::setfill(' ') << convmaxdivPr[i]
                  << std::endl;
    }
}
