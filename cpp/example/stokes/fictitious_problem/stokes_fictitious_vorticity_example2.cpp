#ifdef PROBLEM_UNFITTED_CORIOLIS_STOKES_VORTICITY

namespace Erik_Data_CORIOLIS_STOKES {

R sdBox(const R2 p, const R2 b) {
    R2 absp(abs(p.x), abs(p.y));
    R2 d = absp - b;
    R2 maxd0(max(d.x, 0.0), max(d.y, 0.0));
    return maxd0.norm() + min(max(d.x, d.y), 0.0);
}
R fun_levelSet(const R2 P, const int i) {
    const R2 b(2.0, 1.0); // b.x=width, b.y=height
    const R2 shift(4., 0.);
    return sdBox(P - shift, b) - 0.01;
}

// [Coriolis example]
R fun_y(const R2 P, int i, int dom) { return P.y; }
R fun_div(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    return 0;
}
R fun_rhs(const R2 P, int i, int dom) {
    // R mu=1;
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 0; // 100*(2-x)*(2-x);
    else if (i == 1)
        return 0;
}
R fun_uin(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return y * (2 - y) / 2; // original paper 4*(2-y)*(y-1)
    else
        return 0;
}
R fun_uout(const R2 P, int i, int dom) {
    R x = P.x;
    R y = P.y;
    if (i == 0)
        return 1. / 4 * (2 - y) * (y - 1); // original paper 4*(2-y)*(y-1)
    else
        return 0;
}
} // namespace Erik_Data_CORIOLIS_STOKES

using namespace Erik_Data_CORIOLIS_STOKES;

int main(int argc, char **argv) {
    typedef TestFunction<2> FunTest;
    typedef FunFEM<Mesh2> Fun_h;
    typedef Mesh2 Mesh;
    typedef ActiveMeshT2 CutMesh;
    typedef FESpace2 Space;
    typedef CutFESpaceT2 CutSpace;

    const double cpubegin = CPUtime();
    MPIcf cfMPI(argc, argv);

    int nx = 11;
    int ny = 11;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 4;
    for (int i = 0; i < iters; ++i) { // i<3

        std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, 0., 0., 4., 2.);
        const R hi           = 1. / (nx - 1);
        const R penaltyParam = 400; // 400
        double uPenParam     = 1e0; // 1e-1, 1e1
        double pPenParam     = 1e0; // 1e-1, 1e1

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(4);
        Space VELh_(Kh, FEvelocity);
        Space SCAh_(Kh, DataFE<Mesh>::P2);

        Space Uh_(Kh, DataFE<Mesh>::P2);
        Space Vh_(Kh, DataFE<Mesh2>::RT1);
        Space Wh_(Kh, DataFE<Mesh2>::P1dc);

        ActiveMesh<Mesh> Khi(Kh);
        Khi.truncate(interface, -1);

        MacroElement<Mesh> macro(Khi, 1.);

        CutSpace VELh(Khi, VELh_);
        CutSpace SCAh(Khi, SCAh_);

        CutSpace Uh(Khi, Uh_);
        CutSpace Vh(Khi, Vh_);
        CutSpace Wh(Khi, Wh_);

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h uhin(VELh, fun_uin);
        Fun_h uhout(VELh, fun_uout);
        Fun_h yfunh(VELh, fun_y);

        CutFEM<Mesh2> stokes(Uh);
        stokes.add(Vh);
        stokes.add(Wh);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest w(Uh, 1, 0), tau(Uh, 1, 0), u(Vh, 2, 0), v(Vh, 2, 0), p(Wh, 1, 0), q(Wh, 1, 0);
        FunTest u1(Vh, 1, 0), u2(Vh, 1, 1), v1(Vh, 1, 0), v2(Vh, 1, 1);

        R mu    = 0.01;
        R omega = 100; // 50 // original paper 100;
        {
            // [Bulk]
            stokes.addBilinear( // coriolis
                -innerProduct(2 * omega * u2, v1) + innerProduct(2 * omega * u1, v2), Khi);
            stokes.addBilinear( // w = curl u
                innerProduct(1. / mu * w, tau) - innerProduct(u, rotgrad(tau)), Khi);
            stokes.addBilinear( // mu Delta u + grad p
                innerProduct(rotgrad(w), v) - innerProduct(p, div(v)), Khi);
            stokes.addLinear(innerProduct(fh.expression(2), v), Khi);
            stokes.addBilinear(+innerProduct(div(u), q), Khi);
            stokes.addBilinear(+innerProduct(p, v * n) +
                                   innerProduct(1. / hi * penaltyParam * u * n, v * n) // stability
                               ,
                               interface);
            stokes.addBilinear(+innerProduct(p, v * n) +
                                   innerProduct(1. / hi * penaltyParam * u * n, v * n) // stability
                               ,
                               Khi, INTEGRAL_BOUNDARY);
            stokes.addLinear(+innerProduct(uhout * t, tau) + innerProduct(uhout * n, 1. / hi * penaltyParam * v * n),
                             Khi, INTEGRAL_BOUNDARY, {2});
            stokes.addLinear(+innerProduct(uhin * t, tau) + innerProduct(uhin * n, 1. / hi * penaltyParam * v * n), Khi,
                             INTEGRAL_BOUNDARY, {4});
            // [Sets uniqueness of the pressure]
            // stokes.addLagrangeMultiplier(
            //   innerProduct(1.,p), 0
            //   , Khi
            // );
            // [Sets uniqueness of the pressure in another way such that divu = 0]
            CutFEM<Mesh2> lagr(Uh);
            lagr.add(Vh);
            lagr.add(Wh);
            Rn zero_vec = lagr.rhs_;
            lagr.addLinear(innerProduct(1., p), Khi);
            Rn lag_row(lagr.rhs_);
            lagr.rhs_ = zero_vec;
            lagr.addLinear(innerProduct(1, v * n), interface);
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , Khi, INTEGRAL_BOUNDARY
            // );
            stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);
            // [Stabilization]
            FunTest grad2un = grad(grad(u) * n) * n;
            FunTest grad2wn = grad(grad(w) * n) * n;
            // // FunTest grad2pn = grad(grad(p)*n)*n;
            // // FunTest grad2divun = grad(grad(div(u))*n)*n;
            stokes.addFaceStabilization(
                +innerProduct(uPenParam * pow(hi, 3) * jump(grad(w) * n), jump(grad(tau) * n)) +
                    innerProduct(uPenParam * pow(hi, 5) * jump(grad2wn), jump(grad2wn)) +
                    innerProduct(uPenParam * pow(hi, 1) * jump(u), jump(v)) +
                    innerProduct(uPenParam * pow(hi, 3) * jump(grad(u) * n), jump(grad(v) * n)) +
                    innerProduct(uPenParam * pow(hi, 5) * jump(grad2un), jump(grad2un))

                    // +innerProduct(pPenParam*pow(hi,-1)*jump(rotgrad(w)), jump(v))
                    // -innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(rotgrad(tau)))
                    // +innerProduct(uPenParam*pow(hi,1)*jump(grad(rotgrad(w))*n), jump(grad(v)*n))
                    // -innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(rotgrad(tau))*n))

                    - innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
                    innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(u))), jump(grad(q))) +
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v))))

                // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
                // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))
                ,
                Khi, macro);
        }

        // std::cout << integral(Khi,exactp,0) << std::endl;
        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_vort_dof = Uh.get_nb_dof();
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_wh     = stokes.rhs_(SubArray(nb_vort_dof, 0));
        Rn_ data_uh     = stokes.rhs_(SubArray(
            nb_flux_dof, nb_vort_dof)); // Rn_ data_uh = stokes.rhs_(SubArray(nb_vort_dof+nb_flux_dof,nb_vort_dof));
        Rn_ data_ph     = stokes.rhs_(SubArray(
            Wh.get_nb_dof(),
            nb_vort_dof +
                nb_flux_dof)); // Rn_ data_ph = stokes.rhs_(SubArray(stokes_.get_nb_dof(),nb_vort_dof+nb_flux_dof));
        Fun_h wh(Uh, data_wh);
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Wh, data_ph);

        // [Post process pressure]
        // R meanP = integral(Khi,exactp,0);
        // ExpressionFunFEM<Mesh> fem_p(ph,0,op_id);
        // R meanPfem = integral(Khi,fem_p,0);
        // // std::cout << meanP << std::endl;
        // CutFEM<Mesh2> post(Qh);
        // post.addLinear(
        //   innerProduct(1,q)
        //   , Khi
        // );
        // R area = post.rhs_.sum();
        // ph.v -= meanPfem/area;
        // ph.v += meanP/area;

        ExpressionFunFEM<Mesh> uh_0dx(uh, 0, op_dx);
        ExpressionFunFEM<Mesh> uh_1dy(uh, 1, op_dy);

        // Fun_h test(Vh, rt0_vec);
        // ExpressionFunFEM<Mesh> femSol_0dx(test, 0, op_dx);
        // ExpressionFunFEM<Mesh> femSol_1dy(test, 1, op_dy);

        // [Errors]
        {
            // Fun_h soluErr(Vh, fun_exact_u);
            // Fun_h soluh(Vh, fun_exact_u);
            // soluErr.v -= uh.v;
            // soluErr.v.map(fabs);
            // Fun_h divSolh(Wh, fun_div);
            // ExpressionFunFEM<Mesh> femDiv(divSolh, 0, op_id);

            Paraview<Mesh> writer(Khi, "stokes_" + to_string(i) + ".vtk");
            writer.add(uh, "velocity", 0, 2);
            writer.add(ph, "pressure", 0, 1);
            writer.add(uh_0dx + uh_1dy, "divergence");
            // writer.add(femSol_0dx+femSol_1dy+fflambdah, "divergence");
            // writer.add(soluh, "velocityExact" , 0, 2);
            // writer.add(soluErr, "velocityError" , 0, 2);
            // writer.add(solh, "velocityError" , 0, 2);

            // writer.add(fabs(femDiv, "divergenceError");
        }

        // R errU      = L2normCut(uh,fun_exact_u,0,2);
        // R errP      = L2normCut(ph,fun_exact_p,0,1);
        R errDiv    = L2normCut(uh_0dx + uh_1dy, Khi);
        R maxErrDiv = maxNormCut(uh_0dx + uh_1dy, Khi);
        // R errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // R maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

        // ul2.push_back(errU);
        // pl2.push_back(errP);
        divl2.push_back(errDiv);
        divmax.push_back(maxErrDiv);
        h.push_back(hi);
        // if(i==0) {convu.push_back(0); convp.push_back(0);}
        // else {
        //   convu.push_back( log(ul2[i]/ul2[i-1])/log(h[i]/h[i-1]));
        //   convp.push_back(log(pl2[i]/pl2[i-1])/log(h[i]/h[i-1]));
        // }

        nx = 2 * nx - 1;
        ny = 2 * ny - 1;
    }
    std::cout << "\n"
              << std::left << std::setw(10) << std::setfill(' ')
              << "h"
              // << std::setw(15) << std::setfill(' ') << "err_p"
              // << std::setw(15) << std::setfill(' ') << "conv p"
              // << std::setw(15) << std::setfill(' ') << "err u"
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
    for (int i = 0; i < h.size(); ++i) {
        std::cout << std::left << std::setw(10) << std::setfill(' ')
                  << h[i]
                  // << std::setw(15) << std::setfill(' ') << pl2[i]
                  // << std::setw(15) << std::setfill(' ') << convp[i]
                  // << std::setw(15) << std::setfill(' ') << ul2[i]
                  // << std::setw(15) << std::setfill(' ') << convu[i]
                  << std::setw(15) << std::setfill(' ')
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
#endif
