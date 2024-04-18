
#include "../cutfem.hpp"

using mesh_t     = Mesh2;
using funtest_t  = TestFunction<mesh_t>;
using fct_t      = FunFEM<mesh_t>;
using cutmesh_t  = ActiveMesh<mesh_t>;
using space_t    = GFESpace<mesh_t>;
using cutspace_t = CutFESpace<mesh_t>;

double sdBox(const R2 p, const R2 b) {
    R2 absp(abs(p.x), abs(p.y));
    R2 d = absp - b;
    R2 maxd0(max(d.x, 0.0), max(d.y, 0.0));
    return maxd0.norm() + min(max(d.x, d.y), 0.0);
}
double fun_levelSet(const R2 P, const int i) {
    const R2 b(2.0, 1.0);
    const R2 shift(4., 0.);
    return sdBox(P - shift, b) - 0.01;
}

double fun_y(const R2 P, int i, int dom) { return P.y; }

double fun_div(const R2 P, int i, int dom) {
    double x = P.x;
    double y = P.y;
    return 0;
}

double fun_rhs(const R2 P, int i, int dom) {
    double x = P.x;
    double y = P.y;
    if (i == 0)
        return 0;
    else if (i == 1)
        return 0;
}

double fun_uin(const R2 P, int i, int dom) {
    double x = P.x;
    double y = P.y;
    if (i == 0)
        return y * (2 - y) / 2;
    else
        return 0;
}

double fun_uout(const R2 P, int i, int dom) {
    double x = P.x;
    double y = P.y;
    if (i == 0)
        return 1. / 4 * (2 - y) * (y - 1); // original paper 4*(2-y)*(y-1)
    else
        return 0;
}
// double fun_p(const R2 P, int i, int dom ) {
//   double x = P.x;
//   double y = P.y;
//   return 10*(x*x-y*y)*(x*x-y*y);// - 670.171/(pi*interfaceRad*interfaceRad);
// }
} // namespace Erik_Data_CORIOLIS_STOKESRT

using namespace Erik_Data_CORIOLIS_STOKESRT;

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
    // int d = 2;

    vector<double> ul2, pl2, divmax, divl2, h, convu, convp;

    int iters = 4;
    for (int i = 0; i < iters; ++i) { // i<3

        // std::cout << "\n ------------------------------------- " << std::endl;
        Mesh Kh(nx, ny, 0., 0., 4., 2.);
        // Mesh Kh("../Lmesh"+to_string(i)+".msh");
        // if(i==1) Mesh Kh("../Lmesh1.msh");
        // else if(i==2) Mesh Kh("../Lmesh2.msh");

        const double hi = 1. / (nx - 1);

        Space Lh(Kh, DataFE<Mesh2>::P1);
        Fun_h levelSet(Lh, fun_levelSet);
        InterfaceLevelSet<Mesh> interface(Kh, levelSet);

        Lagrange2 FEvelocity(4);
        Space VELh_(Kh, FEvelocity);
        Space SCAh_(Kh, DataFE<Mesh>::P2);

        Space Wh(Kh, DataFE<Mesh2>::BDM1);
        // Space Wh(Kh, FEvelocity);
        Space Qh(Kh, DataFE<Mesh2>::P0); // FOR MIXEDSPACE

        ActiveMesh<Mesh> Khi(Kh);
        Khi.truncate(interface, -1);

        MacroElement<Mesh> macro(Khi, 1.);

        CutSpace VELh(Khi, VELh_);
        CutSpace SCAh(Khi, SCAh_);

        CutSpace Vh(Khi, Wh);
        CutSpace Ph(Khi, Qh);

        Fun_h fh(VELh, fun_rhs); // interpolates fun_rhs to fh of type Fun_h
        Fun_h uhin(VELh, fun_uin);
        Fun_h uhout(VELh, fun_uout);
        Fun_h yfunh(VELh, fun_y);
        // Fun_h exactph(SCAh, fun_exact_p); ExpressionFunFEM<Mesh> exactp(exactph,0,op_id);

        CutFEM<Mesh2> stokes(Vh);
        stokes.add(Ph);

        Normal n;
        Tangent t;
        /* Syntax:
        FunTest (fem space, #components, place in space)
        */
        FunTest u(Vh, 2, 0), p(Ph, 1, 0), v(Vh, 2, 0), q(Ph, 1, 0);
        FunTest u1(Vh, 1, 0), u2(Vh, 1, 1), v1(Vh, 1, 0), v2(Vh, 1, 1);

        double mu    = 0.01;
        double omega = 100; // 50 // original paper 100;

        // [grad] DOESNT WORK AT ALL!!!?
        {
            stokes.addBilinear(contractProduct(mu * grad(u), grad(v)) - innerProduct(p, div(v)) +
                                   innerProduct(div(u), q)

                                   - innerProduct(2 * omega * u2, v1) + innerProduct(2 * omega * u1, v2),
                               Khi);
            stokes.addLinear(innerProduct(fh.expression(2), v), Khi);
            const double sigma      = 1e0; // 4
            const double pp_tangent = 1e0; // 4=sigma, 400, 400
            const double pp_normal  = 4e2; // 0.8, 4, 40
            double uPenParam        = 1e0; // 4, 1e-1, 1e1
            double pPenParam        = 1e0; // 4, 1e-1, 1e1
            // [NECESSARY FOR RT / BDM ELEMENTS]
            stokes.addBilinear(-innerProduct(average(grad(u * t) * n, 0.5, 0.5), jump(v * t)) +
                                   innerProduct(jump(u * t), average(grad(v * t) * n, 0.5, 0.5)) +
                                   innerProduct(1. / hi * (sigma * jump(u * t)), jump(v * t)),
                               Khi, INTEGRAL_INNER_EDGE_2D);
            // const double pp_13 = 1e2; // 1e2, 4e2
            // const double pp_24 = 1e2; // 1e2, 4e3
            // const double pp_itf = 1e2; // 1e2, 1e2
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v) // natural
                                   + innerProduct(u, mu * grad(v) * n) +
                                   innerProduct(1. / hi * pp_tangent * u * t, v * t) +
                                   innerProduct(1. / hi * pp_normal * u * n, v * n) + innerProduct(p, v * n),
                               interface);
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v) // natural
                                   + innerProduct(u, mu * grad(v) * n) +
                                   innerProduct(1. / hi * pp_tangent * u * t, v * t) +
                                   innerProduct(1. / hi * pp_normal * u * n, v * n) + innerProduct(p, v * n),
                               Khi, INTEGRAL_BOUNDARY, {2, 4});
            stokes.addBilinear(-innerProduct(mu * grad(u) * n, v) // natural
                                   + innerProduct(u, mu * grad(v) * n) +
                                   innerProduct(1. / hi * pp_tangent * u * t, v * t) +
                                   innerProduct(1. / hi * pp_normal * u * n, v * n) + innerProduct(p, v * n),
                               Khi, INTEGRAL_BOUNDARY, {1, 3});
            stokes.addLinear(+innerProduct(uhin.expression(2), mu * grad(v) * n) +
                                 innerProduct(uhin * t, 1. / hi * pp_tangent * v * t) +
                                 innerProduct(uhin * n, 1. / hi * pp_normal * v * n),
                             Khi, INTEGRAL_BOUNDARY, {4});
            stokes.addLinear(+innerProduct(uhout.expression(2), mu * grad(v) * n) +
                                 innerProduct(uhout * t, 1. / hi * pp_tangent * v * t) +
                                 innerProduct(uhout * n, 1. / hi * pp_normal * v * n),
                             Khi, INTEGRAL_BOUNDARY, {2});
            // [Sets uniqueness of the pressure]
            // stokes.addLagrangeMultiplier(
            //   innerProduct(1.,p), 0
            //   , Khi
            // );
            CutFEM<Mesh2> lagr(Vh);
            lagr.add(Ph);
            Rn zero_vec = lagr.rhs_;
            lagr.addLinear(innerProduct(1, p), Khi);
            Rn lag_row(lagr.rhs_);
            lagr.rhs_ = zero_vec;
            lagr.addLinear(innerProduct(1, v * n), interface);
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , Khi, INTEGRAL_BOUNDARY
            // );
            stokes.addLagrangeVecToRowAndCol(lag_row, lagr.rhs_, 0);

            FunTest grad2un = grad(grad(u) * n) * n;
            // FunTest grad2pn = grad(grad(p)*n)*n;
            // FunTest grad2divun = grad(grad(div(u))*n)*n;
            stokes.addFaceStabilization(
                // [2k-1/2k+1 and changing to 2k+1/2k+1 gives the same bad effect!]
                // +innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(v))
                // +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(v)*n))
                // +innerProduct(uPenParam*pow(hi,3)*jump(grad2un), jump(grad2un))
                // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
                // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

                // [Pressure robust!]
                +innerProduct(uPenParam * pow(hi, -1) * jump(u), jump(v)) +
                    innerProduct(uPenParam * pow(hi, 1) * jump(grad(u) * n), jump(grad(v) * n)) +
                    innerProduct(uPenParam * pow(hi, 3) * jump(grad2un), jump(grad2un)) -
                    innerProduct(pPenParam * pow(hi, 1) * jump(p), jump(div(v))) +
                    innerProduct(pPenParam * pow(hi, 1) * jump(div(u)), jump(q)) -
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(p)), jump(grad(div(v)))) +
                    innerProduct(pPenParam * pow(hi, 3) * jump(grad(div(v))), jump(grad(q)))
                // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
                // +innerProduct(pPenParam*pow(hi,5)*jump(grad2divun) , jump(grad2pn))
                ,
                Khi, macro);
        }

        // [grad -> eps]
        {
            // stokes.addBilinear(
            //   contractProduct(2*mu*Eps(u),Eps(v))
            //   - innerProduct(p,div(v))
            //   + innerProduct(div(u),q)

            //   - innerProduct(2*omega*u2,v1)
            //   + innerProduct(2*omega*u1,v2)
            //   , Khi
            // );
            // stokes.addLinear(
            //   innerProduct(fh.expression(2),v)
            //   , Khi
            // );
            // stokes.addBilinear(
            //   - innerProduct(average(2*mu*(Eps(u)*t)*n,0.5,0.5), jump(v*t))
            //   + innerProduct(jump(u*t), average(2*mu*(Eps(v)*t)*n,0.5,0.5))
            //   + innerProduct(1./hi*sigma*(jump(u*t)), jump(v*t))
            //   , Khi
            //   , INTEGRAL_INNER_EDGE_2D
            // );
            // stokes.addBilinear(
            //   // - innerProduct(2*mu*Eps(u)*n, v)  //natural
            //   // + innerProduct(u, 2*mu*Eps(v)*n)
            //   + innerProduct(10*1./hi*penaltyParam*u, v)
            //   , interface
            // );
            // stokes.addBilinear(
            //   // - innerProduct(2*mu*Eps(u)*n, v)  //natural
            //   // + innerProduct(u, 2*mu*Eps(v)*n)
            //   + innerProduct(1./hi*penaltyParam*u, v)
            //   , Khi, INTEGRAL_BOUNDARY, {1,3}
            // );
            // stokes.addBilinear(
            //   - innerProduct(2*mu*Eps(u)*n, v)  //natural
            //   + innerProduct(u, 2*mu*Eps(v)*n)
            //   + innerProduct(1./hi*penaltyParam*u, v)
            //   , Khi, INTEGRAL_BOUNDARY, {2,4}
            // );
            // stokes.addLinear(
            //   + innerProduct(uhout.expression(2), 2*mu*Eps(v)*n)
            //   + innerProduct(uhout.expression(2), 1./hi*penaltyParam*v)

            //   , Khi, INTEGRAL_BOUNDARY, {2}
            // );
            // stokes.addLinear(
            //   + innerProduct(uhin.expression(2), 2*mu*Eps(v)*n)
            //   + innerProduct(uhin.expression(2), 1./hi*penaltyParam*v)

            //   , Khi, INTEGRAL_BOUNDARY, {4}
            // );

            // // [Sets uniqueness of the pressure]
            // // double meanP = integral(Khi,exactp,0);
            // // stokes.addLagrangeMultiplier(
            // //   innerProduct(1.,p), 0
            // //   , Khi
            // // );
            // CutFEM<Mesh2> lagr(Vh); lagr.add(Ph);
            // Rn zero_vec = lagr.rhs_;
            // lagr.addLinear(
            //   innerProduct(1.,p)
            //   , Khi
            // );
            // Rn lag_row(lagr.rhs_); lagr.rhs_ = zero_vec;
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , interface
            // );
            // lagr.addLinear(
            //   innerProduct(1, v*n)
            //   , Khi, INTEGRAL_BOUNDARY
            // );
            // stokes.addLagrangeVecToRowAndCol(lag_row,lagr.rhs_,0);

            // FunTest grad2un = grad(grad(u)*n)*n;
            // FunTest grad2pn = grad(grad(p)*n)*n;
            // // FunTest grad2divun = grad(grad(div(u))*n)*n;
            // stokes.addFaceStabilization( // [h^(2k+1) h^(2k+1)]
            //   // +innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(v))
            //   // +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(v)*n))
            //   // // +innerProduct(uPenParam*pow(hi,5)*jump(grad2un), jump(grad2un))
            //   // +innerProduct(pPenParam*pow(hi,1)*jump(p), jump(q))
            //   // +innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(q)))

            //   +innerProduct(uPenParam*pow(hi,-1)*jump(u), jump(v))
            //   +innerProduct(uPenParam*pow(hi,1)*jump(grad(u)*n), jump(grad(v)*n))
            //   // +innerProduct(uPenParam*pow(hi,3)*jump(grad2un), jump(grad2un))
            //   -innerProduct(pPenParam*pow(hi,1)*jump(p), jump(div(v)))
            //   +innerProduct(pPenParam*pow(hi,1)*jump(div(u)), jump(q))
            //   -innerProduct(pPenParam*pow(hi,3)*jump(grad(p)), jump(grad(div(v))))
            //   +innerProduct(pPenParam*pow(hi,3)*jump(grad(div(v))) , jump(grad(q)))
            //   // -innerProduct(pPenParam*pow(hi,5)*jump(grad2pn), jump(grad2divun))
            //   // +innerProduct(pPenParam*pow(hi,5)*jump(grad2divun) , jump(grad2pn))
            //   , Khi
            //   , macro
            // );
        }

        // std::cout << integral(Khi,exactp,0) << std::endl;
        matlab::Export(stokes.mat_, "mat" + to_string(i) + "Cut.dat");
        stokes.solve();

        // EXTRACT SOLUTION
        int nb_flux_dof = Vh.get_nb_dof();
        Rn_ data_uh     = stokes.rhs_(SubArray(nb_flux_dof, 0));
        Rn_ data_ph     = stokes.rhs_(SubArray(Ph.get_nb_dof(), nb_flux_dof));
        Fun_h uh(Vh, data_uh);
        Fun_h ph(Ph, data_ph);

        // [Post process pressure]
        // double meanP = integral(Khi,p0,0);
        // ExpressionFunFEM<Mesh> fem_p(ph,0,op_id);
        // double meanPfem = integral(Khi,fem_p,0);
        // // std::cout << meanP << std::endl;
        // CutFEM<Mesh2> post(Qh);
        // post.addLinear(
        //   innerProduct(1,q)
        //   , Khi
        // );
        // double area = post.rhs_.sum();
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

        // double errU      = L2normCut(uh,fun_exact_u,0,2);
        // double errP      = L2normCut(ph,fun_exact_p,0,1);
        double errDiv    = L2normCut(uh_0dx + uh_1dy, Khi);
        double maxErrDiv = maxNormCut(uh_0dx + uh_1dy, Khi);
        // double errDiv    = L2normCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);
        // double maxErrDiv = maxNormCut(femSol_0dx+femSol_1dy+fflambdah,fun_div,Khi);

        // // solExVec -= stokesDiv.rhs;
        // for(int i=0;i<solExVec.size();++i){
        //   solExVec(i) = fabs(solExVec(i)-stokesDiv.rhs(i));
        // }
        //
        // Fun_h solEx(mixedSpace, solExVec);
        // writerS.add(solEx,"uh_err", 0, 2);

        // writerS.add(solEx,"uh_ex", 0, 2);
        // writerS.add(solEx,"ph_ex", 2, 1);

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
