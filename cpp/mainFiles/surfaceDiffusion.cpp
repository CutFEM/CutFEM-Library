//
// Created by Sebastian Myrbäck on 2022-01-26.
//

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#  include "cfmpi.hpp"
#include "finiteElement.hpp"
#include "baseCutProblem.hpp"
#include "../num/gnuplot.hpp"
#include "projection.hpp"
#include <typeinfo>

// ok for the interface problem

//// Problem 1

// RHS
//double fun_rhs0(const R2 P, int elementComp) {
//    return (9*P.y*(3*P.x*P.x - P.y*P.y))/(P.x*P.x + P.y*P.y);
//}
//
//// Exact solution
//double fun_uSurface(const R2 P,  int elementComp) {
//    return 3.0*P.x*P.x*P.y - P.y*P.y*P.y;
//}

//// Problem 2
// RHS
double fun_rhs0(const R2 P, int elementComp) {
  return 4*P.x*P.y/(P.x*P.x+P.y*P.y);
    // return -(6*P.x*P.y*(pow(P.x,4) - 4*P.x*P.x*P.y*P.y + pow(P.y,4)))/(pow(P.x*P.x + P.y*P.y,4));
}

// Exact solution
double fun_uSurface(const R2 P,  int elementComp) {
  return P.x*P.y;
    // return pow(P.x,3)*pow(P.y,3)/(pow(P.x*P.x+P.y*P.y,3));
}
double fun_uSurfaceT(const R2 P,  int elementComp, double t) {
  return P.x*P.y;
    // return pow(P.x,3)*pow(P.y,3)/(pow(P.x*P.x+P.y*P.y,3));
}

// Level-set function
double fun_levelSet(const R2 P, const int i) {
    return sqrt((P.x)*(P.x) + (P.y)*(P.y)) - 1.0;
}

static void exactRHSintegration(const FESpace2& Vh, const Interface2& interface, Rn& rhs, R (*f)(const R2, int )){
    typedef typename FESpace2::FElement FElement;
    typedef typename Mesh2::Element Element;
    typedef typename FElement::QFB QFB;
    typedef typename Mesh2::Partition Partition;
    typedef typename QFB::QuadraturePoint QuadraturePoint;


    KNMK<double> fv(Vh[0].NbDoF(),1,1); //  the value for basic fonction
    What_d Fop = Fwhatd(1);
    const QFB &qfb(*QF_Simplex<R1>(5));

    for(int iface=interface.first_element(); iface<interface.last_element(); iface+=interface.next_element()) {


        const int kb = interface.idxElementOfFace(iface);   // idx on
        const typename Interface2::FaceIdx& face = interface[iface];  // the face
        const R meas = interface.computeDx(face).norm();
        const double h = meas;
        const R2 linear_normal(-interface.normal(iface));

        int lastop = getLastop(op_id, op_id);
        const int kv = Vh.idxElementFromBackMesh(kb,0);
        const FElement& FKv(Vh[kv]);
        double measK = FKv.getMeasure();


        for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

            QuadraturePoint ip(qfb[ipq]); // integration point
            const R2 mip = interface.mapToFace(face,(R1)ip);
            const R Cint = meas * ip.getWeight();
            FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element

            double val_fh = f(mip,0);
            double Cst = Cint;

            for(int i = FKv.dfcbegin(0); i < FKv.dfcend(0); ++i) {
                rhs(FKv.loc2glb(i)) +=  Cint * val_fh * fv(i,0,op_id);
            }
        }
    }
}

typedef Mesh2 Mesh;
typedef FESpace2 FESpace;
typedef TestFunction<2> FunTest;
typedef FunFEM<Mesh2> Fun_h;

void errorAnalysis(int argc, char** argv) {
    // INITIALIZE MPI
    // =====================================================
    MPIcf cfMPI(argc, argv);

    const int nms = 40;

    int mesh_sizes[nms];
    double errors_interface[nms];


    for (int i = 0; i < nms; i++) {

        //std::cout << n << ", ";
        mesh_sizes[i] = 10 + 20*i;

        // MESH AND PROBLEM PARAMETERS
        // =====================================================
        int nx = mesh_sizes[i], ny = mesh_sizes[i];
        double lx = 3., ly = 3.;
        double h = lx/nx;

        std::cout << "n = " << nx << std::endl;

        double betaE = 0.25, betaF = 0.25, gamma = 10;
        double cF = 0.25, cGamma = 0.25;

        // CONSTRUCTION OF THE MESH
        // =====================================================
        Mesh Th(nx, ny, -1.5, -1.5, lx, ly);


        // CONSTRUCTION OF THE LEVELSET AND THE INTERFACE
        // =====================================================
        FESpace Lh(Th, DataFE<Mesh>::P1); //dc are the discontinuous elements
        Fun_h levelSet(Lh, fun_levelSet);
        Interface2 interface(Th, levelSet.v);

        // CONSTRUCTION OF THE FE SPACE AND THE CUT SPACE

        //FESpace Vh(Th, DataFE<Mesh>::P1dc);     // background mesh FE Space
        Mesh2 cutTh(interface);                   // Mesh for interface
        FESpace Sh(cutTh, interface, DataFE<Mesh>::P1dc);  // FE Space for the interface mesh
        //Sh.backSpace = &Vh;

        // OBJECTS NEEDED FOR THE PROBLEM
        // =====================================================
        CutFEM<Mesh> surfdiff(Sh);  // define problem on interface space Sh

        Normal n;
        Tangent t;

        KN_<double> data_sh(surfdiff.rhs);

        Fun_h us(Sh, data_sh);  // function representing surface solution
        Fun_h fh0(Sh, fun_rhs0);  // create a function on the cutMesh for the rhs

        FunTest u0(Sh,1), v0(Sh,1);   // Test functions on Omega 0, i.e on the interface.

        //// ----------- ASSEMBLY OF LINEAR SYSTEM ------------- ////

        // Integral on interface for surface variable (K_{h,0})
        surfdiff.addBilinear(
                innerProduct(gradS(u0),gradS(v0))            // (3.12)
                , interface
        );

        /// ---------- Point evaluation on E_{h,0} ------------ ///

//        surfdiff.addBilinear(
//                - innerProduct(average(grad(u0)*t, 0.5,-0.5),jump(v0)) // (3.12)/(3.81)
//                - innerProduct(jump(u0), average(grad(v0)*t,0.5,-0.5))    // (3.13)/(3.82) added for symmetry
//                + innerProduct(betaE/h*jump(u0),jump(v0))                             // (3.13)/(3.82)
//                , interface
//                , nodeEvaluation
//        );

        // ------------- Stabilization of the SURFACE --------------- //

        // Full stabililzation

        // Paper variant
//        surfdiff.addBilinear(
//                innerProduct(betaF/h/h*jump(u0), jump(v0))
//                + innerProduct(gamma*jump(grad(u0))*n, jump(grad(v0))*n)
//                , innerEdge
//        );

        // Sara's variant
        surfdiff.addBilinear(
                innerProduct(cF*h*h*jump(grad(u0)*n), jump(grad(v0)*n))
                , innerEdge
        );
        surfdiff.addBilinear(
                innerProduct(cGamma*h*h*grad(u0), grad(v0))
                , interface
        );

        // Add linear part on faces (rhs)
//        surfdiff.addLinear(
//                innerProduct(fh0.expression(),v0),
//                interface
//        );

        exactRHSintegration(Sh,interface, surfdiff.rhs, fun_rhs0);

        surfdiff.addLagrangeMultiplier(
                innerProduct(1.,v0),
                interface, 0
        );

        // Solve linear system
        surfdiff.solve();

        // Compute L2 error
        R errU0 = L2norm(us, fun_uSurface, 0,1);

        std::cout << " || u0 - u0_ex||_2  = \t" << errU0 << std::endl;

        errors_interface[i] = errU0;
        std::cout << "error surf = " << errU0 << std::endl;
    }

    for (auto mesh:mesh_sizes) {
        std::cout << mesh << ", ";
    }

    // Write errors to binary file
    FILE *fp;
    fp = fopen("../../outputFiles/surfDiff/errorsDiff.bin", "wb");
    fwrite(errors_interface, sizeof(double), nms, fp);
    fclose(fp);
}

void solve(int argc, char** argv) {

    // INITIALIZE MPI
    // =====================================================
    MPIcf cfMPI(argc,argv);

    int nx = 21, ny = 21;
    double lx = 3., ly = 3.;
    double h = lx/nx;
    std::cout << "h = " << h << std::endl;

    // Stabilization parameters
    double betaE = 1000, betaF = 1, gamma = 1;
    double cF = 0.25, cGamma = 0.25;

    // CONSTRUCTION OF THE MESH
    // =====================================================

    Mesh Th(nx, ny, -1.5, -1.5, lx, ly);

    // CONSTRUCTION OF THE LEVELSET AND THE INTERFACE
    // =====================================================
    FESpace Lh(Th, DataFE<Mesh>::P1); //dc are the discontinuous elements
    Fun_h levelSet(Lh, fun_levelSet);
    Interface2 interface(Th, levelSet.v);

    // CONSTRUCTION OF THE FE SPACE AND THE CUT SPACE
    //FESpace Vh(Th, DataFE<Mesh>::P1dc);     // background mesh FE Space
    Mesh2 cutTh(interface);                   // Mesh for interface
    FESpace Sh(cutTh, interface, DataFE<Mesh>::P1);  // FE Space for the interface mesh
    //Sh.backSpace = &Vh;

    /*
    gnuplot::save(Th, "../../outputFiles/surfaceDiffusion/Th.dat");
    gnuplot::save(cutTh, "../../outputFiles/surfaceDiffusion/cutTh.dat");
    gnuplot::save(interface, "../../outputFiles/surfaceDiffusion/interface.dat");
    */

    // OBJECTS NEEDED FOR THE PROBLEM
    // =====================================================
    CutFEM<Mesh> surfdiff(Sh);  // define problem on interface space Sh

    Normal n;
    Tangent t;
    int idx_s0 = Sh.NbDoF(); // where v0 start in the solution array
    KN_<double> data_sh(surfdiff.rhs);

    Fun_h us(Sh, data_sh);  // function representing surface solution
    Fun_h fh0(Sh, fun_rhs0);  // create a function on the cutMesh for the rhs

    FunTest u0(Sh,1), v0(Sh,1);   // Trial and test functions on Omega 0, i.e on the interface.

    //// ----------- ASSEMBLY OF LINEAR SYSTEM ------------- ////

    // Integral on interface for surface variable (K_{h,0})
    surfdiff.addBilinear(
            innerProduct(gradS(u0),gradS(v0))            // (3.12)
            , interface
    );

    /// ---------- Point evaluation on E_{h,0} ------------ ///

//        surfdiff.addBilinear(
//                - innerProduct(average(grad(u0)*t, 0.5,-0.5),jump(v0)) // (3.12)/(3.81)
//                - innerProduct(jump(u0), average(grad(v0)*t,0.5,-0.5))    // (3.13)/(3.82) added for symmetry
//                + innerProduct(betaE/h*jump(u0),jump(v0))                             // (3.13)/(3.82)
//                , interface
//                , nodeEvaluation
//        );

    // ------------- Stabilization of the SURFACE --------------- //

    // Full stabililzation

    // Paper variant
//        surfdiff.addBilinear(
//                innerProduct(betaF/h/h*jump(u0), jump(v0))
//                + innerProduct(gamma*jump(grad(u0))*n, jump(grad(v0))*n)
//                , innerEdge
//        );

    // Sara's variant
    surfdiff.addBilinear(
            innerProduct(1e-2*jump(grad(u0)*n), jump(grad(v0)*n))
            , innerEdge
    );
    surfdiff.addBilinear(
            innerProduct(1e-2* grad(u0)*n, grad(v0)*n)
            , interface
    );

    // Add linear part on faces (rhs)
//    surfdiff.addLinear(
//            innerProduct(fh0.expression(),v0),
//            interface
//    );

    exactRHSintegration(Sh,interface, surfdiff.rhs, fun_rhs0);

    surfdiff.addLagrangeMultiplier(
            innerProduct(1.,v0),
            interface, 0
    );

    // matlab::Export(surfdiff.mat, "matrix1.dat");
    // matlab::Export(surfdiff.rhs, "rhs1.dat");

    // Solve linear system
    surfdiff.solve();

    // Compute L2 error
    R errU0 = L2normSurf(us, fun_uSurfaceT, 0.,0.,1.);
    std::cout << " || u0 - u0_ex||_2  = \t" << errU0 << std::endl;


    // PRINT THE SOLUTION TO PARAVIEW
    // =====================================================
    if(MPIcf::IamMaster()){
        Fun_h uexsurf (Sh, fun_uSurface);
        Paraview2 writerU0(Sh, "surfDiff.vtk");
        writerU0.add(us, "surfDiff", 0, 1);
        writerU0.add(levelSet, "levelSet", 0, 1);
        writerU0.add(uexsurf, "surfDiff_ex", 0, 1);
    }

}

int main(int argc, char** argv ) {

    solve(argc, argv);

    return 1;
}
