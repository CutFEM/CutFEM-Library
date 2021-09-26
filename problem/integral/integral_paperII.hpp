#ifndef INTEGRAL_PAPERII_HPP_
#define INTEGRAL_PAPERII_HPP_


#include "../paperII.hpp"

template<class M>
class IntegralPaperII  {

  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef CutGFESpace<Mesh> CutFESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::Rd Rd;
  typedef typename FElement::RdHat RdHat;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::BorderElement BorderElement;
  typedef typename Mesh::Partition Partition;
  typedef GenericInterface<Mesh> Interface;
  typedef FEMFunction<FESpace> Fun;

  typedef typename TypeOperationMV<Rd::d>::OperationMV MyOp;
  typedef typename FElement::QFB QFB;
  typedef typename FElement::QF QF;
  typedef GQuadraturePoint<R1> QPT;

  const GPaperII<M>& problem;
  const CutFESpace& Vh;
  const FESpace& Sh;
  const QFB &qfb;
  const QF &qf;
  const QuadratureFormular1d &qft;
  const GenericInterface<Mesh>& interface;
  const Fun& mapping;
  const TimeSlab& In;
  const QPT tq;

  KNMK<double> w, wn, ws, wt;
  What_d whatd;
  R Edu = 1e-2, Edp = 1e-2, Eds = 1e-2;

public:

  IntegralPaperII(const GPaperII<M>& pp, const TimeSlab& IIn, const QPT ttq, const int itq = GTime::iter_step_quadrature, const Fun& mapp = DataFunFEM<M>::Id) :
  Vh(*pp.Vh),
  Sh(*pp.Sh),
  problem(pp),
  qfb(*QF_Simplex<RdHatBord>(3)),
  qf(*QF_Simplex<RdHat>(5)),
  qft(*pp.qft),
  interface(*(Vh.gamma[itq])),
  mapping(mapp),
  In(IIn),
  tq(ttq)
  {
    w.init(Vh[0].NbDoF(),Rd::d+1,op_dz+1);
    ws.init(Sh[0].NbDoF(),1,op_dz+1);
    wt.init(In.NbDoF(),1,op_dz+1);

    whatd = Fop_D0|Fop_D1;
    In.BF(tq.x, wt);                // compute time basic funtions

  }

  void initStab(R edu, R edp) {
    Edu = edu; Edp=edp;
    // w.init(nbDoF,op_Dall);
    wn.init(Vh[0].NbDoF(),Rd::d+1,op_dz+1);//,op_Dall);

    // whatd = Fop_D0|Fop_D1|Fop_D2;
  }
  void initStab(R eds) {
    Eds=eds;
    // w.init(nbDoF,op_Dall);
    wn.init(Sh[0].NbDoF(),1,op_dz+1);//,op_Dall);
    // whatd = Fop_D0|Fop_D1|Fop_D2;
  }

  void locNonLinAssemblyFull(const FElement& FK, Rnm& NLM, Rn& FL, const Partition& cutK, const Rn& w0) {

    typedef typename QF::QuadraturePoint QuadraturePoint;

    const int domain = FK.whichDomain();
    ElementSignEnum the_part = cutK.what_part(domain);
    const R rho_ = problem.rho(domain);
    const int nbDoF = FK.NbDoF();
    const int nbDoFTime = In.NbDoF();

    Rn valu0(Rd::d, 0.);       // all componants
    Rnm Dvalu0(Rd::d, Rd::d);  // all componants
    NLM = 0.0; FL = 0.0;

    for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
    it != cutK.element_end(the_part); ++it){

      const R mes = cutK.mesure(it);
      for(int ipq = 0; ipq < qf.n; ++ipq)  {

        QuadraturePoint ip(qf[ipq]);
        Rd mip = cutK.toK(it, ip);
        const R Cint = mes * ip.a * In.T.mesure()* tq.a ;

        FK.BF(FK.T.toKref(mip), w);

        valu0 = 0.; Dvalu0 = 0.;
        // compute u0(xq, tq)
        for(int it=0; it<nbDoFTime; ++it) {
          const R tval = wt(it,0,op_id);
          for(int ci=0;ci<Rd::d;++ci) { // dont need the pressure
            for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
              valu0(ci) += w0(i+it*nbDoF) * w(i,ci,op_id) * tval;
            }
          }
        }
        // compute gradu0(xq, tq)
        for(int it=0; it<nbDoFTime; ++it) {
          const R tval = wt(it,0,op_id);
          for(int ci=0;ci<Rd::d;++ci) { // dont need the pressure
            for(int cj=0;cj<Rd::d;++cj) {
              for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
                Dvalu0(ci, cj) += w0(i+it*nbDoF) * w(i,ci,MyOp::D(cj)) * tval;
              }
            }
          }
        }

        Rd uDivu;
        for(int d=0;d<Rd::d;++d) {
          uDivu[d] = 0;
          for(int ci=0;ci<Rd::d;++ci)
          uDivu[d] += valu0(ci) * Dvalu0(d , ci);
        }


        // (u0.grad)u0,v)
        for(int jt=0; jt<nbDoFTime; ++jt) {
          for(int ci=0;ci<Rd::d;++ci) { // dont need the pressure
            for(int j = FK.dfcbegin(ci); j < FK.dfcend(ci); ++j){
              FL(j+jt*nbDoF) += Cint * uDivu[ci] * w(j,ci,op_id) * wt(jt,0,op_id);
            }
          }
        }

        // (gradw0 . u + w0 div u, r)
        for(int it=0; it<nbDoFTime; ++it) {
          for(int jt=0; jt<nbDoFTime; ++jt) {
            const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
            for(int ci=0;ci<Rd::d;++ci) {
              for(int cj=0;cj<Rd::d;++cj) {

                for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
                  for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){

                    NLM(i+it*nbDoF, j+jt*nbDoF) += Cint *
                    (
                      (ci == cj) * MyOp::innerProduct(w, valu0, ci, j)
                      + w(j,cj,op_id) * Dvalu0( ci, cj)
                    ) * w(i,ci,op_id) * tval;
                  }
                }
              }
            }
          }
        }



      }

    }


  }

  void locFullAssembly(const FElement& FK, Rnm& ML, Rn& VL, Rn& VP, const Partition& cutK, const Rn& w0) {

    typedef typename QF::QuadraturePoint QuadraturePoint;
    ML = 0.0; VL = 0.0; VP = 0.0;

    const int domain = FK.whichDomain();
    ElementSignEnum the_part = cutK.what_part(domain);
    const R mu_ = problem.mu(domain);
    const R rho_ = problem.rho(domain);

    const int nbDoF = FK.NbDoF();
    const int nbDoFTime = In.NbDoF();
    const int cp  = Rd::d;           // the pressure componante

    Rn fxj(nbDoF*nbDoFTime, 0.0);
    problem.init_RHS(FK,fxj, tq.x);

    for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
    it != cutK.element_end(the_part); ++it){

      const R mes = cutK.mesure(it);
      for(int ipq = 0; ipq < qf.n; ++ipq)  {

        QuadraturePoint ip(qf[ipq]);
        Rd mip = cutK.toK(it, ip);
        const R Cint = mes * ip.a * In.T.mesure()* tq.a ;
        const R Cint0 = mes * ip.a;

        FK.BF(FK.T.toKref(mip), w);

        // rho(dt u, v ) + (E(u),E(v))
        for(int it=0; it<nbDoFTime; ++it) {
          for(int jt=0; jt<nbDoFTime; ++jt) {
            const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
            for(int ci=0; ci<Rd::d;++ci) {
              for(int cj=0; cj<Rd::d;++cj) {
                for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
                  for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
                    ML(i+it*nbDoF, j+jt*nbDoF) +=  Cint *
                    (	rho_ * w(i,ci,op_id)*wt(it,0,op_id) * w(j,cj,op_id)*wt(jt,0,op_dx)
                    + mu_ * tval * ((ci==cj)? MyOp::crossTerm(w,ci,i,j)
                    : w(i,ci,MyOp::D(cj)) * w(j,cj,MyOp::D(ci))) );
                  }
                }
              }
            }
          }
        }

        // // (grad p, v) - (u, grad q)
        // // for(int ci=0; ci<Rd::d;++ci) {
        // //   for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
        // //     for(int jp=FK.dfcbegin(cp); jp<FK.dfcend(cp); ++jp){
        // //       ML(i, jp) += Cint * w(jp,MyOp::D(ci)) * w(i,op_id);
        // //       ML(jp, i) -= Cint * w(jp,MyOp::D(ci)) * w(i,op_id);
        // //     }
        // //   }
        // // }

        // -(p, div v) + (div u, q)
        for(int it=0; it<nbDoFTime; ++it) {
          for(int jt=0; jt<nbDoFTime; ++jt) {
            const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
            for(int ci=0; ci<Rd::d;++ci) {
              for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
                for(int jp=FK.dfcbegin(cp); jp<FK.dfcend(cp); ++jp){
                  ML(i+it*nbDoF, jp+jt*nbDoF) -= Cint * w(jp,cp,op_id) * w(i,ci,MyOp::D(ci)) * tval;
                  ML(jp+jt*nbDoF, i+it*nbDoF) += Cint * w(jp,cp,op_id) * w(i,ci,MyOp::D(ci)) * tval;
                }
              }
            }
          }
        }

        // (f,v)
        // for(int it=0; it<nbDoFTime; ++it) {
        for(int jt=0; jt<nbDoFTime; ++jt) {
          const R tval = wt(jt,0,op_id);//*wt(it,op_id);
          for(int ci=0; ci<Rd::d;++ci) {
            for(int j = FK.dfcbegin(ci); j < FK.dfcend(ci); ++j){
              // for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
              VL(j+jt*nbDoF) +=
              Cint  * rho_ * w(j,ci,op_id) * problem.fun_source(mip, tq.x)[ci] * tval;
              //  Cint  * rho_ * w(i,op_id) * w(j,op_id) * fxj(i+it*nbDoF) * tval;
            }
          }
        }
        //	  }
        // }

        //(p, 1)
        if(domain == 0 && GTime::iter_step_quadrature == qft.n-1) {
          for(int jt=0; jt<nbDoFTime; ++jt) {
            for(int jp=FK.dfcbegin(cp); jp<FK.dfcend(cp); ++jp){
              VP(jp + jt*nbDoF) += Cint0 * w(jp,cp,op_id) * wt(jt,0,op_id);
            }
          }
        }
        if(GTime::iter_step_quadrature == 0) {
          for(int ci=0; ci<Rd::d;++ci) {
            for(int j = FK.dfcbegin(ci); j < FK.dfcend(ci); ++j){
              for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
                ML(i, j) += Cint0 * w(i,ci,op_id) * w(j,ci,op_id) * rho_;
                VL(j)    += Cint0 * w(j,ci,op_id) * w(i,ci,op_id) * w0(i) * rho_;

              }
            }
          }
        }


      } // loop quad point
    }

  }




  void locSurfAssembly(const FElement& FK, const FElement& FKs, const int iface, Rnm& ML, Rn& VL, Fun& curvature, const Rn& w0) {

    assert(problem.sigma);

    typedef typename QFB::QuadraturePoint QuadraturePoint;
    const typename Interface::FaceIdx& face = interface[iface];  // the face
    const R mes = interface.computeDx(face).norm();    // integration constant
    Rd normal(interface.normal(iface));
    normal /= -normal.norm();

    const int nbDoF = FK.NbDoF();
    const int nbDoFs = FKs.NbDoF();
    const int nbDoFTime = In.NbDoF();
    const int cp = Rd::d;
    const int idxs0 = 2*(nbDoF*nbDoFTime) ;
    R eps = problem.eps;

    KNM<R> DUn(nbDoF, nbDoF);
    const int jumpSign[2] = {1, -1};

    const R h    = Vh[0].T.lenEdge(0);
    const R invh = 1./h;
    const R invh2 = invh * invh;
    const R gamma = mes * invh;
    const R lambdaG = (problem.kappa1*problem.mu1
      + problem.kappa2*problem.mu2)*invh2*(100 + 10*gamma);

      ML = 0.0; VL = 0.0;
      for(int ipq = 0; ipq < qfb.n; ++ipq)  {             // loop over integration points

        QuadraturePoint ip(qfb[ipq]);
        const Rd mip = interface.mapToFace(face,(RdHatBord)ip);

        const int kcurv = curvature.Vh->idxElementFromBackMesh(Vh.Th(FK.T));
        // Rd kapN = curvature.eval(kcurv, mip); kapN = -kapN;
        Rd kapN = 4*normal;


        const R Cint  = mes * ip.a * In.T.mesure()* tq.a ;
        const R Cint0 = mes * ip.a;

        FK.BF(FK.T.toKref(mip), w);
        FKs.BF(FKs.T.toKref(mip), ws);

        DUn = 0;
        for(int ci=0; ci<Rd::d;++ci) {
          for(int cj=0; cj<Rd::d;++cj) {
            for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
              for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
                DUn(i,j) +=  ((ci==cj)? MyOp::Dxn(w,normal,ci,i,j)
                : w(i,ci,op_id) * w(j,cj,MyOp::D(ci))*normal[cj]);
              }
            }
          }
        }


        // -({mu E(u).n},[[v]]  ) - ([[u]],{mu E(v).n}  ) + lambdaG ([[u]], [[v]])
        for(int di=0;di<2;++di) {
          const int i0 = di*nbDoF*nbDoFTime;
          for(int dj=0;dj<2;++dj) {
            const int j0 = dj*nbDoF*nbDoFTime;

            for(int it=0; it<nbDoFTime; ++it) {
              for(int jt=0; jt<nbDoFTime; ++jt) {
                const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
                for(int ci=0;ci<Rd::d;++ci) {
                  for(int cj=0;cj<Rd::d;++cj) {

                    for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
                      for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
                        ML(i0+i+it*nbDoF, j0+j+jt*nbDoF) += Cint
                        * ( (-1.)
                        * (
                          problem.mu(di) * problem.kappa(di) * DUn(j,i) * jumpSign[dj]  +
                          problem.mu(dj) * problem.kappa(dj) * DUn(i,j) * jumpSign[di]
                        ) +
                        (ci == cj)*lambdaG * jumpSign[di] * w(i,ci,op_id) *
                        jumpSign[dj] * w(j,cj,op_id)
                      )
                      * tval;
                    }
                  }
                }
              }
            }
          }

          //  divergence form
          // ({p}, [[v.n]] ) - ([[u.n]], {q} )
          for(int it=0; it<nbDoFTime; ++it) {
            for(int jt=0; jt<nbDoFTime; ++jt) {
              const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
              for(int ci=0;ci<Rd::d;++ci) {
                for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
                  for(int j = FK.dfcbegin(cp); j < FK.dfcend(cp); ++j){

                    ML(i0+i+it*nbDoF, j0+j+jt*nbDoF) +=
                    Cint * (
                      w(i,ci,op_id) * normal[ci] * jumpSign[di] *
                      problem.kappa(dj) * w(j,cp,op_id)
                    ) * tval;
                    ML(j0+j+jt*nbDoF, i0+i+it*nbDoF) -=
                    Cint  * (
                      w(i,ci,op_id) * normal[ci] * jumpSign[di] *
                      problem.kappa(dj) * w(j,cp,op_id)
                    ) * tval;
                  }
                }
              }
            }
          }


          // 	  // grad form
          // // 	  // (<u.n>, [[q]]) - ([[p]], <v.n>)
          // // 	//   for(int ci=0; ci<Rd::d;++ci) {
          // // 	//     for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
          // // 	//       for(int jp=FK.dfcbegin(cp); jp<FK.dfcend(cp); ++jp){

          // // 	// 	ML(i0+i, j0+jp) -=
          // // 	// 	  Cint * (
          // // 	// 		  w(i,op_id) * normal[ci] * jumpSign[dj] *
          // // 	// 		problem.kappaO(di) * w(jp,op_id)
          // // 	// 		);
          // // 	// 	ML(j0+jp, i0+i) +=
          // // 	// 	  Cint  * (
          // // 	// 		   w(i,op_id) * normal[ci] * jumpSign[dj] *
          // // 	// 		   problem.kappaO(di) * w(jp,op_id)
          // // 	// 		   );


          // //       	//     }
          // //       	//   }
          // //       	// }
        }
      }

      //  surfactant integral
      // (dt w, r ) + eps(gradG w, gradG r )
      for(int it=0; it<nbDoFTime; ++it) {
        for(int jt=0; jt<nbDoFTime; ++jt) {
          const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
          for(int i = FKs.dfcbegin(0); i < FKs.dfcend(0); ++i){
            for(int j = FKs.dfcbegin(0); j < FKs.dfcend(0); ++j){

              ML(idxs0+i+it*nbDoFs, idxs0+j+jt*nbDoFs) += Cint *
              ( ws(i,0,op_id)*wt(it,0,op_id) * ws(j,0,op_id)*wt(jt,0,op_dx)
              + eps * (MyOp::gradSurface(ws,normal,i),  MyOp::gradSurface(ws,normal,j))
              * tval
            );

          }
        }
      }
    }



    // -(sigma(w) K, <v.n>) + (gradG sigma(w), <v>)
    // for(int di=0;di<2;++di) {
    //   const int i0 = di*nbDoF*nbDoFTime;;
    //   for(int it=0; it<nbDoFTime; ++it) {
    //     for(int jt=0; jt<nbDoFTime; ++jt) {
    //       const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
    //       for(int ci=0;ci<Rd::d;++ci) {
    //         for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
    //           for(int j = FKs.dfcbegin(0); j < FKs.dfcend(0); ++j){
    //             // Rd Dsigma;
    //             // for(int ll=0;ll<Rd::d;++ll) Dsigma[ll] = problem.sigma(ws(j,MyOp::D(ll)));
    //             // Rd Gs = MyOp::gradSurface(Dsigma,normal);
    //             Rd Gs = MyOp::gradSurface(ws,normal,j);
    //             Rd Ds;
    //             for(int ll=0;ll<Rd::d;++ll)
    //             Ds[ll]=problem.Dsigma(Gs[ll]);
    //
    //             ML(i0 + i+it*nbDoF, idxs0+j+jt*nbDoFs)  +=
    //             Cint * (
    //               (-1.) * w(i,ci,op_id) * problem.sigma(ws(j,0,op_id)) * kapN[ci]
    //               + Ds[ci] * w(i,ci,op_id)
    //             ) * problem.kappaO(di) * tval;
    //
    //           }
    //         }
    //       }
    //     }
    //   }
    // }


    // // no surfactant case
    // // (sigma K, <v.n>)
    for(int di=0;di<2;++di) {
    	const int i0 = di*nbDoF*nbDoFTime;;
    	for(int it=0; it<nbDoFTime; ++it) {
    	  for(int ci=0;ci<Rd::d;++ci) {
    	    for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
    	      VL[i0 + i+it*nbDoF]  += Cint * w(i,ci,op_id) * problem.sigma(0)
            * kapN[ci] * problem.kappaO(di) * wt(it,0,op_id);
    	    }
    	  }
    	}
    }



    if(GTime::iter_step_quadrature == 0) {
      for(int j = FKs.dfcbegin(0); j < FKs.dfcend(0); ++j){
        for(int i = FKs.dfcbegin(0); i < FKs.dfcend(0); ++i){

          ML(idxs0+i, idxs0+j) += Cint0 * ws(i,0,op_id) * ws(j,0,op_id) ;
          VL(idxs0 + j)        += Cint0 * ws(i,0,op_id) * ws(j,0,op_id) * w0(i) ;
        }
      }
    }



  } // end loop over quadrature point
}




void locBoundAssembly(const FElement& FK, const BorderElement& face, Rnm& ML, Rn& VL, const Rd normal) {
  ML = 0.0; VL = 0.0;
  typedef typename QFB::QuadraturePoint QuadraturePoint;

  const int domain = FK.whichDomain();
  const R mu_ = problem.mu(domain);
  const R rho_ = problem.rho(domain);

  const int nbDoF = FK.NbDoF();
  const int nbDoFTime = In.NbDoF();
  const int cp  = Rd::d;           // the pressure componante

  Rn gxj(nbDoF*nbDoFTime, 0.0);
  problem.init_boundary(FK, gxj, In.Pt(tq.x));


  KNM<R> DUn(nbDoF, nbDoF);
  const R mes = face.mesure();                       // integration constant

  const R h    = Vh[0].T.lenEdge(0);
  const R hyperH = pow(h, Rd::d-1);
  const R gammaK = mes / hyperH;
  const R invh = 1./h;
  const R invh2 = invh * invh;
  const R lambdaG = 10*invh2;
  //const R lambdaG = mu_/h *(10 + 10*gammaK);


  for(int ipq = 0; ipq < qfb.n; ++ipq)  {             // loop over integration points

    QuadraturePoint ip(qfb[ipq]);
    const Rd mip = face((RdHatBord)ip);

    const R Cint = mes * ip.a * In.T.mesure()* tq.a ;
    FK.BF(FK.T.toKref(mip), w);            // compute basic funtions

    DUn = 0;
    for(int ci=0; ci<Rd::d;++ci) {
      for(int cj=0; cj<Rd::d;++cj) {
        for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
          for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
            DUn(i,j) +=  ((ci==cj)? MyOp::Dxn(w,normal,ci,i,j)
            : w(i,ci,op_id) * w(j,cj,MyOp::D(ci))*normal[cj]);
          }
        }
      }
    }

    //ah = -(mu E(u).n, v) - (u, muE(v).n) + lambda (u,v)
    for(int it=0; it<nbDoFTime; ++it) {
      for(int jt=0; jt<nbDoFTime; ++jt) {
        const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
        for(int ci=0;ci<Rd::d;++ci) {
          for(int cj=0;cj<Rd::d;++cj) {
            for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
              for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){

                ML(i+it*nbDoF, j+jt*nbDoF) +=
                Cint * ((-1) * mu_ * ( DUn(j,i) + DUn(i,j) ) +
                (ci==cj) * (lambdaG * w(i,ci,op_id) * w(j,cj,op_id)) )
                *tval;

              }
            }
          }
        }
      }
    }



    // ah = (p, v.n ) - (u.n, q )  // from divergence form
    for(int it=0; it<nbDoFTime; ++it) {
      for(int jt=0; jt<nbDoFTime; ++jt) {
        const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
        for(int ci=0;ci<Rd::d;++ci) {
          for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
            for(int j = FK.dfcbegin(cp); j < FK.dfcend(cp); ++j){

              ML(i+it*nbDoF, j+jt*nbDoF) +=
              Cint * w(i,ci,op_id) * normal[ci] * w(j,cp,op_id) * tval;
              ML(j+jt*nbDoF, i+it*nbDoF) -=
              Cint * w(i,ci,op_id) * normal[ci] * w(j,cp,op_id) * tval;

            }
          }
        }
      }
    }

    // lh = (g, lambda v - muE(v).n)
    for(int it=0; it<nbDoFTime; ++it) {
      const R tval = wt(it,0,op_id);
      for(int ci=0;ci<Rd::d;++ci) {
        for(int cj=0;cj<Rd::d;++cj) {
          for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
            for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
              VL(i+it*nbDoF) += Cint * gxj(j) * tval *
              (lambdaG * w(j,cj,op_id) * w(i,ci,op_id) * (ci == cj)
              - mu_ * DUn(j,i));
            }
          }
        }
      }
    }

    // lh = -(g, nq)
    for(int jt=0; jt<nbDoFTime; ++jt) {
      const R tval = wt(jt,0,op_id);
      for(int ci=0;ci<Rd::d;++ci) {
        for(int j = FK.dfcbegin(cp); j < FK.dfcend(cp); ++j){
          for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
            VL(j+jt*nbDoF) -=
            Cint * gxj(i) * w(i,ci,op_id)*w(j,cp,op_id) * normal[ci] * tval;
          }
        }
      }
    }

  }// end loop quadrature point
}

void locNonLinAssembly(const FElement& FK, const FElement& FKs, const int iface, Rn& FL, Rnm& NLM, const Rn& w0) {

  assert(problem.sigma);

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  const typename Interface::FaceIdx& face = interface[iface];  // the face
  const R mes = interface.computeDx(face).norm();    // integration constant
  Rd normal(interface.normal(iface));
  normal /= normal.norm();

  const int nbDoF = FK.NbDoF();
  const int nbDoFs = FKs.NbDoF();
  const int nbDoFTime = In.NbDoF();
  const int cs = Rd::d+1;
  const int idxs0 = nbDoF*nbDoFTime;

  Rn valu0(Rd::d + 2, 0.);  // all componants
  Rnm Dvalu0(Rd::d + 2, Rd::d);  // all componants

  FL = 0.0;  NLM =0.;
  for(int ipq = 0; ipq < qfb.n; ++ipq)  {             // loop over integration points

    QuadraturePoint ip(qfb[ipq]);
    const Rd mip = interface.mapToFace(face,(RdHatBord)ip);

    const R Cint = mes * ip.a * In.T.mesure()* tq.a ;
    FK.BF(FK.T.toKref(mip), w);
    FKs.BF(FKs.T.toKref(mip), ws);

    valu0 = 0.; Dvalu0 = 0.;
    // compute u0(xq, tq)
    for(int it=0; it<nbDoFTime; ++it) {
      const R tval = wt(it,0,op_id);
      for(int ci=0;ci<Rd::d;++ci) { // dont need the pressure
        for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
          valu0(ci) += w0(i+it*nbDoF) * w(i,ci,op_id) * tval;
        }
      }
      for(int i = FKs.dfcbegin(0); i < FKs.dfcend(0); ++i){
        valu0(cs) += w0(idxs0+i+it*nbDoFs) * ws(i,0,op_id) * tval;
      }
    }
    // compute gradu0(xq, tq)
    for(int it=0; it<nbDoFTime; ++it) {
      const R tval = wt(it,0,op_id);
      for(int ci=0;ci<Rd::d;++ci) { // dont need the pressure
        for(int cj=0;cj<Rd::d;++cj) {
          for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
            Dvalu0(ci, cj) += w0(i+it*nbDoF) * w(i,ci,MyOp::D(cj)) * tval;
          }
        }
      }
      for(int cj=0;cj<Rd::d;++cj) {
        for(int i = FKs.dfcbegin(0); i < FKs.dfcend(0); ++i){
          Dvalu0(cs, cj) += w0(idxs0+i+it*nbDoFs) * ws(i,0,MyOp::D(cj)) * tval;
        }
      }
    }


    // (grad w0.u0 + w0 gradG.u0, r)
    R divGu0 = MyOp::divergenceSurface(Dvalu0, normal);
    R uGradw = 0;
    for(int ci=0;ci<Rd::d;++ci) uGradw += Dvalu0(cs , ci)*valu0(ci);

    for(int jt=0; jt<nbDoFTime; ++jt) {
      for(int j = FKs.dfcbegin(0); j < FKs.dfcend(0); ++j){
        FL(idxs0+j+jt*nbDoFs) += Cint
        * (
          uGradw + valu0(cs) * divGu0
        ) * ws(j, 0,op_id) * wt(jt,0,op_id);
      }
    }


    // (gradw0 . u + w0 divS u, r)
    for(int it=0; it<nbDoFTime; ++it) {
      for(int jt=0; jt<nbDoFTime; ++jt) {
        const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
        for(int cj=0;cj<Rd::d;++cj) {
          for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
            Rd gradS = MyOp::gradSurface(w, normal, j);
            for(int i = FKs.dfcbegin(0); i < FKs.dfcend(0); ++i){
              NLM(idxs0+i+it*nbDoFs, j+jt*nbDoF) += Cint *

              (
                Dvalu0(cs, cj) * w(j,cj,op_id)
                + valu0(cs) * gradS[cj]
              ) * ws(i,0, op_id) * tval;
            }
          }
        }
      }
    }

    // (gradw . u0 + w divS u0, r)
    for(int it=0; it<nbDoFTime; ++it) {
      for(int jt=0; jt<nbDoFTime; ++jt) {
        const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
        for(int j = FKs.dfcbegin(0); j < FKs.dfcend(0); ++j){
          R gradwu0 = 0;
          for(int cj=0;cj<Rd::d;++cj) gradwu0 += ws(j, 0,MyOp::D(cj)) * valu0(cj);

          for(int i = FKs.dfcbegin(0); i < FKs.dfcend(0); ++i){
            NLM(idxs0+i+it*nbDoFs, idxs0+j+jt*nbDoFs) += Cint *
            (
              gradwu0 + ws(j, 0,op_id) * divGu0
            ) * ws(i, 0,op_id) * tval;
          }
        }
      }
    }
  }


}





            // return the index of the neighbor element : -1 if none
            int locStabilization(KNM<R>&  ML, KNM<R>&  ML_M, KNM<R>&  MLN,
              const FElement& FK, const int ifac) {
                typedef typename QFB::QuadraturePoint QuadraturePoint;

                const int nbDoFTime = In.NbDoF();
                const int nbDoF = FK.NbDoF();

                ML = 0.0; ML_M = 0.0; MLN = 0.0;
                const R h = Vh[0].T.lenEdge(0);
                int the_domain = FK.whichDomain();
                int k_back = Vh.Th(FK.T);
                int ifacn = ifac;
                int kn_back = Vh.Th.ElementAdj(k_back,ifacn);
                if(kn_back == -1) return -1;                                 // no neighboor

                int kn = Vh.idxElementFromBackMesh(kn_back, the_domain);   // not in the domain
                if(kn ==-1) return -1;

                const FElement & FKn(Vh[kn]);                     // the neighboor finite elemen
                Rd normal = (FK.number < kn)? FK.T.N(ifac) : FKn.T.N(ifacn);
                const R mes = FK.T.mesureBord(ifac);


                Rn Cst(Rd::d+1);
                Cst = h*Edu*mes*problem.mu(the_domain) * tq.a *In.T.mesure();
                Cst(Rd::d) = h*h*h*Edp*mes*1./(problem.mu(the_domain)) * tq.a *In.T.mesure();;


                for(int ipq = 0; ipq < qfb.n; ++ipq)  {            // loop over integration points

                  QuadraturePoint ip(qfb[ipq]);
                  const Rd mip = FK.T(FK.T.toKref((RdHatBord)ip, ifac));
                  FK.BF (FK.T.toKref(mip) , w);
                  FKn.BF(FKn.T.toKref(mip), wn);

                  Rn Cint(Cst); Cint *= ip.a;


                  for(int it=0; it<nbDoFTime; ++it) {
                    for(int jt=0; jt<nbDoFTime; ++jt) {
                      const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
                      for(int c=0; c<=Rd::d;++c) {
                        for(int j = FK.dfcbegin(c); j < FK.dfcend(c); ++j){
                          for(int i = FK.dfcbegin(c); i < FK.dfcend(c); ++i){
                            ML(i+it*nbDoF,j+jt*nbDoF) += Cint(c) * MyOp::innerProduct(w, normal,c,i)
                            * MyOp::innerProduct(w, normal,c,j) * tval;
                            ML_M(i+it*nbDoF,j+jt*nbDoF) += Cint(c) * MyOp::innerProduct(w, normal,c,i)
                            * MyOp::innerProduct(wn, normal,c,j) * tval;
                            MLN(i+it*nbDoF,j+jt*nbDoF)  += Cint(c) * MyOp::innerProduct(wn, normal,c,i)
                            * MyOp::innerProduct(wn, normal,c,j) * tval;
                          }
                        }
                      }
                    }
                  }
                  // // Need to add whatd and init correctly w and wn
                  //   // R Cint2 = h*h*h*Edu*mes*problem.mu(the_domain)*ip.a;
                  //   // for(int c=0; c<Rd::d;++c) {
                  //   // 	for(int j = FK.dfcbegin(c); j < FK.dfcend(c); ++j){
                  //   // 	  for(int i = FK.dfcbegin(c); i < FK.dfcend(c); ++i){
                  //   // 	    ML(i,j)   += Cint2 * MyOp::normal2Der(w, normal,i)
                  //   // 	      * MyOp::normal2Der(w, normal,j);
                  //   // 	    ML_M(i,j) += Cint2 * MyOp::normal2Der(w, normal,i)
                  //   // 	      * MyOp::normal2Der(wn, normal,j);

                  //   // 	  }
                  //   // 	}
                  //   // }
                } // end loop over quad point in space
                return kn;
              }

              int locFaceStabilization(KNM<R>&  ML, KNM<R>&  ML_M,
                const FElement& FK, const int ifac) {


                  typedef typename QFB::QuadraturePoint QuadraturePoint;

                  const int nbDoFTime = In.NbDoF();
                  const int nbDoF = FK.NbDoF();

                  const R h = Sh[0].T.lenEdge(0);
                  int ifacn = ifac;
                  int kn = Sh.Th.ElementAdj(FK.number,ifacn);
                  if(kn == -1) return -1;                            // no neighboor

                  const FElement & FKn(Sh[kn]);                     // the neighboor finite element
                  Rd normal = (FK.number < kn)? FK.T.N(ifac) : FKn.T.N(ifacn);
                  const R mes = FK.T.mesureBord(ifac);
                  ML=0; ML_M=0;

                  for(int ipq = 0; ipq < qfb.n; ++ipq)  {            // loop over integration points

                    QuadraturePoint ip(qfb[ipq]);
                    const Rd mip = FK.T(FK.T.toKref((RdHatBord)ip, ifac));

                    const double Cint = mes * ip.a * Eds * h * In.T.mesure() * tq.a  ;

                    FK.BF (FK.T.toKref(mip) , ws);
                    FKn.BF(FKn.T.toKref(mip), wn);

                    for(int it=0; it<nbDoFTime; ++it) {
                      for(int jt=0; jt<nbDoFTime; ++jt) {
                        const R tval = wt(it,0,op_id)*wt(jt,0,op_id);
                        for(int j = FK.dfcbegin(0); j < FK.dfcend(0); ++j){
                          for(int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i){

                            ML(i+it*nbDoF,j+jt*nbDoF) += Cint * MyOp::innerProduct(ws, normal,i)
                            * MyOp::innerProduct(ws, normal,j) * tval;
                            ML_M(i+it*nbDoF,j+jt*nbDoF) += Cint * MyOp::innerProduct(ws, normal,i)
                            * MyOp::innerProduct(wn, normal,j) * tval;

                          }
                        }
                      }
                    }

                  }
                  return kn;
                }


                void locL2norm(const FElement& FK, const Partition& cutK,
                  R& errU,  Rd (*fv)(const Rd, const R, const int),
                  R& errP,  R  (*fp)(const Rd, const R, const int),
                  const Rn& w0) {

                    typedef typename QF::QuadraturePoint QuadraturePoint;

                    const int domain = FK.whichDomain();
                    ElementSignEnum the_part = cutK.what_part(domain);
                    const int nbDoF = FK.NbDoF();
                    const int nbDoFTime = In.NbDoF();
                    const int cp  = Rd::d;           // the pressure componante

                    errU = errP = 0;
                    for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
                    it != cutK.element_end(the_part); ++it){

                      const R mes = cutK.mesure(it);
                      for(int ipq = 0; ipq < qf.n; ++ipq)  {             // loop over integration points

                        QuadraturePoint ip(qf[ipq]);
                        Rd mip = cutK.toK(it, ip);
                        const R Cint = mes * ip.a;

                        FK.BF(whatd, FK.T.toKref(mip), w);

                        Rd uh = fv(mip, tq.x, domain);
                        for(int ci=0; ci<Rd::d;++ci) {
                          for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
                            for(int it=0; it<nbDoFTime; ++it) {
                              uh[ci] -= w0(i + it*nbDoF) * w(i,ci, op_id) * wt(it,0,op_id) ;
                            }
                          }
                        }
                        R ph= fp(mip, tq.x, domain);
                        for(int i = FK.dfcbegin(cp); i < FK.dfcend(cp); ++i){
                          for(int it=0; it<nbDoFTime; ++it) {
                            ph -= w0(i + it*nbDoF) * w(i,cp, op_id) * wt(it,0,op_id);
                          }
                        }

                        errU += Cint* (uh, uh);
                        errP += Cint* ph*ph;

                      }

                    }
                  }

                  void locL2norm(const FElement& FK, const int iface,
                    R& errS,  R  (*fs)(const Rd, const R),
                    const Rn& w0) {

                      typedef typename QFB::QuadraturePoint QuadraturePoint;
                      const typename Interface::FaceIdx& face = interface[iface];  // the face
                      const R mes = interface.computeDx(face).norm();    // integration constant

                      const int nbDoF = FK.NbDoF();
                      const int nbDoFTime = In.NbDoF();

                      errS = 0;

                      for(int ipq = 0; ipq < qfb.n; ++ipq)  {             // loop over integration points

                        QuadraturePoint ip(qfb[ipq]);
                        Rd mip = interface.mapToFace(face,(RdHatBord)ip);
                        const R Cint = mes * ip.a;

                        FK.BF(FK.T.toKref(mip), w);

                        R sh= fs(mip, tq.x);
                        for(int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i){
                          for(int it=0; it<nbDoFTime; ++it) {
                            sh -= w0(i + it*nbDoF) * w(i,0, op_id) * wt(it,0,op_id);
                          }
                        }
                        errS += Cint* sh*sh;

                      }
                    }




                  };





                  #endif
