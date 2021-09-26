#ifndef INTEGRAL_ICESHEET_HPP_
#define INTEGRAL_ICESHEET_HPP_


#include "iceSheet.hpp"

class IntegralIceSheet {//}: public Integration {

  typedef Mesh2 Mesh;
  typedef FESpace2 FESpace;
  typedef CutFESpace2 CutFESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::Rd Rd;
  typedef typename FElement::RdHat RdHat;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::BorderElement BorderElement;
  typedef typename Mesh::Partition Partition;
  typedef GenericInterface<Mesh2> Interface;
  typedef Fun2 Fun;

  typedef typename TypeOperationMV<Rd::d>::OperationMV MyOp;
  typedef typename FElement::QFB QFB;
  typedef typename FElement::QF QF;


  enum labelBoundary { surface = 1, floating = 2, bedRock  = 3, symmetry = 4, cliff = 5};

  const IceSheet& problem;
  const CutFESpace& Vh;
  const QFB &qfb;
  const QF &qf;
  const Marker& marker;

  int nbDoF;
  KNMK<double> w, wn;
  What_d whatd;
  R Edu = 1e-2, Edp = 1e-2;



public:
  R vol = 0;
  IntegralIceSheet(const IceSheet& pp) :
  // Integration((*pp.Vh)[0].NbDoF(),
  // pp.Vh->NbDoF(),
  // 2),
  Vh(*pp.Vh),
  problem(pp),
  qfb(*QF_Simplex<RdHatBord>(5)),
  qf(*QF_Simplex<RdHat>(5)),
  marker(*(Vh.marker[0]))
  {
    nbDoF = Vh[0].NbDoF();
    w.init(nbDoF,Rd::d+1,op_dz+1);
    whatd = Fop_D0|Fop_D1;

  }

  void initStab(R edu, R edp) {
    Edu = edu; Edp=edp;
    w.init(nbDoF,Rd::d,op_Dall);
    wn.init(nbDoF,Rd::d,op_Dall);
    whatd = Fop_D0|Fop_D1;//|Fop_D2;
  }



  void locFullAssembly(const FElement& FK, Rnm& ML, Rn& VL, Rn& VP, const Partition& cutK) {

    typedef typename QF::QuadraturePoint QuadraturePoint;
    ML = 0.0; VL = 0.0; VP = 0.0;

    const int domain = 0;//FK.whichDomain();
    ElementSignEnum the_part = cutK.what_part(domain);
    const R mu_ = problem.mu();
    const R rho_ = problem.rho(Ice);
    const int nbDoF = FK.NbDoF();
    const int cp  = Rd::d;           // the pressure componante

    for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
    it != cutK.element_end(the_part); ++it){

      const R mes = cutK.mesure(it);
      for(int ipq = 0; ipq < qf.n; ++ipq)  {

        QuadraturePoint ip(qf[ipq]);
        Rd mip = cutK.toK(it, ip);
        const R Cint = mes * ip.a;
        vol += Cint;
        FK.BF(FK.T.toKref(mip), w);

        // (E(u),E(v))
        for(int ci=0; ci<Rd::d;++ci) {
          for(int cj=0; cj<Rd::d;++cj) {
            for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
              for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){

                ML(i,j) +=  Cint* mu_ *
                ((ci==cj)? MyOp::crossTerm(w,ci,i,j)
                : w(i,ci,MyOp::D(cj)) * w(j,cj,MyOp::D(ci)));
              }
            }
          }
        }

        // -(p, div v) + (div u, q)
        for(int ci=0; ci<Rd::d;++ci) {
          for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
            for(int jp=FK.dfcbegin(cp); jp<FK.dfcend(cp); ++jp){
              ML(i, jp) -= Cint * w(jp,cp,op_id) * w(i,ci,MyOp::D(ci));
              ML(jp, i) += Cint * w(jp,cp,op_id) * w(i,ci,MyOp::D(ci));
            }
          }
        }

        // 	// (f,v)
        for(int ci=0; ci<Rd::d;++ci) {
          for(int j = FK.dfcbegin(ci); j < FK.dfcend(ci); ++j){
            // for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
            VL(j) +=
            Cint  * rho_ * w(j,ci,op_id) * problem.frhs(mip)[ci];//w(j,op_id) * fxj(i);
          }
        }
        //   	}

        //(p, 1)
        for(int jp=FK.dfcbegin(cp); jp<FK.dfcend(cp); ++jp){
          VP(jp) += Cint * w(jp,cp,op_id);
        }

      } // loop quad point
    }

  }



  void locSurfAssembly(const FElement& FK, const int iface, Rnm& ML, Rn& VL) {

    ML = 0.0; VL = 0.0;
    typedef typename QFB::QuadraturePoint QuadraturePoint;
    const FaceMarker& face = marker.getFace(iface);  // the face
    const R mes = face.mesure();
    Rd normal(-face.normal());
    Rd tangente = normal.perp();

    KNM<R> Sn(nbDoF, Rd::d);
    KN<R> Snn(nbDoF), Stn(nbDoF);
    const R mu_ = problem.mu();
    // const int jumpSign[2] = {1, -1};

    const R h    = Max(Vh[0].T.lenEdge(0),Vh[0].T.lenEdge(1),Vh[0].T.lenEdge(2)) ;
    const R invh = 1./h;
    const R invh2 = invh * invh;
    const R gamma = mes * invh;
    const R Beta = problem.C;

    const int cp = Rd::d;

    for(int ipq = 0; ipq < qfb.n; ++ipq)  {             // loop over integration points

      QuadraturePoint ip(qfb[ipq]);
      const Rd mip = face.mapToFace( ip);

      const R Cint = mes * ip.a;
      FK.BF(FK.T.toKref(mip), w);

      if(  face.lab == symmetry || face.lab == bedRock ){

// this is wrong
        for(int ci=0; ci<Rd::d;++ci) {
          for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
            Snn(i) = 2*mu_*(w(i,ci,MyOp::D(0))*normal[0]
                          + w(i,ci,MyOp::D(1))*normal[1])*normal[ci];
            Stn(i) = 2*mu_*(w(i,ci,MyOp::D(0))*normal[0]
                          + w(i,ci,MyOp::D(1))*normal[1])*tangente[ci];

          }
        }
        for(int i = FK.dfcbegin(cp); i < FK.dfcend(cp); ++i){
          Snn(i) = -w(i,cp,op_id);
          Stn(i) = -w(i,cp,op_id);
        }
      }


      // Flotting lower boundary or Front boundary {x=1.6e8}
      //       (pw.n, v)_{\Gamma_{bf}}{cliff}
      if(face.lab == floating || face.lab == cliff) {

        for(int ci=0;ci<Rd::d;++ci) {
          for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
            VL(i)  -= Cint * w(i,ci,op_id) * problem.Pwater(mip) * normal[ci];
          }
        }
      }

      if(  face.lab == symmetry){

        //-(n(sigma.n), v.n) - (u.n, n(sigma.n) )
        for(int ci=0;ci<=Rd::d;++ci) {
          for(int cj=0;cj<Rd::d;++cj) {
            for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
              for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
                ML(i, j) -= Cint * (Snn(i)*w(j,cj,op_id)*normal[cj]);
                ML(j, i) -= Cint * (Snn(i)*w(j,cj,op_id)*normal[cj]
                                  + Stn(i)*w(j,cj,op_id)*tangente[cj]);
              }
            }
          }
        }
        //(B.u,v) + c/h(u.n,v.n)
        for(int ci=0;ci<Rd::d;++ci) {
          for(int j = FK.dfcbegin(ci); j < FK.dfcend(ci); ++j){
            for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
              ML(i, j) += Cint * (Beta*w(i,ci,op_id)*w(j,ci,op_id) +
              1e3*invh*w(i,ci,op_id)*normal[ci]*w(j,ci,op_id)*normal[ci]);
            }
          }
        }


      }


      // BedRock lower boundary
      if(face.lab == bedRock ){//}|| face.lab == symmetry){// || face.lab == floating || face.lab == cliff){

        //-(n(sigma.n), v.n) - (u.n, n(sigma.n) )
        for(int ci=0;ci<=Rd::d;++ci) {
          for(int cj=0;cj<Rd::d;++cj) {
            for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
              for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
                ML(i, j) -= Cint * (Snn(i)*w(j,cj,op_id)*normal[cj]);
                ML(j, i) -= Cint * (Snn(i)*w(j,cj,op_id)*normal[cj]);
              }
            }
          }
        }
        //(B.u,v) + c/h(u.n,v.n)
        for(int ci=0;ci<Rd::d;++ci) {
          for(int cj=0;cj<Rd::d;++cj) {
            for(int j = FK.dfcbegin(cj); j < FK.dfcend(cj); ++j){
              for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
                ML(i, j) += Cint * (Beta*w(i,ci,op_id)*w(j,cj,op_id)*(ci==cj)
                +1e3*invh*invh*w(i,ci,op_id)*normal[ci]*w(j,cj,op_id)*normal[cj]
              );
            }
          }
        }
      }
    }// end if bedRock

    }

  }




  // return the index of the neighbor element : -1 if none
  int locStabilization(KNM<R>&  ML, KNM<R>&  ML_M, KNM<R>&  MLN, const FElement& FK, const int ifac) {
    typedef typename QFB::QuadraturePoint QuadraturePoint;

    ML = 0.0; ML_M = 0.0; MLN = 0.0;
    const R h    = Max(Vh[0].T.lenEdge(0),Vh[0].T.lenEdge(1),Vh[0].T.lenEdge(2)) ;    int the_domain = 1;//FK.whichDomain();
    int k_back = Vh.Th(FK.T);
    int ifacn = ifac;
    int kn_back = Vh.Th.ElementAdj(k_back,ifacn);
    if(kn_back == -1) return -1;                                 // no neighboor

    int kn = Vh.idxElementFromBackMesh(kn_back, 0);//the_domain);   // not in the domain
    if(kn ==-1) return -1;

    const FElement & FKn(Vh[kn]);                     // the neighboor finite elemen
    Rd normal = (FK.number < kn)? FK.T.N(ifac) : FKn.T.N(ifacn);
    const R mes = FK.T.mesureBord(ifac);

    int nbNode = FK.NbNode();
    R Dkn[nbNode], DknN[nbNode];

    Rn Cst(Rd::d+1);
    Cst = h*Edu*mes*problem.mu();
    Cst(Rd::d) = h*h*h*Edp*mes*1./(problem.mu());

    for(int ipq = 0; ipq < qfb.n; ++ipq)  {            // loop over integration points

      QuadraturePoint ip(qfb[ipq]);
      const Rd mip = FK.T(FK.T.toKref((RdHatBord)ip, ifac));
      FK.BF (whatd, FK.T.toKref(mip) , w);
      FKn.BF(whatd, FKn.T.toKref(mip), wn);

      Rn Cint(Cst); Cint *= ip.a;

      for(int c=0; c<=Rd::d;++c) {
        for(int j = FK.dfcbegin(c); j < FK.dfcend(c); ++j){
          for(int i = FK.dfcbegin(c); i < FK.dfcend(c); ++i){
            ML(i,j) += Cint(c) * MyOp::innerProduct(w, normal,i)
            * MyOp::innerProduct(w, normal,j);
            ML_M(i,j) += Cint(c) * MyOp::innerProduct(w, normal,i)
            * MyOp::innerProduct(wn, normal,j);
            MLN(i,j)  += Cint(c) * MyOp::innerProduct(wn, normal,i)
            * MyOp::innerProduct(wn, normal,j);
          }
        }
      }

      // R Cint2 = h*h*h*Edu*mes*Stokes.mu(the_domain)*ip.a;
      // for(int c=0; c<Rd::d;++c) {
      // 	for(int j = FK.dfcbegin(c); j < FK.dfcend(c); ++j){
      // 	  for(int i = FK.dfcbegin(c); i < FK.dfcend(c); ++i){
      // 	    ML(i,j)   += Cint2 * MyOp::normal2Der(w, normal,i)
      // 	      * MyOp::normal2Der(w, normal,j);
      // 	    ML_M(i,j) += Cint2 * MyOp::normal2Der(w, normal,i)
      // 	      * MyOp::normal2Der(wn, normal,j);

      // 	  }
      // 	}
      // }
    } // end loop over quad point in space
    return kn;
  }




  //   void locL2norm(const FElement& FK, const Partition& cutK,
  // 		 R& errU,  Rd (*fv)(const Rd, int),
  // 		 R& errP,  R  (*fp)(const Rd, int) ) {

  //     typedef typename QF::QuadraturePoint QuadraturePoint;

  //     const int domain = FK.whichDomain();
  //     ElementSignEnum the_part = cutK.what_part(domain);
  //     const int nbDoF = FK.NbDoF();
  //     const int cp  = Rd::d;           // the pressure componante
  //     errU = errP = 0;
  //     for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
  //   	it != cutK.element_end(the_part); ++it){

  //       const R mes = cutK.mesure(it);
  //       for(int ipq = 0; ipq < qf.n; ++ipq)  {             // loop over integration points

  //       	QuadraturePoint ip(qf[ipq]);
  //       	Rd mip = cutK.toK(it, ip);
  //       	const R Cint = mes * ip.a;

  //       	FK.BF(whatd, FK.T.toKref(mip), w);

  // 	Rd uh = fv(mip, domain);
  //       	for(int ci=0; ci<Rd::d;++ci) {
  // 	  for(int i = FK.dfcbegin(ci); i < FK.dfcend(ci); ++i){
  // 	    uh[ci] -= Stokes.rhs(FK(i)) * w(i, op_id);
  // 	  }
  // 	}
  // 	R ph= fp(mip, domain);
  // 	for(int i = FK.dfcbegin(cp); i < FK.dfcend(cp); ++i){
  // 	  ph -= Stokes.rhs(FK(i)) * w(i, op_id);
  // 	}

  // 	errU += Cint* (uh, uh);
  // 	errP += Cint* ph*ph;

  //       }

  //     }
  //   }




};





#endif
