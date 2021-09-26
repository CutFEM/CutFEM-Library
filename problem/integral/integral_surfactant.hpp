#ifndef INTEGRAL_SURFACTANT_HPP_
#define INTEGRAL_SURFACTANT_HPP_



template<class M>
class IntegralSurfactant : public TimeIntegration {

  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::Rd Rd;
  typedef typename FElement::RdHatBord RdHatBord;

  typedef GFESpace<Mesh1> TimeSpace;
  typedef typename TimeSpace::FElement TimeSlab;
  typedef GenericInterface<Mesh> Interface;

  typedef typename TypeOperationMV<Rd::d>::OperationMV MyOp;
  typedef typename FElement::QFB QFB;
  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef GQuadraturePoint<R1> QPT;

  typedef FEMVecFunction<FESpace> FunVec ;

  R (*funRHS)(const Rd, const R) = nullptr;

  const QFB &qfb;
  const TimeSlab& In;
  const GenericInterface<Mesh>& interface;
  const FunVec& mapping;
  const QPT tq;

public:

  Rd (*funVELOCITY)(const Rd, const R) = nullptr;
  R (*funDIV)(const Rd, const R) = nullptr;

  IntegralSurfactant(const FESpace& Vh, const Interface& gamma,
		     const TimeSlab& IIn, const QPT ttq,
		     R (*funnRHS)(const Rd, const R),
		     const FunVec& mapp = DataFunFEM<M>::Id
		     ) :
    TimeIntegration(Vh[0].nbDoF(), IIn.nbDoF(),
		Vh.global.NbOfNodes,
		mapp.degre(), IIn.degre()
		),
    qfb(*QF_Simplex<RdHatBord>(15)),//2*deg+1)),
    In(IIn),
    interface(gamma),
    mapping(mapp),
    tq(ttq),
    funRHS(funnRHS)
  {

    w.init(nbDoF,1,op_dz+1);
    wt.init(nbDoFTime,1,op_dz+1);
    whatd = Fop_D0|Fop_D1;

    In.BF(whatd, tq.x, wt);                // compute basic funtions

  };


  void locAssembly(KNM<R>&  ML, KN<R>& VL, KN<R>& LAG,
    const FElement& FK, const int iface,
    const FEMVecFunction<FESpace>& B, const double eps) {

      const typename Interface::FaceIdx& face = interface[iface];  // the face
      const Rd dx = interface.computeDx(face);            // integration constant
      const Rd linear_normal(interface.normal(iface));
      const int kb = interface.idxElementOfFace(iface);  // idx on backMesh
      const double t = In.Pt(R1(tq.x)).x;


      for(int ipq = 0; ipq < qfb.n; ++ipq)  {

        QuadraturePoint ip(qfb[ipq]);
        const Rd mip = interface.mapToFace(face,(RdHatBord)ip);
        const Rd mappip = mapping(FK.number, mip);

        KNM<R> jac = mapping.jacobian(FK.number, mip);
        const KNM<R> invJ(MyOp::inverse(jac).t());
        Rd normal = MyOp::multiply(invJ,linear_normal);

        const R Cint = dx.norm() * In.T.mesure()* tq.a * ip.a
        * MyOp::determinant(jac)*normal.norm();
        const R Cint0 = dx.norm() * ip.a
        * MyOp::determinant(jac)*normal.norm();

        normal /= normal.norm();

        const Rd Beta = (funVELOCITY)? funVELOCITY(mappip, t) : B(kb, mappip);
        const R divB = (funDIV)? funDIV(mappip,t) : B.divergence(kb,mappip,linear_normal);


        FK.BF(whatd, FK.T.toKref(mip), w);                  // compute basic funtions
        KNMK<R> wDef(MyOp::transform_dBF(invJ, w));

        for(int it=0; it<nbDoFTime; ++it) {
          for(int jt=0; jt<nbDoFTime; ++jt) {
            for(int j = 0; j < nbDoF; ++j){
              for(int i = 0; i < nbDoF; ++i){

                ML(i+it*nbDoF, j+jt*nbDoF) += Cint * (w(i,0,op_id)*wt(it,0,op_id))
                * (w(j,0,op_id)*wt(jt,0,op_dx))			
                + Cint * (MyOp::innerProduct(wDef, Beta, j)
                + w(j,0,op_id)*divB) * w(i,0,op_id)
                * wt(it,0,op_id) * wt(jt,0,op_id)
                + Cint * eps * (MyOp::gradSurface(wDef,normal,i),
                MyOp::gradSurface(wDef,normal,j))
                * wt(it,0,op_id) * wt(jt,0,op_id)
                ;
              }
            }
          }
        }

        for(int it=0; it<nbDoFTime; ++it) {
          for(int i = 0; i < nbDoF; ++i){
            VL(i + it*nbDoF) += Cint * w(i,0,op_id)*wt(it,0,op_id) * funRHS(mappip, t);
            LAG(i + it*nbDoF) +=  Cint0  * w(i,0,op_id) * wt(it,0,op_id);

          }
        }
      }
    }


  R initialCondition(KNM<R>&  ML, KN<R>& VL,
		    const FElement& FK, const int iface,
		    FEMFunction<FESpace>& initSol) {

    const int k = FK.number;

    const typename Interface::FaceIdx& face = interface[iface];  // the face
    const Rd linear_normal(interface.normal(iface));
    const Rd dx = interface.computeDx(face);             // integration constant
    const R len = dx.norm();
    R sum = 0;
    for(int ipq = 0; ipq < qfb.n; ++ipq)  {

      QuadraturePoint ip(qfb[ipq]);
      const Rd mip = interface.mapToFace(face,(RdHatBord)ip);
      const Rd mappip = mapping(FK.number, mip);


      KNM<R> jac = mapping.jacobian(FK.number, mip);
      const KNM<R> invJ(MyOp::inverse(jac).t());
      const Rd normal = MyOp::multiply(invJ,linear_normal);

      const R Cint = len * ip.a * MyOp::determinant(jac)*normal.norm();
      sum+=Cint;
      FK.BF(whatd, FK.T.toKref(mip), w);                // compute basic funtions

      for(int j = 0; j < nbDoF; ++j){
  	for(int i = 0; i < nbDoF; ++i){
  	  ML(i, j) += Cint * w(i,0,op_id) * w(j,0,op_id);
  	}
      }

      for(int i = 0; i < nbDoF; ++i){
  	VL(i) +=  Cint  * w(i,0,op_id)* initSol(k, mip);
      }
    }
    return sum;
  }





  /*
     Integral for other formulation
  */

  void locAssembly2(KNM<R>&  ML, KN<R>& VL, KN<R>& LAG,
		   const FElement& FK, const int iface,
		   const FEMVecFunction<FESpace>& B, const double eps) {

    const typename Interface::FaceIdx& face = interface[iface];  // the face
    const Rd dx = interface.computeDx(face);            // integration constant
    const Rd linear_normal(interface.normal(iface));
    const int kb = interface.idxElementOfFace(iface);  // idx on backMesh
    const double t = In.Pt(R1(tq.x)).x;

    for(int ipq = 0; ipq < qfb.n; ++ipq)  {

      QuadraturePoint ip(qfb[ipq]);
      const Rd mip = interface.mapToFace(face,(RdHatBord)ip);
      const Rd mappip = mapping(FK.number, mip);

      KNM<R> jac = mapping.jacobian(FK.number, mip);
      const KNM<R> invJ(MyOp::inverse(jac).t());
      Rd normal = MyOp::multiply(invJ,linear_normal);

      const R Cint = dx.norm() * In.T.mesure()* tq.a * ip.a
	* MyOp::determinant(jac)*normal.norm();


      normal /= normal.norm();

      const Rd Beta = (funVELOCITY)? funVELOCITY(mappip, t) : B(kb, mappip);
      // const R divB = (funDIV)? funDIV(mappip,t) : B.divergence(kb,mappip,linear_normal);

      FK.BF(whatd, FK.T.toKref(mip), w);                  // compute basic funtions
      KNMK<R> wDef(MyOp::transform_dBF(invJ, w));

      for(int it=0; it<nbDoFTime; ++it) {
	for(int jt=0; jt<nbDoFTime; ++jt) {
	  for(int j = 0; j < nbDoF; ++j){
	    for(int i = 0; i < nbDoF; ++i){
	      ML(i+it*nbDoF, j+jt*nbDoF) += Cint * (-1.0) * w(i,0,op_id)*wt(it,0,op_dx)
		* w(j,0,op_id)*wt(jt,0,op_id)
		- Cint * MyOp::innerProduct(wDef, Beta, i)*w(j,0,op_id)
			  * wt(it,0,op_id) * wt(jt,0,op_id)
		+ Cint * eps * (MyOp::gradSurface(wDef,normal,i),
				MyOp::gradSurface(wDef,normal,j))
		* wt(it,0,op_id) * wt(jt,0,op_id)
		;
	    }
	  }
	}
      }

      for(int it=0; it<nbDoFTime; ++it) {
      	for(int i = 0; i < nbDoF; ++i){
      	  VL(i + it*nbDoF) += Cint * w(i,0,op_id)*wt(it,0,op_id) * funRHS(mappip, t);
	  // LAG(i + it*nbDoF) +=  Cint0  * w(i,0,op_id) * wt(it,0,op_id);
      	}
      }
    }
  }


  R initialCondition2(KN<R>& VL,
		      const FElement& FK, const int iface,
		      FEMFunction<FESpace>& initSol) {

    const int k = FK.number;

    const typename Interface::FaceIdx& face = interface[iface];  // the face
    const Rd linear_normal(interface.normal(iface));
    const Rd dx = interface.computeDx(face);             // integration constant
    const R len = dx.norm();
    R sum = 0;
    for(int ipq = 0; ipq < qfb.n; ++ipq)  {

      QuadraturePoint ip(qfb[ipq]);
      const Rd mip = interface.mapToFace(face,(RdHatBord)ip);
      const Rd mappip = mapping(FK.number, mip);

      KNM<R> jac = mapping.jacobian(FK.number, mip);
      const KNM<R> invJ(MyOp::inverse(jac).t());
      const Rd normal = MyOp::multiply(invJ,linear_normal);

      const R Cint = len * ip.a * MyOp::determinant(jac)*normal.norm();
      sum+=Cint;
      FK.BF(whatd, FK.T.toKref(mip), w);                // compute basic funtions

      for(int i = 0; i < nbDoF; ++i){
  	VL(i) +=  Cint  * w(i,0,op_id)* initSol(k, mip);
      }
    }
    return sum;
  }

  R lastCondition2(KNM<R>&  ML,
		   const FElement& FK, const int iface) {
    const int k = FK.number;

    const typename Interface::FaceIdx& face = interface[iface];  // the face
    const Rd linear_normal(interface.normal(iface));
    const Rd dx = interface.computeDx(face);             // integration constant
    const R len = dx.norm();
    R sum = 0;
    for(int ipq = 0; ipq < qfb.n; ++ipq)  {

      QuadraturePoint ip(qfb[ipq]);
      const Rd mip = interface.mapToFace(face,(RdHatBord)ip);
      const Rd mappip = mapping(FK.number, mip);


      KNM<R> jac = mapping.jacobian(FK.number, mip);
      const KNM<R> invJ(MyOp::inverse(jac).t());
      const Rd normal = MyOp::multiply(invJ,linear_normal);

      const R Cint = len * ip.a * MyOp::determinant(jac)*normal.norm();
      sum+=Cint;
      FK.BF(whatd, FK.T.toKref(mip), w);                // compute basic funtions

      for(int it=0; it<nbDoFTime; ++it) {
      	for(int jt=0; jt<nbDoFTime; ++jt) {
	  for(int j = 0; j < nbDoF; ++j){
	    for(int i = 0; i < nbDoF; ++i){
	      ML(i+it*nbDoF, j+jt*nbDoF) += Cint * w(i,0,op_id) * w(j,0,op_id)
		* wt(it,0,op_id) * wt(jt,0,op_id);
	    }
	  }
	}
      }
    }
    return sum;
  }



};



#endif
