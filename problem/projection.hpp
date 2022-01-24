#ifndef PROJETION_HPP_
#define PROJETION_HPP_

#include "baseProblem.hpp"




/*
  apply L2 projection to a vector
- uh is the projected vector
 */
template<typename M>
void projection(FunFEM<M> & fh, FunFEM<M> & ph){
	typedef typename M::Rd Rd;
	typedef TestFunction<Rd::d> FunTest;
	typedef GFESpace<M> FESpace;
	const FESpace& Fh(*fh.Vh);
	const FESpace& Vh(*ph.Vh);
	assert(&Fh.Th == &Vh.Th);

	FEM<M> projection(Vh);

	FunTest u(Vh, Vh.N), v(Vh, Vh.N);
	projection.addBilinear(innerProduct(u,v));
	projection.addLinear(innerProduct(fh.expression(),v));
	projection.solve();
	ph.v = projection.rhs;
}

/*
  projection on FEM in Time
 */
template<typename M>
void projection(FunFEM<M> & fh, FunFEM<M> & ph, const TimeSlab& In){
	typedef typename M::Rd Rd;
	typedef TestFunction<Rd::d> FunTest;
	typedef GFESpace<M> FESpace;

	const FESpace& Fh(*fh.Vh);
	const FESpace& Vh(*ph.Vh);

	assert(&Fh.Th == &Vh.Th);

	FEM<M> projection(Vh, In.Vh);
	FunTest u(Vh, Vh.N), v(Vh, Vh.N);

	projection.addBilinear(innerProduct(u,v), In);
	projection.addLinear(innerProduct(fh.expression(),v), In);

	projection.solve();
	ph.v = projection.rhs;
}




/*
  apply L2 projection to a vector for cutProblem
- uh is the projected vector
 */
template<typename M>
void projection(const FunFEM<M> & fh, FunFEM<M> & ph, const GenericInterface<M>& inter){
	typedef typename M::Rd Rd;
	typedef TestFunction<Rd::d> FunTest;
	typedef GFESpace<M> FESpace;
	typedef GenericInterface<M> Interface;

	const FESpace& Fh(*fh.Vh);
	const FESpace& Vh(*ph.Vh);

	assert(&Fh.Th == &Vh.Th);

	FEM<M> projection(Vh);
	FunTest u(Vh, Vh.N), v(Vh, Vh.N);
	projection.addBilinear(innerProduct(u,v), inter);
	projection.addLinear(innerProduct(fh.expression(),v), inter);

// in case of time mesh, the interface does not go through all the elements at one time
	projection.addDiagonal(1e-14);
	// projection.addLagrangeMultiplier(innerProduct(1.,v), inter, 0);

// matlab::Export(projection.mat, "matProj.dat");
// matlab::Export(projection.rhs, "rhsProj.dat");
//
// getchar();

	projection.solve();
	KN_<double> uh(projection.rhs(SubArray(ph.size(), 0)));
	ph.v = uh;
}


/*
  apply L2 projection to a vector in case of cutFEM
- uh is the projected vector
 */
template<typename M>
void projection(const FunFEM<M> & fh, FunFEM<M> & ph, const GenericInterface<M>& inter, double valLagrange){
	typedef typename M::Rd Rd;
	typedef TestFunction<Rd::d> FunTest;
	typedef GFESpace<M> FESpace;
	typedef GenericInterface<M> Interface;

	const FESpace& Fh(*fh.Vh);
	const FESpace& Vh(*ph.Vh);

	assert(&Fh.Th == &Vh.Th);

	FEM<M> projection(Vh);
	FunTest u(Vh, Vh.N), v(Vh, Vh.N);
	projection.addBilinear(innerProduct(u,v), inter);
	projection.addLinear(innerProduct(fh.expression(),v), inter);

// in case of time mesh, the interface does not go through all the elements at one time
	projection.addDiagonal(1e-14);


	projection.addLagrangeMultiplier(innerProduct(1.,v), inter, valLagrange);


	projection.solve();
	KN_<double> uh(projection.rhs(SubArray(ph.size(), 0)));
	ph.v = uh;
}




/*
  apply L2 projection to a vector in case of Space time cutFEM
- uh is the projected vector
 */
template<typename M>
void projection(const FunFEM<M> & fh, FunFEM<M> & ph, const TimeSlab& In, const TimeInterface<M>& interface, double valLagrange){
	typedef typename M::Rd Rd;
	typedef TestFunction<Rd::d> FunTest;
	typedef GFESpace<M> FESpace;
	typedef GenericInterface<M> Interface;

	const FESpace& Fh(*fh.Vh);
	const FESpace& Vh(*ph.Vh);

	assert(&Fh.Th == &Vh.Th);
	int exact =  exactLobatto_nPt(interface.size());
	FEM<M> projection(Vh, In.Vh, 5, exact);
	std::cout << interface.size() << std::endl;
	std::cout << projection.qTime.n << std::endl;

	assert(interface.size() == projection.qTime.n);
	FunTest u(Vh, Vh.N), v(Vh, Vh.N);
	const CutFEM_Parameter& h(Parameter::h);
	Normal n;

	projection.addBilinear(innerProduct(u,v), interface, In);
	projection.addLinear(innerProduct(fh.expression(),v), interface, In);

	projection.addEdgeIntegral(
		innerProduct(h*jump(grad(u).t()*n), 1e-2*h*jump(grad(v).t()*n)),
		In
	);
 // in case of time mesh, the interface does not go through all the elements at one time
 // projection.addDiagonal(1e-14);
 projection.addLagrangeMultiplier(innerProduct(1.,v), *interface[0], In, projection.qTime[0].x, valLagrange);
 // projection.addLagrangeMultiplier(innerProduct(1.,v), *interface[1], In, projection.qTime[1].x, valLagrange);
 // projection.addLagrangeMultiplier(innerProduct(1.,v), *interface[2], In, projection.qTime[2].x, valLagrange);
 projection.addLagrangeMultiplier(innerProduct(1.,v), interface, In, valLagrange);
 //
 // matlab::Export(projection.mat, "matProjLap.dat");
 // matlab::Export(projection.rhs, "rhsProjLap.dat");

	projection.solve();
	KN_<double> uh(projection.rhs(SubArray(ph.size(), 0)));
	ph.v = uh;
}



template<typename M>
void projection(FunFEM<M> & ph, double (*f)(const typename M::Rd, int i)){

	typedef typename M::Rd Rd;
	typedef TestFunction<Rd::d> FunTest;
	typedef GFESpace<M> FESpace;
	const FESpace& Vh(*ph.Vh);
	typedef typename FESpace::FElement::QF::QuadraturePoint QuadraturePoint;
	typedef typename FESpace::FElement FElement;

	FEM<M> projection(Vh);
	const CutFEM_Parameter& h(Parameter::h);
	FunTest u(Vh, Vh.N), v(Vh, Vh.N);
	projection.addBilinear(innerProduct(u,v));

	What_d Fop = Fwhatd(1);
	for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
		const FElement& FK(Vh[k]);
		const R meas = FK.getMeasure();
		KNMK<double> fv(FK.NbDoF(),1,1); //  the value for basic fonction

		for(int ipq = 0; ipq < projection.qf.getNbrOfQuads(); ++ipq)  {
			QuadraturePoint ip(projection.qf[ipq]); // integration point
			Rd mip = FK.map(ip); // quadrature point in global geometry
			const R Cint = meas * ip.getWeight();

			FK.BF(Fop,ip, fv);
			R val_fh = f(mip, 0);

			for(int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
        projection(FK.loc2glb(i)) +=  Cint *  fv(i,0,op_id) * val_fh ;
      }
		}

  }

	projection.solve();
	ph.v = projection.rhs;
}


template<typename M>
void projection(FunFEM<M> & ph, double (*f)(const typename M::Rd, int i, int dd)){

	typedef typename M::Rd Rd;
	typedef TestFunction<Rd::d> FunTest;
	typedef GFESpace<M> FESpace;
	const FESpace& Vh(*ph.Vh);
	typedef typename FESpace::FElement::QF::QuadraturePoint QuadraturePoint;
	typedef typename FESpace::FElement FElement;

	CutFEM<M> projection(Vh);
	const CutFEM_Parameter& h(Parameter::h);
	FunTest u(Vh, Vh.N), v(Vh, Vh.N);
	projection.addBilinear(innerProduct(u,v));
	projection.addFaceStabilization(
		  innerProduct(jump(u), 1e-2*jump(v))
		+ innerProduct((h^2)*jump(grad(u)), 1e-2*jump(grad(v)))
	);

	What_d Fop = Fwhatd(1);
	for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {
		const FElement& FK(Vh[k]);
		const R meas = FK.getMeasure();
		KNMK<double> fv(FK.NbDoF(),1,1); //  the value for basic fonction

		for(int ipq = 0; ipq < projection.qf.getNbrOfQuads(); ++ipq)  {
			QuadraturePoint ip(projection.qf[ipq]); // integration point
			Rd mip = FK.map(ip); // quadrature point in global geometry
			const R Cint = meas * ip.getWeight();

			FK.BF(Fop,ip, fv);
			R val_fh = f(mip, 0, 0);

			for(int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
        projection(FK.loc2glb(i)) +=  Cint *  fv(i,0,op_id) * val_fh ;
      }
		}

  }

	projection.solve();
	ph.v = projection.rhs;
}




#endif
