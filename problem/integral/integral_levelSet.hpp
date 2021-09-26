#ifndef INTEGRAL_LAPLACEBELTRAMI_HPP_
#define INTEGRAL_LAPLACEBELTRAMI_HPP_


template<class M>
class IntegralLevelSet : public Integration {

  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::Rd Rd;
  typedef typename FElement::RdHat RdHat;

  typedef typename TypeOperationMV<Rd::d>::OperationMV MyOp;
  typedef typename FElement::QF QF;
  typedef typename QF::QuadraturePoint QuadraturePoint;
  
  typedef FEMFunction<FESpace> Fun;
  typedef FEMVecFunction<FESpace> FunVec;

  const QF &qf;
  const FunVec& Betap, Beta;
  const Fun& Up;
  const R dt;
public:
  
  IntegralLevelSet(const FESpace& Vh,
		   const FunVec& velp, const FunVec& vel,
		   const Fun& up, const R ddt) :
    Integration(Vh[0].NbDoF(),Vh.global.NbOfNodes, Vh[0].degre()),
    qf(*QF_Simplex<RdHat>(2*deg+1)),
    Betap(velp),
    Beta(vel),
    Up(up),
    dt(ddt)
  {
    
    w.init(nbDoF,1,op_dz+1);
    whatd = Fop_D0|Fop_D1;
  }


  void locAssembly(KNM<R>&  ML, KN<R>& VL,		   
		   const FElement& FK) {
    
    const R mes = FK.T.mesure();
    const R invdt = 1./dt;
    const R dt2 = invdt*invdt;
    const R h = FK.T.lenEdge(0);

    for(int ipq = 0; ipq < qf.n; ++ipq)  {             // loop over integration points
    
    QuadraturePoint ip(qf[ipq]);
    Rd mip = FK.T(ip);
    
    const R Cint = mes * ip.a;
    FK.BF(whatd, ip, w);                           // compute basic funtions

    Rd Beta = beta(FK.number, mip);
    Rd Betap = betap(FK.number, mip);
    R  up = Up(FK.number, mip);
    Rd Dup = Up.D(FK.number, mip);

    const R Tsd = 2. / (sqrt(dt2 + Norme2_2(Beta) * 1. / h / h));
    
    for(int i=0;i<nbDoF;++i){
      for(int j=0;j<nbDoF;++j){
      	ML(i,j) += Cint * ( invdt * w(j,0,op_id) + 0.5 * Op::innerProduct(w,Beta,j) )
	  * ( w(i,0,op_id) + Tsd * Op::innerProduct(w,Beta,i) );
      }
      VL(i) +=  Cint * (invdt * up - 0.5 * (Betap, Dup) )    
	* ( w(i,0,op_id) + Tsd * Op::innerProduct(w,Beta,i) );
    }
  }
    return sum;
  }


  
  
};






#endif
