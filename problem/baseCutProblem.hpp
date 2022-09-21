#ifndef BASE_CUTPROBLEM_HPP
#define BASE_CUTPROBLEM_HPP


// #include "baseProblem.hpp"


template<typename Mesh>
class BaseCutFEM : public BaseFEM<Mesh> {

  typedef ActiveMesh<Mesh> CutMesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::Rd Rd;
  typedef typename FElement::QF QF;
  typedef typename FElement::QFB QFB;
  typedef typename Mesh::Element Element;
  typedef typename Mesh::BorderElement BorderElement;


public:
  BaseCutFEM(const QuadratureFormular1d& qt, const ProblemOption& option) : BaseFEM<Mesh>(qt, option) {}
  // BaseCutFEM(const list<FESpace*>& vh,const ProblemOption& option) : BaseFEM<Mesh>(vh, option) {}
  BaseCutFEM(const FESpace& vh,const ProblemOption& option) : BaseFEM<Mesh>(vh, option) {}


  // Integral on K
  void addBilinear(  const ListItemVF<Rd::d>&, const CutMesh&) ;
  void addBilinear(  const ListItemVF<Rd::d>& VF, const CutMesh&, int itq, const TimeSlab& In) ;
  void addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh& Th, const TimeSlab& In);
  void addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh& Th, const TimeSlab& In, int itq);
  void addLinear(  const ListItemVF<Rd::d>& VF, const CutMesh&) ;
  void addLinear(  const ListItemVF<Rd::d>& VF, const CutMesh&, int itq, const TimeSlab& In) ;
  void addLinear(const ListItemVF<Rd::d>& VF, const CutMesh& Th, const TimeSlab& In);
  void addLinear(const ListItemVF<Rd::d>& VF, const CutMesh& Th, const TimeSlab& In, int itq);
  void addElementContribution(const ListItemVF<Rd::d>& VF, const int k, const TimeSlab* In, int itq, double cst_time);

  void addBilinear(const ListItemVF<Rd::d>&, const CutMesh&, const CExtension&, const int);
  void addLinear(const ListItemVF<Rd::d>&, const CutMesh&, const CExtension&, const int);
  void addElementContributionOtherSide(const ListItemVF<Rd::d>&, const int, const TimeSlab*, int, double);

  // integral on innerFace
  void addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh&, const CFacet& b);
  void addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh&, const CFacet& b, const TimeSlab& In);
  void addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh&, const CFacet& b, const TimeSlab& In, int itq);
  void addLinear  (const ListItemVF<Rd::d>& VF, const CutMesh&, const CFacet& b);
  void addLinear(const ListItemVF<Rd::d>& VF, const CutMesh&, const CFacet& b, const TimeSlab& In);
  void addLinear(const ListItemVF<Rd::d>& VF, const CutMesh&, const CFacet& b, const TimeSlab& In, int itq);
  void addFaceContribution(const ListItemVF<Rd::d>& VF, const std::pair<int,int>& e1, const std::pair<int,int>& e2, const TimeSlab* In, int itq, double cst_time);


  // integral on boundary
  void addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh&, const CBorder& b,list<int> label = {});
  void addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh&, const CBorder& b, const TimeSlab& In,list<int> label = {});
  void addBilinear(const ListItemVF<Rd::d>& VF, const CutMesh&, const CBorder& b, const TimeSlab& In, int itq,list<int> label = {});

  void addLinear  (const ListItemVF<Rd::d>& VF, const CutMesh&, const CBorder& b,list<int> label = {});
  void addLinear(const ListItemVF<Rd::d>& VF, const CutMesh&, const CBorder& b, const TimeSlab& In,list<int> label = {});
  void addLinear(const ListItemVF<Rd::d>& VF, const CutMesh&, const CBorder& b, const TimeSlab& In, int itq,list<int> label = {});
  void addBorderContribution(const ListItemVF<Rd::d>& VF, const Element& K,const BorderElement& BE, int ifac, const TimeSlab* In, int itq, double cst_time);

  void setDirichlet(const FunFEM<Mesh>& gh, const CutMesh& Th, list<int> label = {});

  // integral on interface
  void addBilinear(const ListItemVF<Rd::d>& VF, const Interface<Mesh>& gamma,list<int> label = {}) {return BaseFEM<Mesh>::addBilinear(VF, gamma, label);}
  void addBilinear(const ListItemVF<Rd::d>& VF, const Interface<Mesh>& gamma, const TimeSlab& In, int itq, list<int> label = {}) {return BaseFEM<Mesh>::addBilinear(VF, gamma, In, itq,label);}
  void addBilinear(const ListItemVF<Rd::d>& VF, const TimeInterface<Mesh>& gamma, const TimeSlab& In, list<int> label = {}){return BaseFEM<Mesh>::addBilinear(VF, gamma, In, label);}
  void addBilinear(const ListItemVF<Rd::d>& VF, const TimeInterface<Mesh>& gamma, const TimeSlab& In, int itq, list<int> label = {}) {return BaseFEM<Mesh>::addBilinear(VF, gamma, In, itq, label);}

  void addLinear  (const ListItemVF<Rd::d>& VF, const Interface<Mesh>& gamma,list<int> label = {}) {return BaseFEM<Mesh>::addLinear(VF, gamma, label);}
  void addLinear(const ListItemVF<Rd::d>& VF, const Interface<Mesh>& gamma, const TimeSlab& In, int itq, list<int> label = {}) {return BaseFEM<Mesh>::addLinear(VF, gamma, In, itq,label);}
  void addLinear(const ListItemVF<Rd::d>& VF, const TimeInterface<Mesh>& gamma, const TimeSlab& In, list<int> label = {}) {return BaseFEM<Mesh>::addLinear(VF, gamma, In,label);}
  void addLinear(const ListItemVF<Rd::d>& VF, const TimeInterface<Mesh>& gamma, const TimeSlab& In, int itq, list<int> label = {}) {return BaseFEM<Mesh>::addLinear(VF, gamma, In, itq, label);}
  // void addInterfaceContribution(const ListItemVF<Rd::d>& VF, const Interface<Mesh>& gamma, int ifac, double tid, const TimeSlab* In, double cst_time);

  // integral on inner Ridge / intersction with interface
  void addBilinear(const ListItemVF<Rd::d>& VF, const Interface<Mesh>& gamma, const CRidge& innerRidge,list<int> label = {});
  void addBilinear(const ListItemVF<Rd::d>& VF, const TimeInterface<Mesh>& gamma, const CRidge& innerRidge, const TimeSlab& In, list<int> label = {});
  void addBilinear(const ListItemVF<Rd::d>& VF, const TimeInterface<Mesh>& gamma, const CRidge& innerRidge, const TimeSlab& In, int itq,list<int> label = {});

  void addLinear  (const ListItemVF<Rd::d>& VF, const Interface<Mesh>& gamma, const CRidge& innerRidge,list<int> label = {});
  void addLinear(const ListItemVF<Rd::d>& VF, const TimeInterface<Mesh>& gamma, const CRidge& innerRidge, const TimeSlab& In, list<int> label = {});
  void addLinear(const ListItemVF<Rd::d>& VF, const TimeInterface<Mesh>& gamma, const CRidge& innerRidge, const TimeSlab& In, int itq,list<int> label = {});
  void addInterfaceRidgeContribution(const ListItemVF<Rd::d>& VF, const Interface<Mesh>& interface, int ifac, const TimeSlab* In, int itq, double cst_time);


  // Face stabilization
  void addFaceStabilization(const ListItemVF<Rd::d>& VF, const CutMesh&);
  void addFaceStabilization(const ListItemVF<Rd::d>& VF, const CutMesh&,const TimeSlab& In);
  void addFaceStabilization(const ListItemVF<Rd::d>& VF, const CutMesh&,const TimeSlab& In, int itq);
  void addFaceStabilization(const ListItemVF<Rd::d>& VF, const CutMesh&, const MacroElement<Mesh>& );

  // Lagrange multiplier
  void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, double val, const CutMesh&);
  void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, double val, const CutMesh&, const int k);
  void addLagrangeContribution(const ListItemVF<Rd::d>& VF, const int k);
  void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, double val, const CutMesh&,const CBorder& b,list<int> label = {});
  void addLagrangeBorderContribution(const ListItemVF<Rd::d>& VF, const Element& K,const BorderElement& BE, int ifac, const TimeSlab* In, int itq, double cst_time);
  void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, double val, const CutMesh& Th, const CExtension& ext, const int epsE);
  void addLagrangeContributionOtherSide(const ListItemVF<Rd::d>& VF, const int k, const int epsE);


  void removeDofForHansbo(const FESpace& Vh) ;

  // For time problem
  void saveSolution(const Rn&);
  void initialSolution(Rn&);


};

template<typename Mesh>
class CutFEM : public BaseCutFEM<Mesh>, public Solver{
  typedef GFESpace<Mesh> FESpace;

public:

  CutFEM(const QuadratureFormular1d& qt, const ProblemOption& option = defaultProblemOption) : BaseCutFEM<Mesh>(qt, option), Solver(option) {}
  CutFEM(const FESpace& vh        , const ProblemOption& option = defaultProblemOption) : BaseCutFEM<Mesh>(vh, option), Solver(option) {}
  void solve() {
    double tt0 = MPIcf::Wtime();
    Solver::solve(this->mat_, this->rhs_);
    if(this->verbose_>0)std::cout << " Real Time Solver \t \t " << MPIcf::Wtime() - tt0 << std::endl;
  }
  void solve(string solverName) {
    double tt0 = MPIcf::Wtime();
    this->solver_name_ = solverName;
    Solver::solve(this->mat_, this->rhs_);
    if(this->verbose_>0)std::cout << " Real Time Solver \t \t " << MPIcf::Wtime() - tt0 << std::endl;
  }
  void solve(std::map<std::pair<int,int>,R> & A, Rn & b, string solverName = "mumps") {
    R tt0 = MPIcf::Wtime();
    this->solver_name_ = solverName;
    Solver::solve(A, b);
    if(this->verbose_>0)std::cout << " Real Time Solver \t \t " << MPIcf::Wtime() - tt0 << std::endl;
  }

};


#include "baseCutProblem_Function.hpp"




#ifdef OLD_PROBLEM

template<typename M>
class BaseCutProblem : public BaseProblem<M> {
public:
  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::Rd Rd;
  typedef typename FElement::QF QF;
  typedef typename FElement::QFB QFB;
  typedef typename TypeCutData<Rd::d>::CutData CutData;
  typedef GenericInterface<Mesh> Interface;
  typedef GenericMapping<Mesh> Mapping;


  // typedef ListItemVF<Rd::d> ListItemVF;
  // typedef ItemVF<Rd::d> ItemVF;

  using BaseProblem<M>::Vh;
  using BaseProblem<M>::qf;
  using BaseProblem<M>::qfb;
  using BaseProblem<M>::qTime;


  BaseCutProblem(const FESpace& vh) : BaseProblem<M>(vh)
  {
    if(Vh->gamma.size() == 0){
      std::cout << " There is no interface linked to your FESpace " << std::endl;
      assert(0);
    }
  }

  BaseCutProblem(const list<FESpace*>& vh) : BaseProblem<M>(vh) {

  }


  BaseCutProblem(const QuadratureFormular1d& qT, int orderSpace=5)
  :BaseProblem<M>(qT,orderSpace) {}

public:
  void saveSolution(const Rn&);
  void initialSolution(Rn&);

public :
// STANDARD CUTFEM
// integral on element
  void addBilinear(const ListItemVF<Rd::d>& VF);
  void addBilinear(const ListItemVF<Rd::d>& VF, const MacroElement&);
  void addLinear(const ListItemVF<Rd::d>& VF);
  void addBilinearFormExtDomain(const ListItemVF<Rd::d>& VF,  const R epsE = 0);
  void addBilinearFormExtDomain(const ListItemVF<Rd::d>& VF, const MacroElement& macro, const R epsE);
  void addLinearFormExtDomain(const ListItemVF<Rd::d>& VF, const R epsE = 0);
  void addLinearFormExtDomain(const ListItemVF<Rd::d>& VF, const MacroElement& macro, const R epsE = 0);

  // integral on the boundary
  void addBilinear(const ListItemVF<Rd::d>& VF, const CBorder& b, list<int> label = {}){ return BaseProblem<M>::addBilinear(VF, b, label); };
  void addLinear(const ListItemVF<Rd::d>& VF, const CBorder& b, list<int> label = {}){ return BaseProblem<M>::addLinear(VF, b, label); };
  // integral on inner Edge/Face
  void addBilinear(const ListItemVF<Rd::d>& VF, const CFacet& b){ return BaseProblem<M>::addBilinear(VF,b); };
  void addBilinear(const ListItemVF<Rd::d>& VF, const CFacet& b, const GMacro& macro){ return BaseProblem<M>::addBilinear(VF, b,macro); };
  void addLinear  (const ListItemVF<Rd::d>& VF, const CFacet& b){ return BaseProblem<M>::addLinear(VF, b); };
  // face stabilization
  void addFaceStabilization(const ListItemVF<Rd::d>& VF);
  void addFaceStabilization(const ListItemVF<Rd::d>& VF, const MacroElement& bigMac);
  // lagrange multiplier
  void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, double val,const int itq = 0);

  protected :

  void addElementMat         (const ListItemVF<Rd::d>& VF , const int k,const bool extend=false, const R epsE = 0);
  void addElementRHS         (const ListItemVF<Rd::d>& VF, const int k,const bool extend=false, const R epsE = 0);
  void addElementLagrange    (const ListItemVF<Rd::d>& VF , const int k, const int itq);
  void addElementMatEdge     (const ListItemVF<Rd::d>& VF, const int k, const int ifac);
  void addElementRHSEdge     (const ListItemVF<Rd::d>& VF, const int k, const int ifac);
  void addElementMatEdge     (KNMK<double>& basisFunMat,const ListItemVF<Rd::d>& VF, const int k, int ifac);
  void addElementMatBorder   (const ListItemVF<Rd::d>& VF, const int ifac, int dom = 0);
  void addElementRHSBorder   (const ListItemVF<Rd::d>& VF, const int ifac, int dom = 0);


public:
  // TIME CUTFEM
  // on element
  void addBilinear   (const ListItemVF<Rd::d>& VF, const TimeSlab& In);
  void addLinear     (const ListItemVF<Rd::d>& VF, const TimeSlab& In);
  // stabilization
  void addFaceStabilization    (const ListItemVF<Rd::d>& VF, const TimeSlab& In) ;
  // lagrange multiplier
  void addLagrangeMultiplier   (const ListItemVF<Rd::d>& VF, double val, const TimeSlab& In) ;
  void addLagrangeMultiplier   (const ListItemVF<Rd::d>& VF, double val, int itq, const TimeSlab& In) ;

private:
  void addElementMat      (const ListItemVF<Rd::d>& VF, const int k, const TimeSlab& In);
  void addElementRHS      (const ListItemVF<Rd::d>& VF, const int k, const TimeSlab& In);
  void addElementMatEdge  (const ListItemVF<Rd::d>& VF, const int k, const int ifac, const TimeSlab& In);
  void addElementLagrange (const ListItemVF<Rd::d>& VF, const int k, const TimeSlab& In);


public:
  // FEM on surface
  void addBilinear(const ListItemVF<Rd::d>& VF, const Interface& gamma, list<int> label = {},const Mapping& mapping = DataMapping<Mesh>::Id){return BaseProblem<M>::addBilinear(VF, gamma, label, mapping);};
  void addBilinear(const ListItemVF<Rd::d>& VF, const Interface& gamma, const GMacro& macro, list<int> label = {},const Mapping& mapping = DataMapping<Mesh>::Id){return BaseProblem<M>::addBilinear(VF, gamma, macro, label, mapping);};
  void addBilinear(const ListItemVF<Rd::d>& VF, const Interface& gamma, const CFacet& b){return BaseProblem<M>::addBilinear(VF, gamma, b);};
  void addLinear  (const ListItemVF<Rd::d>& VF, const Interface& gamma, list<int> label = {}, const Mapping& mapping = DataMapping<Mesh>::Id){return BaseProblem<M>::addLinear(VF, gamma, label, mapping);};
  void addLinear  (const ListItemVF<Rd::d>& VF, const Interface& gamma, const CBorder& b, list<int> label = {}) {return BaseProblem<M>::addLinear(VF, gamma, b, label);};


  void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, const Interface& gamma, double val, const Mapping& mapping = DataMapping<Mesh>::Id){return BaseProblem<M>::addLagrangeMultiplier(VF, gamma, val, mapping);};
  void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, const Interface& gamma,  const TimeSlab& In, double tq, double val, const Mapping& mapping = DataMapping<Mesh>::Id){return BaseProblem<M>::addLagrangeMultiplier(VF, gamma, In, tq, val, mapping);};

  // SpaceTime FEM on surface
  void addBilinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, const TimeSlab& In, KN<const Mapping*> mapping = {}){return BaseProblem<M>::addBilinear(VF, gamma,In, mapping);};
  void addBilinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, int itq, const TimeSlab& In, KN<const Mapping*> mapping = {}){return BaseProblem<M>::addBilinear(VF, gamma,itq, In, mapping);};
  void addLinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, const TimeSlab& In, KN<const Mapping*> mapping = {}){return BaseProblem<M>::addLinear(VF, gamma, In, mapping);};
  void addLinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, int itq, const TimeSlab& In, KN<const Mapping*> mapping = {}){return BaseProblem<M>::addLinear(VF, gamma, itq, In, mapping);};
  void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, const TimeSlab& In, double val, KN<const Mapping*> mapping = {}){return BaseProblem<M>::addLagrangeMultiplier(VF,gamma,In,val,mapping);};

public:

  R computeCoef(const ItemVF<Rd::d>&,int, const typename Mesh::Partition& cutK ) const;

};



template<typename M>
R BaseCutProblem<M>::computeCoef(const ItemVF<Rd::d>& item, int domain, const typename Mesh::Partition& cutK) const {
  R val = 1;
  double h = cutK.getEdgeLength();
  double meas = cutK.mesure(domain);
  double measK = cutK.T.mesure();
  for(int l=0;l<2;++l) {
    const vector<string>& listCoef = (l==0)?item.coefu : item.coefv;
    int domCoef = domain;
    for(int i=0;i<listCoef.size();++i) {
      string coef = listCoef[i];

      if(this->parameterList.find(coef)) {
        CutFEM_Parameter& p(*this->parameterList.listParameter[coef]);
        val *= p(domCoef, h, meas, measK);
      }
    }
  }
  return val;
}


/*
          Add Bilinear Form
*/
template<typename M>
void BaseCutProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF) {
  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

    addElementMat(VF, k);

  }
}

template<typename M>
void BaseCutProblem<M>::addBilinearFormExtDomain(const ListItemVF<Rd::d>& VF, const R epsE) {
  bool extend = true;
  // [loops through ALL elements, addElementMat checks if cut]
  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {
    addElementMat(VF, k); // [add one standard contribution]
    if(Vh->isCut(k)) { // [makes sure we extend only cut elements]
      addElementMat(VF, k, extend, epsE); // [add epsilon worth of other side, the supplied arguments limits to only cut elements]
    }

  }
}

template<typename M>
void BaseCutProblem<M>::addBilinearFormExtDomain(const ListItemVF<Rd::d>& VF, const MacroElement& macro, const R epsE) {
  bool extend = true;
  for(auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {
    for(int i=0;i<it->second.idx_element.size();++i){
      int k = it->second.idx_element[i];
      // addElementMat(VF, k); // [add one standard contribution]
      if(Vh->isCut(k)) { // [makes sure we extend only cut elements]
        addElementMat(VF, k, extend, epsE); // [add epsilon worth of other side, the supplied arguments limits to only cut elements]
      }
    }
  }
}

template<typename M>
void BaseCutProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF, const MacroElement& macro) {
  for(auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {
    for(int i=0;i<it->second.idx_element.size();++i){
      int k = it->second.idx_element[i];
      addElementMat(VF, k);

    }
  }
}

template<typename M>
void BaseCutProblem<M>::addElementMat(const ListItemVF<Rd::d>& VF, const int k, const bool extend, const R epsE) {
  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;

  const FElement& FK((*Vh)[k]);
  const int domain = FK.whichDomain();

  const int kb = Vh->Th(FK.T);
  CutData cutData(Vh->getInterface(0).getCutData(kb));     // get the cut data
  const Partition& cutK =  Partition(FK.T, cutData);  // build the cut

  ElementSignEnum domPart = cutK.what_part(domain);
  if(extend) {
    assert(domain != -1); // [if for some reason gets called in standard no cut pb]
    int otherDom = (domain == 0); // [gets the other domain]
    domPart = cutK.what_part(otherDom);
  }
  const R h = FK.T.lenEdge(0);
  for(typename Partition::const_element_iterator it = cutK.element_begin(domPart);
  it != cutK.element_end(domPart); ++it){

    const R meas = cutK.mesure(it);

    for(int l=0; l<VF.size();++l) {
      if(!VF[l].on(domain)) continue;

      int lastop = getLastop(VF[l].du, VF[l].dv);
      const int ku = (VF[l].domu != -1)? VF[l].fespaceU->idxElementFromBackMesh(kb, VF[l].domu) : k;
      const int kv = (VF[l].domv != -1)? VF[l].fespaceV->idxElementFromBackMesh(kb, VF[l].domv) : k;
      const FElement& FKu((*VF[l].fespaceU)[ku]);
      const FElement& FKv((*VF[l].fespaceV)[kv]);
      bool same = (VF[l].fespaceU==VF[l].fespaceV);
      RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
      RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
      What_d Fop = Fwhatd(lastop);

      R coef = computeCoef(VF[l],domain,cutK);

      this->initIndex(FKu, FKv);
      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = cutK.toK(it, ip);
        Rd cut_ip = FKv.T.toKref(mip); // in Kref
        const R Cint = (extend)?epsE * meas * ip.getWeight() : meas * ip.getWeight();

        FKu.BF(Fop, cut_ip, fu); // need point in local reference element
        if(!same) FKv.BF(Fop,cut_ip, fv);

        double val = VF[l].fxu_backMesh(kb, domain, mip);

        for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
          for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
            // (*this)(FKv.loc2glb(i),FKu.loc2glb(j)) +=  Cint * coef * VF[l].c * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du)*val;
            this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j)) +=  Cint * coef * VF[l].c * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du)*val;
          }
        }
      }
    }
  }

  this->resetIndex();
  this->addLocalContribution();
  // getchar();
}

template<typename M>
void BaseCutProblem<M>::addElementMatEdge(const ListItemVF<Rd::d>& VF, const int k, const int ifac) {
  typedef typename Mesh::Partition Partition;
  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;
  assert(VF[0].fespaceU && VF[0].fespaceV);
  const FESpace& Vhu = *VF[0].fespaceU;
  const FESpace& Vhv = *VF[0].fespaceV;
  const FElement& FK(Vhu[k]);
  const int kb = Vhu.idxElementInBackMesh(k);
  CutData cutData(Vh->getInterface(0).getCutData(kb));

  if(!Vhu.isCutSpace() || cutData.face_isnot_cut(ifac)){
  // if(cutData.edge_isnot_cut(ifac)){
    BaseProblem<M>::addElementMatEdge(VF, k, ifac);
    return;
  }
  double ss = 0.;
  // for(int the_domain=0;the_domain<2;++the_domain){
  int the_domain = FK.whichDomain();
  int ifacn = ifac;
  int kn = Vhu.getNeighborElement(k, ifacn, the_domain);


  if(kn == -1) return;         // border edge
  if(k > kn) return;           // only compute one time

  const Partition& cutK =  Partition(FK.T, cutData);  // build the cut

  Rd normal = FK.T.N(ifac);
  const R meas_cut = cutK.mesure(the_domain);
  // const R measK = cutK.mesure(the_domain);

  const R h = FK.T.lenEdge(ifac);
  const R meas = cutK.mesureEdge(ifac, the_domain);
  ss += meas;

  const FElement& FKn(Vhu[kn]);
  const int kbn = Vhu.idxElementInBackMesh(kn);
  CutData cutDatan(Vh->getInterface(0).getCutData(kbn));
  const Partition& cutKn =  Partition(FKn.T, cutDatan);  // build the cut

  double measK = cutK.mesure(the_domain) + cutKn.mesure(the_domain);


  for(int l=0; l<VF.size();++l) {

    assert(VF[l].fespaceU == VF[l].fespaceV);
    int lastop = getLastop(VF[l].du, VF[l].dv);

    // R coef = BaseProblem<M>::computeCoef(VF[l],h,meas, measK,the_domain) * VF[l].getCoef(normal);
    R coef = BaseCutProblem<M>::computeCoefEdge(VF[l],h,meas, measK,meas_cut, the_domain) * VF[l].getCoef(normal);

    const int ku = (VF[l].domu == 0)? k : kn;
    const int kv = (VF[l].domv == 0)? k : kn;
    const FElement& FKu(Vhu[ku]);
    const FElement& FKv(Vhv[kv]);
    this->initIndex(FKu, FKv);
    int kuback = Vh->idxElementInBackMesh(ku);
    int kvback = Vh->idxElementInBackMesh(kv);

    bool same = (ku == kv );
    What_d Fop = Fwhatd(lastop);
    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction


    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = cutK.toEdge(ifac, (RdHatBord)ip, the_domain); // mip is in the global edge
      const R Cint = meas * ip.getWeight();

      FKu.BF(Fop,FKu.T.toKref(mip), fu); // need point in local reference element
      if(!same)FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element

      double val = VF[l].fx_backMesh_U(kuback,the_domain, mip, normal)
                  *VF[l].fx_backMesh_V(kvback,the_domain, mip, normal);
      double Cst = Cint * VF[l].c  * val * coef;

      for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
        for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
          this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j))  +=  Cst * fv(i,VF[l].cv,VF[l].dv) * fu(j,VF[l].cu,VF[l].du);
        }
      }
    }
  }
// }
  this->resetIndex();
  this->addLocalContribution();

}


/*
Add Extended Linear Form
*/
template<typename M>
void BaseCutProblem<M>::addLinear(const ListItemVF<Rd::d>& VF) {
  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {
    // std::cout << k << std::endl;
    addElementRHS( VF, k);
  }
}

template<typename M>
void BaseCutProblem<M>::addLinearFormExtDomain(const ListItemVF<Rd::d>& VF, const R epsE) {
  bool extend = true;

  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {
    addElementRHS(VF, k);
    if(Vh->isCut(k)) { // [makes sure we extend only cut elements]
      addElementRHS(VF, k, extend, epsE); // [add epsilon worth of other side, the supplied arguments limits to only cut elements]
    }
  }
}

template<typename M>
void BaseCutProblem<M>::addLinearFormExtDomain(const ListItemVF<Rd::d>& VF, const MacroElement& macro, const R epsE) {
  bool extend = true;

  for(auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {
    for(int i=0;i<it->second.idx_element.size();++i){
      int k = it->second.idx_element[i];
      // addElementRHS(VF, k);
      if(Vh->isCut(k)) { // [makes sure we extend only cut elements]
        addElementRHS(VF, k, extend, epsE); // [add epsilon worth of other side, the supplied arguments limits to only cut elements]
      }
    }
  }
}

template<typename M>
void BaseCutProblem<M>::addElementRHS(const ListItemVF<Rd::d>& VF, const int k, const bool extend, const R epsE) {
  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;

  const FElement& FK((*Vh)[k]);
  const int domain = FK.whichDomain();
  const R h = FK.T.lenEdge(0);

  const int kb = Vh->Th(FK.T);
  CutData cutData(Vh->getInterface(0).getCutData(kb));     // get the cut data
  const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
  ElementSignEnum domPart = cutK.what_part(domain);

  if(extend) {
    assert(domain != -1); // [if for some reason gets called in standard no cut pb]
    int otherDom = (domain == 0); // [gets the other domain]
    domPart = cutK.what_part(otherDom);
  }

  for(typename Partition::const_element_iterator it = cutK.element_begin(domPart);
  it != cutK.element_end(domPart); ++it){

    const R meas = cutK.mesure(it);

    for(int l=0; l<VF.size();++l) {
      if(!VF[l].on(domain)) continue;

      int lastop = getLastop(0, VF[l].dv);
      What_d Fop = Fwhatd(lastop);
      const int kv = (VF[l].domu != -1)? VF[l].fespaceV->idxElementFromBackMesh(kb, VF[l].domv) : k;
      const FElement& FKv((*VF[l].fespaceV)[kv]);
      RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction

      R coef = computeCoef(VF[l],domain,cutK);
      this->initIndex( VF[l].fespaceV);

      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = cutK.toK(it, ip);
        Rd cut_ip = FKv.T.toKref(mip); // in Kref
        const R Cint = (extend)? epsE * meas * ip.getWeight() :meas * ip.getWeight();

        FKv.BF(Fop,cut_ip, fv);

        R val_fh = VF[l].fxu_backMesh(kb, domain, mip);

        for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
          (*this)(FKv.loc2glb(i)) +=  Cint * coef * val_fh * VF[l].c * fv(i,VF[l].cv,VF[l].dv);
        }
      }
    }
  }
  this->resetIndex();
}



template<typename M>
void BaseCutProblem<M>::addElementMatBorder(const ListItemVF<Rd::d>& VF, const int ifac, int dddd) {

  typedef typename Mesh::Partition Partition;
  typedef typename Mesh::BorderElement BorderElement;
  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;


  if(!Vh->isCutSpace()){
    BaseProblem<M>::addElementMatBorder(VF, ifac, dddd);
    return;
  }

  const BorderElement & BE(Vh->Th.be(ifac)); // The face
  int e; // index of face of triangle corresp to edge (0,1,2)
  const int kb = Vh->Th.BoundaryElement(ifac, e); // index of element (triangle), ifaceK gets modified inside
  if(!Vh->containBackElement(kb)) return;

  CutData cutData(Vh->getInterface(0).getCutData(kb));     // get the cut data
  if(cutData.edge_isnot_cut(e) ){
    int domain = (cutData.sign[0] == -1)? 1 : 0;
    BaseProblem<M>::addElementMatBorder(VF, ifac, domain);
    return;
  }
  return;


  const Partition& cutK =  Partition(Vh->Th[kb], cutData);  // build the cut

  // integration on both domain
  for(int dom=0; dom<2; ++dom) {
    const int k = Vh->idxElementFromBackMesh(kb, dom);
    const FElement& FK((*Vh)[k]);

    // std::cout << " hey " << ifac << "\t" << k << std::endl;

    Rd normal = FK.T.N(e);
    const R measK = cutK.mesure(dom);
    const R h = FK.T.lenEdge(e);
    const R meas = cutK.mesureEdge(e, dom);


    for(int l=0; l<VF.size();++l) {

      int lastop = getLastop(VF[l].du, VF[l].dv);
      What_d Fop = Fwhatd(lastop);

      const int ku = VF[l].fespaceU->idxElementFromBackMesh(kb,dom);
      const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb,dom);
      const FElement& FKu((*VF[l].fespaceU)[ku]);
      const FElement& FKv((*VF[l].fespaceV)[kv]);
      bool same = (VF[l].fespaceU==VF[l].fespaceV);
      RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
      RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction

      R coef = BaseProblem<M>::computeCoef(VF[l],h,meas,measK,dom)*VF[l].getCoef(normal);

      // this->initIndex(VF[l].fespaceU, VF[l].fespaceV);
      this->initIndex(FKu, FKv);

      for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qfb[ipq]); // integration point

        const Rd mip = cutK.toEdge(e, (RdHatBord)ip, dom); // mip is in the global edge


        // std::cout << ip.x << "\t" << mip << std::endl;
        const R Cint = meas * ip.getWeight();

        FKu.BF(Fop, FKu.T.toKref(mip), fu); // need point in local reference element
        if(!same) FKv.BF(Fop,FKv.T.toKref(mip), fv);

        double val = VF[l].fxu(ku, mip);
        double Cst = Cint * coef * VF[l].c  * val;

        for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
          for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
            // (*this)(FKv.loc2glb(i),FKu.loc2glb(j)) +=  Cst * fv(i,VF[l].cv,VF[l].dv) * fu(j,VF[l].cu,VF[l].du);
            this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j)) +=  Cst * fv(i,VF[l].cv,VF[l].dv) * fu(j,VF[l].cu,VF[l].du);
          }
        }
      }
    }

  }

  this->resetIndex();
  this->addLocalContribution();
}

template<typename M>
void BaseCutProblem<M>::addElementRHSBorder(const ListItemVF<Rd::d>& VF, const int ifac, int dom) {
  typedef typename Mesh::Partition Partition;
  typedef typename Mesh::BorderElement BorderElement;
  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;

  if(!Vh->isCutSpace()){
    BaseProblem<M>::addElementRHSBorder(VF, ifac, dom);
    return;
  }


  const BorderElement & BE(Vh->Th.be(ifac)); // The face
  int e; // index of face of triangle corresp to edge (0,1,2)
  const int kb = Vh->Th.BoundaryElement(ifac, e); // index of element (triangle), ifaceK gets modified inside
  if(!Vh->containBackElement(kb)) return;

  CutData cutData(Vh->getInterface(0).getCutData(kb));     // get the cut data
  if(cutData.edge_isnot_cut(e) ){
    int domain = (cutData.sign[(e+1)%3] == -1)? 1 : 0;
    if(domain >= Vh->nbDomain()) return;
    BaseProblem<M>::addElementRHSBorder(VF, ifac, domain);
    return;
  }



  const Partition& cutK =  Partition(Vh->Th[kb], cutData);  // build the cut
  double val = 0.;
  // integration on both domain
  for(int dom=0; dom<Vh->nbDomain(); ++dom) {
    const int k = Vh->idxElementFromBackMesh(kb, dom);
    const FElement& FK((*Vh)[k]);
    Rd normal = FK.T.N(e);
    const R measK = cutK.mesure(dom);
    const R h = FK.T.lenEdge(e);
    const R meas = cutK.mesureEdge(e, dom);

    for(int l=0; l<VF.size();++l) {
      int lastop = getLastop(0, VF[l].dv);
      const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb, dom);
      const FElement& FKv((*VF[l].fespaceV)[kv]);

      RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
      What_d Fop = Fwhatd(lastop);

      R coef = BaseProblem<M>::computeCoef(VF[l],h,meas,measK,dom);
      this->initIndex(VF[l].fespaceV);

      for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qfb[ipq]); // integration point
        const Rd mip = cutK.toEdge(e, (RdHatBord)ip, dom); // mip is in the global edge
        const R Cint = meas * ip.getWeight();
// val += Cint;
        FKv.BF(Fop,FK.T.toKref(mip), fv);

        double cst_normal = VF[l].getCoef(normal);
        double val_fh = VF[l].fxu_backMesh(kb, dom, mip, normal);

        for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
          (*this)(FKv.loc2glb(i)) +=  Cint * coef * cst_normal * val_fh * VF[l].c * fv(i,VF[l].cv,VF[l].dv);
        }
      }
    }
  }
   this->resetIndex();
}


/*
Add Stabilization
*/
template<typename M>
void BaseCutProblem<M>::addFaceStabilization(const ListItemVF<Rd::d>& VF) {
typedef typename Mesh::Element Element;
const FESpace& Sh =(VF[0].fespaceU)? *VF[0].fespaceU : *Vh;
  for(int k=Sh.first_element(); k<Sh.last_element(); k+= Sh.next_element()) {

    if(!Sh.isCut(k)) continue;
    for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces
      BaseProblem<M>::addElementMatEdge(VF, k, ifac);
    }
  }
}

template<typename M>
void BaseCutProblem<M>::addFaceStabilization(const ListItemVF<Rd::d>& VF, const MacroElement& bigMac) {

  for(auto me=bigMac.macro_element.begin(); me!=bigMac.macro_element.end();++me) {
    for(auto it=me->second.inner_edge.begin(); it!=me->second.inner_edge.end();++it){

      int k = it->first;
      int ifac = it->second;
      BaseProblem<M>::addElementMatEdge(VF, k, ifac);
    }
  }
}


template<typename M>
void BaseCutProblem<M>::addElementRHSEdge(const ListItemVF<Rd::d>& VF, const int k, const int ifac) {

  typedef typename Mesh::Partition Partition;
  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;
  assert(VF[0].fespaceV);
  const FESpace& Vhv = *VF[0].fespaceV;
  const FElement& FK(Vhv[k]);
  const int kb = Vhv.idxElementInBackMesh(k);

  CutData cutData(Vh->getInterface(0).getCutData(kb));

  if(!Vhv.isCutSpace() || cutData.edge_isnot_cut(ifac)){
    BaseProblem<M>::addElementRHSEdge(VF, k, ifac);
    return;
  }
  int the_domain = FK.whichDomain();
  int ifacn = ifac;
  int kn = Vhv.getNeighborElement(k, ifacn, the_domain);

  if(kn == -1) return;         // border edge
  if(k > kn) return;           // only compute one time

  const Partition& cutK =  Partition(FK.T, cutData);  // build the cut

  Rd normal = FK.T.N(ifac);
  const R meas_cut = cutK.mesure(the_domain);
  const R h = FK.T.lenEdge(ifac);
  const R meas = cutK.mesureEdge(ifac, the_domain);


  const FElement& FKn(Vhv[kn]);
  const int kbn = Vhv.idxElementInBackMesh(kn);
  CutData cutDatan(Vh->getInterface(0).getCutData(kbn));
  const Partition& cutKn =  Partition(FKn.T, cutDatan);  // build the cut

  double measK = cutK.mesure(the_domain) + cutKn.mesure(the_domain);




  for(int l=0; l<VF.size();++l) {

    assert(VF[l].fespaceV);
    int lastop = getLastop(VF[l].du, VF[l].dv);
    // R coef = BaseProblem<M>::computeCoef(VF[l],h,meas, measK,the_domain) * VF[l].getCoef(normal);
    R coef = BaseCutProblem<M>::computeCoefEdge(VF[l],h,meas, measK,meas_cut, the_domain) * VF[l].getCoef(normal);

    const int ku = (VF[l].domu == 0)? k : kn;
    const int kv = (VF[l].domv == 0)? k : kn;
    const FElement& FKv(Vhv[kv]);
    this->initIndex(FKv);
    int kuback = Vh->idxElementInBackMesh(ku);
    int kvback = Vh->idxElementInBackMesh(kv);

    bool same = (ku == kv );
    What_d Fop = Fwhatd(lastop);
    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction


    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = cutK.toEdge(ifac, (RdHatBord)ip, the_domain); // mip is in the global edge
      const R Cint = meas * ip.getWeight();

      FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element

      double val = VF[l].fx_backMesh_U(kuback,the_domain, mip, normal)
                  *VF[l].fx_backMesh_V(kvback,the_domain, mip, normal);
      double Cst = Cint * VF[l].c  * val * coef;

      for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
          (*this)(FKv.loc2glb(i)) += Cst*fv(i,VF[l].cv,VF[l].dv);
      }
    }
  }
  this->resetIndex();
}

/*
Add Lagrange multiplier
*/

template<typename M>
void BaseCutProblem<M>::addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, double val, const int itq) {

  int ndf = this->rhs.size();
  this->rhs.resize(ndf+1);
  this->rhs(ndf) = val;
  // this->nDoF += 1;
  this->resetIndex();
  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

    addElementLagrange(VF, k, itq);

  }
}

template<typename M>
void BaseCutProblem<M>::addElementLagrange(const ListItemVF<Rd::d>& VF , const int k, const int itq){
  typedef typename Mesh::Partition Partition;
  typedef typename QF::QuadraturePoint QuadraturePoint;
  int nend = this->rhs.size()-1;

  for(int l=0; l<VF.size();++l) {

    const FESpace& Wh(*VF[l].fespaceV);
    const FElement& FK(Wh[k]);
    const int domain = FK.whichDomain();
    if(!VF[l].on(domain)) continue;

    const int kb = Wh.Th(FK.T);
    CutData cutData(Wh.getInterface(itq).getCutData(kb));     // get the cut data
    const Partition& cutK =  Partition(FK.T, cutData);  // build the cut
    ElementSignEnum the_part = cutK.what_part(domain);

    const R h = FK.T.lenEdge(0);

    int lastop = getLastop(VF[l].dv, 0);
    What_d Fop = Fwhatd(lastop);
    RNMK_ fv(this->databf,FK.NbDoF(),FK.N,lastop);

    for(typename Partition::const_element_iterator it = cutK.element_begin(the_part);
    it != cutK.element_end(the_part); ++it){

      const R meas = cutK.mesure(it);

      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = cutK.toK(it, ip);
        Rd cut_ip = FK.T.toKref(mip); // in Kref

        const R Cint = meas * ip.getWeight();

        FK.BF(Fop,cut_ip, fv); // need point in local reference element

        for(int j = FK.dfcbegin(VF[l].cv); j < FK.dfcend(VF[l].cv); ++j) {
          (*this)(nend, FK.loc2glb(j)+this->mapIdx0[&Wh]) +=  Cint * VF[l].c * fv(j,VF[l].cv,VF[l].dv);
          (*this)(FK.loc2glb(j)+this->mapIdx0[&Wh], nend) +=  Cint * VF[l].c * fv(j,VF[l].cv,VF[l].dv);
        }
      }
    }
  }
}



inline int getIdDomain(int id0, int domain) {
  if(domain == -1) return id0;
  return id0 + domain;
}

template<typename M>
void BaseCutProblem<M>::saveSolution(const Rn& sol) {

  typedef typename Mesh::Partition Partition;

  this->mapU0.clear();
  int nbTime = (*this->Ih)[0].NbDoF();
  int id_domain_0 = 0;

  for (typename std::map<const FESpace*,int>::const_iterator q=this->mapIdx0.begin(); q != this->mapIdx0.end(); ++q) {
    const FESpace& Wh     = *q->first;
    const int n0          =  q->second;
    const FESpace& backWh = Wh.getBackSpace();

    KN_<double> solS(sol(SubArray(Wh.NbDoF()*nbTime, n0)));

    for(int k=0;k<Wh.NbElement();++k) {
      const FElement& FK(Wh[k]);
      const int domain = FK.whichDomain();
      int id_domain = getIdDomain(id_domain_0, domain);

      int kback = Wh.idxElementInBackMesh(k);
      const FElement& FKback(backWh[kback]);

      for(int ic=0; ic<Wh.N;++ic) {
        for(int i=FK.dfcbegin(ic); i<FK.dfcend(ic);++i) {
          R val = 0.;
          for(int it=0;it<nbTime;++it) {
            val += solS(FK.loc2glb(i,it));

          }
          this->mapU0[std::make_pair(id_domain, FKback(i))] = val;
        }
      }
    }
    // std::cout << id_domain_0 << "\t" << Wh.getNumberOfSubDomain() << std::endl;

    id_domain_0 += Wh.getNumberOfSubDomain();
  }
  // std::cout <<  " test \t" <<  1 +	6.7582e-13 << std::endl;
}

template<typename M>
void BaseCutProblem<M>::initialSolution(Rn& u0) {

  this->DF.clear(); this->NL.clear();
  typedef typename Mesh::Partition Partition;
  int nbTime = (*this->Ih)[0].NbDoF();
  u0.init(this->nDoF);
  if(this->mapU0.size() == 0) {
    std::cout << " Default Initial solution " << std::endl;
    return;
  }

  int id_domain_0 = 0;

  for (typename std::map<const FESpace*,int>::const_iterator q=this->mapIdx0.begin(); q != this->mapIdx0.end(); ++q) {

    const FESpace& Wh     = *q->first;
    const int n0          = q->second;
    const FESpace& backWh = Wh.getBackSpace();

    KN_<double> u0S(u0(SubArray(Wh.NbDoF()*nbTime, n0)));

    for(int k=0;k<Wh.NbElement();++k) {

      const FElement& FK(Wh[k]);
      const int domain = FK.whichDomain();
      int id_domain = getIdDomain(id_domain_0, domain);

      int kb = Wh.idxElementInBackMesh(k);

      const FElement& FKback(backWh[kb]);

      CutData cutData(Wh.getInterface(0).getCutData(kb));     // get the cut data
      const Partition& cutK =  Partition(FK.T, cutData);
      ElementSignEnum the_part = cutK.what_part(domain);

      if( the_part != NoElement) {
        for(int ic=0; ic<Wh.N;++ic) {   // ESSAYER VH->N
          for(int i=FK.dfcbegin(ic); i<FK.dfcend(ic);++i) {
            u0S(FK(i)) = this->mapU0[std::make_pair(id_domain, FKback(i))];
          }
        }
      }
    }
    id_domain_0 += Wh.getNumberOfSubDomain();
  }

  this->mapU0.clear();
}


#include "baseCutProblemTime.hpp"
#endif






#endif
