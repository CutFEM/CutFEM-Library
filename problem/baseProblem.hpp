#ifndef BASE_PROBLEM_HPP
#define BASE_PROBLEM_HPP

#include "../FESpace/FESpace.hpp"
#include "problem.hpp"
#include "itemVF.hpp"
#include "GenericMapping.hpp"
// #include <omp.h>




template<typename M>
class BaseProblem : public ShapeOfNonLinProblem , public Solver {
public:
  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename FElement::Rd Rd;
  typedef typename FElement::QF QF;
  typedef typename FElement::QFB QFB;
  typedef TestFunction<Rd::d> FunTest;
  typedef GenericInterface<Mesh> Interface;
  typedef GenericMapping<Mesh> Mapping;
  typedef FunFEM<Mesh> Fun_h;

  CutFEM_ParameterList parameterList;

  const FESpace * Vh;
  const FESpace1 * Ih = nullptr;
  const QF &qf;
  const QFB &qfb;
  const QuadratureFormular1d& qTime;
  std::map<const FESpace*,int> mapIdx0;
  std::map<const FESpace*,int> mapIdx0_K;


  double* databf = nullptr;
  double tgv = 1e9;

  BaseProblem(const FESpace& vh, int orderSpace=5) :
  ShapeOfNonLinProblem(vh.NbDoF()),
  Solver(),
  Vh(&vh),
  qfb(*QF_Simplex<typename FElement::RdHatBord>(orderSpace)),
  qf(*QF_Simplex<typename FElement::RdHat>(orderSpace)),
  qTime(*Lobatto(1))
   {
    databf = new double[20*vh[0].NbDoF()*vh.N*4]; //4 <=> op_dz
    mapIdx0[Vh] = 0;
    mapIdx0_K[Vh] = 0;
    this->nt = Vh->NbElement();
  }

  BaseProblem(const list<FESpace*>& vh) :
    ShapeOfNonLinProblem(0),
    Solver(),
    Vh(*vh.begin()),
    qfb(*QF_Simplex<typename FElement::RdHatBord>(5)),
    qf(*QF_Simplex<typename FElement::RdHat>(5)),
    qTime(*Lobatto(1))
 {
   this->mapIdx0.clear();
   this->nDoF = 0;
   int NN = 0; int df = 0;
   this->nt=0;
   for(auto it = vh.begin();it!=vh.end();++it) {
     this->mapIdx0[*it] = this->nDoF;
     mapIdx0_K[Vh] = nt;
     this->nDoF += (*it)->NbDoF();
     NN += (*it)->N;
     df += (**it)[0].NbDoF();
     this->nt += (*it)->NbElement();
   }
   this->rhs.init(this->nDoF); rhs=0.;
   databf = new double[20*df*NN*op_DDall];
 }

  BaseProblem(const FESpace& vh, const FESpace1& ih, int orderSpace=5, int orderTime=3) :
  ShapeOfNonLinProblem(vh.NbDoF()*ih[0].NbDoF()),
  Solver(),
  Vh(&vh),
  Ih(&ih),
  qfb(*QF_Simplex<typename FElement::RdHatBord>(orderSpace)),
  qf(*QF_Simplex<typename FElement::RdHat>(orderSpace)),
  qTime(*Lobatto(orderTime))
 {
    databf = new double[20*vh[0].NbDoF()*vh.N*4]; //4 <=> op_dz
    mapIdx0[Vh] = 0;
    mapIdx0_K[Vh] = 0;
    this->nt=0;

  }


  BaseProblem(const QuadratureFormular1d& qT, int orderSpace=5) :
  ShapeOfNonLinProblem(),
  Solver(),
  Vh(nullptr),
  Ih(nullptr),
  qfb(*QF_Simplex<typename FElement::RdHatBord>(orderSpace)),
  qf(*QF_Simplex<typename FElement::RdHat>(orderSpace)),
  qTime(qT) {}


  void initSpace(const GFESpace<Mesh>& vh, const TimeSlab& In) {
    Vh = &vh;
    this->Ih = &In.Vh;
    this->nDoF = Vh->NbDoF() * In.NbDoF();
    this->rhs.resize(this->nDoF); this->rhs = 0.0;
    if(!this->databf) this->databf = new double[10*vh[0].NbDoF()*vh.N*4]; //4 <=> op_dz

    this->mapIdx0.clear();
    this->mapIdx0[Vh] = 0;
    this->mapIdx0_K[Vh] = 0;
    this->nt=0;

  }

  void add(const FESpace& Qh, const TimeSlab& IIn) {
    this->mapIdx0[&Qh] = this->nDoF;
    this->nDoF += Qh.NbDoF() * IIn.NbDoF();
    this->rhs.resize(this->nDoF); this->rhs = 0.0;
    this->mapIdx0_K[&Qh] = this->nt;
    this->nt += Qh.NbElement();
  }
  void add(const FESpace& Qh) {
    this->mapIdx0[&Qh] = this->nDoF;
    this->nDoF += Qh.NbDoF();
    this->mapIdx0_K[&Qh] = this->nt;
    this->nt += Qh.NbElement();
    this->rhs.resize(this->nDoF); this->rhs = 0.0;
  }


protected:
  void initIndex(const FESpace* Vhu, const FESpace* Vhv) {
    this->index_i0 = mapIdx0[Vhv];
    this->index_j0 = mapIdx0[Vhu];
  }
  void initIndex(const FESpace* Vhv) {
    this->index_i0 = mapIdx0[Vhv];
  }
  void initIndex(const FElement& FKu, const FElement& FKv) {
    this->index_i0 = mapIdx0[&FKv.Vh];
    this->index_j0 = mapIdx0[&FKu.Vh];
  }
  void initIndex( const FElement& FKv) {
    this->index_i0 = mapIdx0[&FKv.Vh];
  }

  virtual void solve();
  virtual void assembly() {assert(0);};
public :

void addDiagonal(double epsilon_machine);
void addDiagonal(const FESpace& Qh, double );
void setDiagonal(const FESpace& Qh, double );

// STANDARD FEM
// -----------------------------------------------------------------------------------
// on element
void addBilinear(const ListItemVF<Rd::d>& VF);
void addLinear(const ListItemVF<Rd::d>& VF);
// on boundary
void addBilinear(const ListItemVF<Rd::d>& VF, const CBorder& b, list<int> label = {});
void addLinear(const ListItemVF<Rd::d>& VF, const CBorder& b  , list<int> label = {});
void addStrongBC(const ExpressionVirtual& gh, list<int> label = {});
void addStrongBC(std::list<ExpressionFunFEM<typename typeMesh<Rd::d>::Mesh>> gh, list<int> label={});
// on innerEdge
void addBilinear(const ListItemVF<Rd::d>& VF, const CHyperFace& b);
void addLinear(const ListItemVF<Rd::d>& VF, const CHyperFace& b);
void addBilinear(const ListItemVF<Rd::d>& VF, const CHyperFace& b, const GMacro& macro);


// lagrange multiplier
void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, double val) ;


protected:
  virtual void addElementMat(const ListItemVF<Rd::d>& VF , const int k);
  virtual void addElementRHS(const ListItemVF<Rd::d>& VF, const int k);
  virtual void addElementMatBorder(const ListItemVF<Rd::d>& VF, const int ifac, int dom = 0);
  virtual void addElementRHSBorder(const ListItemVF<Rd::d>& VF, const int k, int dom = 0);
  void setElementStrongBC(int ifac, const ExpressionVirtual& gh);
  virtual void addElementMatEdge(const ListItemVF<Rd::d>& VF, const int k, const int ifac);
  virtual void addElementRHSEdge(const ListItemVF<Rd::d>& VF, const int k, const int ifac);

  virtual void addElementLagrange(const ListItemVF<Rd::d>& VF , const int k);


public:
  // SURFACE FEM
  // ---------------------------------------------------
  // on element
  void addBilinear(const ListItemVF<Rd::d>& VF, const Interface& gamma,list<int> label = {}, const Mapping& mapping = DataMapping<Mesh>::Id);
  void addBilinear(const ListItemVF<Rd::d>& VF, const Interface& gamma,const GMacro& macro, list<int> label = {}, const Mapping& mapping = DataMapping<Mesh>::Id);
  void addLinear(const ListItemVF<Rd::d>& VF, const Interface& gamma, list<int> label = {},const Mapping& mapping = DataMapping<Mesh>::Id);
  // on boundary
  void addLinear(const ListItemVF<Rd::d>& VF, const Interface& gamma, const CBorder& b, list<int> label = {});
  // on inner edges => node eval in 2D
  void addBilinear(const ListItemVF<Rd::d>& VF, const Interface& gamma, const CHyperFace& b);
  // lagrange multiplier
  void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, const Interface& gamma, double val, const Mapping& mapping = DataMapping<Mesh>::Id) ;
  void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, const Interface& gamma,  const TimeSlab& In, double tq, double val, const Mapping& mapping = DataMapping<Mesh>::Id) ;



private:
  void addElementMat(const ListItemVF<Rd::d>& VF, const Interface& gamma, const int ifac,const Mapping& mapping);
  void addElementRHS(const ListItemVF<Rd::d>& VF, const Interface& gamma, const int ifac,const Mapping& mapping);
  void addElementMatEdge(const ListItemVF<Rd::d>& VF,const Interface& interface, const int iface);
  void addElementRHSBorder(const ListItemVF<Rd::d>& VF,const Interface& interface, const int iface);
  void addElementLagrange(const ListItemVF<Rd::d>& VF , const Interface& gamma, const int k,const Mapping& mapping);
  void addElementLagrange(const ListItemVF<Rd::d>& VF , const Interface& gamma, const int k,const TimeSlab& In, double tq,const Mapping& mapping);


// STANDART TIME FEM
public:
  // on element
  void addBilinear(const ListItemVF<Rd::d>& VF, const TimeSlab& In);
  void addLinear(const ListItemVF<Rd::d>& VF, const TimeSlab& In);
  // on boundary
  void addBilinear(const ListItemVF<Rd::d>& VF, const TimeSlab& In, const CBorder& b, list<int> label = {});
  void addLinear(const ListItemVF<Rd::d>& VF  , const TimeSlab& In, const CBorder& b, list<int> label = {});
  void addStrongBC(const ExpressionVirtual& gh, const TimeSlab& In, list<int> label = {});
  void addStrongBC(std::list<ExpressionFunFEM<typename typeMesh<Rd::d>::Mesh>> gh, const TimeSlab& In, list<int> label={});
  // on inner edges
  void addBilinear(const ListItemVF<Rd::d>& VF, const TimeSlab& In, const CHyperFace& b);

protected:
  virtual void addElementMat(const ListItemVF<Rd::d>& VF, const int k, const TimeSlab& In);
  virtual void addElementRHS(const ListItemVF<Rd::d>& VF, const int k, const TimeSlab& In);
  void addElementMatBorder(const ListItemVF<Rd::d>& VF, const int ifac, const TimeSlab& In);
  void addElementRHSBorder(const ListItemVF<Rd::d>& VF, const int k, const TimeSlab& In);
  void setElementStrongBC(int ifac, const TimeSlab& In, const ExpressionVirtual& gh);
  virtual void addElementMatEdge(const ListItemVF<Rd::d>& VF, const int k, const int ifac, const TimeSlab& In);


public:
  // SURFACE TIME FEM
  // on elements
  void addBilinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, const TimeSlab& In, KN<const Mapping*> mapping = {});
  void addBilinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, int itq, const TimeSlab& In, KN<const Mapping*> mapping = {});
  void addLinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, const TimeSlab& In, KN<const Mapping*> mapping = {});
  void addLinear(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, int itq, const TimeSlab& In, KN<const Mapping*> mapping = {});

  // lagrange multiplier
  void addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, const TimeInterface<M>& gamma, const TimeSlab& In, double val, KN<const Mapping*> mapping = {}) ;


private:
  void addElementMat(const ListItemVF<Rd::d>& VF, const TimeInterface<M>&  interface, const int iface, const TimeSlab& In, const Mapping& mapping);
  void addElementRHS(const ListItemVF<Rd::d>& VF, const TimeInterface<M>&  interface, const int iface, const TimeSlab& In, const Mapping& mapping);
  void addElementLagrange(const ListItemVF<Rd::d>& VF , const TimeInterface<M>& gamma, const int k, const TimeSlab& In, const Mapping& mapping);


public:
  virtual R computeCoef(const ItemVF<Rd::d>&, double, double, double, int d = 0) const;
  R computeCoefInterface(const ItemVF<Rd::d>& item, double h, double meas, double measK) const;
  R computeCoefEdge(const ItemVF<Rd::d>&, double h, double meas, double measK, double meas_Cut, int d = 0) const;



  ~BaseProblem(){
    if(databf) delete databf;
  }
};




template<typename M>
class FEM : public BaseProblem<M>{
  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef GenericMapping<Mesh> Mapping;

public:
  FEM(const FESpace& vh) : BaseProblem<M>(vh, 5) {}
  FEM(const FESpace& vh, int orderSpace) : BaseProblem<M>(vh, orderSpace) {}
  FEM(const FESpace& vh, const FESpace1& Ih, int orderSpace=5, int orderTime=3) : BaseProblem<M>(vh, Ih,orderSpace, orderTime) {}

  FEM(const list<FESpace*>& vh) : BaseProblem<M>(vh) {}


  void solve() {
    // matlab::Export(this->mat, "matRT.dat");
    // matlab::Export(this->rhs, "rhs.dat");
    R t0 = CPUtime();
    Solver::solve(this->mat, this->rhs);
    R t1 = CPUtime();
    // std::cout << " Time Solver \t \t " << t1 - t0 << std::endl;
  }

  void solve(std::map<std::pair<int,int>,R> & A, Rn & b) {
    R tt0 = MPIcf::Wtime();
    Solver::solve(A, b);
    std::cout << " Real Time Solver \t \t " << MPIcf::Wtime() - tt0 << std::endl;

  }

};

template<typename M>
void BaseProblem<M>::addDiagonal(double epsilon_machine) {

  for(int i=0;i<this->nDoF;++i){
    this->mat[std::make_pair(i,i)] += epsilon_machine;
  }
}
template<typename M>
void BaseProblem<M>::addDiagonal(const FESpace& Qh, double val){
  int idxBegin = mapIdx0[&Qh];
  int idxEnd   = Qh.nbDoF;
  for(int i=idxBegin;i<idxEnd ;++i){
    this->mat[std::make_pair(i,i)] += val;
  }
}
template<typename M>
void BaseProblem<M>::setDiagonal(const FESpace& Qh, double val){
  int idxBegin = mapIdx0[&Qh];
  int idxEnd   = Qh.nbDoF;
  for(int i=idxBegin;i<idxEnd ;++i){
    this->mat[std::make_pair(i,i)] = val;
  }
}

template<typename M>
R BaseProblem<M>::computeCoef(const ItemVF<Rd::d>& item, double h, double meas, double measK, int domain) const {
  R val = 1;
  for(int l=0;l<2;++l) {
    const vector<string>& listCoef = (l==0)?item.coefu : item.coefv;
    for(int i=0;i<listCoef.size();++i) {
      string coef = listCoef[i];

      if(parameterList.find(coef)) {
        CutFEM_Parameter& p(*parameterList.listParameter[coef]);
        val *= p(domain, h, meas, measK);
      }
    }
  }
  return val;
}

template<typename M>
R BaseProblem<M>::computeCoefInterface(const ItemVF<Rd::d>& item, double h, double meas, double measK) const {
  R val = 1;
  for(int l=0;l<2;++l) {
    const vector<string>& listCoef = (l==0)?item.coefu : item.coefv;
    int domCoef = (l==0)?item.domu : item.domv;
    if(domCoef == -1) domCoef = 0;
    assert(domCoef == 0 || domCoef == 1);
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

template<typename M>
R BaseProblem<M>::computeCoefEdge(const ItemVF<Rd::d>& item, double h, double meas, double measK, double meas_Cut, int domain) const {
  R val = 1;
  for(int l=0;l<2;++l) {
    const vector<string>& listCoef = (l==0)?item.coefu : item.coefv;
    for(int i=0;i<listCoef.size();++i) {
      string coef = listCoef[i];

      if(this->parameterList.find(coef)) {
        CutFEM_Parameter& p(*this->parameterList.listParameter[coef]);
        val *= p(domain, h, meas, measK, meas_Cut);
      }
    }
  }
  return val;
}



template<typename M>
void BaseProblem<M>::solve() {
  R t0 = CPUtime();
  this->assembly();
  // std::cout << " Time assembly \t \t " << CPUtime()- t0 << std::endl;

  t0 = CPUtime();
  Solver::solve(mat, rhs);

  R t1 = CPUtime();
  // std::cout << " Time Solver \t \t " << t1 - t0 << std::endl;

}



/*
       ADD BILINEAR FORM DOMAIN
*/
template<typename M>
void BaseProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF) {

//   int iam = 0, np = 1;
// // //#pragma omp parallel default(shared) private(iam, np)
//   np = omp_get_num_threads();
//   iam = omp_get_thread_num();
//   printf("Hello from thread %d out of %d from process %d out of %d \n",
//   iam, np, MPIcf::my_rank(), MPIcf::size());
// #pragma omp  for private(index_i0, index_j0)
//default(shared)
  // private(index_i0, index_j0), firstprivate(localContributionMatrix)
  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

    addElementMat(VF, k);

//    #pragma omp critical
    this->addLocalContribution();

  }
}


template<typename M>
void BaseProblem<M>::addElementMat(const ListItemVF<Rd::d>& VF, const int k) {
  // std::map<std::pair<int,int>,R> localCont;
  typedef typename QF::QuadraturePoint QuadraturePoint;

  const FElement& FK((*Vh)[k]);
  const R meas = FK.getMeasure();
  const R h = FK.T.lenEdge(0);
  // const int kb = Vh->Th(FK.T);

  for(int l=0; l<VF.size();++l) {
    int lastop = getLastop(VF[l].du, VF[l].dv);

    const int ku = VF[l].fespaceU->idxElementFromBackMesh(k, VF[l].domu);
    const int kv = VF[l].fespaceV->idxElementFromBackMesh(k, VF[l].domv);
    // const int ku = (VF[l].domu != -1)? VF[l].fespaceU->idxElementFromBackMesh(kb, VF[l].domu) : k;
    // const int kv = (VF[l].domv != -1)? VF[l].fespaceV->idxElementFromBackMesh(kb, VF[l].domv) : k;
    const FElement& FKu((*VF[l].fespaceU)[ku]);
    const FElement& FKv((*VF[l].fespaceV)[kv]);
    bool same = (VF[l].fespaceU==VF[l].fespaceV);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);


    R coef = computeCoef(VF[l],h,meas,meas);
    this->initIndex(VF[l].fespaceU, VF[l].fespaceV);

    for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qf[ipq]); // integration point
      const Rd mip = FK.map(ip);
      const R Cint = meas * ip.getWeight();

      FKu.BF(Fop, ip, fu); // need point in local reference element
      if(!same) FKv.BF(Fop,ip, fv);
      VF[l].applyFunNL(fu,fv);

      double val = VF[l].fxu(ku, mip);
      for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
        for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
          // (*this)(FKv.loc2glb(i),FKu.loc2glb(j)) +=  Cint * coef * val * VF[l].c * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du);
          this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j)) +=  Cint * coef * val * VF[l].c * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du);
          // localCont[std::make_pair(FKv.loc2glb(i),FKu.loc2glb(j))] +=  Cint * coef * val * VF[l].c * fv(i,VF[l].cv,VF[l].dv)*fu(j,VF[l].cu,VF[l].du);
        }
      }
    }

  }
  this->resetIndex();
  // #pragma omp critical
  // for (auto q=localCont.begin(); q != localCont.end(); ++q) {
  //   (*this)(q->first.first,q->first.second) += q->second;
  // }
  // localCont.clear();
}




/*
        ADD BILINEAR FORM BOUNDARY
*/

template<typename M>
void BaseProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF, const CBorder& b, list<int> label) {
  typedef typename Mesh::BorderElement BorderElement;
  bool all_label = (label.size() == 0);

  for( int ifac = Vh->first_boundary_element(); ifac < Vh->last_boundary_element(); ifac+=Vh->next_boundary_element()) {
    const BorderElement & face(Vh->Th.be(ifac)); // The face
    if(util::contain(label, face.lab) || all_label) {

      addElementMatBorder(VF, ifac);
    }
  }
}

template<typename M>
void BaseProblem<M>::addElementMatBorder(const ListItemVF<Rd::d>& VF, const int ifac, int dom) {

  typedef typename Mesh::BorderElement BorderElement;
  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;

  const BorderElement & BE(Vh->Th.be(ifac)); // The face
  int ifaceK; // index of face of triangle corresp to edge (0,1,2)
  const int kb = Vh->Th.BoundaryElement(ifac, ifaceK); // index of element (triangle), ifaceK gets modified inside
  const int k = Vh->idxElementFromBackMesh(kb, dom);
  const FElement& FK((*Vh)[k]);

  Rd normal = FK.T.N(ifaceK);
  const R meas = BE.mesure();
  const R measK = FK.getMeasure();
  const R h = FK.T.lenEdge(0);

  int nb_face_onB = 0;
  for(int i=0;i<Element::nea;++i){
    int ib = i;
    if(Vh->Th.ElementAdj(kb,ib) == -1) nb_face_onB += 1;
  }
  assert(nb_face_onB > 0);
  double measOnB = nb_face_onB*meas;

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

    R coef = BaseProblem<M>::computeCoef(VF[l],h,measOnB,measK,dom);
    this->initIndex(FKu, FKv);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = BE((RdHatBord)ip); // mip is in the global edge
      const R Cint = meas * ip.getWeight();

      FKu.BF(Fop, FKu.T.toKref(mip), fu); // need point in local reference element
      if(!same) FKv.BF(Fop,FKv.T.toKref(mip), fv);

      double cst_normal = VF[l].getCoef(normal);
      double val = VF[l].fxu(ku, mip);
      double Cst = Cint * coef * VF[l].c * cst_normal * val;

      for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
        for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
          // (*this)(FKv.loc2glb(i),FKu.loc2glb(j)) +=  Cst * fv(i,VF[l].cv,VF[l].dv) * fu(j,VF[l].cu,VF[l].du);
          this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j)) +=  Cst * fv(i,VF[l].cv,VF[l].dv) * fu(j,VF[l].cu,VF[l].du);
        }
      }
    }
  }
  this->resetIndex();
  this->addLocalContribution();
}



/*
     ADD LINEAR FORM DOMAIN
*/
template<typename M>
void BaseProblem<M>::addLinear(const ListItemVF<Rd::d>& VF) {

  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {
    addElementRHS(VF, k);
  }
}

template<typename M>
void BaseProblem<M>::addElementRHS(const ListItemVF<Rd::d>& VF, const int k) {

  typedef typename QF::QuadraturePoint QuadraturePoint;
  const FElement& FK((*Vh)[k]);
  const R meas = FK.getMeasure();
  const R h = FK.T.lenEdge(0);
  for(int l=0; l<VF.size();++l) {
    int lastop = getLastop(0, VF[l].dv);

    // const int ku = VF[l].fespaceU->idxElementFromBackMesh(k, VF[l].domu);
    const int kv = VF[l].fespaceV->idxElementFromBackMesh(k, VF[l].domv);
    // const FElement& FKu((*VF[l].fespaceU)[ku]);
    const FElement& FKv((*VF[l].fespaceV)[kv]);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    R coef = computeCoef(VF[l],h,meas,meas);
    this->initIndex(VF[l].fespaceV);

    for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qf[ipq]); // integration point
      Rd mip = FK.map(ip); // quadrature point in global geometry
      const R Cint = meas * ip.getWeight();

      FKv.BF(Fop,ip, fv);

      R val_fh = VF[l].fxu(k, mip);

      for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
        (*this)(FKv.loc2glb(i)) +=  Cint * val_fh * coef * VF[l].c * fv(i,VF[l].cv,VF[l].dv);
      }
    }
  }
  this->resetIndex();
}



/*
     ADD LINEAR FORM BORDER
*/
template<typename M>
void BaseProblem<M>::addLinear(const ListItemVF<Rd::d>& VF, const CBorder& b, list<int> label) {
  typedef typename Mesh::BorderElement BorderElement;
  bool all_label = (label.size() == 0);

  for( int ifac = Vh->first_boundary_element(); ifac < Vh->last_boundary_element(); ifac+=Vh->next_boundary_element()) {        // loop over boundary faces

    const BorderElement & face(Vh->Th.be(ifac)); // The face
    if(util::contain(label, face.lab) || all_label) {

      addElementRHSBorder(VF, ifac);
    }
  }
}

template<typename M>
void BaseProblem<M>::addElementRHSBorder(const ListItemVF<Rd::d>& VF, const int ifac, int dom) {
  typedef typename Mesh::BorderElement BorderElement;
  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;

  const BorderElement & BE(Vh->Th.be(ifac)); // The face
  int ifaceK; // index of face of triangle corresp to edge (0,1,2)
  const int kb = Vh->Th.BoundaryElement(ifac, ifaceK); // index of element (triangle), ifaceK gets modified inside
  const int k = Vh->idxElementFromBackMesh(kb, dom);

  const FElement& FK((*Vh)[k]);
  Rd normal = FK.T.N(ifaceK);
  const R meas = BE.mesure();
  const R measK = FK.getMeasure();
  const R h = FK.T.lenEdge(0);

  int nb_face_onB = 0;
  for(int i=0;i<Element::nea;++i){
    int ib = i;
    if(Vh->Th.ElementAdj(kb,ib) == -1) nb_face_onB += 1;
  }
  assert(nb_face_onB > 0);
  double measOnB = nb_face_onB*meas;
  double val = 0.;
  for(int l=0; l<VF.size();++l) {
    int lastop = getLastop(0, VF[l].dv);
    const int kv = VF[l].fespaceV->idxElementFromBackMesh(kb, dom);
    const FElement& FKv((*VF[l].fespaceV)[kv]);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    R coef = BaseProblem<M>::computeCoef(VF[l],h,measOnB,measK,dom);
    this->initIndex(VF[l].fespaceV);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = BE((RdHatBord)ip); // mip is in the global edge
      const R Cint = meas * ip.getWeight();

      FKv.BF(Fop,FKv.T.toKref(mip), fv);

      double cst_normal = VF[l].getCoef(normal);
      double val_fh = VF[l].fxu_backMesh(kb, dom, mip, normal);

      val += Cint;
      for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
        (*this)(FKv.loc2glb(i)) +=  Cint * coef * cst_normal * val_fh * VF[l].c * fv(i,VF[l].cv,VF[l].dv);
      }
    }
  }
   this->resetIndex();
}




/*
    ADD LAGRANGE MULTIPLIER
*/
template<typename M>
void BaseProblem<M>::addLagrangeMultiplier(const ListItemVF<Rd::d>& VF, double val) {

  int ndf = rhs.size();
  rhs.resize(ndf+1);
  rhs(ndf) = MPIcf::IamMaster()*val;
  // nDoF += 1;

  for(int k=Vh->first_element(); k<Vh->last_element(); k+= Vh->next_element()) {

    addElementLagrange(VF, k);

  }
}

template<typename M>
void BaseProblem<M>::addElementLagrange(const ListItemVF<Rd::d>& VF , const int k){

  typedef typename QF::QuadraturePoint QuadraturePoint;
  int nend = rhs.size()-1;

  for(int l=0; l<VF.size();++l) {
    const FESpace& Wh(*VF[l].fespaceV);
    const FElement& FK(Wh[k]);
    const R meas = FK.getMeasure();


    int lastop = getLastop(VF[l].dv, 0);
    What_d Fop = Fwhatd(lastop);
    RNMK_ fv(this->databf,FK.NbDoF(),FK.N,lastop); //  the value for basic fonction

    for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qf[ipq]); // integration point
      const R Cint = meas * ip.getWeight();

      FK.BF(Fop,(Rd)ip, fv); // need point in local reference element

      for(int j = FK.dfcbegin(VF[l].cv); j < FK.dfcend(VF[l].cv); ++j) {
        (*this)(nend, FK.loc2glb(j)+mapIdx0[&Wh]) +=  Cint * VF[l].c * fv(j,VF[l].cv,VF[l].dv);
        (*this)(FK.loc2glb(j)+mapIdx0[&Wh], nend) +=  Cint * VF[l].c * fv(j,VF[l].cv,VF[l].dv);
      }
    }
  }
}




/*
     ADD EDGE STABILIZATION
*/
template<typename M>
void BaseProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF,  const CHyperFace& b) {
  typedef typename Mesh::Element Element;
  const FESpace& Sh =(VF[0].fespaceU)? *VF[0].fespaceU : *Vh;

  for(int k=Sh.first_element(); k<Sh.last_element(); k+= Sh.next_element()) {
    for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces
      addElementMatEdge(VF, k, ifac);

    }
  }
}

template<typename M>
void BaseProblem<M>::addBilinear(const ListItemVF<Rd::d>& VF,  const CHyperFace& b, const GMacro& macro) {

  for(auto it = macro.macro_element.begin(); it != macro.macro_element.end(); ++it) {

    for(int i=0;i<it->second.inner_edge.size();++i){
      int k = it->second.inner_edge[i].first;
      int ie = it->second.inner_edge[i].second;
      addElementMatEdge(VF, k, ie);
    }
  }
}

template<typename M>
void BaseProblem<M>::addElementMatEdge(const ListItemVF<Rd::d>& VF, const int k, const int ifac) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;

  assert(VF[0].fespaceU && VF[0].fespaceV);
  const FESpace& Vhu = *VF[0].fespaceU;
  const FESpace& Vhv = *VF[0].fespaceV;
  const FElement& FK(Vhu[k]);
  const Element & K = FK.T;                  // the triangle
  int the_domain = FK.whichDomain();

  int ifacn = ifac;
  int kn = Vhu.getNeighborElement(k, ifacn, the_domain);
  if(kn == -1) return;         // border edge
  if(k > kn) return;           // only compute one time

  const Rd normal = FK.T.N(ifac);
  const R meas    = FK.T.mesureBord(ifac);
  const R measK   = FK.getMeasure();
  const R h       = FK.T.lenEdge(ifac);//meas;

  const FElement& FKn(Vhu[kn]);
  const R meas_macroK  = measK + FKn.getMeasure();

  for(int l=0; l<VF.size();++l) {

    assert(VF[l].fespaceU == VF[l].fespaceV);
    int lastop = getLastop(VF[l].du, VF[l].dv);
    double coef = BaseProblem<M>::computeCoefEdge(VF[l],h,meas,meas_macroK,measK,the_domain) * VF[l].getCoef(normal) ;

    const int ku = (VF[l].domu == 0)? k : kn;
    const int kv = (VF[l].domv == 0)? k : kn;
    bool same = (ku == kv);

    const FElement& FKu(Vhu[ku]);
    const FElement& FKv(Vhv[kv]);
    this->initIndex(FKu, FKv);
    int kuback = Vh->idxElementInBackMesh(ku);
    int kvback = Vh->idxElementInBackMesh(kv);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = K(K.toKref((RdHatBord)ip, ifac));
      const R Cint = meas * ip.getWeight();

      FKu.BF(Fop,FKu.T.toKref(mip), fu); // need point in local reference element
      if(!same)FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element
      double val = VF[l].fx_backMesh_U(kuback,the_domain, mip, normal)
                  *VF[l].fx_backMesh_V(kvback,the_domain, mip, normal);
      double Cst = Cint * VF[l].c  * val * coef;


      for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
        for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
          this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j)) += Cst* fu(j,VF[l].cu,VF[l].du)*fv(i,VF[l].cv,VF[l].dv);
        }
      }
    }
  }
  this->resetIndex();
  this->addLocalContribution();
}

template<typename M>
void BaseProblem<M>::addLinear(const ListItemVF<Rd::d>& VF, const CHyperFace& b) {
  typedef typename Mesh::Element Element;
  const FESpace& Sh =(VF[0].fespaceU)? *VF[0].fespaceU : *Vh;

  for(int k=Sh.first_element(); k<Sh.last_element(); k+= Sh.next_element()) {
    for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces
      addElementRHSEdge(VF, k, ifac);
    }
  }
}

template<typename M>
void BaseProblem<M>::addElementRHSEdge(const ListItemVF<Rd::d>& VF, const int k, const int ifac) {

  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;

  assert(VF[0].fespaceV);
  const FESpace& Vhv = *VF[0].fespaceV;

  const FElement& FK(Vhv[k]);
  const Element & K = FK.T;                  // the triangle
  int the_domain = FK.whichDomain();

  int ifacn = ifac;
  int kn = Vhv.getNeighborElement(k, ifacn, the_domain);
  if(kn == -1) return;         // border edge
  if(k > kn) return;           // only compute one time

  const Rd normal = FK.T.N(ifac);
  const R meas    = FK.T.mesureBord(ifac);
  const R measK   = FK.getMeasure();
  const R h       = FK.T.lenEdge(ifac);//meas;

  for(int l=0; l<VF.size();++l) {
    int lastop = getLastop(0,VF[l].dv);
    double coef = BaseProblem<M>::computeCoef(VF[l],h,meas,measK,the_domain) * VF[l].getCoef(normal) ;

    const int kf = (VF[l].domu == 0)? k : kn;
    const int kv = (VF[l].domv == 0)? k : kn;

    const FElement& FKv(Vhv[kv]);
    this->initIndex(FKv);
    int kback  = Vh->idxElementInBackMesh(kv);
    int kfback = Vh->idxElementInBackMesh(kf);

    RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
    What_d Fop = Fwhatd(lastop);

    for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {

      QuadraturePoint ip(qfb[ipq]); // integration point
      const Rd mip = K(K.toKref((RdHatBord)ip, ifac));
      const R Cint = meas * ip.getWeight();

      FKv.BF(Fop,FKv.T.toKref(mip), fv); // need point in local reference element
      double val =
                  VF[l].fx_backMesh_U(kfback,the_domain, mip, normal)
                  *
                  VF[l].fx_backMesh_V(kback ,the_domain, mip, normal);
      double Cst = Cint * VF[l].c  * val * coef;
      for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
        (*this)(FKv.loc2glb(i)) += Cst*fv(i,VF[l].cv,VF[l].dv);
      }
    }
  }
  this->resetIndex();
}

/*
     ADD STRONG BC
*/
template<typename M>
void BaseProblem<M>::addStrongBC(std::list<ExpressionFunFEM<typename typeMesh<Rd::d>::Mesh>> gh, list<int> label) {

  typedef typename Mesh::BorderElement BorderElement;
  bool all_label = (label.size() == 0);

  for( int ifac = Vh->first_boundary_element(); ifac < Vh->last_boundary_element(); ifac+=Vh->next_boundary_element()) {
    const BorderElement & face(Vh->Th.be(ifac)); // The face
    if(util::contain(label, face.lab) || all_label) {
      for(auto it = gh.begin(); it != gh.end(); ++it) {
        setElementStrongBC(ifac, *it);
      }
    }
  }
}

template<typename M>
void BaseProblem<M>::addStrongBC(const ExpressionVirtual& gh, list<int> label) {

  typedef typename Mesh::BorderElement BorderElement;
  bool all_label = (label.size() == 0);

  for( int ifac = Vh->first_boundary_element(); ifac < Vh->last_boundary_element(); ifac+=Vh->next_boundary_element()) {
    const BorderElement & face(Vh->Th.be(ifac)); // The face
    if(util::contain(label, face.lab) || all_label) {
      setElementStrongBC(ifac, gh);
    }
  }
}

template<typename M>
void BaseProblem<M>::setElementStrongBC(int ifac, const ExpressionVirtual& gh) {

  typedef typename Mesh::BorderElement BorderElement;
  typedef typename QFB::QuadraturePoint QuadraturePoint;
  typedef typename FElement::RdHatBord RdHatBord;
  typedef typename Mesh::Element Element;


  const BorderElement & BE(Vh->Th.be(ifac)); // The face
  int ifaceK; // index of face of triangle corresp to edge (0,1,2)
  const int kb = Vh->Th.BoundaryElement(ifac, ifaceK); // index of element (triangle), ifaceK gets modified inside
  const int k = Vh->idxElementFromBackMesh(kb, 0);
  const FElement& FK((*Vh)[k]);


  for(int i = FK.dfcbegin(gh.cu); i < FK.dfcend(gh.cu); ++i) {

    // std::cout << " dof id local " << i << std::endl;
    // std::cout << " df on What " << FK.DFOnWhat(i) << std::endl;
    // std::cout << " id face loc" << ifaceK << std::endl;

    if(Element::onWhatBorder[ifaceK][FK.DFOnWhat(i)]) {
      long df = FK.loc2glb(i);

      (*this)(df,df) = tgv         ;
      (*this)(df)    = tgv * gh(df);

    }
  }
}



#include "baseProblemTime.hpp"
#include "baseSurfProblem.hpp"





//
// template<typename M>
// void BaseProblem<M>::addElementMatEdge(const ListItemVF<Rd::d>& VF, const int k) {
//
//   typedef typename QFB::QuadraturePoint QuadraturePoint;
//   typedef typename FElement::RdHatBord RdHatBord;
//   typedef typename Mesh::Element Element;
//
//   const FESpace& Sh =(VF[0].fespaceU)? *VF[0].fespaceU : *Vh;
//   this->initIndex(VF[0].fespaceU, VF[0].fespaceV);
//
//   const FElement& FK(Sh[k]);
//   const Element & K = FK.T;                  // the triangle
//   double measK = FK.getMeasure();
//   for(int ifac = 0; ifac < Element::nea; ++ifac) {    //loop over the edges / faces
//     int ifacn = ifac;
//     int kn = Sh.Th.ElementAdj(k,ifacn);
//     if(kn == -1) continue;
//     if(k > kn) continue;      // only compute one time
//
//     const R h = FK.T.lenEdge(ifac);
//
//     const FElement & FKn(Sh[kn]);                     // the neighboor finite element
//     Rd normal = FK.T.N(ifac);
//
//     const R meas = FK.T.mesureBord( ifac);
//
//     for(int l=0; l<VF.size();++l) {
//
//       assert(VF[l].fespaceU == VF[l].fespaceV);
//       int lastop = getLastop(VF[l].du, VF[l].dv);
//
//       double coef = computeCoef(VF[l],h,meas,measK) * VF[l].getCoef(normal) ;
//
//       const int ku = (VF[l].domu == 0)? k : kn;
//       const int kv = (VF[l].domv == 0)? k : kn;
//       bool same = (ku == kv && VF[l].fespaceU==VF[l].fespaceV);
//       const FElement& FKu(Sh[ku]);
//       const FElement& FKv(Sh[kv]);
//       this->initIndex(FKu, FKv);
//
//
//       RNMK_ fv(this->databf,FKv.NbDoF(),FKv.N,lastop); //  the value for basic fonction
//       RNMK_ fu(this->databf+ (same ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value for basic fonction
//       What_d Fop = Fwhatd(lastop);
//
//       for(int ipq = 0; ipq < qfb.getNbrOfQuads(); ++ipq)  {
//
//         QuadraturePoint ip(qfb[ipq]); // integration point
//         const Rd ip_Kref = K.toKref((RdHatBord)ip, ifac);
//         const R Cint = meas * ip.getWeight();
//
//         FKu.BF(Fop,ip_Kref, fu); // need point in local reference element
//         if(!same) FKv.BF(Fop,ip_Kref, fv);
//
//         for(int i = FKv.dfcbegin(VF[l].cv); i < FKv.dfcend(VF[l].cv); ++i) {
//           for(int j = FKu.dfcbegin(VF[l].cu); j < FKu.dfcend(VF[l].cu); ++j) {
//             // (*this)(FKv.loc2glb(i),FKu.loc2glb(j))  +=  Cint * coef * VF[l].c * fu(j,VF[l].cu,VF[l].du)*fv(i,VF[l].cv,VF[l].dv);
//             this->addToLocalContribution(FKv.loc2glb(i),FKu.loc2glb(j)) += Cint * coef * VF[l].c * fu(j,VF[l].cu,VF[l].du)*fv(i,VF[l].cv,VF[l].dv);
//           }
//         }
//       }
//     }
//   }
//   this->resetIndex();
//   this->addLocalContribution();
// }
//




#endif
