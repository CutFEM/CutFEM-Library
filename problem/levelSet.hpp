#ifndef LEVEL_SET_HPP
#define LEVEL_SET_HPP

#include "baseProblem.hpp"


template<typename M>
class GLevelSet : public BaseProblem<M> {//public ShapeOfLinProblem , public Solver  {
public:

  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename Mesh::Element Element;
  typedef typename FElement::Rd Rd;

  typedef GenericInterface<Mesh> Interface;
  typedef GLevelSet<M> LevelSet;
  typedef FunFEM<Mesh> Fun;

  list<int> label_strongBC;


  GLevelSet(const FESpace &);
  GLevelSet(const Fun& up, const Fun& Betap, const Fun& Beta, double dt);

  void solve(const Fun& up, const Fun& Betap, const Fun& Beta, double dt);
  virtual void assembly(const Fun&, const Fun&, const Fun&, double);

  void solve(const Fun& up, const Fun& Betap, const Fun& Beta, double dt, const Interface& interface);

  void setBoundary( const Interface& interface);
  void setStrongBC(  std::list<int> lab) {
    label_strongBC = lab;
  }

  template<typename T>
  void test(T&& lambda){
    lambda (R2(1,2));
  }

};


template<typename M>
GLevelSet<M>::GLevelSet(const FESpace & vh) : BaseProblem<M>(vh) {};

template<typename M>
GLevelSet<M>::GLevelSet(const Fun& up, const Fun& Betap, const Fun& Beta, double dt) :
BaseProblem<M>(*up.Vh) {};

template<typename M>
void GLevelSet<M>::solve(const Fun& up, const Fun& Betap, const Fun& Beta, double dt) {
  this->rhs = 0.;
  assembly(up, Betap, Beta, dt);
  R t0 = CPUtime();
  Solver::solve(this->mat, this->rhs);
  R t1 = CPUtime();
}

template<typename M>
void GLevelSet<M>::solve(const Fun& up, const Fun& Betap, const Fun& Beta, double dt, const Interface& interface) {
  this->rhs = 0.;
  assembly(up, Betap, Beta, dt);
  setBoundary(interface);
  R t0 = CPUtime();
  Solver::solve(this->mat, this->rhs);
  R t1 = CPUtime();
}


template<typename M>
void GLevelSet<M>::setBoundary( const Interface& interface) {
  typedef typename Mesh::BorderElement BorderElement;

  for( int ifac = this->Vh->Th.first_element(); ifac < this->Vh->Th.last_boundary_element(); ifac+=this->Vh->Th.next_element()) {
    const BorderElement & face(this->Vh->Th.be(ifac)); // The face
    if(contain(label_strongBC, face.lab)) {

      int ifaceK; // index of face of triangle corresp to edge (0,1,2)
      const int kb = this->Vh->Th.BoundaryElement(ifac, ifaceK); // index of element (triangle), ifaceK gets modified inside
      const int k = this->Vh->idxElementFromBackMesh(kb, 0);
      const FElement& FK((*(this->Vh))[k]);

      for(int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
        if(Element::onWhatBorder[ifaceK][FK.DFOnWhat(i)]) {

          Rd P = FK.Pt(i);
          double mind = interface.distance(P);
          long df = FK.loc2glb(i);

          (*this)(df,df) = this->tgv       ;
          (*this)(df)    = this->tgv * mind;
        }
      }
    }
  }
}


class LevelSet2 : public GLevelSet<Mesh2> {
public :
  LevelSet2(const FESpace & vh) : LevelSet(vh) {};
  LevelSet2(  const Fun& up, const Fun& Betap, const Fun& Beta, double dt) :
  LevelSet(up, Betap, Beta, dt) {
    solve(up, Betap, Beta, dt);
  };
  void assembly(const Fun&, const Fun&, const Fun&, double);

};


class LevelSet3 : public GLevelSet<Mesh3> {
public :
  LevelSet3(const FESpace & vh)  : LevelSet(vh) {};
  LevelSet3(  const Fun& up, const Fun& Betap, const Fun& Beta, double dt) :
  LevelSet(up, Betap, Beta, dt) {
    solve(up, Betap, Beta, dt);
  };
  void assembly(const Fun&, const Fun&, const Fun&, double);
};
// 
// namespace LevelSet{
//
//   template<typename Fun>
//   KN<double> move_2D(const Fun&, const Fun&, const Fun&, double);
//   template<typename Fun>
//   KN<double> move_3D(const Fun&, const Fun&, const Fun&, double);
//
//
//   KN<double> move(const FunFEM<MeshT2> &, const FunFEM<MeshT2> &, const FunFEM<MeshT2> &, double);
//   KN<double> move(const FunFEM<MeshQ2> &, const FunFEM<MeshQ2> &, const FunFEM<MeshQ2> &, double);
//   KN<double> move(const FunFEM<MeshT3> &, const FunFEM<MeshT3> &, const FunFEM<MeshT3> &, double);
//   KN<double> move(const FunFEM<MeshQ3> &, const FunFEM<MeshQ3> &, const FunFEM<MeshQ3> &, double);
//
// };




#endif
