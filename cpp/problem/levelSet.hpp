#ifndef LEVEL_SET_HPP
#define LEVEL_SET_HPP

#include "baseProblem.hpp"


// template<typename M>
// class GLevelSet : public BaseProblem<M> {//public ShapeOfLinProblem , public Solver  {
// public:
//
//   typedef M Mesh;
//   typedef GFESpace<Mesh> FESpace;
//   typedef typename FESpace::FElement FElement;
//   typedef typename Mesh::Element Element;
//   typedef typename FElement::Rd Rd;
//
//   // typedef GenericInterface<Mesh> Interface;
//   // typedef GLevelSet<M> LevelSet;
//   typedef FunFEM<Mesh> Fun;
//
//   list<int> label_strongBC;
//
//
//   GLevelSet(const FESpace &);
//   GLevelSet(const Fun& up, const Fun& Betap, const Fun& Beta, double dt);
//
//   void solve(const Fun& up, const Fun& Betap, const Fun& Beta, double dt);
//   virtual void assembly(const Fun&, const Fun&, const Fun&, double);
//
//   void solve(const Fun& up, const Fun& Betap, const Fun& Beta, double dt, const Interface& interface);
//
//   void setBoundary( const Interface& interface);
//   void setStrongBC(  std::list<int> lab) {
//     label_strongBC = lab;
//   }
//
//   template<typename T>
//   void test(T&& lambda){
//     lambda (R2(1,2));
//   }
//
// };
//
//
// template<typename M>
// GLevelSet<M>::GLevelSet(const FESpace & vh) : BaseProblem<M>(vh) {};
//
// template<typename M>
// GLevelSet<M>::GLevelSet(const Fun& up, const Fun& Betap, const Fun& Beta, double dt) :
// BaseProblem<M>(*up.Vh) {};
//
// template<typename M>
// void GLevelSet<M>::solve(const Fun& up, const Fun& Betap, const Fun& Beta, double dt) {
//   this->rhs = 0.;
//   assembly(up, Betap, Beta, dt);
//   R t0 = CPUtime();
//   Solver::solve(this->mat, this->rhs);
//   R t1 = CPUtime();
// }
//
// template<typename M>
// void GLevelSet<M>::solve(const Fun& up, const Fun& Betap, const Fun& Beta, double dt, const Interface& interface) {
//   this->rhs = 0.;
//   assembly(up, Betap, Beta, dt);
//   setBoundary(interface);
//   R t0 = CPUtime();
//   Solver::solve(this->mat, this->rhs);
//   R t1 = CPUtime();
// }
//
//
// template<typename M>
// void GLevelSet<M>::setBoundary( const Interface& interface) {
//   typedef typename Mesh::BorderElement BorderElement;
//
//   for( int ifac = this->Vh->Th.first_element(); ifac < this->Vh->Th.last_boundary_element(); ifac+=this->Vh->Th.next_element()) {
//     const BorderElement & face(this->Vh->Th.be(ifac)); // The face
//     if(contain(label_strongBC, face.lab)) {
//
//       int ifaceK; // index of face of triangle corresp to edge (0,1,2)
//       const int kb = this->Vh->Th.BoundaryElement(ifac, ifaceK); // index of element (triangle), ifaceK gets modified inside
//       const int k = this->Vh->idxElementFromBackMesh(kb, 0);
//       const FElement& FK((*(this->Vh))[k]);
//
//       for(int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
//         if(Element::onWhatBorder[ifaceK][FK.DFOnWhat(i)]) {
//
//           Rd P = FK.Pt(i);
//           double mind = interface.distance(P);
//           long df = FK.loc2glb(i);
//
//           (*this)(df,df) = this->tgv       ;
//           (*this)(df)    = this->tgv * mind;
//         }
//       }
//     }
//   }
// }


// class LevelSet2 : public GLevelSet<Mesh2> {
// public :
//   LevelSet2(const FESpace & vh) : LevelSet(vh) {};
//   LevelSet2(  const Fun& up, const Fun& Betap, const Fun& Beta, double dt) :
//   LevelSet(up, Betap, Beta, dt) {
//     solve(up, Betap, Beta, dt);
//   };
//   void assembly(const Fun&, const Fun&, const Fun&, double);
//
// };
//
//
// class LevelSet3 : public GLevelSet<Mesh3> {
// public :
//   LevelSet3(const FESpace & vh)  : LevelSet(vh) {};
//   LevelSet3(  const Fun& up, const Fun& Betap, const Fun& Beta, double dt) :
//   LevelSet(up, Betap, Beta, dt) {
//     solve(up, Betap, Beta, dt);
//   };
//   void assembly(const Fun&, const Fun&, const Fun&, double);
// };
//



namespace LevelSet{

  template<typename Mesh>
  void move_2D(const FunFEM<Mesh>& up, const FunFEM<Mesh>& Betap, const FunFEM<Mesh>& Beta, double dt, FunFEM<Mesh>& ls) {
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FESpace::Rd Rd;
    typedef typename FElement::QF QF;
    typedef typename QF::QuadraturePoint QuadraturePoint;

    const FESpace&  Vh(*up.Vh);

    FEM<Mesh> levelSet(Vh);
    KNMK<double> fu(Vh[0].NbDoF(),1,op_Dall); //  the value for basic fonction

    const QF& qf(levelSet.get_quadrature_formular_K());

    for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {

      const FElement& FK(Vh[k]);
      double h = FK.T.hElement();
      double meas = FK.getMeasure();

      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = FK.map(ip);
        double Cint = meas * ip.getWeight();

        FK.BF(ip, fu); // need point in local reference element
        double Bx = Beta.eval(k, mip, 0, op_id);
        double By = Beta.eval(k, mip, 1, op_id);
        double Bpx = Betap.eval(k, mip, 0, op_id);
        double Bpy = Betap.eval(k, mip, 1, op_id);
        double Up = up.eval(k, mip, 0, op_id);
        double dxup = up.eval(k, mip, 0, op_dx);
        double dyup = up.eval(k, mip, 0, op_dy);
        double normB = Bx*Bx + By*By;
        double Tsd =  2. / (sqrt(1./dt/dt + normB * 1. / h / h));


        for(int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
          for(int j = FK.dfcbegin(0); j < FK.dfcend(0); ++j) {
            levelSet(FK.loc2glb(i),FK.loc2glb(j)) +=  Cint * (
              ((1./dt)*fu(j,0,op_id) + 0.5*(Bx*fu(j,0,op_dx)+ By*fu(j,0,op_dy)))
              *
              (fu(i,0,op_id) + Tsd * (Bx*fu(i,0,op_dx) + By*fu(i,0,op_dy)))
            );
          }
          levelSet(FK.loc2glb(i)) +=
          Cint*(
            ((1./dt)*Up - 0.5*(Bpx*dxup + Bpy*dyup))
            *
            (fu(i,0,op_id) + Tsd * (Bx*fu(i,0,op_dx) + By*fu(i,0,op_dy)))
          );
        }


        //     if(label_strongBC.size() > 0){
        //       ExpressionFunFEM<Mesh2> Up(up, 0, op_id);
        //       this->addStrongBC(Up, label_strongBC);
        //     }
      }
    }


    levelSet.solve();
    ls.init(levelSet.rhs_);
  }


  template<typename Mesh>
  void move_3D(const FunFEM<Mesh>& up, const FunFEM<Mesh>& Betap, const FunFEM<Mesh>& Beta, double dt, FunFEM<Mesh>& ls) {
    typedef GFESpace<Mesh> FESpace;
    typedef typename FESpace::FElement FElement;
    typedef typename FESpace::Rd Rd;
    typedef typename FElement::QF QF;
    typedef typename QF::QuadraturePoint QuadraturePoint;

    const FESpace&  Vh(*up.Vh);

    FEM<Mesh> levelSet(Vh);
    KNMK<double> fu(Vh[0].NbDoF(),1,op_Dall); //  the value for basic fonction

    const QF& qf(levelSet.get_quadrature_formular_K());

    for(int k=Vh.first_element(); k<Vh.last_element(); k+= Vh.next_element()) {

      const FElement& FK(Vh[k]);
      double h = FK.T.hElement();
      double meas = FK.getMeasure();

      for(int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq)  {

        QuadraturePoint ip(qf[ipq]); // integration point
        Rd mip = FK.map(ip);
        double Cint = meas * ip.getWeight();

        FK.BF(ip, fu); // need point in local reference element
        double Bx = Beta.eval(k, mip, 0, op_id);
        double By = Beta.eval(k, mip, 1, op_id);
        double Bz = Beta.eval(k, mip, 2, op_id);
        double Bpx = Betap.eval(k, mip, 0, op_id);
        double Bpy = Betap.eval(k, mip, 1, op_id);
        double Bpz = Betap.eval(k, mip, 2, op_id);
        double Up = up.eval(k, mip, 0, op_id);
        double dxup = up.eval(k, mip, 0, op_dx);
        double dyup = up.eval(k, mip, 0, op_dy);
        double dzup = up.eval(k, mip, 0, op_dz);
        double normB = Bx*Bx + By*By + Bz*Bz;
        double Tsd =  2. / (sqrt(1./dt/dt + normB * 1. / h / h));


        for(int i = FK.dfcbegin(0); i < FK.dfcend(0); ++i) {
          for(int j = FK.dfcbegin(0); j < FK.dfcend(0); ++j) {
            levelSet(FK.loc2glb(i),FK.loc2glb(j)) +=  Cint * (
              ((1./dt)*fu(j,0,op_id) + 0.5*(Bx*fu(j,0,op_dx)+ By*fu(j,0,op_dy)+ Bz*fu(j,0,op_dz)))
              *
              (fu(i,0,op_id) + Tsd * (Bx*fu(i,0,op_dx) + By*fu(i,0,op_dy)+ Bz*fu(i,0,op_dz)))
            );
          }
          levelSet(FK.loc2glb(i)) +=
          Cint*(
            ((1./dt)*Up - 0.5*(Bpx*dxup + Bpy*dyup + Bpz*dzup))
            *
            (fu(i,0,op_id) + Tsd * (Bx*fu(i,0,op_dx) + By*fu(i,0,op_dy) + Bz*fu(i,0,op_dz)))
          );
        }


        //     if(label_strongBC.size() > 0){
        //       ExpressionFunFEM<Mesh2> Up(up, 0, op_id);
        //       this->addStrongBC(Up, label_strongBC);
        //     }
      }
    }


    levelSet.solve();
    ls.init(levelSet.rhs_);
  }

};

namespace LevelSet{

  // template<typename Mesh>
  // KN<double> move_2D(const FunFEM<Mesh>&, const FunFEM<Mesh>&, const FunFEM<Mesh>&, double, FunFEM<Mesh>&);
  // template<typename Fun>
  // KN<double> move_3D(const Fun&, const Fun&, const Fun&, double);


  void move(const FunFEM<Mesh2> &up, const FunFEM<Mesh2> &betap, const FunFEM<Mesh2> &beta, double dt, FunFEM<Mesh2> & ls) {
    move_2D<Mesh2>(up, betap, beta, dt, ls);
  }
  // KN<double> move(const FunFEM<MeshQ2> &, const FunFEM<MeshQ2> &, const FunFEM<MeshQ2> &, double);
  void move(const FunFEM<Mesh3> &up, const FunFEM<Mesh3> &betap, const FunFEM<Mesh3> &beta, double dt, FunFEM<Mesh3> & ls) {
    move_3D<Mesh3>(up, betap, beta, dt, ls);
  }
  void move(const FunFEM<MeshHexa> &up, const FunFEM<MeshHexa> &betap, const FunFEM<MeshHexa> &beta, double dt, FunFEM<MeshHexa> & ls) {
    move_3D<MeshHexa>(up, betap, beta, dt, ls);
  }
  // KN<double> move(const FunFEM<MeshQ3> &, const FunFEM<MeshQ3> &, const FunFEM<MeshQ3> &, double);

};


#endif
