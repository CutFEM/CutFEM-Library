#ifndef GENERIC_MAPPING_HPP_
#define GENERIC_MAPPING_HPP_

// #include "problem.hpp"


template<typename Mesh>
class Mapping {
public:
  typedef CutFESpace<Mesh> CutSpace;
  typedef typename CutSpace::FElement FElement;
  typedef typename Mesh::Element Element;
  typedef typename FElement::Rd Rd;
  typedef typename FElement::RdHat RdHat;
  typedef typename FElement::QF QF;
  // typedef GenericInterface<Mesh> Interface;
  typedef FunFEM<Mesh> Fun_h;
  // typedef typename Mesh::Partition Partition;


public:
  const CutSpace&            Vh;
  const ActiveMesh<Mesh>&    Th;

  KN<double> deformation_;
  KN<Ubyte>  ni_;

  const Fun_h& levelSet_Pk;
  const Fun_h& levelSet_P1;

  // const int itq;

// protected:
//   const double current_time_ = 0;

public:
  double max_deform_;

  // Mapping() {}
  Mapping(const CutSpace& VVh, const Fun_h& levelSetk, const Fun_h& levelSet1)
  : Vh(VVh), Th(Vh.cutTh),
  levelSet_Pk(levelSetk), levelSet_P1(levelSet1) {
    assert(Vh.N == Rd::d);
    assert(&(levelSetk.Vh->Th) == &(Vh.Th));

    deformation_.resize(Vh.NbDoF()); deformation_ = 0.;
    ni_.resize(Vh.nbNode); ni_ = static_cast<Ubyte>(0);

    make_mapping();
  };


private:

  void make_mapping();
  void searchCorrespondingPoint(const int k,
    const Rd & init_point, const R goal_val,
    const Rd & search_dir, Rd & final_point) const;

public:
  void computeInverseJacobian(int kb, Rd x, KNM_<double>& invJ) const {
    assert(invJ.N() == Rd::d && invJ.M() == Rd::d);
    assert(invJ.N()==2);
    invJ = 0.;
    int D[3] = {op_dx, op_dy, op_dz};
    int k = Vh.idxElementFromBackMesh(kb, 0);
    const FElement& FK(Vh[k]);
    KNMK<double> bf(FK.NbDoF(),Rd::d,op_dz+1);
    FK.BF(FK.T.toKref(x), bf);

    for(int d=0, dv=Rd::d-1;d<Rd::d;++d, dv--) {
      for( int ci = 0, civ = Rd::d-1; ci<Rd::d; ++ci, civ--) {
        R s = (civ == dv)? 1 : -1;
        for(int i=0, j=FK.dfcbegin(ci);j<FK.dfcend(ci);++j,++i) {
          R val = (deformation_[FK(j)] + FK.Pt(i)[ci]);
          invJ(civ, dv) += val * bf(j,ci,D[d]) * s;
        }
      }
    }
    R det = invJ(0,0)*invJ(1,1) - invJ(0,1)*invJ(1,0);

    invJ /= det;
  }
  void transform(const FElement& FK, KNMK_<double>& bf, const KNM_<double>& invJ) const {
    if( bf.K() == 1 ) return;
    int op[3] = {op_dx, op_dy, op_dz};
    KNMK<R> x(bf.N(), bf.M(), bf.K()); x = bf;
    for( int ci = 0; ci<bf.M(); ++ci) {
      for(int df=FK.dfcbegin(ci);df<FK.dfcend(ci);++df) {
        for(int i=0;i<invJ.N();++i){
          bf(df,ci,op[i]) = 0;
          for(int j=0;j<invJ.M();++j){
            bf(df,ci,op[i]) += invJ(i,j)*x(df,ci,op[j]);
          }
        }
      }
    }
  }
  Rd map(const int kb, const Rd x) const {

    int k = Vh.idxElementFromBackMesh(kb, 0);
    const FElement& FK(Vh[k]);
    KNMK<double> bf(FK.NbDoF(),Rd::d,op_id+1);
    FK.BF(Fop_D0,FK.T.toKref(x), bf );

    Rd Pt = x;
    for( int ci = 0; ci<Rd::d; ++ci) {
      for(int i=0, j=FK.dfcbegin(ci);j<FK.dfcend(ci);++j,++i) {
        Pt[ci] += (deformation_[FK(j)]) * bf(j,ci,op_id) ;
      }
    }
    return Pt;
  }

};


template<typename M>
void Mapping<M>::make_mapping() {
  int D[3] = {op_dx, op_dy, op_dz};

  const int nbNode = Vh.TFE[0]->NbPtforInterpolation;
  for(int k=0;k<Vh.NbElement();++k) {

    const FElement& FK(Vh[k]);
    max_deform_ = 0.1 * FK.T.lenEdge(0);                           // the triangle

    int kb = Vh.idxElementInBackMesh(k);
    if (!Th.isCut(k, 0)) continue;

    // The nodes are not deform, only edges
    for(int l=Rd::d+1; l<nbNode; ++l) {

      Rd mip = FK.Pt((int) l);
      Rd normal, final_point;
      R1 goal_val;

      for(int i=0;i<Rd::d;++i) normal[i] = levelSet_Pk.eval(kb, mip, 0, D[i]);
      goal_val = levelSet_P1.eval(kb, mip, 0, op_id);    // the goal value

      searchCorrespondingPoint(k, mip, goal_val, normal, final_point);
      Rd ref_dist = final_point - mip;
      const double ref_dist_size = ref_dist.norm();

      if ((max_deform_ >= 0.0) && (ref_dist_size > max_deform_))
      {
        ref_dist = ref_dist * (max_deform_ / ref_dist_size);
      }

      for(int d = 0; d < Rd::d; ++d) {
        int id = FK.dfcbegin(d) + l;
        deformation_[FK(id)] += ref_dist[d] ;
      }
      ni_[FK[l]] += 1;
    }
  }
  for(int i=0;i<ni_.size();++i) {
    if(ni_[i] == 0) continue;
    for(int nd = 0, j = Vh.FirstDFOfNode(i); nd < Rd::d; ++nd, ++j) {
      deformation_[j] /= ni_[i];
    }
  }
}
template<typename M>
void Mapping<M>::searchCorrespondingPoint(
  const int k,
  const Rd & init_point, const R goal_val,
  const Rd & search_dir, Rd & final_point) const {

    int D[3] = {op_dx, op_dy, op_dz};
    Rd curr_ip = init_point;
    int it = 0;
    const R tol = 1e-14;
    int iterMax = 20;
    const FElement &FK(Vh[k]);
    int kb = Vh.idxElementInBackMesh(k);
    for ( it = 0 ; it < iterMax ; ++it) {
      R curr_val;
      Rd curr_grad;

      curr_val   = levelSet_Pk.eval(kb, curr_ip, 0, op_id);
      for(int i=0;i<Rd::d;++i) curr_grad[i] = levelSet_Pk.eval(kb, curr_ip, 0, D[i]);

      const R curr_def = goal_val - curr_val;
      const R dphin = (curr_grad, search_dir);
      if( fabs(curr_def) < tol ) break;
      curr_ip += (curr_def / dphin) * search_dir;

    }
    if( it >= iterMax-1) {
      final_point = init_point;
    }
    else {
      final_point = curr_ip;
    }

  }



  typedef Mapping<Mesh2> Mapping2;
  typedef Mapping<Mesh3> Mapping3;

  // class Mapping2 : public Mapping< Mesh2> {
  //
  //   typedef Mapping< Mesh2> GMapping;
  // public:
  //
  //   Mapping2() : GMapping() {}
  //   Mapping2(const FESpace& fh, Fun_h& lset) : GMapping(fh, lset) {}
  //
  //   R errorInfty(const Interface2& interface, R (*fun)(const R2, int i));
  //
  // };
  // class Mapping3 : public Mapping< Mesh3 > {
  //
  //   typedef Mapping< Mesh3> GMapping;
  // public:
  //   Mapping3() : GMapping() {}
  //   Mapping3(const FESpace& fh, Fun_h& lset) : GMapping(fh, lset) {}
  //
  // };
  //
  //
  // template<class mesh>
  // struct DataMapping
  // {
  //   static  Mapping<mesh> & Id;
  // };


  #endif
