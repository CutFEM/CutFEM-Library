#ifndef GENERIC_MAPPING_HPP_
#define GENERIC_MAPPING_HPP_

#include "problem.hpp"


template<typename M>
class GenericMapping {
public:
  typedef M Mesh;
  typedef GFESpace<Mesh> FESpace;
  typedef typename FESpace::FElement FElement;
  typedef typename Mesh::Element Element;
  typedef typename FElement::Rd Rd;
  typedef typename FElement::RdHat RdHat;
  typedef typename FElement::QF  QF;
  typedef GenericInterface<Mesh> Interface;
  typedef FunFEM<Mesh> Fun_h;
  typedef typename Mesh::Partition Partition;


public:
  KN<double> deformation_;
  KN<Ubyte>  ni_;

  FESpace const * Vh = nullptr;

protected:
  const R current_time = 0;

public:
  R max_deform;

  GenericMapping() {}
  GenericMapping(const FESpace& fh, const Fun_h& levelSet)
  : Vh(&fh) {
    assert(Vh->N == Rd::d);
    max_deform = 0.1 * fh[0].T.lenEdge(0);
    make_mapping(levelSet);
  };

  void init(const FESpace& fh, const Fun_h& levelSet) {
    Vh = &fh;
    assert(Vh->N == Rd::d);
    max_deform = 0.1 * fh[0].T.lenEdge(0);
    make_mapping(levelSet);
  }
private:

  void make_mapping(const  Fun_h& levelSet );
  void searchCorrespondingPoint(
    const Fun_h& levelSet, const int k,
    const Rd & init_point, const R goal_val,
    const Rd & search_dir, Rd & final_point,
    const Fun_h& ui) const;

  R evalP1 (const int k, const Rd& x, R* loc_ls) {
    KNMK<double> bf(Element::nv,1,1);
    const FElement& FK((*Vh)[k]);
    DataFE<Mesh>::P1.FB(Fop_D0,FK.T, FK.T.toKref(x), bf);

    R val = 0.;
    for(int j=0;j<Element::nv;++j) val += loc_ls[j]*bf(j,0,op_id);
    return val;
  }

public:
  virtual int idxElementFromBackMesh(int kb) const {
    return Vh->idxElementFromBackMesh(kb);
  }


  virtual void computeInverseJacobian(int k, Rd x, KNM_<double>& invJ) const {
    assert(invJ.N() == Rd::d && invJ.M() == Rd::d);
    assert(invJ.N()==2);
    invJ = 0.;
    int D[3] = {op_dx, op_dy, op_dz};
    const FElement& FK((*Vh)[k]);
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

  virtual void transform(const FElement& FK, KNMK_<double>& bf, const KNM_<double>& invJ) const {
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

  virtual Rd transform(const int k, const Rd x) const {

    const FElement& FK((*Vh)[k]);
    KNMK<double> bf(FK.NbDoF(),Rd::d,op_id+1);
    FK.BF(Fop_D0,FK.T.toKref(x), bf );

    Rd Pt;
    for( int ci = 0; ci<Rd::d; ++ci) {
      for(int i=0, j=FK.dfcbegin(ci);j<FK.dfcend(ci);++j,++i) {
        R val = (deformation_[FK(j)] + FK.Pt(i)[ci]);
        Pt[ci] += val * bf(j,ci,op_id);
      }
    }
    return Pt;
  }

};


  template<typename M>
  void GenericMapping<M>::make_mapping(const Fun_h& levelSet) {
    int D[3] = {op_dx, op_dy, op_dz};

    const int nbNode = (*Vh)[0].tfe->NbPtforInterpolation;
    const int nbDoF = (*Vh)[0].NbDoF();
    deformation_.resize(Vh->NbDoF()); deformation_ = 0.;
    ni_.resize(Vh->nbNode); ni_ = static_cast<Ubyte>(0);

    R loc_lset[levelSet.size(0)];//nbDoF];
    for(int k=0;k<Vh->NbElement();++k) {

      const FElement &FK((*Vh)[k]);
      const Element & K = FK.T;                            // the triangle

      int kb = Vh->idxElementInBackMesh(k);
      int kls = levelSet.idxElementFromBackMesh(kb);
      levelSet.eval(loc_lset, kls);
      const Partition& cutK =  Partition(FK.T, loc_lset);
      if (cutK.isnot_cut()) continue;

      // The nodes are not deform, only edges
      // for(int l=Rd::d+1; l<nbNode; ++l) {
      for(int l=0; l<nbNode; ++l) {

        Rd mip = FK.Pt((int) l);
        Rd normal, final_point;
        R1 goal_val;


        for(int i=0;i<Rd::d;++i) normal[i] = levelSet.eval(kls, mip, 0, D[i]);
        goal_val = this->evalP1(k, mip, loc_lset);    // the goal value

        searchCorrespondingPoint(levelSet, k, mip, goal_val, normal, final_point, levelSet);
        Rd ref_dist = final_point - mip;
        const double ref_dist_size = ref_dist.norm();


        if ((max_deform >= 0.0) && (ref_dist_size > max_deform))
        {
          ref_dist = ref_dist * (max_deform / ref_dist_size);
        }

        for(int d = 0; d < Rd::d; ++d) {
          int id = FK.dfcbegin(d) + l;
          deformation_(FK(id)) += ref_dist[d] ;
        }
        ni_( FK[l] ) += 1;
      }
    }
    for(int i=0;i<Vh->nbNode;++i) {
      if(ni_(i) == 0) continue;
      for(int nd = 0, j = Vh->FirstDFOfNode(i); nd < Rd::d; ++nd, ++j) {
        deformation_(j) /= ni_(i);
      }
    }
  }

  template<typename M>
  void GenericMapping<M>::searchCorrespondingPoint(
    const Fun_h& levelSet, const int k,
    const Rd & init_point, const R goal_val,
    const Rd & search_dir, Rd & final_point,
    const Fun_h& ui) const {

    int D[3] = {op_dx, op_dy, op_dz};
    Rd curr_ip = init_point;
    int it = 0;
    const R tol = 1e-14;
    int iterMax = 20;

    const FElement &FK((*Vh)[k]);
    int kb = Vh->idxElementInBackMesh(k);
    int kls = levelSet.idxElementFromBackMesh(kb);
    for ( it = 0 ; it < iterMax ; ++it) {
      R curr_val;
      Rd curr_grad;

      curr_val   = levelSet.eval(kls, curr_ip);
      for(int i=0;i<Rd::d;++i) curr_grad[i] = levelSet.eval(kls, curr_ip, 0, D[i]);

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




  class Mapping2 : public GenericMapping< Mesh2> {

    typedef GenericMapping< Mesh2> GMapping;
  public:

    Mapping2() : GMapping() {}
    Mapping2(const FESpace& fh, Fun_h& lset) : GMapping(fh, lset) {}

    R errorInfty(const Interface2& interface, R (*fun)(const R2, int i));

  };




  class Mapping3 : public GenericMapping< Mesh3 > {

    typedef GenericMapping< Mesh3> GMapping;
  public:
    Mapping3() : GMapping() {}
    Mapping3(const FESpace& fh, Fun_h& lset) : GMapping(fh, lset) {}

  };



  template<class mesh>
  struct DataMapping
  {
    static  GenericMapping<mesh> & Id;
  };


  #endif
